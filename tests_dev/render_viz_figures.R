# Figure runner for the visualization migration.
#
# Renders each registered visualization function on each dataset and saves one
# PDF per (dataset, function) under tests_dev/figures/<phase>/<dataset>/.
# A failing figure does not abort the run; the summary decides the exit code.
#
# Usage:
#   Rscript tests_dev/render_viz_figures.R                     # all datasets, all registered figures
#   Rscript tests_dev/render_viz_figures.R --dataset=slideseq  # one dataset
#   Rscript tests_dev/render_viz_figures.R --phase=P0          # one phase
#   Rscript tests_dev/render_viz_figures.R --only=spatialC2CDistPlot

setwd(local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp); library(ggplot2)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/visualization.R")
source("R/visualization_widget.R")
source("R/modeling.R")
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
DATASET.filter <- get_arg("dataset", "all")
PHASE.filter <- get_arg("phase", "all")
ONLY.filter <- get_arg("only", NULL)

for (f in list.files("tests_dev/fixtures", pattern = "[.]R$", full.names = TRUE))
  source(f)

FIG.root <- "tests_dev/figures"
dir.create(FIG.root, showWarnings = FALSE, recursive = TRUE)

# dataset chat builders (cached; built on first use)
load_chat <- function(dataset) {
  switch(dataset,
         simulated = make_simulated_chat(),
         slideseq = make_slideseq_chat(),   # full 29,536-cell puck
         visium = make_visium_chat(),
         stop("unknown dataset: ", dataset, call. = FALSE))
}
DATASETS <- c("simulated", "slideseq", "visium")

# figure registry: phase, dataset, figure name, expression (evaluated with `chat` bound)
FIGURES <- list(
  list(phase = "P0", dataset = "simulated", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", enriched.only = TRUE))),
  list(phase = "P0", dataset = "simulated", name = "spatialC2CDistPlot-geometry",
       expr = quote(spatialC2CDistPlot(chat, enriched.only = FALSE))),
  list(phase = "P0", dataset = "simulated", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat, image = FALSE))),
  list(phase = "P0", dataset = "slideseq", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", enriched.only = TRUE))),
  list(phase = "P0", dataset = "slideseq", name = "spatialC2CDistPlot-geometry",
       expr = quote(spatialC2CDistPlot(chat, enriched.only = FALSE))),
  list(phase = "P0", dataset = "slideseq", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat, image = FALSE))),
  list(phase = "P0", dataset = "visium", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", enriched.only = TRUE))),
  list(phase = "P0", dataset = "visium", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat)))
)

results <- list()
for (spec in FIGURES) {
  if (!identical(DATASET.filter, "all") && !identical(DATASET.filter, spec$dataset)) next
  if (!identical(PHASE.filter, "all") && !identical(PHASE.filter, spec$phase)) next
  if (!is.null(ONLY.filter) && !identical(ONLY.filter, spec$name)) next

  tag <- paste(spec$phase, spec$dataset, spec$name, sep = "/")
  out.dir <- file.path(FIG.root, spec$phase, spec$dataset)
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  out.pdf <- file.path(out.dir, paste0(spec$name, ".pdf"))
  cat("== ", tag, "\n", sep = "")

  status <- "PASS"
  detail <- ""
  t0 <- Sys.time()
  chat <- tryCatch(load_chat(spec$dataset),
                   error = function(e) { status <<- "FIXTURE-FAIL"; detail <<- conditionMessage(e); NULL })
  p <- NULL
  if (!is.null(chat)) {
    p <- tryCatch({
      env <- list(chat = chat)
      eval(spec$expr, envir = env)
    }, error = function(e) {
      status <<- "FAIL"
      detail <<- conditionMessage(e)
      NULL
    })
  }
  if (inherits(p, "patchwork")) p <- patchwork::wrap_plots(p, nrow = 1)
  if (inherits(p, "ggplot")) {
    tryCatch({
      ggplot2::ggsave(out.pdf, p, width = 7, height = 5.5, dpi = 150)
      cat("   saved:", out.pdf, "\n")
    }, error = function(e) {
      status <<- "SAVE-FAIL"; detail <<- conditionMessage(e)
    })
  } else if (status == "PASS") {
    status <- "FAIL"; detail <- "expression did not return a ggplot/patchwork object"
  }
  cat("   -> ", status, if (nzchar(detail)) paste0(" | ", detail) else "", " (",
      round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s)\n", sep = "")
  results[[tag]] <- list(status = status, detail = detail, pdf = out.pdf)
}

df <- do.call(rbind, lapply(names(results), function(k) {
  data.frame(figure = k, status = results[[k]]$status, detail = results[[k]]$detail)
}))
cat("\n==== summary ====\n")
print(df, right = FALSE)
fails <- sum(df$status != "PASS")
cat("total:", nrow(df), "| pass:", nrow(df) - fails, "| fail:", fails, "\n")
if (fails > 0) quit(status = 1)
