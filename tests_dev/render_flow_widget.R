# Render communication-flow widgets for full migration fixtures.
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
  library(methods); library(Matrix); library(cli); library(dplyr)
  library(ggplot2); library(colorspace); library(Rcpp); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
source("R/visualization_widget.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
for (f in list.files("tests_dev/fixtures", pattern = "[.]R$", full.names = TRUE)) source(f)
args <- commandArgs(TRUE)
dataset <- if (length(args)) args[[1L]] else "simulated"
chat <- switch(dataset, simulated = make_simulated_chat(), visium = make_visium_chat(),
               slideseq = make_slideseq_chat(), stop("unknown dataset"))
signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
out.dir <- file.path("tests_dev/figures/P3", dataset)
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
file <- file.path(out.dir, "netVisual_CommunFlowWidget.html")
result <- tryCatch({
  chat <- computeCommunField(chat, signaling.name = signal)
  widget <- netVisual_CommunFlowWidget(chat, signaling = signal)
  if (!inherits(widget, "htmlwidget")) stop("result is not an htmlwidget")
  htmlwidgets::saveWidget(widget, file, selfcontained = FALSE,
                          libdir = "netVisual_CommunFlowWidget_files")
  if (!file.exists(file) || file.info(file)$size <= 1000L) stop("HTML missing or empty")
  ""
}, error = function(e) conditionMessage(e))
status <- if (nzchar(result)) "FAIL" else "HTML-WRITTEN; BROWSER-NOT-TESTED"
summary <- data.frame(dataset = dataset, cells = nrow(chat@meta), signal = signal,
                      status = status, detail = result, file = file)
write.csv(summary, file.path(out.dir, "netVisual_CommunFlowWidget-summary.csv"), row.names = FALSE)
print(summary)
if (nzchar(result)) quit(status = 1L)
