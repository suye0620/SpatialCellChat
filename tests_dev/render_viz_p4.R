# Render P4 htmlwidgets on real migration fixtures; browser validation is separate.
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
  library(ggplot2); library(colorspace); library(Rcpp)
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
if (identical(dataset, "slideseq")) options(future.globals.maxSize = 20 * 1024^3)
output.root <- Sys.getenv("SC_P4_OUT_ROOT", unset = "tests_dev/figures/P4")
if (!nzchar(output.root)) stop("SC_P4_OUT_ROOT must not be empty", call. = FALSE)
out.dir <- file.path(output.root, dataset)
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
gene <- rownames(assay(chat, "norm"))[1L]
stack.genes <- head(rownames(assay(chat, "norm")), 3L)
stack.layers <- t(as.matrix(assay(chat, "norm")[stack.genes, , drop = FALSE]))
rownames(stack.layers) <- colnames(assay(chat, "norm"))
colnames(stack.layers) <- stack.genes
ensure_gi <- function() {
  if (is.null(chat@netP$centr.cell)) {
    signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
    chat <<- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = FALSE,
                                           signaling.name = signal)
  }
  dimnames(chat@netP$centr.cell)[[3L]][1L]
}
specs <- list(
  plotly_spatialDimPlot_2d = function() plotly_spatialDimPlot(chat, method = "2d"),
  plotly_spatialDimPlot_3d = function() plotly_spatialDimPlot(chat, method = "3d"),
  plotly_spatialStackPlot = function() plotly_spatialStackPlot(
    chat, layers = stack.layers, title.name = "Spatial feature stack"
  ),
  plotly_spatialLRpairPlot = function() plotly_spatialLRpairPlot(chat, features = gene)[[1L]],
  spatialGiPlot_3d = function() {
    signal <- ensure_gi()
    spatialGiPlot(chat, signaling.name = signal, measure = "indeg", do.binary = FALSE,
                  method = "3d", return.object = FALSE)
  }
)
results <- lapply(names(specs), function(name) {
  cat("==", dataset, name, "\n")
  file <- file.path(out.dir, paste0(name, ".html"))
  message <- tryCatch({
    p <- plotly::plotly_build(specs[[name]]())
    ids <- unlist(lapply(p$x$data, function(x) x$customdata), use.names = FALSE)
    if (name == "spatialGiPlot_3d") {
      count <- sum(vapply(p$x$data[seq_len(min(2L, length(p$x$data)))],
                          function(x) length(x$x), integer(1)))
      if (count != nrow(chat@meta)) stop("Gi trace point count differs from input cell count")
    } else if (name == "plotly_spatialStackPlot") {
      data.traces <- Filter(function(x) !identical(x$hoverinfo, "none"), p$x$data)
      count <- sum(vapply(data.traces, function(x) length(x$x), integer(1)))
      if (length(data.traces) != 1L || count != nrow(chat@meta) * length(stack.genes))
        stop("stack trace count or point coverage differs from input cells and layers")
    } else if (!setequal(ids, rownames(chat@meta))) {
      stop("plot trace does not contain all input cells")
    }
    # The isolated recheck can share one dependency directory per dataset.
    libdir <- if (identical(output.root, "tests_dev/figures/P4")) paste0(name, "_files") else "lib"
    htmlwidgets::saveWidget(p, file, selfcontained = FALSE, libdir = libdir)
    if (!file.exists(file) || file.info(file)$size <= 1000L) stop("HTML missing or empty")
    ""
  }, error = function(e) conditionMessage(e))
  status <- if (nzchar(message)) "FAIL" else "HTML-WRITTEN"
  cat(status, message, "\n")
  data.frame(dataset = dataset, figure = name, cells = nrow(chat@meta),
             feature = if (name == "plotly_spatialLRpairPlot") gene else "",
             status = status, detail = message, file = file)
})
results <- do.call(rbind, results)
# htmlwidgets copies the JS assets, but not the accompanying license files.
license.paths <- c(
  "plotly-latest.min.js" = system.file("htmlwidgets", "lib", "plotlyjs", "LICENSE", package = "plotly"),
  "typedarray.min.js" = system.file("htmlwidgets", "lib", "typedarray", "LICENSE", package = "plotly")
)
if (any(!nzchar(license.paths))) stop("Plotly dependency license files are unavailable", call. = FALSE)
assets <- list.files(out.dir, recursive = TRUE, full.names = TRUE)
for (asset in assets[basename(assets) %in% names(license.paths)]) {
  license.file <- file.path(dirname(asset), "LICENSE")
  if (!file.copy(license.paths[[basename(asset)]], license.file, overwrite = TRUE))
    stop("Unable to copy dependency license to ", license.file, call. = FALSE)
}
write.csv(results, file.path(out.dir, "summary.csv"), row.names = FALSE)
print(results)
if (any(results$status == "FAIL")) quit(status = 1L)
