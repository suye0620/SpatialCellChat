# Re-run cell-level inference from the full Slide-seq pre-inference checkpoint.
setwd(local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
options(future.globals.maxSize = 20 * 1024^3)
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp)
  library(ggplot2); library(colorspace); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
out.dir <- "tests_dev/figures/P0/slideseq"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out.dir, "spatialC2CDistPlot.pdf")
status <- "PASS"
detail <- tryCatch({
  checkpoint <- "tests_dev/test_data/slideseq_p0_fresh_pre_filter.rds"
  if (file.exists(checkpoint)) {
    chat <- readRDS(checkpoint)
    .sc_c2c_inference_cache(chat)
    cat("Using verified fresh pre-filter checkpoint\n")
  } else {
    obj <- readRDS("tests_dev/test_data/MBM05Chat_pre_inference.rds")
    chat <- createSpatialCellChat(
      object = obj@assay$norm, meta = obj@meta, group.by = "rctd_cell_type",
      input.assay = "norm", datatype = "spatial", coordinates = obj@images$coordinates,
      spatial.factors = obj@images$spatial.factors
    )
    chat@DB <- obj@DB
    chat@LR <- obj@LR
    assay(chat, "signaling") <- obj@assay$signaling
    rm(obj)
    gc()
    chat <- computeCommunProb(chat, n = 1, interaction.range = 250,
                              scale.distance = 9, nthreads = 4L, verbose = TRUE)
    .sc_c2c_inference_cache(chat)
    saveRDS(chat, checkpoint, compress = FALSE)
    cat("Saved fresh pre-filter checkpoint:", checkpoint, "\n")
  }
  chat <- filterProbability(chat)
  plot <- spatialC2CDistPlot(chat, split.by = "annotation", mode = "communication")
  ggplot2::ggsave(out, plot, width = 7, height = 5.5, dpi = 150)
  if (!file.exists(out) || file.info(out)$size <= 1000L) stop("PDF missing or empty")
  paste0("fresh inference; cells=", nrow(chat@meta), "; LR layers=",
         dim(chat@net$cell$prob)[3L], "; n=1; interaction.range=250; scale.distance=9")
}, error = function(e) { status <<- "FAIL"; conditionMessage(e) })
result <- data.frame(figure = "P0/slideseq/spatialC2CDistPlot", status = status,
                     detail = detail, file = out)
write.csv(result, "tests_dev/figures/summary-slideseq-P0-spatialC2CDistPlot-fresh.csv", row.names = FALSE)
print(result)
if (status != "PASS") quit(status = 1L)
