# Render the full Visium chord with explicitly adjusted presentation parameters.
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
  library(methods); library(Matrix); library(cli); library(dplyr); library(ggplot2)
  library(colorspace); library(circlize); library(ComplexHeatmap); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
source("tests_dev/fixtures/make_visium_pp1.R")
chat <- make_visium_chat()
out <- "tests_dev/figures/P2/visium/netVisual_chord_gene-adjusted.pdf"
grDevices::pdf(out, width = 20, height = 20)
tryCatch(netVisual_chord_gene(chat, signaling = chat@netP$pathways[1],
                             show.legend = FALSE, lab.cex = 0.4, small.gap = 0.25),
         finally = grDevices::dev.off())
if (!file.exists(out) || file.info(out)$size <= 1000L) stop("PDF missing or empty")
cat("Rendered full pathway", chat@netP$pathways[1], "to", out,
    "with 20x20 inch page, lab.cex=0.4, small.gap=0.25; visual review required.\n")
