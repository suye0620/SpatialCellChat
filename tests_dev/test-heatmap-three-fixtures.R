# Re-render the single-object heatmap on all three canonical fixtures.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(ggplot2)
  library(circlize); library(ComplexHeatmap); library(igraph); library(colorspace); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R"); source("R/database.R"); source("R/utilities.R")
source("R/spatial.R"); source("R/modeling.R"); source("R/analysis.R"); source("R/visualization.R")
source("tests_dev/fixtures/make_simulated_chat.R"); source("tests_dev/fixtures/make_visium_pp1.R")
source("tests_dev/fixtures/make_slideseq_window.R")
results <- lapply(c("simulated", "visium", "slideseq"), function(dataset) {
  chat <- switch(dataset, simulated = make_simulated_chat(), visium = make_visium_chat(), slideseq = make_slideseq_chat())
  expected <- as.matrix(chat@net$group$count)
  out.dir <- file.path("tests_dev/figures/P2", dataset, "heatmap-cancellation-check")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(out.dir, "netVisual_heatmap.pdf")
  grDevices::pdf(file, width = 8, height = 6)
  result <- tryCatch(netVisual_heatmap(chat, measure = "count"), finally = grDevices::dev.off())
  plotted <- result@ht_list[[1L]]@matrix
  stopifnot(identical(dimnames(plotted), dimnames(expected)),
            all(plotted[expected != 0] == expected[expected != 0]),
            file.exists(file), file.info(file)$size > 1000)
  data.frame(dataset = dataset, cells = nrow(chat@meta), nonzero = sum(expected != 0),
             status = "PDF-WRITTEN", file = file)
})
summary <- do.call(rbind, results)
write.csv(summary, "tests_dev/figures/P2/heatmap-cancellation-regression-2026-11-10.csv", row.names = FALSE)
print(summary)
