setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(tidyr)
  library(purrr); library(sf); library(ggplot2); library(colorspace); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
source("tests_dev/fixtures/make_simulated_chat.R")
source("tests_dev/fixtures/make_visium_pp1.R")
results <- lapply(c("simulated", "visium"), function(dataset) {
  chat <- switch(dataset, simulated = make_simulated_chat(), visium = make_visium_chat())
  signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
  chat <- computeCommunField(chat, signaling.name = signal)
  plot <- netVisual_CommunFieldGrid(chat, signaling = signal, pattern = "outgoing", arrow.skip = 0)
  current <- ggplot_build(plot)$data[[2L]]
  stopifnot(all(current$alpha[current$mag == 0] == 0),
            all(current$alpha[current$mag > 0] == 0.9))
  out.dir <- file.path("tests_dev/figures/P3", dataset, "zero-vector-display-check")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(out.dir, "netVisual_CommunFieldGrid.pdf")
  ggsave(file, plot, width = 7, height = 5.5)
  plot$layers[[2L]]$mapping$alpha <- NULL
  plot$layers[[2L]]$aes_params$alpha <- 0.9
  original <- ggplot_build(plot)$data[[2L]]
  shared <- c("x", "y", "dx", "dy", "mag", "norm_mag", "size")
  stopifnot(identical(current[, shared], original[, shared]))
  data.frame(dataset = dataset, cells = nrow(chat@meta), signal = signal,
             zero_transparent = sum(current$mag == 0), positive_visible = sum(current$mag > 0),
             geometry_unchanged = TRUE, output = file)
})
summary <- do.call(rbind, results)
write.csv(summary, "tests_dev/figures/P3/zero-vector-display-fixtures-2026-11-10.csv", row.names = FALSE)
print(summary)
