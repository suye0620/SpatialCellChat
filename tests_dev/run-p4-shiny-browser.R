setwd(local({
  args <- commandArgs(FALSE)
  file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file)) dirname(dirname(normalizePath(file))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr)
  library(ggplot2); library(colorspace); library(Rcpp); library(sf)
  library(tidyr); library(purrr); library(spatstat.sparse)
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
source("tests_dev/fixtures/make_slideseq_window.R")

args <- commandArgs(TRUE)
port <- if (length(args)) as.integer(args[[1L]]) else 18765L
dataset <- if (length(args) >= 2L) args[[2L]] else "simulated"
chat <- switch(dataset,
  simulated = make_simulated_chat(),
  visium = make_visium_chat(),
  slideseq = make_slideseq_chat(),
  stop("dataset must be simulated, visium, or slideseq", call. = FALSE)
)
norm <- assay(chat, "norm")
expressed <- which(Matrix::rowSums(norm != 0) > 0)
if (!length(expressed)) stop("selected fixture has no expressed genes", call. = FALSE)
feature <- rownames(norm)[expressed[[1L]]]
app <- plotly_spatialLRpairPlot_shiny(
  chat,
  features = feature,
  port = port
)
shiny::runApp(app, host = "127.0.0.1", port = port,
              launch.browser = FALSE, quiet = TRUE)
