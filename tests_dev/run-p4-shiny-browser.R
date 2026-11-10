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
  library(ggplot2); library(colorspace); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
source("tests_dev/fixtures/make_simulated_chat.R")

args <- commandArgs(TRUE)
port <- if (length(args)) as.integer(args[[1L]]) else 18765L
chat <- make_simulated_chat()
app <- plotly_spatialLRpairPlot_shiny(
  chat,
  features = rownames(assay(chat, "norm"))[[1L]],
  port = port
)
shiny::runApp(app, host = "127.0.0.1", port = port,
              launch.browser = FALSE, quiet = TRUE)
