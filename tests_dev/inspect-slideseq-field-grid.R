# Inspect Slide-seq field display sampling without changing the computed field.
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
source("tests_dev/fixtures/make_slideseq_window.R")

chat <- make_slideseq_chat()
signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
chat <- computeCommunField(chat, signaling.name = signal)
field <- as.matrix(chat@netP$cell$field$outgoing[[signal]])
field.nonzero <- sum(rowSums(field^2) > 0)
timing <- system.time(
  plot <- netVisual_CommunFieldGrid(chat, signaling = signal, pattern = "outgoing",
                                   point.size = 0.22, image.alpha = 0.4)
)
uv <- plot$layers[[2L]]$data
base <- ggplot2::ggplot_build(plot)$data[[2L]]
default.skip <- plot$layers[[2L]]$stat_params$skip.x
denser <- plot
denser$layers[[2L]]$stat_params$skip.x <- 0L
denser$layers[[2L]]$stat_params$skip.y <- 0L
no.skip <- ggplot2::ggplot_build(denser)$data[[2L]]
out.dir <- Sys.getenv("SC_GRID_DIAGNOSTIC_DIR", "tests_dev/figures/P3/slideseq")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out.dir, "netVisual_CommunFieldGrid-no-skip.pdf")
ggplot2::ggsave(out, denser, width = 7, height = 5.5)
summary <- data.frame(dataset = "slideseq", signal = signal, cells = nrow(field),
                      nonzero_field_cells = field.nonzero, grid_points = nrow(uv),
                      nonzero_grid_points = sum(uv$u^2 + uv$v^2 > 0, na.rm = TRUE),
                      plot_constructor_seconds = unname(timing[["elapsed"]]),
                      default_skip = default.skip,
                      default_sampled_rows = nrow(base), no_skip_sampled_rows = nrow(no.skip),
                      output = out)
write.csv(summary, file.path(out.dir, "grid-sampling-diagnostics.csv"), row.names = FALSE)
print(summary)
