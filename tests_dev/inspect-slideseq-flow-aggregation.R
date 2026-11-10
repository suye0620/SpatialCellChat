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
source("tests_dev/fixtures/make_slideseq_window.R")
chat <- make_slideseq_chat()
signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
chat <- computeCommunField(chat, signaling.name = signal)
timing <- system.time(p <- netVisual_CommunFlow(chat, signaling = signal, pattern = "outgoing", point.size = 0.22, image.alpha = 0.4, res = 8))
out.dir <- "tests_dev/figures/P3/slideseq/aggregation-check"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(file.path(out.dir, "netVisual_CommunFlow-adjusted.pdf"), p, width = 7, height = 5.5, dpi = 150, limitsize = FALSE)
write.csv(data.frame(dataset = "slideseq", signal = signal, cells = nrow(chat@images$coordinates),
                     plot_constructor_seconds = unname(timing[["elapsed"]]),
                     parameters = "outgoing; point.size=0.22; image.alpha=0.4; res=8"),
          file.path(out.dir, "flow-aggregation-diagnostics.csv"), row.names = FALSE)
print(timing)
