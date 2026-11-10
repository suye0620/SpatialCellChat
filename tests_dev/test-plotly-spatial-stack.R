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

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
expect_error <- function(label, expr) {
  check(label, inherits(try(force(expr), silent = TRUE), "try-error"))
}

chat <- make_simulated_chat()
cells <- rownames(chat@images$coordinates)
ncell <- length(cells)
original <- chat
layers <- cbind(
  `Layer<&` = seq_len(ncell),
  reverse = -seq_len(ncell),
  half = seq_len(ncell) / 2
)
rownames(layers) <- rev(cells)
colnames(layers) <- c("Layer<&", "reverse", "half")
aligned.layers <- layers[cells, , drop = FALSE]
p <- plotly_spatialStackPlot(chat, layers, z.axis.space = 2)
traces <- Filter(function(trace) !inherits(trace, "plotly_colorbar"), p$x$data)
trace <- traces[[1L]]
check("builds one scatter3d trace for all layers",
      length(traces) == 1L && identical(trace$type, "scatter3d") &&
        length(trace$x) == ncell * ncol(layers))
check("aligns layer rows to coordinates by cell name",
      identical(as.character(trace$customdata), rep(cells, ncol(layers))) &&
        identical(as.numeric(trace$z), rep(c(0, 2, 4), each = ncell)) &&
        isTRUE(all.equal(as.numeric(trace$marker$color), as.numeric(aligned.layers))))
check("uses one colorbar and one global mixed-sign color domain",
      isTRUE(trace$marker$showscale) &&
        sum(vapply(p$x$data, inherits, logical(1), "plotly_colorbar")) == 1L &&
        trace$marker$cmin == -ncell &&
        trace$marker$cmax == ncell)
zaxis <- p$x$layout$scene$zaxis
check("z-axis ticks identify the supplied layers",
      identical(as.numeric(zaxis$tickvals), c(0, 2, 4)) &&
        identical(as.character(zaxis$ticktext), colnames(layers)))
check("hover text escapes layer names and preserves the numeric value",
      any(grepl("Layer&lt;&amp;", trace$text, fixed = TRUE)) &&
        any(grepl("cell:", trace$text, fixed = TRUE)) &&
        any(grepl("value:", trace$text, fixed = TRUE)))
check("does not mutate the source object", identical(chat, original))

partial.na <- layers
partial.na[1L, 1L] <- NA_real_
partial.plot <- plotly_spatialStackPlot(chat, partial.na)
check("omits partial NA values without dropping other layer points",
      length(partial.plot$x$data[[1L]]$x) == ncell * ncol(layers) - 1L)

single <- layers[, "half", drop = FALSE]
single.plot <- plotly_spatialStackPlot(chat, single)
single.traces <- Filter(function(trace) !inherits(trace, "plotly_colorbar"), single.plot$x$data)
check("supports a single named layer", length(single.traces) == 1L &&
        length(single.traces[[1L]]$x) == ncell &&
        identical(as.character(single.plot$x$layout$scene$zaxis$ticktext), "half"))

zeros <- layers[, 1L, drop = FALSE]
zeros[,] <- 0
zero.plot <- plotly_spatialStackPlot(chat, zeros)
check("keeps an all-zero layer on a finite color domain",
      zero.plot$x$data[[1L]]$marker$cmin == -1e-6 &&
        zero.plot$x$data[[1L]]$marker$cmax == 1e-6)

all.na <- layers
all.na[, 1L] <- NA_real_
expect_error("rejects an all-NA layer", plotly_spatialStackPlot(chat, all.na))
infinite <- layers
infinite[1L, 1L] <- Inf
expect_error("rejects infinite values", plotly_spatialStackPlot(chat, infinite))
nan <- layers
nan[1L, 1L] <- NaN
expect_error("rejects NaN values", plotly_spatialStackPlot(chat, nan))
bad.names <- layers
rownames(bad.names)[1L] <- "unknown-cell"
expect_error("rejects mismatched cell-name sets", plotly_spatialStackPlot(chat, bad.names))
duplicate.cells <- layers
rownames(duplicate.cells)[1L] <- rownames(duplicate.cells)[2L]
expect_error("rejects duplicate cell names", plotly_spatialStackPlot(chat, duplicate.cells))
duplicate.layers <- layers
colnames(duplicate.layers) <- c("same", "same", "third")
expect_error("rejects duplicate layer names", plotly_spatialStackPlot(chat, duplicate.layers))
expect_error("rejects invalid z-axis spacing",
             plotly_spatialStackPlot(chat, layers, z.axis.space = 0))
expect_error("rejects invalid marker size",
             plotly_spatialStackPlot(chat, layers, point.size = Inf))

cat("All spatial stack checks passed.\n")
