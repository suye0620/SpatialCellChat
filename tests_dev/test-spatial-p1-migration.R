# Focused regression checks for the P1 visualization migration.
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
  library(future.apply); library(ggplot2); library(colorspace); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")

check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
chat <- make_simulated_chat()
chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = FALSE)
signal <- dimnames(chat@netP$centr.cell)[[3]][1]

single <- spatialGiPlot(chat, signaling.name = signal, measure = "indeg",
                        do.binary = FALSE, n = 4, return.object = FALSE)
check("single Gi measure returns a plot", inherits(single, "ggplot"))
gi.chat <- chat
set.seed(123)
gi.chat@images$coordinates <- gi.chat@images$coordinates +
  matrix(runif(length(gi.chat@images$coordinates), -0.01, 0.01), ncol = 2L)
original.gi <- spatialGiPlot(gi.chat, signaling.name = signal, measure = "indeg",
                             do.binary = FALSE, n = 4, return.object = TRUE, show.plot = FALSE)
reordered.gi <- gi.chat
reordered.gi@images$coordinates <- reordered.gi@images$coordinates[rev(seq_len(nrow(gi.chat@meta))), , drop = FALSE]
reordered.gi <- spatialGiPlot(reordered.gi, signaling.name = signal, measure = "indeg",
                              do.binary = FALSE, n = 4, return.object = TRUE, show.plot = FALSE)
gi.columns <- setdiff(names(original.gi@meta), names(chat@meta))
check("Gi results follow cell names after coordinates are reordered",
      isTRUE(all.equal(original.gi@meta[, gi.columns], reordered.gi@meta[, gi.columns],
                       check.attributes = FALSE)))
g <- plotStatistics_Gi(data.frame(x = seq_len(20), y = seq_len(20)^2,
                                  score = rep(1, 20)),
                       group.by = chat@idents, group.highlight = rep(FALSE, 20),
                       n = 4, method = "computation")
check("constant Gi remains undefined", all(is.nan(g$statistic_G)))
constant <- spatialGiPlot(chat, signaling.name = signal, measure = "indeg",
                          do.binary = TRUE, n = 4, return.object = FALSE)
check("undefined Gi has a renderable plot", inherits(constant, "ggplot"))

lee.data <- suppressWarnings(spatialLeePlot(chat, group.by1 = chat@idents,
                                            group.by2 = chat@idents, type = "data"))
lee.heatmap <- suppressWarnings(spatialLeePlot(chat, group.by1 = chat@idents,
                                               group.by2 = chat@idents, type = "heatmap"))
check("Lee default cutoff retains original values", isTRUE(all.equal(unname(diag(as.matrix(lee.data))), rep(1, nrow(lee.data)))) &&
        all(is.na(as.matrix(lee.data)[row(lee.data) != col(lee.data)])))
pdf.file <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf.file)
suppressWarnings(ComplexHeatmap::draw(lee.heatmap))
grDevices::dev.off()
check("single-valued Lee heatmap is actually drawn", file.info(pdf.file)$size > 0)
unlink(pdf.file)

chat <- identifyCellTopics(chat, slot.name = "netP", pattern = "outgoing",
                           topic.k = 2, verbose = FALSE)
topic <- chat@netP$topic$outgoing$cell
check("topic row axis is the original cell axis", identical(rownames(topic), rownames(chat@meta)))
p <- spatialTopicPlot(chat, slot.name = "netP", pattern = "outgoing",
                      topics = "Topic_1", combine = FALSE)[[1]]
expected <- as.numeric(topic[rownames(chat@images$coordinates), 1])
check("topic plot uses original loading", isTRUE(all.equal(p$data$feature.data, expected)))
reordered <- chat
reordered@images$coordinates <- reordered@images$coordinates[rev(seq_len(nrow(topic))), , drop = FALSE]
p2 <- spatialTopicPlot(reordered, slot.name = "netP", pattern = "outgoing",
                       topics = "Topic_1", combine = FALSE)[[1]]
check("topic values follow cell names", isTRUE(all.equal(p2$data$feature.data, rev(expected))))
custom <- spatialTopicPlot(chat, slot.name = "netP", pattern = "outgoing",
                           topics.key = "Custom_", topics = "Custom_1", combine = FALSE)
check("topic prefix retains historical behavior", inherits(custom[[1]], "ggplot"))
cat("All P1 migration checks passed.\n")
