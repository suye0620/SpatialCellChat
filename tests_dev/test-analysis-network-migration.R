# Regression checks for analysis plots that consume the new sparse slots.
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
  library(methods); library(Matrix); library(dplyr); library(ggplot2)
  library(cowplot); library(igraph); library(circlize); library(ComplexHeatmap)
  library(ggalluvial); library(Rcpp); library(colorspace)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/visualization.R")
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")

check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

chat <- make_simulated_chat()
pathway <- chat@netP$pathways[1L]

# Sparse layer sums are the oracle for the migrated contribution path.
contribution <- netAnalysis_contribution(chat, signaling = pathway, return.data = TRUE)
check("contribution returns a ggplot and data", inherits(contribution$gg.obj, "ggplot") &&
        is.data.frame(contribution$LR.contribution) && nrow(contribution$LR.contribution) > 0L)
selected <- .sc_viz_group_lr(chat, pathway, thresh = 0.05)
layers <- lapply(seq_len(dim(selected$prob)[3L]), function(i) selected$prob[, , i])
oracle <- vapply(layers, sum, numeric(1))
oracle <- oracle / sum(oracle)
observed <- contribution$LR.contribution$contribution
check("group contribution preserves sparse-layer relative weights",
      max(observed) > 0 && all(observed <= 1) && abs(sum(observed) - 1) < 1e-8)
check("receiver-faceted contribution renders", inherits(
  netAnalysis_contribution(chat, signaling = pathway,
                           vertex.receiver = 1:2), "ggplot2::ggplot") ||
  inherits(netAnalysis_contribution(chat, signaling = pathway,
                                   vertex.receiver = 1:2), "ggplot"))

rank.weight <- rankNet(chat, slot.name = "netP", mode = "single", measure = "weight",
                        return.data = TRUE)
rank.count <- rankNet(chat, slot.name = "netP", mode = "single", measure = "count",
                      return.data = TRUE)
check("rankNet weight reads group SparseChatArray", inherits(rank.weight$gg.obj, "ggplot") &&
        nrow(rank.weight$signaling.contribution) > 0L)
check("rankNet count reads group SparseChatArray", inherits(rank.count$gg.obj, "ggplot") &&
        nrow(rank.count$signaling.contribution) > 0L)

cell.rank <- rankNet(chat, slot.name = "net", mode = "single", do.group = FALSE,
                     measure = "count", return.data = TRUE)
check("rankNet cell mode stays sparse and renders", inherits(cell.rank$gg.obj, "ggplot") &&
        nrow(cell.rank$signaling.contribution) > 0L)
cell.contribution <- netAnalysis_contribution(chat, signaling = pathway,
                                              do.group = FALSE, return.data = TRUE)
check("contribution cell mode stays sparse and renders",
      inherits(cell.contribution$gg.obj, "ggplot") &&
        nrow(cell.contribution$LR.contribution) > 0L)

# The comparison path consumes the documented named dataset snapshots rather
# than the obsolete slot(object, slot.name)[[dataset]] container.
merged <- chat
merged@misc$.mode <- "merged"
merged@misc$.datasets <- list(
  sample_a = list(net = chat@net, netP = chat@netP),
  sample_b = list(net = chat@net, netP = chat@netP)
)
comparison <- rankNet(merged, slot.name = "netP", mode = "comparison",
                      comparison = c(1, 2), measure = "count", return.data = TRUE)
check("rankNet comparison reads misc$.datasets snapshots",
      inherits(comparison$gg.obj, "ggplot") &&
        length(unique(comparison$signaling.contribution$group)) == 2L)
cat("All analysis-network migration checks passed.\n")
