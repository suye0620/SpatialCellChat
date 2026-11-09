# Regression tests for plotGeneExpression with the current 11-slot schema.

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
  library(methods)
  library(Matrix)
  library(dplyr)
  library(future.apply)
  library(cli)
  library(ggplot2)
  library(Seurat)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/database.R")
source("R/analysis.R")
source("R/visualization.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

cells <- paste0("cell", 1:6)
genes <- c("L1", "R1", "L2", "R2", "L3", "R3")
expr <- Matrix::Matrix(
  matrix(seq_len(length(genes) * length(cells)), nrow = length(genes),
         dimnames = list(genes, cells)),
  sparse = TRUE
)
meta <- data.frame(label = factor(rep(c("A", "B"), each = 3L)), row.names = cells)
chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")
chat@assay$signaling <- expr

lr <- data.frame(
  interaction_name = c("LR1", "LR2"),
  ligand = c("L1", "L2"),
  receptor = c("R1", "R2"),
  pathway_name = c("P1", "P1"),
  stringsAsFactors = FALSE,
  row.names = c("LR1", "LR2")
)
prob <- array(1, dim = c(2L, 2L, 2L),
              dimnames = list(c("A", "B"), c("A", "B"), rownames(lr)))
chat@LR <- list(LRsig = lr)
chat@DB <- list(
  interaction = lr,
  complex = data.frame(subunit_1 = character(), subunit_2 = character(),
                       row.names = character()),
  geneInfo = data.frame(Symbol = genes, stringsAsFactors = FALSE)
)
chat@net$group <- list(prob = prob, LR.sig = "LR2")

plot <- plotGeneExpression(chat, signaling = "P1", type = "dot", scale = FALSE)
check("pathway expression plot reads assay$signaling and new LR cache structure",
      inherits(plot, "ggplot"))

reordered <- chat
reordered@assay$signaling <- expr[, rev(cells), drop = FALSE]
reordered.plot <- plotGeneExpression(reordered, signaling = "P1", type = "dot",
                                     scale = FALSE)
check("metadata and idents are aligned to signaling columns by cell name",
      identical(plot$data, reordered.plot$data))

missing.signaling <- chat
missing.signaling@assay$signaling <- NULL
missing.error <- tryCatch({
  plotGeneExpression(missing.signaling, features = "L1", type = "dot")
  ""
}, error = conditionMessage)
check("missing signaling assay reports the current layer prerequisite",
      grepl("assay\\$signaling is empty", missing.error))
