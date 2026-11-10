# Regression tests for extractEnrichedLR on the current object schema.

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
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/database.R")
source("R/merging.R")
source("R/analysis.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

cells <- paste0("cell", 1:4)
genes <- c("L1", "R1", "L2", "R2", "L3", "R3")
expr <- Matrix::Matrix(1, nrow = length(genes), ncol = length(cells), sparse = TRUE,
                       dimnames = list(genes, cells))
meta <- data.frame(label = factor(c("A", "A", "B", "B")), row.names = cells)
chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")

lr <- data.frame(
  interaction_name = c("LR1", "LR2", "LR3"),
  ligand = c("L1", "L2", "L3"),
  receptor = c("R1", "R2", "R3"),
  pathway_name = c("P1", "P1", "P2"),
  stringsAsFactors = FALSE,
  row.names = c("LR1", "LR2", "LR3")
)
prob <- array(1, dim = c(2L, 2L, 3L),
              dimnames = list(c("A", "B"), c("A", "B"), rownames(lr)))
chat@LR <- list(LRsig = lr)
chat@DB <- list(
  interaction = lr,
  complex = data.frame(subunit_1 = character(), subunit_2 = character(),
                       row.names = character()),
  geneInfo = data.frame(Symbol = genes, stringsAsFactors = FALSE)
)
chat@net$group <- list(prob = prob, LR.sig = "LR2")
chat@net$cell <- list(prob = prob, LR.sig = "LR1")

group <- extractEnrichedLR(chat, "P1", geneLR.return = TRUE, thresh = 1e-10)
check("group cache selects LR names and matching genes; thresh does not re-test",
      identical(group$pairLR$interaction_name, "LR2") &&
        identical(group$geneLR, c("L2", "R2")))

cell <- extractEnrichedLR(chat, "P1", geneLR.return = TRUE, do.group = FALSE)
check("cell cache is read from net$cell$LR.sig",
      identical(cell$pairLR$interaction_name, "LR1") &&
        identical(cell$geneLR, c("L1", "R1")))

no.network <- chat
no.network@net <- list(group = list(), cell = list())
all.candidates <- extractEnrichedLR(no.network, "P1", geneLR.return = TRUE,
                                    enriched.only = FALSE)
check("enriched.only=FALSE returns every LR candidate in the pathway",
      identical(all.candidates$pairLR$interaction_name, c("LR1", "LR2")) &&
        identical(all.candidates$geneLR, c("L1", "L2", "R1", "R2")))

no.hits <- chat
no.hits@net$group$LR.sig <- character()
empty <- extractEnrichedLR(no.hits, "P1", geneLR.return = TRUE)
check("empty cache intersection returns no LR rows and no genes",
      nrow(empty$pairLR) == 0L && identical(empty$geneLR, character()))

merged <- chat
merged@misc$.mode <- "merged"
merged@misc$.datasets <- list(
  sample1 = list(LR = chat@LR,
                 net = list(group = list(prob = prob, LR.sig = "LR1"))),
  sample2 = list(LR = chat@LR,
                 net = list(group = list(prob = prob, LR.sig = "LR2")))
)
merged.result <- extractEnrichedLR(merged, "P1", geneLR.return = TRUE)
check("merged mode unions cached LR results from misc$.datasets",
      identical(merged.result$pairLR$interaction_name, c("LR1", "LR2")) &&
        setequal(merged.result$geneLR, c("L1", "R1", "L2", "R2")))

merged.local.lr <- merged
merged.local.lr@misc$.datasets$sample1$LR$LRsig <- lr["LR1", , drop = FALSE]
merged.local.lr@misc$.datasets$sample2$LR$LRsig <- lr["LR3", , drop = FALSE]
local.result <- extractEnrichedLR(merged.local.lr, "P1", geneLR.return = TRUE,
                                  enriched.only = FALSE)
check("merged candidate extraction uses each snapshot LR metadata",
      identical(local.result$pairLR$interaction_name, "LR1") &&
        identical(local.result$geneLR, c("L1", "R1")))
