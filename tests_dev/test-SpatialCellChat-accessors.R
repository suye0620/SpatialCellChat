# Phase 2 accessor tests for the final SpatialCellChat schema

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
})
source("R/SpatialCellChat_class.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

expect_error <- function(label, expr) {
  result <- try(force(expr), silent = TRUE)
  check(label, inherits(result, "try-error"))
}

expr <- Matrix::rsparsematrix(4, 3, density = 0.5)
rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
meta_data <- data.frame(
  label = factor(c("A", "B", "A")),
  sample = c("s1", "s1", "s1"),
  row.names = colnames(expr)
)
chat <- createSpatialCellChat(
  expr,
  meta = meta_data,
  group.by = "label",
  datatype = "RNA"
)

check("assay() returns the complete assay list", identical(assay(chat), chat@assay))
check("assay() selects a named layer", identical(assay(chat, "norm"), chat@assay$norm))
expect_error("assay() rejects unknown layers", assay(chat, "project"))

new_scale <- Matrix::Matrix(matrix(seq_len(length(expr)), nrow = nrow(expr),
                                    ncol = ncol(expr)), sparse = TRUE)
dimnames(new_scale) <- dimnames(expr)
chat2 <- chat
assay(chat2, "scale") <- new_scale
check("assay<- writes a named layer", identical(assay(chat2, "scale"), new_scale))

check("meta() returns metadata", identical(meta(chat), chat@meta))
check("meta() selects columns", identical(meta(chat, "label"), chat@meta[, "label", drop = FALSE]))
expect_error("meta() rejects unknown columns", meta(chat, "missing"))

meta2 <- meta(chat)
meta2$batch <- factor(c("b1", "b1", "b2"))
chat2 <- chat
meta(chat2) <- meta2
check("meta<- writes metadata", identical(meta(chat2), meta2))
expect_error("meta<- rejects cell-name mismatch", {
  bad_meta <- meta2
  rownames(bad_meta)[1] <- "other"
  meta(chat2) <- bad_meta
})

check("idents() returns a factor", is.factor(idents(chat)) && identical(idents(chat), chat@idents))
check("images() returns image state", identical(images(chat), chat@images))
expect_error("images() rejects unknown keys", images(chat, "missing"))

coords <- matrix(seq_len(6), ncol = 2L, byrow = TRUE)
rownames(coords) <- colnames(expr)
chat2 <- chat
images(chat2, "coordinates") <- coords
check("images<- writes coordinates", identical(images(chat2, "coordinates"), coords))

lr <- list(LRsig = data.frame(ligand = "L", receptor = "R"))
chat2 <- chat
LR(chat2) <- lr
check("LR<- replaces LR metadata", identical(LR(chat2), lr))
check("DB() returns the database", identical(DB(chat), chat@DB))

db <- list(interaction = data.frame(ligand = "L", receptor = "R"))
chat2 <- chat
DB(chat2) <- db
check("DB<- replaces the database", identical(DB(chat2), db))

chat2 <- chat
params(chat2, "interaction.range") <- 250
check("params<- writes a parameter", identical(params(chat2, "interaction.range"), 250))
check("params() returns the complete parameter list", is.list(params(chat2)))

chat2 <- chat
count_mat <- Matrix::Diagonal(ncol(expr))
dimnames(count_mat) <- list(colnames(expr), colnames(expr))
communication(chat2, "net", "cell", "count") <- count_mat
check("communication<- writes cell count", inherits(
  communication(chat2, "net", "cell", "count"), "diagonalMatrix"))
expect_error("communication() rejects unknown measures", communication(chat2, "net", "cell", "missing"))
expect_error("params() rejects unknown parameters", params(chat, "missing"))

check("updated object remains valid", isTRUE(validateSpatialCellChat(chat2)))
cat("All SpatialCellChat accessor checks passed.\n")
