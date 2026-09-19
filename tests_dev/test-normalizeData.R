# Regression tests for explicit raw/norm expression-layer handling.

setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

raw <- matrix(c(2, 0, 1, 3, 4, 2), nrow = 3L, byrow = FALSE,
              dimnames = list(paste0("gene", 1:3), paste0("cell", 1:2)))
expected <- log1p(sweep(raw, 2L, colSums(raw), "/") * 10000)

norm.dense <- normalizeData(raw, verbose = FALSE)
check("dense normalization follows library-size log1p formula",
      is.matrix(norm.dense) && max(abs(norm.dense - expected)) < 1e-12)

raw.sparse <- methods::as(raw, "dgCMatrix")
norm.sparse <- normalizeData(raw.sparse, verbose = FALSE)
check("sparse normalization preserves dgCMatrix and dimnames",
      inherits(norm.sparse, "dgCMatrix") && identical(dimnames(norm.sparse), dimnames(raw)))
check("sparse and dense normalization agree",
      max(abs(as.matrix(norm.sparse) - expected)) < 1e-12)

meta <- data.frame(label = factor(c("A", "B")), row.names = colnames(raw))
chat <- createSpatialCellChat(raw.sparse, meta = meta, group.by = "label",
                              input.assay = "raw")
check("raw constructor stores raw and normalized layers",
      inherits(chat@assay[["raw"]], "dgCMatrix") &&
        inherits(chat@assay[["norm"]], "dgCMatrix"))
check("raw constructor records normalization parameters",
      identical(params(chat, "input.assay"), "raw") &&
        isTRUE(params(chat, "normalize")))

assay(chat, "scale") <- chat@assay[["norm"]]
chat <- normalizeData(chat, verbose = FALSE)
check("object normalization refreshes norm and clears derived layers",
      is.null(chat@assay[["scale"]]) && is.null(chat@assay[["smooth"]]) &&
        is.null(chat@assay[["signaling"]]) &&
        max(abs(as.matrix(chat@assay[["norm"]]) - expected)) < 1e-12)

zero.raw <- raw.sparse
zero.raw[, 1L] <- 0
check("zero-library cells are rejected",
      inherits(try(normalizeData(zero.raw, verbose = FALSE), silent = TRUE), "try-error"))
check("raw input cannot bypass normalization",
      inherits(try(createSpatialCellChat(raw.sparse, meta = meta, group.by = "label",
                                         input.assay = "raw", normalize = FALSE),
                   silent = TRUE), "try-error"))
check("normalized input remains an explicit no-raw path",
      is.null(createSpatialCellChat(norm.sparse, meta = meta, group.by = "label",
                                    input.assay = "norm")@assay[["raw"]]))

cat("All normalization checks passed.\n")
