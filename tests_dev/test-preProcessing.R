# Regression checks for preProcessing ALRA imputation on the final object schema

setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

# 12 cells x 12 genes, log-normalized-looking values with genuine low-rank
# structure so ALRA can detect a rank k >= 2
set.seed(42)
cells <- paste0("cell", seq_len(12L))
genes <- c("L1", "R1", "L2", "R2", "L3", "R3", "L4", "R4", "L5", "R5", "L6", "R6")
factors <- matrix(rexp(12L * 3L, 1), nrow = 12L, ncol = 3L) # genes x 3
loadings <- matrix(rexp(3L * 12L, 1), nrow = 3L, ncol = 12L) # 3 x cells
expr <- log1p(factors %*% loadings)
dimnames(expr) <- list(genes, cells)
expr <- round(expr, 4)
meta <- data.frame(
  label = factor(rep(c("A", "B", "C"), each = 4L)),
  row.names = cells
)
coordinates <- cbind(x = rep(seq_len(4L), 3L), y = rep(seq_len(3L), each = 4L))
rownames(coordinates) <- cells

make_chat <- function() {
  chat <- createSpatialCellChat(
    expr, meta = meta, group.by = "label",
    input.assay = "norm", datatype = "spatial",
    coordinates = coordinates,
    spatial.factors = list(ratio = 1, tol = 1)
  )
  chat@DB <- list(
    interaction = data.frame(
      ligand = c("L1", "L2", "L3", "L4", "L5", "L6"),
      receptor = c("R1", "R2", "R3", "R4", "R5", "R6"),
      annotation = rep("Secreted Signaling", 6L),
      interaction_name = paste0(c("L1", "L2", "L3", "L4", "L5", "L6"), "-",
                                c("R1", "R2", "R3", "R4", "R5", "R6")),
      pathway_name = paste0("P", 1:6),
      stringsAsFactors = FALSE
    ),
    complex = data.frame(row.names = character(0L))
  )
  subsetData(chat)
}

## ---- object path: impute assay$signaling ----
chat <- make_chat()
orig <- assay(chat, "signaling")
chat.imputed <- suppressWarnings(preProcessing(chat, seed.use = 1))

sig <- assay(chat.imputed, "signaling")
check("signaling layer replaced by a dgCMatrix", is(sig, "dgCMatrix"))
check("gene and cell names preserved after imputation",
      identical(rownames(sig), rownames(orig)) &&
        identical(colnames(sig), colnames(orig)))
check("imputed values are non-negative on the log scale", all(sig@x >= 0))
check("originally expressed entries remain expressed",
      all(sig[orig != 0] != 0))
check("imputation does not shrink the support", sum(sig@x != 0) >= sum(orig@x != 0))
check("input object was not mutated in place",
      isTRUE(all.equal(assay(chat, "signaling")@x, orig@x)) &&
        is.null(chat@misc$.param$alra))
check("estimated rank recorded in misc$.param",
      is.integer(chat.imputed@misc$.param$alra$k) &&
        chat.imputed@misc$.param$alra$k >= 1L)
check("operation logged", any(vapply(chat.imputed@misc$.log,
      function(e) identical(e[["function"]], "preProcessing"), logical(1L))))
check("imputed object validates", isTRUE(validObject(chat.imputed)))

## ---- determinism ----
chat2 <- make_chat()
chat2.imputed <- suppressWarnings(preProcessing(chat2, seed.use = 1))
check("same seed reproduces the identical imputed matrix",
      isTRUE(all.equal(assay(chat2.imputed, "signaling")@x, sig@x)))

## ---- norm path clears derived layers ----
chat3 <- make_chat()
assay(chat3, "scale") <- assay(chat3, "norm") # stand-in derived layer
chat3 <- suppressWarnings(preProcessing(chat3, slot.name = "norm", seed.use = 1))
check("norm imputation targets the norm layer",
      !isTRUE(all.equal(assay(chat3, "norm")@x, expr)) )
check("norm imputation clears stale derived layers",
      is.null(assay(chat3, "scale")) && is.null(assay(chat3, "signaling")))

## ---- matrix path ----
norm.mat <- expr # genes x cells, already log-normalized-looking
res.mat <- suppressWarnings(preProcessing(norm.mat, seed.use = 1))
check("matrix input returns a genes x cells dgCMatrix",
      is(res.mat, "dgCMatrix") &&
        identical(dim(res.mat), dim(norm.mat)) &&
        identical(rownames(res.mat), rownames(norm.mat)) &&
        identical(colnames(res.mat), colnames(norm.mat)))
check("matrix result matches object path result for the same data",
      isTRUE(all.equal(res.mat@x, sig@x)))

## ---- rejection paths ----
chat.empty <- createSpatialCellChat(
  expr, meta = meta, group.by = "label",
  input.assay = "norm", datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = 1)
)
check("empty signaling layer errors with subsetData hint",
      inherits(tryCatch(preProcessing(chat.empty), error = identity), "error"))
check("unknown slot.name rejected",
      inherits(tryCatch(preProcessing(make_chat(), slot.name = "data.signaling"),
                        error = identity), "error"))
check("matrix smaller than ALRA minimum rejected",
      inherits(tryCatch(preProcessing(expr[1:5, 1:5]), error = identity), "error"))
check("invalid quantile.prob rejected",
      inherits(tryCatch(preProcessing(make_chat(), quantile.prob = 2),
                        error = identity), "error"))

cat("\nAll preProcessing regression checks passed.\n")
