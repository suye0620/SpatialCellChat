# Visium producer regression: final 11-slot input/output and mixture semantics
setwd(local({
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg)) dirname(dirname(normalizePath(file_arg))) else getwd()
}))
source(".Rprofile")

suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(dplyr)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
close_enough <- function(a, b, tol = 1e-12) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol, check.attributes = FALSE))
}

cells <- paste0("cell", 1:4)
group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
names(group) <- cells
expr <- matrix(
  1,
  nrow = 2,
  ncol = length(cells),
  dimnames = list(c("ligand", "receptor"), cells)
)
coords <- matrix(
  c(0, 0, 3, 1, 1, 5, 8, 2),
  nrow = length(cells),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(cells, c("x", "y"))
)
meta <- data.frame(label = group, row.names = cells)
LRsig <- data.frame(
  ligand = c("ligand", "ligand"),
  receptor = c("receptor", "receptor"),
  pathway_name = c("P1", "P1"),
  row.names = c("LR.edge", "LR.zero")
)

make_chat <- function() {
  chat <- createSpatialCellChat(
    expr,
    meta = meta,
    group.by = "label",
    input.assay = "norm",
    normalize = FALSE,
    datatype = "spatial",
    coordinates = coords,
    coordinate.system = list(unit = "arbitrary", calibrated = FALSE),
    spatial.factors = list(ratio = NULL, tol = NULL)
  )
  chat@assay$signaling <- as(expr, "dgCMatrix")
  chat@LR$LRsig <- LRsig
  chat
}

# P is intentionally asymmetric and contains a stored explicit zero. D is
# non-one-hot, and the second layer is a named all-zero sparse layer.
P <- sparseMatrix(
  i = c(1L, 2L, 1L, 3L, 4L, 2L),
  j = c(2L, 1L, 4L, 1L, 3L, 4L),
  x = c(2, 1, 3, 4, 5, 6),
  dims = c(4L, 4L),
  giveCsparse = TRUE
)
P@x[1L] <- 0
P0 <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(4L, 4L), giveCsparse = TRUE)
prob <- SparseChatArray(
  list(LR.edge = P, LR.zero = P0),
  dimnames = list(cells, cells, c("LR.edge", "LR.zero"))
)
D <- matrix(
  c(0.8, 0.2,
    0.4, 0.6,
    0.3, 0.7,
    0.1, 0.9),
  nrow = 4,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(cells, c("A", "B"))
)
B <- 1 * (D > 0)
E <- P
E@x <- rep.int(1, length(E@x))

oracle_layer <- function(prob, decomposition, avg.type) {
  numerator <- as.matrix(crossprod(decomposition, prob %*% decomposition))
  if (identical(avg.type, "sum")) return(numerator)
  denominator <- as.matrix(crossprod(1 * (decomposition > 0), E %*% (1 * (decomposition > 0))))
  result <- numerator / denominator
  result[is.nan(result)] <- 0
  result
}

permutation_oracle <- function(prob, decomposition, avg.type, seed, nboot) {
  observed <- oracle_layer(prob, decomposition, avg.type)
  set.seed(seed)
  permutation <- matrix(
    replicate(nboot, sample.int(ncol(decomposition), size = ncol(decomposition))),
    nrow = ncol(decomposition),
    ncol = nboot
  )
  boot <- vapply(seq_len(nboot), function(b) {
    decomposition.boot <- decomposition[, permutation[, b], drop = FALSE]
    as.vector(oracle_layer(prob, decomposition.boot, avg.type))
  }, numeric(nrow(observed) * ncol(observed)))
  pval <- matrix(
    rowSums(boot - as.vector(observed) > 0) / nboot,
    nrow = nrow(observed),
    ncol = ncol(observed),
    byrow = FALSE
  )
  pval[observed == 0] <- 1
  pval
}

spot_row_permutation_oracle <- function(prob, decomposition, avg.type, seed, nboot) {
  observed <- oracle_layer(prob, decomposition, avg.type)
  set.seed(seed)
  permutation <- matrix(
    replicate(nboot, sample.int(nrow(decomposition), size = nrow(decomposition))),
    nrow = nrow(decomposition),
    ncol = nboot
  )
  boot <- vapply(seq_len(nboot), function(b) {
    decomposition.boot <- decomposition[permutation[, b], , drop = FALSE]
    as.vector(oracle_layer(prob, decomposition.boot, avg.type))
  }, numeric(nrow(observed) * ncol(observed)))
  pval <- matrix(
    rowSums(boot - as.vector(observed) > 0) / nboot,
    nrow = nrow(observed),
    ncol = ncol(observed),
    byrow = FALSE
  )
  pval[observed == 0] <- 1
  pval
}

base <- make_chat()
base@net$cell$prob <- prob
seed.use <- 37L
nboot <- 19L

for (avg.type in c("sum", "avg")) {
  out <- computeAvgCommunProb_Visium(
    base,
    cell.type.decomposition = D,
    avg.type = avg.type,
    do.permutation = TRUE,
    nboot = nboot,
    seed.use = seed.use
  )
  check(paste(avg.type, "group probability uses SparseChatArray"),
        inherits(out@net$group$prob, "SparseChatArray"))
  check(paste(avg.type, "group and LR names preserve order"),
        identical(dimnames(out@net$group$prob), list(c("A", "B"), c("A", "B"), c("LR.edge", "LR.zero"))))
  check(paste(avg.type, "active LR matches independent D'PD oracle"),
        close_enough(out@net$group$prob[["LR.edge"]], oracle_layer(P, D, avg.type)))
  check(paste(avg.type, "zero LR layer is retained"),
        all(as.matrix(out@net$group$prob[["LR.zero"]]) == 0))
  check(paste(avg.type, "p-value matches group-column permutation oracle"),
        close_enough(out@net$group$pval[["LR.edge"]],
                     permutation_oracle(P, D, avg.type, seed.use, nboot)))
  check(paste(avg.type, "zero LR p-values are one"),
        all(as.matrix(out@net$group$pval[["LR.zero"]]) == 1))
  check(paste(avg.type, "validObject passes"), isTRUE(validObject(out)))
  check(paste(avg.type, "validateSpatialCellChat passes"), isTRUE(validateSpatialCellChat(out)))
  check(paste(avg.type, "group-column and spot-row permutation differ"), {
    wrong <- spot_row_permutation_oracle(P, D, avg.type, seed.use, nboot)
    max(abs(as.matrix(out@net$group$pval[["LR.edge"]]) - wrong)) > 0
  })
  check(paste(avg.type, "runtime cache keys are retained"), {
    required.tmp <- c("LRsig.CCC.counts", "LRsig.GGC.counts", "LRsig.use.idx", "LRsig.use.idx.perm")
    all(required.tmp %in% names(out@net$tmp)) &&
      is.null(out@net$tmp$prob.cell) && is.null(out@net$tmp$cell.type.decomposition)
  })
  check(paste(avg.type, "decomposition is persisted"),
        identical(out@misc$.param$averaging$cell.type.decomposition, D))
  check(paste(avg.type, "cell counts are persisted"),
        close_enough(out@misc$.param$averaging$LRsig.CCC.counts, c(LR.edge = 6, LR.zero = 0)))
  check(paste(avg.type, "active indexes are persisted"),
        identical(as.integer(out@misc$.param$averaging$LRsig.use.idx), 1L)
)
  check(paste(avg.type, "key parameters and operation log are persisted"), {
    identical(out@misc$.param$averaging$do.permutation, TRUE) &&
      identical(out@misc$.param$averaging$nLR, 2L) &&
      length(out@misc$.log) >= 2L &&
      identical(out@misc$.log[[length(out@misc$.log)]][["function"]], "computeAvgCommunProb_Visium")
  })
}

no_perm <- computeAvgCommunProb_Visium(
  base,
  cell.type.decomposition = D,
  avg.type = "avg",
  do.permutation = FALSE,
  nboot = 0,
  seed.use = seed.use
)
check("do.permutation=FALSE writes no group p-value layer", is.null(no_perm@net$group$pval))
check("do.permutation=FALSE is persisted with counts and log", {
  identical(no_perm@misc$.param$averaging$do.permutation, FALSE) &&
    close_enough(no_perm@net$tmp$LRsig.CCC.counts, c(LR.edge = 6, LR.zero = 0)) &&
    identical(no_perm@net$tmp$LRsig.GGC.counts, c(LR.edge = 4, LR.zero = 0)) &&
    length(no_perm@misc$.log) >= 2L
})

# Named metadata labels are reordered to the final cell axis, not silently used
# in their metadata vector order.
named <- base
named@meta$named.group <- setNames(as.character(group)[c(3L, 1L, 4L, 2L)],
                                   c(cells[3L], cells[1L], cells[4L], cells[2L]))
named.out <- computeAvgCommunProb_Visium(
  named, cell.type.decomposition = D, group.by = "named.group",
  avg.type = "sum", do.permutation = FALSE
)
check("named group labels align to the probability cell axis", {
  close_enough(named.out@net$group$prob[["LR.edge"]], oracle_layer(P, D, "sum"))
})

# All validation failures below occur before an output group result or log entry
# can be created.
must_fail_before_compute <- function(label, object, decomposition = D,
                                     do.permutation = FALSE, ...) {
  err <- tryCatch(
    computeAvgCommunProb_Visium(object, cell.type.decomposition = decomposition,
                                do.permutation = do.permutation, ...),
    error = function(e) e
  )
  check(label, inherits(err, "error") && is.null(object@net$group$prob) && length(object@misc$.log) == 1L)
}
missing <- base
missing@net$cell$prob <- NULL
missing@net$prob.cell <- prob
must_fail_before_compute("missing net$cell$prob is rejected without old-slot fallback", missing)
old_only <- base
old_only@net$cell$prob <- NULL
old_only@net$prob.cell <- prob
must_fail_before_compute("old net$prob.cell alone is rejected", old_only)
bad_rows <- D
rownames(bad_rows) <- rev(rownames(bad_rows))
must_fail_before_compute("decomposition row order mismatch is rejected", base, bad_rows)
bad_cols <- D[, c(2L, 1L), drop = FALSE]
colnames(bad_cols) <- c("B", "A")
must_fail_before_compute("decomposition group column order mismatch is rejected", base, bad_cols)
bad_negative <- D
bad_negative[1L, 1L] <- -1
must_fail_before_compute("negative decomposition values are rejected", base, bad_negative)
bad_nonfinite <- D
bad_nonfinite[1L, 1L] <- Inf
must_fail_before_compute("non-finite decomposition values are rejected", base, bad_nonfinite)
no_cell_names <- base
no_cell_names@net$cell$prob <- SparseChatArray(
  list(P, P0), dimnames = list(NULL, NULL, c("LR.edge", "LR.zero"))
)
must_fail_before_compute("unnamed cell axes are rejected", no_cell_names)
non_square <- base
non_square@net$cell$prob <- SparseChatArray(
  list(as(sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(4L, 3L)), "dgCMatrix")),
  dimnames = list(cells, paste0("other", 1:3), "LR.edge")
)
must_fail_before_compute("non-square cell input is rejected", non_square,
                         matrix(1, nrow = 4, ncol = 2,
                                dimnames = list(cells, c("A", "B"))))
must_fail_before_compute("nboot zero is rejected with permutation", base, D,
                         do.permutation = TRUE, nboot = 0)

# relabelSpatialCellChat dispatches to the Visium producer when the cached
# decomposition exists. The chain uses min.cells=1 so this four-cell fixture
# remains eligible for filtering, pathway, and aggregate summaries.
relabel.visium <- make_chat()
relabel.visium@net$cell$prob <- prob
relabel.visium@misc$.param$averaging <- list(cell.type.decomposition = D)
relabel.visium <- relabelSpatialCellChat(
  relabel.visium,
  labelSet = c(B = "AB"),
  nboot = 5L,
  do.permutation = FALSE,
  min.cells = 1L,
  min.cells.sr = 1L,
  min.percent = 0
)
D.new <- cbind(D[, "A", drop = FALSE], AB = D[, "B"])
colnames(D.new) <- c("A", "AB")
check("relabel Visium dispatch keeps merged decomposition oracle", {
  close_enough(relabel.visium@net$group$prob[["LR.edge"]], oracle_layer(P, D.new, "avg")) &&
    identical(relabel.visium@misc$.param$averaging$cell.type.decomposition, D.new) &&
    identical(dimnames(relabel.visium@net$group$prob)[[1L]], c("A", "AB")) &&
    isTRUE(validateSpatialCellChat(relabel.visium)) &&
    !is.null(relabel.visium@netP$group$prob)
})

# Without the cache, relabel uses ordinary hard-label averaging. With positive
# expression and min.cells.sr=1, its avg oracle is the one-hot D'PD / B'EB form.
relabel.ordinary <- make_chat()
relabel.ordinary@net$cell$prob <- prob
relabel.ordinary <- relabelSpatialCellChat(
  relabel.ordinary,
  labelSet = c(B = "AB"),
  nboot = 5L,
  do.permutation = FALSE,
  min.cells = 1L,
  min.cells.sr = 1L,
  min.percent = 0
)
H <- matrix(0, nrow = 4L, ncol = 2L,
            dimnames = list(cells, c("A", "AB")))
H[cbind(seq_len(4L), match(c("A", "A", "AB", "AB"), colnames(H)))] <- 1
E.hard <- P
E.hard@x <- rep.int(1, length(E.hard@x))
hard.oracle <- as.matrix(crossprod(H, P %*% H)) /
  as.matrix(crossprod(H, E.hard %*% H))
hard.oracle[is.nan(hard.oracle)] <- 0
check("relabel without decomposition uses ordinary hard-label averaging", {
  is.null(relabel.ordinary@misc$.param$averaging$cell.type.decomposition) &&
    close_enough(relabel.ordinary@net$group$prob[["LR.edge"]], hard.oracle) &&
    isTRUE(validateSpatialCellChat(relabel.ordinary))
})

cat("All Visium producer migration checks passed.\n")
