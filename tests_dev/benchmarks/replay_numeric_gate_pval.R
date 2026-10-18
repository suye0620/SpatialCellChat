#!/usr/bin/env Rscript
# Replay numeric-gate group averaging on a post-filter checkpoint.
#
# Compares the R/Matrix numeric reference with the C++ CSC kernel using the
# same cell-level input, group labels, seed, permutation stream, and
# `signif(mean, 1) >= min.percent` gate. This isolates residual p-value drift
# caused by aggregation order after the pinned-upstream format gate is removed.
#
# Usage:
#   Rscript tests_dev/benchmarks/replay_numeric_gate_pval.R [dump_dir] [n_layers] [nboot]

args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
if (!nzchar(file_arg)) stop("Run this script with Rscript", call. = FALSE)
root <- dirname(dirname(dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))))
setwd(root)

dump_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], winslash = "/") else
  file.path(root, "tests_dev", "benchmarks", "out", "mbm05_serial_tensor_comparison")
n_layers <- if (length(args) >= 2L) as.integer(args[[2L]]) else 30L
nboot <- if (length(args) >= 3L) as.integer(args[[3L]]) else 100L
if (!dir.exists(dump_dir)) stop("dump directory missing: ", dump_dir)
if (!is.finite(n_layers) || n_layers < 1L || !is.finite(nboot) || nboot < 1L)
  stop("n_layers and nboot must be positive integers")

Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/reference-computeAvgCommunProb-numeric.R")

obj_path <- file.path(dump_dir, "new_postfilter_object.rds")
if (!file.exists(obj_path)) stop("missing checkpoint: ", obj_path)
cat("loading checkpoint:", obj_path, "\n")
obj <- readRDS(obj_path)
prob.cell_ <- unclass(obj@net$cell$prob)
group <- obj@idents
sig <- assay(obj, "signaling")
complex_subunits <- .sc_complex_subunits(obj@DB$complex)
cache <- new.env(parent = emptyenv())
LRsig <- dimnames(obj@net$cell$prob)[[3L]]
if (is.null(LRsig)) LRsig <- names(prob.cell_)
if (is.null(LRsig)) LRsig <- paste0("LR", seq_along(prob.cell_))
LRtab <- obj@LR$LRsig[LRsig, , drop = FALSE]
geneL <- as.character(LRtab$ligand)
geneR <- as.character(LRtab$receptor)
nC <- ncol(sig)
numCluster <- nlevels(group)

nnz <- vapply(prob.cell_, function(m) length(m@x), numeric(1))
cand <- which(nnz > 0)
if (!length(cand)) stop("checkpoint has no non-empty cell-level layers")
set.seed(42L)
sample.idx <- sort(sample(cand, min(n_layers, length(cand))))
cat("layers:", length(sample.idx), "of", length(prob.cell_), "nboot:", nboot,
    "groups:", numCluster, "cells:", nC, "\n")

thr <- vapply(tabulate(as.integer(group), nbins = numCluster), .sc_gate_threshold,
              integer(1), min.percent = 0.1)
gv <- as.integer(group) - 1L
layer_rows <- vector("list", length(sample.idx))
diff_rows <- vector("list", length(sample.idx))
for (z in seq_along(sample.idx)) {
  i <- sample.idx[[z]]
  rowsL <- .sc_expr_row(geneL[i], sig, complex_subunits, cache)
  rowsR <- .sc_expr_row(geneR[i], sig, complex_subunits, cache)
  ref <- computeAvgCommunProb_numeric_ref(
    prob.cell_ = prob.cell_[i], dataLavg = matrix(rowsL, nrow = 1L),
    dataRavg = matrix(rowsR, nrow = 1L), group = group,
    avg.type = "avg", nboot = nboot, seed.use = 1L)

  m <- prob.cell_[[i]]
  suppL <- as.numeric(rowsL > 0); suppR <- as.numeric(rowsR > 0)
  sr_out <- tabulate(m@i + 1L, nbins = nC)
  sr_in <- tabulate(rep.int(seq_len(nC), diff(m@p)), nbins = nC)
  Pnew <- cpp_group_avg_obs(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                            gv, numCluster, thr, 5, 0L, nC)
  Pnew <- matrix(Pnew, nrow = numCluster, ncol = numCluster)
  Pnull <- as.vector(Pnew)
  set.seed(1L)
  permutation <- matrix(unlist(lapply(seq_len(nboot), function(...) sample.int(nC, size = nC)),
                             use.names = FALSE), nrow = nC, ncol = nboot)
  Pboot <- cpp_group_avg_perm(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                              gv, numCluster, thr, 5, 0L, permutation, 1L)
  pnew <- rowSums(Pboot - Pnull > 0) / nboot
  Pnew.p <- matrix(pnew, nrow = numCluster, ncol = numCluster)
  Pnew.p[Pnew == 0] <- 1

  Pref <- ref$Prob[, , 1L]
  pref <- ref$Pval[, , 1L]
  prob_delta <- Pnew - Pref
  p_delta <- Pnew.p - pref
  support_mismatch <- sum((Pnew > 0) != (Pref > 0))
  p_mismatch <- p_delta != 0
  common <- (Pnew > 0) & (Pref > 0)
  sig_flip <- (Pnew.p < 0.05) != (pref < 0.05)
  layer_rows[[z]] <- data.frame(
    layer = LRsig[i], layer_index = i,
    support_mismatch = support_mismatch,
    prob_common = sum(common),
    prob_max_abs = max(abs(prob_delta)),
    prob_max_rel = max(abs(prob_delta) / pmax(abs(Pref), 1e-300)),
    pval_mismatch = sum(p_mismatch),
    pval_max_abs = max(abs(p_delta)),
    pval_max_boots = max(abs(p_delta)) * nboot,
    pval_significance_flips = sum(sig_flip),
    stringsAsFactors = FALSE
  )
  if (any(p_mismatch)) {
    ij <- which(p_mismatch, arr.ind = TRUE)
    diff_rows[[z]] <- data.frame(
      layer = LRsig[i], layer_index = i,
      group_i = ij[, 1L], group_j = ij[, 2L],
      prob_numeric_ref = Pref[p_mismatch], prob_kernel = Pnew[p_mismatch],
      pval_numeric_ref = pref[p_mismatch], pval_kernel = Pnew.p[p_mismatch],
      pval_delta = p_delta[p_mismatch],
      significance_flip = sig_flip[p_mismatch],
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("[%d/%d] %s support_mismatch=%d pval_mismatch=%d max_boots=%.0f sig_flips=%d\n",
              z, length(sample.idx), LRsig[i], support_mismatch, sum(p_mismatch),
              max(abs(p_delta)) * nboot, sum(sig_flip)))
  rm(ref, m, Pboot, permutation, Pref, pref, Pnew, Pnew.p, p_delta, prob_delta)
  invisible(gc())
}

layers_out <- do.call(rbind, layer_rows)
diffs_out <- diff_rows[vapply(diff_rows, Negate(is.null), logical(1L))]
diffs_out <- if (length(diffs_out)) do.call(rbind, diffs_out) else
  data.frame(layer = character(), layer_index = integer(), group_i = integer(), group_j = integer(),
             prob_numeric_ref = numeric(), prob_kernel = numeric(),
             pval_numeric_ref = numeric(), pval_kernel = numeric(),
             pval_delta = numeric(), significance_flip = logical(),
             stringsAsFactors = FALSE)
summary <- data.frame(
  metric = c("n_layers", "nboot", "support_mismatch_total", "max_prob_abs", "max_prob_rel",
             "pval_mismatch_total", "max_pval_abs", "max_pval_boots",
             "significance_flips", "layers_with_pval_mismatch"),
  value = c(nrow(layers_out), nboot, sum(layers_out$support_mismatch),
            max(layers_out$prob_max_abs), max(layers_out$prob_max_rel),
            sum(layers_out$pval_mismatch), max(layers_out$pval_max_abs),
            max(layers_out$pval_max_boots), sum(layers_out$pval_significance_flips),
            sum(layers_out$pval_mismatch > 0)),
  stringsAsFactors = FALSE
)
write.csv(layers_out, file.path(dump_dir, "numeric_gate_replay_layers.csv"), row.names = FALSE)
write.csv(diffs_out, file.path(dump_dir, "numeric_gate_replay_pval_diff.csv"), row.names = FALSE)
write.csv(summary, file.path(dump_dir, "numeric_gate_replay_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
cat("DONE numeric gate replay\n")
