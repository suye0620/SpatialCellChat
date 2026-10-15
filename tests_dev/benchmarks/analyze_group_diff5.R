# analyze_group_diff5.R -- FULL-SPACE pval comparison with a coordinate-array
# read that is consistent with the stored triple convention of compare_serial_tensors.
#   old pval: dense 3d array 8x8x1563 (class array confirmed)
#   new pval: SparseChatArray layers, absent entry == TRUE pval 0 (explicit zeros dropped)
# We read the dense array by cbind(i,j,k) indexing at the SAME i/j/k frames used
# for the stored triples, so both sides share one key space.

prefix <- "tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/"
d1 <- 8L; d2 <- 8L; d3 <- 1563L

flat_of <- function(i, j, k) (i - 1) * (d2 * d3) + (j - 1) * d3 + k
unflat <- function(flat) {
  i <- floor((flat - 1) / (d2 * d3)) + 1
  rem <- (flat - 1) - (i - 1) * (d2 * d3)
  j <- floor(rem / d3) + 1
  k <- rem - (j - 1) * d3 + 1
  cbind(i = i, j = j, k = k)
}

raw_pval_old <- readRDS(paste0(prefix, "old_group_pval.rds"))
cat("old pval class:", class(raw_pval_old)[1L],
    " dims:", paste(dim(raw_pval_old), collapse = "x"), "\n")
stopifnot(is.array(raw_pval_old))

raw_pval_new <- readRDS(paste0(prefix, "new_group_pval.rds"))
layers <- if (inherits(raw_pval_new, "SparseChatArray")) unclass(raw_pval_new) else raw_pval_new
stopifnot(is.list(layers), length(layers) == d3)

## full-space value vectors in the canonical key space
p_old <- numeric(d1 * d2 * d3)
p_new <- numeric(d1 * d2 * d3)

## dense old: for every canonical flat key, (i,j,k) = unflat(flat); value at
## raw_pval_old[i, j, k] (R column-major array indexing).
allflat <- seq_len(d1 * d2 * d3)
C <- unflat(allflat)
p_old <- raw_pval_old[cbind(C[, "i"], C[, "j"], C[, "k"])]
cat("old dense zeros:", sum(p_old == 0), "\n")

## sparse new: stored entries; absent -> 0
for (idx in seq_along(layers)) {
  dg <- methods::as(layers[[idx]], "dgCMatrix")
  cc <- diff(dg@p)
  ii <- as.numeric(dg@i) + 1
  jj <- as.numeric(rep.int(seq_along(cc), cc))
  flatv <- flat_of(ii, jj, rep.int(idx, length(ii)))
  p_new[flatv] <- as.numeric(dg@x)
}
cat("new stored nonzero:", sum(p_new != 0),
    " zeros(absent -> 0):", sum(p_new == 0), "\n")

dpv <- p_new - p_old
cat("\n|dpv| != 0 count:", sum(dpv != 0))
cat(" max:", max(abs(dpv)), "\n")

## zero-crossing deltas
cat("old>0 -> new==0:", sum(p_old != 0 & p_new == 0), "\n")
cat("old==0 -> new>0:", sum(p_old == 0 & p_new != 0), "\n")
cat("both nonzero, diff:", sum(p_old != 0 & p_new != 0 & dpv != 0), "\n")
cat("  both-nonzero |diff|==0.01:",
    sum(p_old != 0 & p_new != 0 & dpv != 0 & abs(dpv) == 0.01), "\n")
cat("  both-nonzero |diff|>0.01:",
    sum(p_old != 0 & p_new != 0 & dpv != 0 & abs(dpv) > 0.01), "\n")

## significance flips in full space (pval 0 or <0.05 -> significant)
sig_o <- p_old < 0.05 & p_old >= 0
sig_n <- p_new < 0.05 & p_new >= 0
cat("\nsignificant(<0.05): old", sum(sig_o), " new", sum(sig_n),
    " both", sum(sig_o & sig_n),
    " old-only", sum(sig_o & !sig_n),
    " new-only", sum(!sig_o & sig_n), "\n")

## top flips by |dpv|
ord <- order(abs(dpv) * (dpv != 0), decreasing = TRUE)[seq_len(10)]
fl <- unflat(ord)
cat("\ntop 10 pval flips:\n")
print(data.frame(lr_k = fl[, "k"], i = fl[, "i"], j = fl[, "j"],
                 pval_old = p_old[ord], pval_new = p_new[ord], diff = dpv[ord]),
      row.names = FALSE)

## name check: old side lr name lookup is dimnames [[3]] of old prob (not pval)
og <- readRDS(paste0(prefix, "old_group_prob.rds"))
dn <- if (inherits(og, "sparse3Darray")) unclass(og)$dimnames else dimnames(og)
if (!is.null(dn) && !is.null(dn[[3]])) {
  cat("\nLR names for layers named (first 12):\n")
  print(paste(as.character(dn[[3]])[fl[, "k"]], collapse = ", "))
}
