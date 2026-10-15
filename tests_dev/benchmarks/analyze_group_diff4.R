# analyze_group_diff4.R -- resolve pval comparison over the FULL 8x8x1563 space
# using each side's storage semantics:
#   old rep: whatever class the upstream dump carries (report it). If dense,
#     zeros are REAL pval values. If sparse3Darray, absent coords are unknown,
#     and we report both zero-assumed and zero-excluded variants.
#   new rep: SparseChatArray layers built via my_as_sparse3Darray(...), which
#     drops explicit zeros; absent coords are therefore TRUE pval values of 0.
#
# outputs the full-space pvalue diff stats + the exact top flips.

prefix <- "tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/"
rawp <- readRDS(paste0(prefix, "old_group_pval.rds"))
cat("old_group_pval storage class:", class(rawp)[1L], "\n")
d1 <- 8L; d2 <- 8L; d3 <- 1563L
total <- d1 * d2 * d3
cat("coords total:", total, "\n")

## full-space pval arrays, values at exact coords
val_at <- function(prob) {
  if (inherits(prob, "sparse3Darray")) {
    u <- unclass(prob)
    idx <- (u$i - 1) * (d2 * d3) + (u$j - 1) * d3 + u$k
    v <- numeric(total)
    v[idx] <- u$x
    return(v)
  }
  if (length(dim(prob)) == 3L) return(as.numeric(prob))
  stop("unexpected class")
}
val_layers <- function(prob) {
  layers <- if (inherits(prob, "SparseChatArray")) unclass(prob) else prob
  stopifnot(is.list(layers), length(layers) == d3)
  v <- numeric(total)
  for (idx in seq_along(layers)) {
    m <- methods::as(layers[[idx]], "dgCMatrix")
    cc <- diff(m@p)
    ii <- as.numeric(m@i) + 1
    jj <- as.numeric(rep.int(seq_along(cc), cc))
    flat <- (ii - 1) * (d2 * d3) + (jj - 1) * d3 + idx
    v[flat] <- as.numeric(m@x)
  }
  v
}

po <- val_at(rawp)
cat("old pval zeros in dense rep:", sum(po == 0), "\n")
np <- readRDS(paste0(prefix, "new_group_pval.rds"))
pn <- val_layers(np)
cat("new pval stored:", sum(pn != 0), " zeros(absent->0):", sum(pn == 0), "\n")

dpv <- pn - po
dz <- dpv != 0
cat("\ncommon pval values (nonzero both sides):",
    sum(po != 0 & pn != 0), " differing:", sum(dz), "\n")
cat("|dpv| == 0.01 on common nonzero:", sum(dz & abs(dpv) == 0.01 & po != 0 & pn != 0),
    " / differing with other magnitude:",
    sum(dz & abs(dpv) != 0.01 & po != 0 & pn != 0), "\n")
cat("max |pdiff| on common nonzero:", max(abs(dpv)), " at",
    which.max(abs(dpv)), "\n")

## zero-boundary stats: entries where a nonzero pval collapses to zero or rises
cat("\npo==0, pn!=0:", sum(po == 0 & pn != 0),
    " po!=0, pn==0:", sum(po != 0 & pn == 0), "\n")
cat("max pn among po==0:", max(pn[po == 0]), "\n")
cat("max po among pn==0:", max(po[pn == 0]), "\n")

## significance flips in the FULL space (pval<0.05 OR ==0 means significant)
sig_o <- po < 0.05
sig_n <- pn < 0.05
cat("\nsignificant(<0.05): old", sum(sig_o), " new", sum(sig_n),
    " both", sum(sig_o & sig_n),
    " old-only", sum(sig_o & !sig_n),
    " new-only", sum(!sig_o & sig_n), "\n")

## top flips by pdiff
top <- order(abs(dpv) * (dpv != 0), decreasing = TRUE)[seq_len(10)]
cat("\ntop 10 pval flips:\n")
k <- 1
tgt <- data.frame(flat = top)
tgt$i <- floor((tgt$flat - 1) / (d2 * d3)) + 1
rem <- (tgt$flat - 1) - (tgt$i - 1) * (d2 * d3)
tgt$j <- floor(rem / d3) + 1
tgt$k <- rem - (tgt$j - 1) * d3 + 1
tgt$pval_old <- po[top]; tgt$pval_new <- pn[top]; tgt$diff <- dpv[top]
print(tgt, row.names = FALSE)

## also verify prob-side zero semantics with the same full-space rule
rawgO <- readRDS(paste0(prefix, "old_group_prob.rds"))
vo <- val_at(rawgO)
rawgN <- readRDS(paste0(prefix, "new_group_prob.rds"))
vn <- val_layers(rawgN)
cat("\nprob full-space: old nonzero", sum(vo != 0), " new nonzero", sum(vn != 0),
    " common", sum(vo != 0 & vn != 0),
    " only_old", sum(vo != 0 & vn == 0),
    " only_new", sum(vo == 0 & vn != 0), "\n")
