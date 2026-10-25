# analyze_group_diff6.R -- final doc numbers: counts by direction of pval
# flips, zero-crossings, and which LR layers carry complete 0<->1 flips.
# also: exact relative magnitudes of the common-support prob value diffs and
# the refactored-vs-upstream significance table per direction.

prefix <- "tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/"
d1 <- 8L; d2 <- 8L; d3 <- 1563L
flat_of <- function(i, j, k) (i - 1) * (d2 * d3) + (j - 1) * d3 + k

raw_po <- readRDS(paste0(prefix, "old_group_pval.rds"))
raw_pn_ <- readRDS(paste0(prefix, "new_group_pval.rds"))
stopifnot(is.array(raw_po))
allflat <- seq_len(d1 * d2 * d3)
i_v <- floor((allflat - 1) / (d2 * d3)) + 1
rem <- (allflat - 1) - (i_v - 1) * (d2 * d3)
j_v <- floor(rem / d3) + 1
k_v <- rem - (j_v - 1) * d3 + 1
p_old <- raw_po[cbind(i_v, j_v, k_v)]

layers <- if (inherits(raw_pn_, "SparseChatArray")) unclass(raw_pn_) else raw_pn_
p_new <- numeric(d1 * d2 * d3)
for (idx in seq_along(layers)) {
  dg <- methods::as(layers[[idx]], "dgCMatrix")
  cc <- diff(dg@p)
  ii <- as.numeric(dg@i) + 1
  jj <- as.numeric(rep.int(seq_along(cc), cc))
  p_new[flat_of(ii, jj, rep.int(idx, length(ii)))] <- as.numeric(dg@x)
}
dpv <- p_new - p_old

cat("total coords:", d1 * d2 * d3, "\n")
cat("pval diffs total:", sum(dpv != 0), "\n")
cat("  old==0 -> new==1 (complete flip):", sum(p_old == 0 & p_new == 1), "\n")
cat("  old==1 -> new==0 (complete flip):", sum(p_old == 1 & p_new == 0), "\n")
cat("  old==0 -> new in (0,1):", sum(p_old == 0 & p_new > 0 & p_new < 1), "\n")
cat("  old==1 -> new in (0,1):", sum(p_old == 1 & p_new > 0 & p_new < 1), "\n")
cat("  both-valid quantized shifts:", sum(p_old > 0 & p_old < 1 & p_new > 0 & p_new < 1 & dpv != 0), "\n")
cat("  of those, |diff|>0.5:", sum(p_old > 0 & p_old < 1 & p_new > 0 & p_new < 1 & abs(dpv) > 0.5), "\n")

cat("\nzero-crossing |pval|<=0.05 significance flips in FULL space:\n")
sig_o <- p_old >= 0 & p_old < 0.05
sig_n <- p_new >= 0 & p_new < 0.05
cat("  significant old:", sum(sig_o), " significant new:", sum(sig_n), "\n")
cat("  sig crossed old->insig:", sum(sig_o & !sig_n), " insig->sig:", sum(!sig_o & sig_n), "\n")

## complete-flip layer names
flip_k <- unique(k_v[p_old == 0 & p_new == 1])
cat("\ncomplete-flip (0->1) layer count:", length(flip_k), "\n")

og <- readRDS(paste0(prefix, "old_group_prob.rds"))
dn <- if (inherits(og, "sparse3Darray")) unclass(og)$dimnames else dimnames(og)
lr_names <- if (!is.null(dn) && !is.null(dn[[3]])) as.character(dn[[3]]) else as.character(seq_len(d3))
cat("first 20 complete-flip LR layers:\n")
print(paste(head(lr_names[flip_k], 20), collapse = ", "))

## relative magnitude check on common-support prob value diffs
raw_go <- readRDS(paste0(prefix, "old_group_prob.rds"))
v_go <- if (inherits(raw_go, "sparse3Darray")) unclass(raw_go) else raw_go
raw_gn <- readRDS(paste0(prefix, "new_group_prob.rds"))
v_gn <- NULL; x_gn <- NULL
lay_n <- if (inherits(raw_gn, "SparseChatArray")) unclass(raw_gn) else raw_gn
x_all <- numeric(0); i_all <- numeric(0); j_all <- numeric(0); k_all <- numeric(0)
for (idx in seq_along(lay_n)) {
  dg <- methods::as(lay_n[[idx]], "dgCMatrix")
  cc <- diff(dg@p)
  i_all <- c(i_all, as.numeric(dg@i) + 1)
  j_all <- c(j_all, as.numeric(rep.int(seq_along(cc), cc)))
  x_all <- c(x_all, as.numeric(dg@x))
  k_all <- c(k_all, rep.int(as.numeric(idx), length(dg@i)))
}
if (inherits(raw_go, "sparse3Darray")) {
  u <- unclass(raw_go)
  key_o <- paste(u$i, u$j, u$k, sep = "_")
  x_o <- u$x
} else {
  w <- which(raw_go != 0, arr.ind = TRUE)
  key_o <- paste(w[, 1], w[, 2], w[, 3], sep = "_")
  x_o <- raw_go[w]
}
key_n <- paste(i_all, j_all, k_all, sep = "_")
common_k <- intersect(key_o, key_n)
pos_o <- match(common_k, key_o); pos_n2 <- match(common_k, key_n)
x_n2 <- x_all[pos_n2]
xdiff <- x_o[pos_o] - x_n2
scale_ref <- pmax(abs(x_all[pos_n2]), .Machine$double.xmin)
cat("\ncommon prob support:", length(common_k),
    " differing:", sum(xdiff != 0),
    " max |dx|:", format(max(abs(xdiff)), digits = 6),
    " max |relative|:", format(max(abs(xdiff) / scale_ref), digits = 6), "\n")

cat("\nmax |dx| relative on the top-10 largest common entries:\n")
big10 <- order(abs(x_o[pos_o]), decreasing = TRUE)[seq_len(10)]
print(data.frame(key = common_k[big10],
                 x_old = x_o[pos_o][big10], x_new = x_all[pos_n2][big10],
                 dx = xdiff[big10],
                 rel = xdiff[big10] / scale_ref[big10]), row.names = FALSE)
