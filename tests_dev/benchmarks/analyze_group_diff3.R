# analyze_group_diff3.R -- exact significance accounting + biggest pval flip,
# including pval==0 (most significant) which the nonzero-only lookup used by
# group_edge_diff.csv reports as NA there.

prefix <- "tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/"
load3 <- function(name) readRDS(paste0(prefix, name, ".rds"))

to_trips <- function(prob) {
  if (inherits(prob, "sparse3Darray")) {
    u <- unclass(prob)
    d <- as.integer(u$dim)
    list(i = as.numeric(u$i), j = as.numeric(u$j), k = as.numeric(u$k),
         x = as.numeric(u$x), dims = d)
  } else if (length(dim(prob)) == 3L) {
    w <- which(prob != 0, arr.ind = TRUE)
    list(i = as.numeric(w[, 1L]), j = as.numeric(w[, 2L]), k = as.numeric(w[, 3L]),
         x = as.numeric(prob[w]), dims = as.integer(dim(prob)))
  } else {
    stop("unsupported group pval storage: ", class(prob)[1L])
  }
}
to_layers <- function(prob) {
  layers <- if (inherits(prob, "SparseChatArray")) unclass(prob) else prob
  stopifnot(is.list(layers))
  k_all <- NULL
  x_all <- NULL
  i_all <- NULL
  j_all <- NULL
  for (idx in seq_along(layers)) {
    dg <- methods::as(layers[[idx]], "dgCMatrix")
    cc <- diff(dg@p)
    i_all <- c(i_all, as.numeric(dg@i) + 1)
    j_all <- c(j_all, as.numeric(rep.int(seq_along(cc), cc)))
    x_all <- c(x_all, as.numeric(dg@x))
    k_all <- c(k_all, rep.int(as.numeric(idx), length(dg@i)))
  }
  list(i = i_all, j = j_all, k = k_all, x = x_all,
       dims = as.integer(c(dim(layers[[1L]]), length(layers))))
}

old_g <- to_trips(load3("old_group_prob"))
new_g <- to_layers(load3("new_group_prob"))
old_p <- to_trips(load3("old_group_pval"))
new_p <- to_layers(load3("new_group_pval"))

flat <- function(t_) {
  d2 <- t_$dims[2L]; d3 <- t_$dims[3L]
  (t_$i - 1) * d2 * d3 + (t_$j - 1) * d3 + t_$k
}
## per-side full pval-value key (including pval==0): (flat, pval) from STORED pval
key_p_old <- setNames(old_p$x, as.character(flat(old_p)))
key_p_new <- setNames(new_p$x, as.character(flat(new_p)))
lookup0 <- function(keys, flat_vec) {
  out <- numeric(length(flat_vec))
  hit <- as.character(flat_vec) %in% names(keys)
  out[hit] <- keys[as.character(flat_vec)[hit]]
  out[!hit] <- NA
  out
}

e <- read.csv(paste0(prefix, "group_edge_diff.csv"))
oo <- e$side == "only_old"
n2 <- e$side == "only_new"
cm <- e$side == "common_diff"
flat_of_row <- function(rows) (rows$group_i - 1) * 8 * 1563 +
  (rows$group_j - 1) * 1563 + rows$layer

f_oo <- flat_of_row(e[oo, ])
f_n2 <- flat_of_row(e[n2, ])
f_cm <- flat_of_row(e[cm, ])

cat("== significance counting with pval in [0,1], 0 = most significant ==\n")
p_oo <- lookup0(key_p_old, f_oo)
p_n2 <- lookup0(key_p_new, f_oo)
cat("only_old edges (old-side pval): NA(=0? check)",
    sum(is.na(p_oo)), " p==0:", sum(!is.na(p_oo) & p_oo == 0),
    " p<0.05:", sum(!is.na(p_oo) & p_oo < 0.05),
    " p==1:", sum(!is.na(p_oo) & p_oo == 1), "\n")
cat("only_old edges (new-side pval exists):", sum(!is.na(lookup0(key_p_new, f_oo))), "\n")

p_n2 <- lookup0(key_p_new, f_n2)
cat("only_new edges (new-side pval): NA", sum(is.na(p_n2)),
    " p==0:", sum(!is.na(p_n2) & p_n2 == 0),
    " p<0.05:", sum(!is.na(p_n2) & p_n2 < 0.05),
    " p==1:", sum(!is.na(p_n2) & p_n2 == 1), "\n")
cat("only_new edges (old-side pval exists):", sum(!is.na(lookup0(key_p_old, f_n2))), "\n")

po <- lookup0(key_p_old, f_cm)
pn_ <- lookup0(key_p_new, f_cm)
cat("common_diff edges: pval_old NA count:", sum(is.na(po)),
    " pval_new NA:", sum(is.na(pn_)),
    " diff(old<0.05, new>=0.05):", sum(!is.na(po) & !is.na(pn_) & po < 0.05 & pn_ >= 0.05),
    " diff(new<0.05, old>=0.05):", sum(!is.na(po) & !is.na(pn_) & pn_ < 0.05 & po >= 0.05), "\n")
cat("common_diff pdiff==+-0.01 count:",
    sum(!is.na(po) & !is.na(pn_) & abs(pn_ - po) == 0.01), " of ", sum(cm), "\n")
cat("common_diff |pdiff| != 0.01 count:",
    sum(!is.na(po) & !is.na(pn_) & abs(pn_ - po) != 0.01 & pn_ != po), "\n")
wbig <- which(!is.na(po) & !is.na(pn_) & abs(pn_ - po) >= 0.5)
if (length(wbig)) {
  cat("big |pdiff| >= 0.5 rows on common support:\n")
  print(cbind(lr = e$lr[cm][wbig], i = e$group_i[cm][wbig], j = e$group_j[cm][wbig],
              x_old = e$x_old[cm][wbig], x_new = e$x_new[cm][wbig],
              pval_old = po[wbig], pval_new = pn_[wbig], pdiff = pn_[wbig] - po[wbig]), row.names = FALSE)
}
