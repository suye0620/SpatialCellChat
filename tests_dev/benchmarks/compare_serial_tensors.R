#!/usr/bin/env Rscript
# compare_serial_tensors.R -- compare old (pinned upstream) vs new (refactored)
# SERIAL pipeline tensor dumps produced by benchmark_mbm05_upstream_vs_refactored.R
# with --dump-dir. Reading ~6 GB per side; peak RAM ~25 GB on this fixture.
#
# Usage: Rscript compare_serial_tensors.R <dump_dir>
# Writes inside <dump_dir>:
#   prefilter_cell_comparison.csv   aggregate stats for pre-filter tensors
#   postfilter_cell_comparison.csv  aggregate stats for post-filter tensors
#   group_summary.csv               group tensor aggregate stats
#   group_edge_diff.csv             every differing group support entry (detail)
# and prints the same to stdout.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1L])
if (!nzchar(file_arg)) stop("Run this file with Rscript", call. = FALSE)
root <- dirname(dirname(dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))))
setwd(root)
user_args <- commandArgs(trailingOnly = TRUE)
dump_dir <- if (length(user_args)) normalizePath(user_args[[1L]], winslash = "/") else
  file.path(root, "tests_dev", "benchmarks", "out", "mbm05_serial_tensor_comparison")
stopifnot(dir.exists(dump_dir))
cat("dump_dir:", dump_dir, "\n")

Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
})

## ---- canonicalization -------------------------------------------------

to_i_j_k_x <- function(prob, label) {
  if (inherits(prob, "sparse3Darray")) {
    u <- unclass(prob)
    dname <- u$dimnames
    dimnames_out <- lapply(seq_along(u$dim), function(d)
      if (!is.null(dname) && !is.null(dname[[d]])) as.character(dname[[d]])
      else as.character(seq_len(u$dim[d])))
    return(list(label = label, i = as.numeric(u$i), j = as.numeric(u$j),
                k = as.numeric(u$k), x = as.numeric(u$x),
                dims = as.integer(u$dim), dimnames = dimnames_out))
  }
  if (is.list(prob)) {
    ## SparseChatArray or a bare per-layer list: stored entries per layer.
    layers <- if (inherits(prob, "SparseChatArray")) unclass(prob) else prob
    stopifnot(is.list(layers), length(layers) >= 1L)
    dims <- NULL
    parts_i <- vector("list", length(layers))
    parts_j <- vector("list", length(layers))
    parts_k <- vector("list", length(layers))
    parts_x <- vector("list", length(layers))
    for (idx in seq_along(layers)) {
      m <- layers[[idx]]
      dm <- as.integer(dim(m))
      if (is.null(dims)) dims <- c(dm, as.integer(length(layers)))
      if (!identical(dm, dims[seq_along(dm)])) stop(label, ": layer dims differ")
      if (inherits(m, "sparseMatrix")) {
        dg <- methods::as(m, "dgCMatrix")
        col_counts <- diff(dg@p)
        n_idx <- length(dg@i)
        if (!identical(n_idx, as.integer(dg@p[length(dg@p)])))
          stop(label, ": dgCMatrix length mismatch")
        parts_i[[idx]] <- as.numeric(dg@i) + 1
        parts_j[[idx]] <- as.numeric(rep.int(seq_along(col_counts), col_counts))
        parts_k[[idx]] <- rep.int(as.numeric(idx), n_idx)
        parts_x[[idx]] <- as.numeric(dg@x)
        rm(dg)
      } else {
        stop(label, ": unsupported layer storage class ", class(m)[1L])
      }
    }
    i_all <- unlist(parts_i, use.names = FALSE)
    rm(parts_i); invisible(gc())
    j_all <- unlist(parts_j, use.names = FALSE)
    rm(parts_j); invisible(gc())
    k_all <- unlist(parts_k, use.names = FALSE)
    rm(parts_k); invisible(gc())
    x_all <- unlist(parts_x, use.names = FALSE)
    rm(parts_x); invisible(gc())
    dn <- dimnames(prob)
    dimnames_out <- lapply(seq_along(dims), function(d)
      if (!is.null(dn) && !is.null(dn[[d]])) as.character(dn[[d]])
      else as.character(seq_len(dims[d])))
    return(list(label = label, i = i_all, j = j_all, k = k_all, x = x_all,
                dims = dims, dimnames = dimnames_out))
  }
  if (length(dim(prob)) == 3L) {
    w <- which(prob != 0, arr.ind = TRUE)
    dname <- dimnames(prob)
    dimnames_out <- lapply(seq_along(dim(prob)), function(d)
      if (!is.null(dname) && !is.null(dname[[d]])) as.character(dname[[d]])
      else as.character(seq_len(dim(prob)[d])))
    return(list(label = label, i = as.numeric(w[, 1L]), j = as.numeric(w[, 2L]),
                k = as.numeric(w[, 3L]), x = as.numeric(prob[w]),
                dims = as.integer(dim(prob)), dimnames = dimnames_out))
  }
  stop(label, ": unsupported tensor class ", class(prob)[1L])
}

flat_index <- function(t_) {
  d2 <- t_$dims[2L]; d3 <- t_$dims[3L]
  (t_$i - 1) * d2 * d3 + (t_$j - 1) * d3 + t_$k
}

load_tensor <- function(name) {
  path <- file.path(dump_dir, paste0(name, ".rds"))
  if (!file.exists(path)) stop("missing dump: ", path)
  cat("loading", name, "(", round(file.info(path)$size / 1024^2, 0), "MB )\n")
  to_i_j_k_x(readRDS(path), name)
}

sum_fmt <- function(v) format(sum(v), digits = 15, scientific = TRUE)

num_diff_table <- function(idx_o, x_o, idx_n, x_n) {
  hit <- match(idx_n, idx_o)
  common_n <- sum(!is.na(hit))
  dx <- x_o[hit[!is.na(hit)]] - x_n[!is.na(hit)]
  list(common = common_n,
       only_old = length(idx_o) - common_n,
       only_new = sum(is.na(hit)),
       dx = dx,
       max_abs = if (length(dx)) max(abs(dx)) else NA_real_,
       sum_abs = if (length(dx)) sum(abs(dx)) else NA_real_,
       n_diff_gt_0 = if (length(dx)) sum(dx != 0) else 0L,
       n_diff_gt_1e12 = if (length(dx)) sum(abs(dx) > 1e-12) else 0L)
}

compare_cell_stage <- function(old_name, new_name, out_name) {
  cat("\n===", out_name, ":", old_name, "vs", new_name, "===\n")
  old <- load_tensor(old_name)
  new <- load_tensor(new_name)
  old_stored <- length(old$x); new_stored <- length(new$x)
  old_nz_n <- sum(old$x != 0); new_nz_n <- sum(new$x != 0)
  old_sum <- sum_fmt(old$x); new_sum <- sum_fmt(new$x)
  dims_equal <- isTRUE(identical(old$dims, new$dims))
  old_nz <- list(idx = flat_index(old), x = old$x); old_nz$idx <- old_nz$idx[old$x != 0]
  old_nz$x <- old$x[old$x != 0]
  rm(old); invisible(gc())
  new_nz <- list(idx = flat_index(new), x = new$x); new_nz$idx <- new_nz$idx[new$x != 0]
  new_nz$x <- new$x[new$x != 0]
  rm(new); invisible(gc())
  cat("sorting/aligning...\n")
  ord_o <- order(old_nz$idx)
  idx_o <- old_nz$idx[ord_o]; x_o <- old_nz$x[ord_o]
  rm(old_nz, ord_o); invisible(gc())
  ord_n <- order(new_nz$idx)
  idx_n <- new_nz$idx[ord_n]; x_n <- new_nz$x[ord_n]
  rm(new_nz, ord_n); invisible(gc())
  cmp <- num_diff_table(idx_o, x_o, idx_n, x_n)
  stats <- data.frame(
    metric = c("dims_equal", "old_stored_entries", "new_stored_entries",
               "old_nonzero_entries", "new_nonzero_entries",
               "old_stored_sum", "new_stored_sum",
               "common_nonzero_entries", "only_old_entries", "only_new_entries",
               "max_abs_value_diff_on_common", "sum_abs_value_diff",
               "n_common_diff_gt_0", "n_common_diff_gt_1e-12"),
    value = c(dims_equal, old_stored, new_stored, old_nz_n, new_nz_n,
              old_sum, new_sum, cmp$common, cmp$only_old, cmp$only_new,
              format(cmp$max_abs, digits = 15), format(cmp$sum_abs, digits = 15),
              cmp$n_diff_gt_0, cmp$n_diff_gt_1e12),
    stringsAsFactors = FALSE)
  write.csv(stats, file.path(dump_dir, paste0(out_name, "_comparison.csv")),
            row.names = FALSE)
  print(stats, row.names = FALSE)
  rm(idx_o, x_o, idx_n, x_n, cmp, stats); invisible(gc())
  invisible(NULL)
}

cat("=== cell level ===\n")
compare_cell_stage("old_prefilter_cell_prob", "new_prefilter_cell_prob", "prefilter_cell")
compare_cell_stage("old_postfilter_cell_prob", "new_postfilter_cell_prob", "postfilter_cell")

## ---- group level (full detail) ------------------------------------------

cat("\n=== group level ===\n")
old_g <- load_tensor("old_group_prob")
new_g <- load_tensor("new_group_prob")
old_p <- load_tensor("old_group_pval")
new_p <- load_tensor("new_group_pval")
stopifnot(identical(old_g$dims, new_g$dims),
          identical(old_p$dims, old_g$dims),
          identical(new_p$dims, old_g$dims))
g_names_old <- old_g$dimnames[[3L]]
g_names_new <- new_g$dimnames[[3L]]
if (!identical(g_names_old, g_names_new)) {
  cat("WARNING: group dimnames[[3]] differ; using old-side LR names\n")
}
g_names <- g_names_old

g_o <- nonzero_view_safe <- function(t_) {
  keep <- t_$x != 0
  list(idx = flat_index(t_)[keep], x = t_$x[keep])
}(old_g)
g_n <- g_o <- NULL; g_o <- NULL
g_o_nz <- function(t_) g_o  # placeholder removed below

d1 <- old_g$dims[1L]; d2 <- old_g$dims[2L]; d3 <- old_g$dims[3L]
view_nz <- function(t_) {
  keep <- t_$x != 0
  list(idx = flat_index(t_)[keep], x = t_$x[keep])
}
g_vo <- view_nz(old_g); g_vn <- view_nz(new_g)
ord_go <- order(g_vo$idx); idx_go <- g_vo$idx[ord_go]; x_go <- g_vo$x[ord_go]
ord_gn <- order(g_vn$idx); idx_gn <- g_vn$idx[ord_gn]; x_gn <- g_vn$x[ord_gn]
g_cmp <- num_diff_table(idx_go, x_go, idx_gn, x_gn)

decode_rows <- function(flat, side, x_old, x_new) {
  i <- floor((flat - 1) / (d2 * d3)) + 1
  rem <- (flat - 1) - (i - 1) * d2 * d3
  j <- floor(rem / d3) + 1
  k <- rem - (j - 1) * d3 + 1
  data.frame(side = side, lr = g_names[as.integer(k)],
             group_i = as.integer(i), group_j = as.integer(j),
             layer = as.integer(k),
             x_old = x_old, x_new = x_new,
             dx = if (is.na(x_new[1L]) && is.na(x_old[1L])) NA_real_ else x_old - x_new,
             stringsAsFactors = FALSE)
}

common_idx <- idx_go[!is.na(match(idx_go, idx_gn))]
common_pos_n <- match(common_idx, idx_gn)
common_pos_o <- match(common_idx, idx_go)
common_rows <- decode_rows(common_idx, "common_diff",
                           x_go[common_pos_o], x_gn[common_pos_n])
only_old_flat <- setdiff(idx_go, idx_gn)
only_new_flat <- setdiff(idx_gn, idx_go)
only_old_rows <- decode_rows(only_old_flat, "only_old",
                             x_go[match(only_old_flat, idx_go)],
                             rep(NA_real_, length(only_old_flat)))
only_new_rows <- decode_rows(only_new_flat, "only_new",
                             rep(NA_real_, length(only_new_flat)),
                             x_gn[match(only_new_flat, idx_gn)])
group_rows <- rbind(common_rows[common_rows$dx != 0, ],
                    only_old_rows, only_new_rows)
group_rows <- group_rows[order(group_rows$side, group_rows$layer,
                               group_rows$group_i, group_rows$group_j), ]

## pval support/value comparison
p_vo <- view_nz(old_p); p_vn <- view_nz(new_p)
ord_po <- order(p_vo$idx); idx_po <- p_vo$idx[ord_po]; x_po <- p_vo$x[ord_po]
ord_pn <- order(p_vn$idx); idx_pn <- p_vn$idx[ord_pn]; x_pn <- p_vn$x[ord_pn]
p_cmp <- num_diff_table(idx_po, x_po, idx_pn, x_pn)

group_rows$flat_ <- (group_rows$group_i - 1) * d2 * d3 +
  (group_rows$group_j - 1) * d3 + group_rows$layer
group_rows$pval_old <- x_po[match(group_rows$flat_, idx_po)]
group_rows$pval_new <- x_pn[match(group_rows$flat_, idx_pn)]
group_rows$flat_ <- NULL

group_summary <- data.frame(
  metric = c("dims", "group_lr_names_equal",
             "group_nonzero_old", "group_nonzero_new",
             "group_stored_sum_old", "group_stored_sum_new",
             "group_common_entries", "group_only_old", "group_only_new",
             "max_abs_value_diff_on_common", "n_common_diff_gt_0",
             "n_common_diff_gt_1e-12",
             "pval_nonzero_old", "pval_nonzero_new",
             "pval_common", "pval_only_old", "pval_only_new",
             "n_pval_common_diff_gt_0", "max_abs_pval_diff"),
  value = c(paste(old_g$dims, collapse = "x"),
            isTRUE(identical(g_names_old, g_names_new)),
            length(idx_go), length(idx_gn),
            sum_fmt(g_o_nz_x <- g_vo$x), sum_fmt(g_vn$x),
            g_cmp$common, g_cmp$only_old, g_cmp$only_new,
            format(g_cmp$max_abs, digits = 15), g_cmp$n_diff_gt_0, g_cmp$n_diff_gt_1e12,
            length(idx_po), length(idx_pn),
            p_cmp$common, p_cmp$only_old, p_cmp$only_new,
            p_cmp$n_diff_gt_0, format(p_cmp$max_abs, digits = 15)),
  stringsAsFactors = FALSE)

write.csv(group_rows, file.path(dump_dir, "group_edge_diff.csv"), row.names = FALSE)
write.csv(group_summary, file.path(dump_dir, "group_summary.csv"), row.names = FALSE)

cat("\n--- group summary ---\n")
print(group_summary, row.names = FALSE)
cat("\n--- differing group edges (", nrow(group_rows), " rows) ---\n", sep = "")
print(group_rows, row.names = FALSE)
cat("\nDONE compare_serial_tensors\n")
