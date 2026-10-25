e <- read.csv("tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/group_edge_diff.csv")

## 1. sum gap decomposition
co <- e$side == "common_diff"
oo <- e$side == "only_old"
on2 <- e$side == "only_new"
cat("stored_sum_old - stored_sum_new =", 2.64022099022303 - 2.63671247907007, "\n")
cat("sum only_old x =", sum(e$x_old[oo]), "\n")
cat("sum only_new x =", sum(e$x_new[on2]), "\n")
cat("gap explained by support:", sum(e$x_old[oo]) - sum(e$x_new[on2]),
    " of which old-side:", sum(e$x_old[oo]), "\n\n")

pv <- data.frame(lr = e$lr[co], group_i = e$group_i[co], group_j = e$group_j[co],
                 x_old = e$x_old[co], x_new = e$x_new[co],
                 pval_old = e$pval_old[co], pval_new = e$pval_new[co],
                 pdiff = e$pval_new[co] - e$pval_old[co])
big <- order(abs(pv$pdiff), decreasing = TRUE)[seq_len(12)]
print(pv[big, ], row.names = FALSE)

## 3. significance direction changes on common support
sig_old <- pv$pval_old < 0.05
sig_new <- pv$pval_new < 0.05
cat("\nsignificant (=pval<0.05) on common support: old", sum(sig_old),
    " new", sum(sig_new), " both", sum(sig_old & sig_new),
    " old-only", sum(sig_old & !sig_new), " new-only", sum(!sig_old & sig_new), "\n")

## 4. only_old / only_new pval values
pv_only_old <- data.frame(lr = e$lr[oo], i = e$group_i[oo], j = e$group_j[oo],
                          x = e$x_old[oo], pval = e$pval_old[oo])
pv_only_new <- data.frame(lr = e$lr[on2], i = e$group_i[on2], j = e$group_j[on2],
                          x = e$x_new[on2], pval = e$pval_new[on2])
cat("\nonly_old pval summary: <0.05:", sum(pv_only_old$pval < 0.05, na.rm = TRUE),
    " [==1:", sum(pv_only_old$pval == 1, na.rm = TRUE),
    " of ", sum(!is.na(pv_only_old$pval)), "]\n")
print(head(pv_only_old[order(-pv_only_old$x), ], 10), row.names = FALSE)
cat("\nonly_new pval summary: <0.05:", sum(pv_only_new$pval < 0.05, na.rm = TRUE), "\n")
print(head(pv_only_new[order(-pv_only_new$x), ], 10), row.names = FALSE)

## 5. groups involved
cat("\ngroup pairs (i,j) coverage of ALL diff rows:\n")
ij <- paste(e$group_i, e$group_j, sep = "->")
print(sort(table(ij), decreasing = TRUE)[seq_len(min(12, length(unique(ij))))])
cat("distinct (i,j):", length(unique(ij)), "of 64\n")
