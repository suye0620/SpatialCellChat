gs <- read.csv("tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/group_summary.csv")
cat("--- group_summary.csv ---\n")
print(gs, row.names = FALSE)

e <- read.csv("tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/group_edge_diff.csv")
cat("\nrows:", nrow(e), "\n")
print(table(e$side))

co <- e$side == "common_diff"
oo <- e$side == "only_old"
on2 <- e$side == "only_new"

dxco <- abs(e$dx[co])
cat("\n[common_diff] n=", sum(co), " |dx|>0=", sum(dxco > 0),
    " max=", format(max(dxco), digits = 6),
    " median=", format(median(dxco), digits = 6),
    " sum=", format(sum(dxco), digits = 6), "\n", sep = "")
print(quantile(dxco, c(0.5, 0.9, 0.99, 1)))

cat("\n[only_old] n=", sum(oo),
    " x_old>0.05(pval): ", sum(e$pval_old[oo] < 0.05, na.rm = TRUE),
    " x_old max=", format(max(e$x_old[oo]), digits = 6),
    " x_old min=", format(min(e$x_old[oo]), digits = 6), "\n", sep = "")
cat("[only_new] n=", sum(on2),
    " x_new>0.05(pval): ", sum(e$pval_new[on2] < 0.05, na.rm = TRUE),
    " x_new max=", format(max(e$x_new[on2], na.rm = TRUE), digits = 6), "\n", sep = "")

cat("\n[common_diff] pval_old<0.05:", sum(e$pval_old[co] < 0.05, na.rm = TRUE),
    " pval_new<0.05:", sum(e$pval_new[co] < 0.05, na.rm = TRUE), "\n")
cat("[LR families involved in diff rows]:\n")
lr_tab <- sort(table(e$lr), decreasing = TRUE)
print(head(lr_tab, 25))
cat("n distinct LR layers in diff:", length(unique(e$lr)), "\n")

cat("\n[top 12 |dx| common_diff rows]\n")
top <- order(dxco, decreasing = TRUE)[seq_len(min(12, sum(co)))]
idx_co <- which(co)[top]
print(e[idx_co, c("side", "lr", "group_i", "group_j", "x_old", "x_new", "dx", "pval_old", "pval_new")],
      row.names = FALSE)

cat("\n[top 12 only_old rows by x_old]\n")
ord <- order(e$x_old[oo], decreasing = TRUE)[seq_len(min(12, sum(oo)))]
idx_oo <- which(oo)[ord]
print(e[idx_oo, c("side", "lr", "group_i", "group_j", "x_old", "x_new", "pval_old")],
      row.names = FALSE)
