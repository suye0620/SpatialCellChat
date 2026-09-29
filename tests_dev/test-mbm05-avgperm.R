# MBM05 端到端验证：computeCommunProb(v2) -> filterProbability(v2) -> computeAvgCommunProb(v2 Plan A)
# 1) 全量跑通 + 计时；2) 冻结基线参考在层子样本上对拍（Prob ≤1 ulp、Pval 零差异）；3) 加速比
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({ library(methods); library(Matrix); library(cli) })
source("R/SpatialCellChat_class.R"); source("R/utilities.R"); source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/reference-computeAvgCommunProb-baseline.R")
check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

ckpt <- "tests_dev/test_data/MBM05Chat_pre_inference.rds"
if (!file.exists(ckpt)) stop("checkpoint missing: ", ckpt)
chat <- readRDS(ckpt)
nC <- ncol(assay(chat, "signaling"))
cat("checkpoint: nC =", nC, " idents levels =", nlevels(chat@idents), "\n"); flush.console()

ckpt1 <- ".tmp_mbm05_post_ccp.rds"
if (file.exists(ckpt1)) {
  chat <- readRDS(ckpt1); t1 <- c(elapsed = NA_real_)
  cat("computeCommunProb: restored from checkpoint\n")
} else {
  t1 <- system.time(chat <- computeCommunProb(chat, scale.distance = 9, verbose = FALSE))
  saveRDS(chat, ckpt1, compress = FALSE)
}
cat("computeCommunProb:", t1["elapsed"], "s\n"); flush.console()
# Wave-1 note 已记录的 globals 误报（缓存 env 引用图三重计数）——脚本侧解除上限
options(future.globals.maxSize = Inf)
t2 <- system.time(chat <- filterProbability(chat))
cat("filterProbability:", t2["elapsed"], "s\n"); flush.console()

nLR <- dim(chat@net$cell$prob)[3]
cat("layers:", nLR, " total nnz:", sum(vapply(unclass(chat@net$cell$prob), function(m) length(m@x), numeric(1))), "\n")
t3 <- system.time(chat <- computeAvgCommunProb(chat, nboot = 100, seed.use = 1, verbose = FALSE))
cat("computeAvgCommunProb (Plan A, full):", t3["elapsed"], "s\n"); flush.console()
check("object validates after full pipeline", isTRUE(validateSpatialCellChat(chat)))
prob.group <- chat@net$group$prob
check("group-level output shape", identical(dim(prob.group), c(nlevels(chat@idents), nlevels(chat@idents), nLR)))

# ---- 冻结基线对拍（层子样本）----
set.seed(42)
prob.cell_ <- unclass(chat@net$cell$prob)
nnz.vec <- vapply(prob.cell_, function(m) length(m@x), numeric(1))
cand <- which(nnz.vec > 0)
sample.idx <- sort(sample(cand, min(30L, length(cand))))
group <- chat@idents
LRsig <- dimnames(chat@net$cell$prob)[[3]]
pairLRsig <- chat@LR$LRsig[LRsig, , drop = FALSE]
sig <- assay(chat, "signaling")
complex_subunits <- .sc_complex_subunits(chat@DB$complex)
cache <- new.env(parent = emptyenv())
geneL <- as.character(pairLRsig$ligand); geneR <- as.character(pairLRsig$receptor)

ref.elapsed <- 0
v2.elapsed <- 0
max.rel <- 0
pval.diff <- 0
for (i in sample.idx) {
  rowsL <- .sc_expr_row(geneL[i], sig, complex_subunits, cache)
  rowsR <- .sc_expr_row(geneR[i], sig, complex_subunits, cache)
  tr <- system.time(ref <- computeAvgCommunProb_baseline_ref(
    prob.cell_ = prob.cell_[i], dataLavg = matrix(rowsL, nrow = 1),
    dataRavg = matrix(rowsR, nrow = 1), group = group,
    avg.type = "avg", nboot = 100, seed.use = 1))
  ref.elapsed <- ref.elapsed + tr["elapsed"]
  # v2 单层复算（与主函数同一调用形态；阈值同 .sc_gate_threshold）
  thr <- vapply(tabulate(as.integer(group), nbins = nlevels(group)), .sc_gate_threshold,
                integer(1), min.percent = 0.1)
  gv <- as.integer(group) - 1L
  suppL <- as.numeric(rowsL > 0); suppR <- as.numeric(rowsR > 0)
  m <- prob.cell_[[i]]
  sr_out <- tabulate(m@i + 1L, nbins = nC)
  sr_in <- tabulate(rep.int(seq_len(nC), diff(m@p)), nbins = nC)
  set.seed(1)  # 与主函数同 RNG 序（配对 permutation）
  tv <- system.time({
    v <- cpp_group_avg_obs(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in, gv,
                           nlevels(group), thr, 5, 0L, nC)
    Pnull <- as.vector(prob.group[, , i])
    if (any(Pnull > 0)) {
      permutation <- replicate(100, sample.int(nC, size = nC))
      Pboot <- cpp_group_avg_perm(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in, gv,
                                  nlevels(group), thr, 5, 0L, permutation, 1L)
      p <- rowSums(Pboot - Pnull > 0) / 100
    } else p <- rep(1, length(Pnull))
  })
  v2.elapsed <- v2.elapsed + tv["elapsed"]
  Prob.i <- matrix(v, nlevels(group), nlevels(group))
  rel <- max(abs(Prob.i - ref$Prob[, , 1]) / pmax(abs(ref$Prob[, , 1]), 1e-300))
  max.rel <- max(max.rel, rel)
  Pval.i <- matrix(p, nlevels(group), nlevels(group))
  Pval.i[Prob.i == 0] <- 1
  pval.diff <- pval.diff + sum(Pval.i != ref$Pval[, , 1])
}
check(sprintf("MBM05 %d-layer Prob within 1 ulp of frozen baseline (max rel = %.2e)",
              length(sample.idx), max.rel), max.rel <= 1e-12)
check(sprintf("MBM05 %d-layer Pval zero mismatches", length(sample.idx)), pval.diff == 0)

n.boot.layers <- length(sample.idx)
cat(sprintf("\n==== 加速比（%d 层 × nboot=100） ====\n", n.boot.layers))
cat(sprintf("baseline ref : %8.2f s  (%.1f ms/层)\n", ref.elapsed, 1000 * ref.elapsed / n.boot.layers))
cat(sprintf("全量 Plan A（%d 层）: %.2f s；基线外推全量: %.1f s（%.0f×）\n",
            nLR, t3["elapsed"], ref.elapsed / n.boot.layers * nLR,
            ref.elapsed / n.boot.layers * nLR / max(t3["elapsed"], 1e-9)))
cat("All MBM05 Plan A checks passed.\n")
