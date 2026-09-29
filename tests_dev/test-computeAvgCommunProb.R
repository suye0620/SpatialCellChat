# Plan A: computeAvgCommunProb v2 —— kernel / R 语义规范 / 冻结基线 对拍
# 等价性契约见 agent note 2026-09-28 与 docs/PERMUTATION_TEST_AUDIT.md（统计语义审计）。
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/reference-computeAvgCommunProb-baseline.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
same <- function(a, b) identical(as.numeric(a), as.numeric(b))

## ---- fixture helpers ----
make_layer <- function(nC, n.store, seed, n.explicit.zeros = 0L, values = NULL) {
  set.seed(seed)
  idx <- sample.int(nC * nC, n.store)
  j <- (idx - 1L) %/% nC + 1L
  i <- (idx - 1L) %% nC + 1L
  x <- if (is.null(values)) runif(n.store, 0.01, 5) else values
  m <- sparseMatrix(i = i, j = j, x = x, dims = c(nC, nC), giveCsparse = TRUE)
  if (n.explicit.zeros > 0) {
    zi <- sample.int(n.store, n.explicit.zeros)
    m@x[zi] <- 0   # 注入存储型显式零：基线 binarize 计数含它，kernel 同
  }
  m
}

make_group <- function(sizes, labels) factor(rep(labels, times = sizes))

## R 语义规范：与 kernel 相同的遍历序（CSC 列升序、列内行升序）与乘加分组，
## 逐字标量循环累积 -> 可做 identical 逐位对拍。门控用整数阈值（与 kernel 同）。
.spec_boot <- function(m, suppL, suppR, sr_out, sr_in, g, K, thr, min.cells.sr, do_avg) {
  nC <- ncol(m)
  cntL <- numeric(K); cntR <- numeric(K); Sout <- numeric(K); Sin <- numeric(K)
  for (c in seq_len(nC)) {
    a <- g[c] + 1L
    cntL[a] <- cntL[a] + suppL[c]
    cntR[a] <- cntR[a] + suppR[c]
    Sout[a] <- Sout[a] + sr_out[c]
    Sin[a]  <- Sin[a]  + sr_in[c]
  }
  gL <- as.numeric(cntL >= thr); gR <- as.numeric(cntR >= thr)
  gS <- as.numeric(Sout >= min.cells.sr); gI <- as.numeric(Sin >= min.cells.sr)
  num <- numeric(K * K); den <- numeric(K * K)
  ii <- m@i                                   # 0-based row
  jj0 <- rep.int(seq_len(nC) - 1L, diff(m@p)) # 0-based col (CSC 序)
  for (t in seq_along(m@x)) {
    f <- g[ii[t] + 1L] + g[jj0[t] + 1L] * K + 1L
    num[f] <- num[f] + m@x[t]
    den[f] <- den[f] + 1
  }
  val <- if (do_avg) num / den else num
  val[is.nan(val)] <- 0
  gg <- as.numeric((gL * gS) %o% (gR * gI))   # f = a + b*K 布局（sender 最快）
  val * gg
}

## ---- 1. kernel vs R 语义规范：obs + perm，avg/sum，门控退化分支，nthreads ----
combos <- 0L
for (K in c(1L, 3L, 6L)) {
  for (nC in c(7L, 24L)) {
    sizes <- if (K == 1L) nC else {
      s <- rep(1L, K); s[1] <- nC - (K - 1L); s
    }
    group <- make_group(sizes, paste0("g", seq_len(K)))
    g <- as.integer(group) - 1L
    thr <- vapply(tabulate(as.integer(group), nbins = K), .sc_gate_threshold,
                  integer(1), min.percent = 0.1)
    for (avg_sum in 0:1) {
      for (spec in list(
        list(n.store = 30, zeros = 0L),
        list(n.store = 24, zeros = 6L),   # 存储型显式零
        list(n.store = 0,  zeros = 0L)    # 空层
      )) {
        m <- make_layer(nC, spec$n.store, seed = combos + 1L, n.explicit.zeros = spec$zeros)
        suppL <- as.numeric(runif(nC) > 0.4)
        suppR <- as.numeric(runif(nC) > 0.4)
        sr_out <- tabulate(m@i + 1L, nbins = nC)
        sr_in <- tabulate(rep.int(seq_len(nC), diff(m@p)), nbins = nC)
        # 观测标签
        v.k <- cpp_group_avg_obs(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                                 g, K, thr, 5, avg_sum, nC)
        v.s <- .spec_boot(m, suppL, suppR, sr_out, sr_in, g, K, thr, 5, avg_sum == 0)
        combos <- combos + 1L
        check(sprintf("obs kernel==spec (K=%d nC=%d avg=%d zeros=%d nnz=%d)",
                      K, nC, avg_sum, spec$zeros, length(m@x)),
              same(v.k, v.s))
        # 置换 kernel：与 spec 逐 boot 对拍 + nthreads 一致性
        set.seed(99)
        perm <- replicate(4L, sample.int(nC, nC))
        P.k <- cpp_group_avg_perm(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                                  g, K, thr, 5, avg_sum, perm, 1L)
        P.k4 <- cpp_group_avg_perm(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                                   g, K, thr, 5, avg_sum, perm, 4L)
        check(sprintf("perm nthreads identical (K=%d nC=%d avg=%d)", K, nC, avg_sum),
              same(P.k, P.k4))
        ok <- TRUE
        for (b in 1:4) {
          v.pb <- .spec_boot(m, suppL, suppR, sr_out, sr_in, g[perm[, b]],
                             K, thr, 5, avg_sum == 0)
          if (!same(P.k[, b], v.pb)) ok <- FALSE
        }
        check(sprintf("perm kernel==spec per boot (K=%d nC=%d avg=%d)", K, nC, avg_sum), ok)
        # 恒等置换列 == 观测结果
        v.id <- cpp_group_avg_perm(m@x, m@i, m@p, suppL, suppR, sr_out, sr_in,
                                   g, K, thr, 5, avg_sum, matrix(seq_len(nC), ncol = 1L), 1L)
        check(sprintf("perm(identity)==obs (K=%d nC=%d avg=%d)", K, nC, avg_sum),
              same(v.id[, 1], v.k))
      }
    }
  }
}

## ---- 2. .sc_gate_threshold：对 fixture 组大小穷举验证单调性 + 阈值正确性 ----
for (n in 1:80) {
  for (mp in c(0.05, 0.1, 0.2, 0.5)) {
    dec <- vapply(0:n, function(cnt)
      isTRUE(signif(mean(rep(c(1, 0), c(cnt, n - cnt))), 1) >= mp), logical(1))
    if (!all(diff(as.numeric(dec)) >= 0))
      stop("[FAIL] threshold monotonicity broken: n=", n, " mp=", mp, call. = FALSE)
    thr <- .sc_gate_threshold(n, mp)
    expect <- if (any(dec)) which(dec)[1] - 1L else n + 1L
    if (thr != expect)
      stop("[FAIL] threshold mismatch: n=", n, " mp=", mp, " thr=", thr,
           " expect=", expect, call. = FALSE)
  }
}
check("gate threshold exhaustive n=1..80 x mp {0.05,0.1,0.2,0.5}", TRUE)
# 大 n 抽样验证（分层 + 阈值邻域全测）
for (n in c(1000L, 29536L)) {
  for (mp in c(0.05, 0.1, 0.2)) {
    thr <- .sc_gate_threshold(n, mp)
    pts <- sort(unique(c(0, n, thr + (-3:3),
                         round(seq(0, n, length.out = 120)))))
    pts <- pts[pts >= 0 & pts <= n]
    dec <- vapply(pts, function(cnt)
      isTRUE(signif(mean(rep(c(1, 0), c(cnt, n - cnt))), 1) >= mp), logical(1))
    ok <- all(dec == (pts >= thr))
    if (!ok) stop("[FAIL] large-n threshold inconsistent: n=", n, " mp=", mp, call. = FALSE)
  }
}
check("gate threshold large-n stratified probes (n=1000, 29536)", TRUE)

## ---- 3. 冻结基线端到端对拍（最小 11-slot 对象）----
set.seed(2026)
nC <- 40L
k <- 3L
nLR <- 3L
group <- make_group(c(15L, 12L, 13L), c("A", "B", "C"))
expr <- matrix(runif(6 * nC) * (runif(6 * nC) > 0.5), nrow = 6, ncol = nC,
               dimnames = list(paste0("gene", 1:6), paste0("cell", seq_len(nC))))
expr <- as(expr, "dgCMatrix")
meta <- data.frame(label = group, row.names = colnames(expr))
chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")
sig <- expr
chat@assay$signaling <- sig
LRtab <- data.frame(ligand = c("gene1", "gene5", "C1"),
                    receptor = c("gene2", "gene6", "gene3"),
                    row.names = c("LR1", "LR2", "LR3"))
chat@LR$LRsig <- LRtab
chat@DB$complex <- data.frame(subunit_1 = "gene5", subunit_2 = "gene6",
                              row.names = "C1")
layers <- list(
  LR1 = make_layer(nC, 200, seed = 1L, n.explicit.zeros = 12L),
  LR2 = make_layer(nC, 80,  seed = 2L, n.explicit.zeros = 5L),
  LR3 = make_layer(nC, 0,   seed = 3L)   # 空层：应被 prob.sum>0 剔除
)
prob.array <- SparseChatArray(layers)
dimnames(prob.array) <- list(colnames(expr), colnames(expr), names(layers))
chat@net$cell$prob <- prob.array

# 参考输入：expr 行即旧 tmp$Lavg/Ravg 的 0/1 等价替身（门控内 `1*(x>0)` 逐位一致）
complex_subunits <- .sc_complex_subunits(chat@DB$complex)
cache <- new.env(parent = emptyenv())
rowsL <- t(vapply(seq_len(nLR), function(i)
  .sc_expr_row(as.character(LRtab$ligand[i]), sig, complex_subunits, cache), numeric(nC)))
rowsR <- t(vapply(seq_len(nLR), function(i)
  .sc_expr_row(as.character(LRtab$receptor[i]), sig, complex_subunits, cache), numeric(nC)))

for (avg.type in c("avg", "sum")) {
  ref <- computeAvgCommunProb_baseline_ref(
    prob.cell_ = layers, dataLavg = rowsL, dataRavg = rowsR, group = group,
    avg.type = avg.type, nboot = 20, seed.use = 1)
  chat2 <- computeAvgCommunProb(chat, avg.type = avg.type, nboot = 20, seed.use = 1,
                                verbose = FALSE)
  Prob.new <- array(unlist(lapply(unclass(chat2@net$group$prob), as.matrix)),
                    dim = c(k, k, nLR))
  Pval.new <- array(unlist(lapply(unclass(chat2@net$group$pval), as.matrix)),
                    dim = c(k, k, nLR))
  # 等价契约：门控/计数/模式逐位一致；分块和受基线 BLAS dgemm 归约序不可移植性
    # 约束，允许 ≤1 ulp 级浮点差（实测 ~3.6e-16 相对差）；p 值为离散计数，必须零差异。
  rel <- max(abs(Prob.new - ref$Prob) / pmax(abs(ref$Prob), 1e-300))
  check(sprintf("end-to-end Prob within 1 ulp of frozen baseline (%s): %.2e", avg.type, rel),
        rel <= 1e-12)
  check(sprintf("end-to-end Pval identical vs frozen baseline (%s)", avg.type),
        same(drop(Pval.new), drop(ref$Pval)))
  check(sprintf("object validates (%s)", avg.type),
        isTRUE(validateSpatialCellChat(chat2)))
  check(sprintf("params recorded (%s)", avg.type),
        identical(chat2@misc$.param$averaging$avg.type, avg.type) &&
          !is.null(chat2@misc$.param$averaging$LRsig.GGC.counts))
}

# SparseChatArray 逐层拍平（按存储序拼接 @x；空层贡献零长度段）
flat.s <- function(arr) as.numeric(unlist(lapply(unclass(arr), function(m) m@x)))
# do.permutation = FALSE：pval 槽为 NULL，prob 仍写入
chat3 <- computeAvgCommunProb(chat, nboot = 20, seed.use = 1,
                              do.permutation = FALSE, verbose = FALSE)
check("no-permutation leaves group$pval NULL",
      is.null(chat3@net$group$pval) && !is.null(chat3@net$group$prob))

chat4 <- computeAvgCommunProb(chat, avg.type = "avg", nboot = 20, seed.use = 1,
                              verbose = FALSE)
chat5 <- computeAvgCommunProb(chat, avg.type = "avg", nboot = 20, seed.use = 1,
                              nthreads = 4L, verbose = FALSE)
check("deterministic re-run + nthreads=4 bitwise identical",
      same(flat.s(chat4@net$group$prob), flat.s(chat5@net$group$prob)) &&
        same(flat.s(chat4@net$group$pval), flat.s(chat5@net$group$pval)))

# 空层被正确剔除（LR3 无存储条目 -> GGC 计数 0 -> pval 恒 1）
check("empty layer gets pval 1 and prob 0",
      all(chat4@net$group$pval[[3]]@x == 1) && length(chat4@net$group$prob[[3]]@x) == 0L)

## ---- 4. 拒绝路径 ----
chat.bad <- chat; chat.bad@net$cell$prob <- NULL
check("missing net$cell$prob rejected",
      inherits(try(computeAvgCommunProb(chat.bad, verbose = FALSE), silent = TRUE),
               "try-error"))
check("bad group.by rejected",
      inherits(try(computeAvgCommunProb(chat, group.by = "nope", verbose = FALSE),
                   silent = TRUE), "try-error"))
check("nboot=0 with permutation rejected",
      inherits(try(computeAvgCommunProb(chat, nboot = 0, verbose = FALSE), silent = TRUE),
               "try-error"))

cat("All computeAvgCommunProb Plan A checks passed.\n")
