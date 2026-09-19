## profile_filter_components.R —— filterProbability 组件级剖析（定位瓶颈用）
## 依赖：test_data/MBM05Chat_pre_inference.rds + bench_out/mbm05_new_pre.rds（可选，重层取自情形 A 产物）
## 运行：Rscript tests_dev/benchmarks/profile_filter_components.R
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else getwd() }))
source(file.path("tests_dev", "benchmarks", "_common.R"))

obj <- readRDS("tests_dev/test_data/MBM05Chat_pre_inference.rds")
d <- obj@images$.distance$d.spatial
nC <- ncol(d)
say("[D] d.spatial nnz = %d, nC = %d", length(d@x), nC)

t1 <- system.time(rs <- Matrix::rowSums(d != 0))
say("[T1] rowSums(d != 0)（v1 每层重复）: %.3fs", t1["elapsed"])
t1b <- system.time(rs2 <- tabulate(d@i[d@x != 0] + 1L, nbins = nC))
say("[T1b] tabulate 一次性替代: %.4fs（计数一致: %s）",
    t1b["elapsed"], identical(as.numeric(rs), as.numeric(rs2)))

csr.of <- function(m) {
  ord <- order(m@i, method = "radix")
  rp <- as.integer(c(0L, cumsum(tabulate(m@i + 1L, nbins = nrow(m)))))
  jc <- rep.int(seq_len(ncol(m)), diff(m@p))
  list(ord = ord, rp = rp, jc = jc)
}

pre.path <- file.path(bench_out, "mbm05_new_pre.rds")
if (file.exists(pre.path)) {
  layers <- readRDS(pre.path)
  nnz <- vapply(layers, function(m) length(m@x), numeric(1))
  k <- which.max(nnz)
  P <- layers[[k]]
  say("[D] 最重层 #%d nnz = %d", k, nnz[k])
  cells <- sample(1:nC, 100)

  t2 <- system.time(sx <- P[cells, , drop = FALSE]@x)
  say("[T2] dgC[100 行子集]@x: %.3fs", t2["elapsed"])

  t3a <- system.time(idxs <- lapply(cells, function(ci) which(d[ci, , drop = TRUE] > 0)))
  say("[T3a] which(d[ci,]>0) ×100（v1 d 侧单行提取）: %.3fs", t3a["elapsed"])
  t3b <- system.time(pl <- unlist(lapply(seq_along(cells), function(jj)
    P[cells[jj], idxs[[jj]], drop = TRUE]), use.names = FALSE))
  say("[T3b] P[ci, idx] ×100（v1 P 侧单行提取）: %.3fs", t3b["elapsed"])

  t4a <- system.time(cP <- csr.of(P))
  t4b <- system.time(cD <- csr.of(d))
  pats <- new.env(hash = TRUE)
  t4c <- system.time(pl2 <- unlist(lapply(seq_along(cells), function(jj) {
    ci <- cells[jj]
    idx <- pats[[as.character(ci)]]
    if (is.null(idx)) {
      s0 <- cD$rp[ci]; e0 <- cD$rp[ci + 1L]
      idx <- cD$jc[cD$ord[(s0 + 1L):e0]]
      pats[[as.character(ci)]] <- idx
    }
    vals <- numeric(length(idx))
    s <- cP$rp[ci]; e <- cP$rp[ci + 1L]
    if (e > s) {
      seg <- cP$ord[(s + 1L):e]
      hit <- match(cP$jc[seg], idx)
      keep <- !is.na(hit)
      vals[hit[keep]] <- P@x[seg][keep]
    }
    vals
  }), use.names = FALSE))
  say("[T4] CSR 构建: P %.3fs + d %.3fs | 缓存版 gather ×100: %.3fs",
      t4a["elapsed"], t4b["elapsed"], t4c["elapsed"])
  say("[T4b] 两条采样向量 identical: %s（逐位等价证据）", identical(pl, pl2))
} else {
  say("[SKIP] 未找到情形 A 产物（先跑 bench_mbm05_full.R 可剖析重层）")
}
writeLines(LOG, file.path(bench_out, "profile_log.txt"))
cat("PROFILE DONE\n")
