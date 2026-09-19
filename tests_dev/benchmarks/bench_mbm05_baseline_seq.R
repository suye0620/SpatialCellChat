## bench_mbm05_baseline_seq.R —— 情形 B：冻结基线顺序执行 + 过滤复刻 + 与情形 A 产物逐位对比
## 依赖：bench_mbm05_full.R 的产物（bench_out/mbm05_new_{pre,post,meta}.rds）
## 运行：Rscript tests_dev/benchmarks/bench_mbm05_baseline_seq.R（任意 cwd；长任务，建议 nohup/tmux）
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else getwd() }))
source(file.path("tests_dev", "benchmarks", "_common.R"))
source("tests_dev/reference-computeCommunProb-baseline.R")

say("== 情形 B：基线顺序执行 ==")
obj <- readRDS("tests_dev/test_data/MBM05Chat_pre_inference.rds")
sig <- assay(obj, "signaling")
LRsig <- obj@LR$LRsig
nLR <- nrow(LRsig)
meta <- obj@meta

## 垫补：DB cofactor 引用但不在 signaling 子集的基因补零行
## （∏(1+0)=×1 逐位恒等；使冻结参考 cof_factor（无 presence 过滤）与新实现 presence-drop 语义逐位一致）
cof.cols <- grep("cofactor", colnames(obj@DB$cofactor))
key.nms <- unique(c(LRsig$co_A_receptor, LRsig$co_I_receptor, LRsig$agonist, LRsig$antagonist))
key.nms <- key.nms[!is.na(key.nms) & nzchar(key.nms)]
ref.genes <- unique(as.character(unlist(obj@DB$cofactor[key.nms, cof.cols, drop = FALSE],
                                        use.names = FALSE)))
ref.genes <- ref.genes[ref.genes != ""]
missing.cof <- setdiff(ref.genes, rownames(sig))
say("[B] DB cofactor 引用基因 %d | 缺失于 signaling 子集 %d -> 垫零行",
    length(ref.genes), length(missing.cof))
if (length(missing.cof)) {
  pad <- Matrix::sparseMatrix(i = seq_along(missing.cof), j = rep(1L, length(missing.cof)),
                              x = rep(0, length(missing.cof)),
                              dims = c(length(missing.cof), ncol(sig)))
  rownames(pad) <- missing.cof; colnames(pad) <- colnames(sig)
  sig <- rbind(sig, pad)
  say("[B] 垫补后 signaling: %d x %d", nrow(sig), ncol(sig))
}

say("[B] 冻结基线 baseline_reference_commun_prob（scale.distance=9, use.AGAN=TRUE）...")
gc()
t0 <- Sys.time()
oracle <- baseline_reference_commun_prob(
  data.signaling = sig, LR.use = LRsig, coordinates = obj@images$coordinates,
  spatial.factors = obj@images$spatial.factors,
  complex_input = obj@DB$complex, cofactor_input = obj@DB$cofactor, group = obj@idents,
  Kh = 0.5, n = 1, distance.use = TRUE, interaction.range = 250, scale.distance = 9,
  use.AGAN = TRUE, contact.dependent = TRUE, contact.range = 10,
  contact.dependent.forced = FALSE)
t.base <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
say("[B] 基线 run.time = %.1fs (自报 %.1fs) | nLR1 = %d",
    t.base, oracle$run.time, oracle$parameter$nLR1)
gct <- gc(); say("[GC] Vcells max used %.0f MB", gct[2, 3] * 8 / 1024^2)
base.layers <- oracle$net$tmp$prob.cell
rm(oracle$net); invisible(gc())
nnz.base <- vapply(base.layers, function(m) length(m@x), numeric(1))
say("[B] 基线层数 %d | nnz 总 %d", length(base.layers), sum(nnz.base))

meta.new <- readRDS(file.path(bench_out, "mbm05_new_meta.rds"))
stopifnot(identical(meta.new$layer.names, names(base.layers)),
          meta.new$param$nLR1 == oracle$parameter$nLR1)
nC <- ncol(assay(obj, "signaling"))
nLR1 <- oracle$parameter$nLR1
d.spatial <- oracle$res$d.spatial; Matrix::diag(d.spatial) <- 1
adj.contact <- oracle$res$adj.contact
rm(oracle); invisible(gc())

## 基线 filterProbability 复刻（逐位忠实：dgR 等价 CSR + 模式缓存，采样向量按置换原序重建）
say("[B] 构建模式行索引缓存（CSR 行提取）...")
csr.of <- function(m) {
  ord <- order(m@i, method = "radix")
  rp <- as.integer(c(0L, cumsum(tabulate(m@i + 1L, nbins = nrow(m)))))
  jc <- rep.int(seq_len(ncol(m)), diff(m@p))
  list(ord = ord, rp = rp, jc = jc)
}
dR <- csr.of(d.spatial); aR <- csr.of(adj.contact)
pats <- new.env(hash = TRUE)
pat.row <- function(csr, m, ci) {
  key <- as.character(ci)
  v <- pats[[key]]
  if (is.null(v)) {
    s <- csr$rp[ci]; e <- csr$rp[ci + 1L]
    seg <- csr$ord[(s + 1L):e]
    v <- csr$jc[seg][m@x[seg] > 0]
    pats[[key]] <- v
  }
  v
}
set.seed(666L)
permutation <- replicate(nLR, base::sample(x = 1:nC, size = 100, replace = F))

say("[B] 基线过滤复刻循环...")
gc()
t0 <- Sys.time()
base.filt <- lapply(seq_len(nLR), function(i) {
  csr <- if (i <= nLR1) dR else aR
  mat <- if (i <= nLR1) d.spatial else adj.contact
  Pi <- base.layers[[i]]
  Pc <- csr.of(Pi)
  cells <- permutation[, i, drop = TRUE]
  pooled <- unlist(lapply(cells, function(ci) {
    idx <- pat.row(csr, mat, ci)
    vals <- numeric(length(idx))
    s <- Pc$rp[ci]; e <- Pc$rp[ci + 1L]
    if (e > s) {
      seg <- Pc$ord[(s + 1L):e]
      hit <- match(Pc$jc[seg], idx)
      keep <- !is.na(hit)
      vals[hit[keep]] <- Pi@x[seg][keep]
    }
    vals
  }), use.names = FALSE)
  q <- quantile(pooled, probs = 0.95)
  if (q == 0) return(Pi)
  scMatrixTruncation(Pi, cutoff = q, remain.cutoff.v = TRUE, repr = "C")
})
t.filt.base <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
gct <- gc(); say("[GC] Vcells max used %.0f MB", gct[2, 3] * 8 / 1024^2)
nnz.bf <- vapply(base.filt, function(m) length(m@x), numeric(1))
say("[B] 基线过滤复刻: %.1fs | post 总 nnz = %d (保留率 %.1f%%)",
    t.filt.base, sum(nnz.bf), 100 * sum(nnz.bf) / sum(nnz.base))

nC.dbl <- as.double(nC)
layer.diff <- function(a, b) {
  ok <- identical(as.numeric(a@x), as.numeric(b@x)) && identical(a@i, b@i) && identical(a@p, b@p)
  if (ok) return(list(ok = TRUE, maxdiff = 0))
  at <- methods::as(a, "dgTMatrix"); bt <- methods::as(b, "dgTMatrix")
  u <- union(at@i + at@j * nC.dbl, bt@i + bt@j * nC.dbl)
  va <- numeric(length(u)); va[match(at@i + at@j * nC.dbl, u)] <- at@x
  vb <- numeric(length(u)); vb[match(bt@i + bt@j * nC.dbl, u)] <- bt@x
  list(ok = FALSE, maxdiff = max(abs(va - vb)))
}

say("[B] 加载新 post-layers 并逐层对比...")
new.post <- readRDS(file.path(bench_out, "mbm05_new_post.rds"))
post.res <- data.frame(layer = names(base.layers),
  nnz_base_post = nnz.bf, nnz_new_post = vapply(new.post, function(m) length(m@x), numeric(1)),
  identical = FALSE, maxdiff = 0.0, stringsAsFactors = FALSE)
for (k in seq_along(new.post)) {
  r <- layer.diff(base.filt[[k]], new.post[[k]])
  post.res$identical[k] <- r$ok; post.res$maxdiff[k] <- r$maxdiff
}
rm(new.post, base.filt); invisible(gc())
say("[B] POST-FILTER 逐位一致: %d / %d | maxdiff max = %.3g",
    sum(post.res$identical), nrow(post.res), max(post.res$maxdiff))
write.csv(post.res, file.path(bench_out, "mbm05_post_compare.csv"), row.names = FALSE)

say("[B] 加载新 pre-layers 并逐层对比...")
new.pre <- readRDS(file.path(bench_out, "mbm05_new_pre.rds"))
pre.res <- data.frame(layer = names(base.layers),
  nnz_base = nnz.base, nnz_new = vapply(new.pre, function(m) length(m@x), numeric(1)),
  identical = FALSE, maxdiff = 0.0, stringsAsFactors = FALSE)
for (k in seq_along(new.pre)) {
  r <- layer.diff(base.layers[[k]], new.pre[[k]])
  pre.res$identical[k] <- r$ok; pre.res$maxdiff[k] <- r$maxdiff
}
rm(new.pre, base.layers); invisible(gc())
say("[B] PRE-FILTER 逐位一致: %d / %d | maxdiff max = %.3g",
    sum(pre.res$identical), nrow(pre.res), max(pre.res$maxdiff))
write.csv(pre.res, file.path(bench_out, "mbm05_pre_compare.csv"), row.names = FALSE)
gct <- gc(); say("[GC] Vcells max used %.0f MB", gct[2, 3] * 8 / 1024^2)
writeLines(LOG, file.path(bench_out, "mbm05_B_log.txt"))
cat("SCENARIO B DONE\n")
