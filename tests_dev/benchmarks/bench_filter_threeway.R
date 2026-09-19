## bench_filter_threeway.R —— filterProbability 三方逐位对比（v2 优化版 / v1 冻结版 / 基线复刻）
## 数据：tests_dev/simulated_data/SpatialChat_2.rds（2000 细胞 × 35 LR）
## 运行：Rscript tests_dev/benchmarks/bench_filter_threeway.R（任意 cwd）
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else getwd() }))
source(file.path("tests_dev", "benchmarks", "_common.R"))
source(file.path("tests_dev", "benchmarks", "filterProbability_v1_body.R"))

obj <- readRDS("tests_dev/simulated_data/SpatialChat_2.rds")
nC <- ncol(obj@data.signaling)
LRsig <- obj@LR$LRsig
coord <- obj@images$coordinates
coord <- coord[colnames(obj@data.signaling), , drop = FALSE]
meta <- obj@meta
rownames(meta) <- as.character(meta$cell_id)
meta <- meta[colnames(obj@data.signaling), , drop = FALSE]
chat <- createSpatialCellChat(object = obj@data.signaling, meta = meta, group.by = "cell_type",
                              input.assay = "norm", datatype = "spatial",
                              coordinates = coord, spatial.factors = obj@images$spatial.factors)
chat@DB <- obj@DB
chat@LR$LRsig <- LRsig
assay(chat, "signaling") <- assay(chat, "norm")
chat.cc <- computeCommunProb(chat, LR.use = LRsig, raw.use = TRUE, Kh = 0.5, n = 1,
                             distance.use = TRUE, interaction.range = 250, scale.distance = 1,
                             use.AGAN = TRUE, contact.dependent = TRUE, contact.range = 10,
                             contact.dependent.forced = FALSE, nthreads = 1L, verbose = FALSE)
prob.in <- unclass(chat.cc@net$cell$prob)
nLR <- nrow(LRsig)
say("[DATA] nC=%d nLR=%d | pre nnz 总 %d", nC, nLR,
    sum(vapply(prob.in, function(m) length(m@x), numeric(1))))

t1 <- system.time(chat.f1 <- filterProbability_v1(chat.cc, nboot = 100, seed.use = 666L, thresh = 0.05))
l1 <- unclass(chat.f1@net$cell$prob)
say("[v1] %.2fs | post nnz %d", t1["elapsed"], sum(vapply(l1, function(m) length(m@x), numeric(1))))

t2 <- system.time(chat.f2 <- filterProbability(chat.cc, nboot = 100, seed.use = 666L, thresh = 0.05))
l2 <- unclass(chat.f2@net$cell$prob)
say("[v2] %.2fs | post nnz %d", t2["elapsed"], sum(vapply(l2, function(m) length(m@x), numeric(1))))
say("[CHECK] 输入对象未被修改: %s",
    identical(lapply(prob.in, function(m) as.numeric(m@x)),
              lapply(unclass(chat.cc@net$cell$prob), function(m) as.numeric(m@x))))

quantile.prob <- 0.95
set.seed(666L)
d.spatial <- chat.cc@images$.distance$d.spatial; Matrix::diag(d.spatial) <- 1
adj.contact <- chat.cc@images$.distance$adj.contact
nLR1 <- chat.cc@misc$.param$communication$nLR1
permutation <- replicate(nLR, base::sample(x = 1:nC, size = 100, replace = F))
t3 <- system.time(ref <- lapply(seq_len(nLR), function(i) {
  d_spatial <- if (i <= nLR1) d.spatial else adj.contact
  Prob.cell.i <- prob.in[[i]]
  sv <- purrr::map(.x = permutation[, i, drop = TRUE], .f = function(ci) {
    idx <- which(d_spatial[ci, , drop = TRUE] > 0)
    Prob.cell.i[ci, idx, drop = TRUE]
  }) %>% unlist()
  q <- quantile(sv, probs = quantile.prob)
  if (q == 0) return(Prob.cell.i)
  scMatrixTruncation(Prob.cell.i, cutoff = q, remain.cutoff.v = TRUE, repr = "C")
}))
say("[REF] 基线复刻 %.2fs", t3["elapsed"])

cmp <- function(a, b) identical(a@x, b@x) && identical(a@i, b@i) && identical(a@p, b@p)
r21 <- vapply(seq_len(nLR), function(k) cmp(l2[[k]], l1[[k]]), logical(1))
r2r <- vapply(seq_len(nLR), function(k) cmp(l2[[k]], ref[[k]]), logical(1))
r1r <- vapply(seq_len(nLR), function(k) cmp(l1[[k]], ref[[k]]), logical(1))
say("[RESULT] v2 vs v1   逐位一致: %d/%d", sum(r21), nLR)
say("[RESULT] v2 vs 基线  逐位一致: %d/%d", sum(r2r), nLR)
say("[RESULT] v1 vs 基线  逐位一致: %d/%d", sum(r1r), nLR)
say("[RESULT] 计时: v1 %.2fs | v2 %.2fs | 提升 %.1fx",
    t1["elapsed"], t2["elapsed"], t1["elapsed"] / t2["elapsed"])
writeLines(LOG, file.path(bench_out, "filter_threeway_log.txt"))
cat("THREE-WAY DONE\n")
