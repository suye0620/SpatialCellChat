## bench_mbm05_full.R —— 情形 A：新代码顺序执行全量（29,536 细胞 × 1563 LR）
## computeCommunProb（内存模式）→ 保存 pre-layers → filterProbability → 保存 post-layers
## 运行：Rscript tests_dev/benchmarks/bench_mbm05_full.R（任意 cwd；需 test_data/MBM05Chat_pre_inference.rds）
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else getwd() }))
source(file.path("tests_dev", "benchmarks", "_common.R"))

say("== 情形 A：新 computeCommunProb 顺序执行 ==")
obj <- readRDS("tests_dev/test_data/MBM05Chat_pre_inference.rds")
sig <- assay(obj, "signaling")
say("[DATA] signaling: %d x %d (nnz %d) | LRsig %d 对",
    nrow(sig), ncol(sig), length(sig@x), nrow(obj@LR$LRsig))
dpar <- obj@images$.distance$parameters
if (!is.null(dpar)) say("[DATA] 距离缓存参数: %s",
                        paste(names(dpar), dpar, sep = "=", collapse = " "))

gc()
t0 <- Sys.time()
chat.cc <- computeCommunProb(obj, LR.use = obj@LR$LRsig, raw.use = TRUE, Kh = 0.5, n = 1,
                             distance.use = TRUE, interaction.range = 250, scale.distance = 9,
                             use.AGAN = TRUE, contact.dependent = TRUE, contact.range = 10,
                             contact.dependent.forced = FALSE, nthreads = 1L, verbose = TRUE)
t.comp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
gctab <- function(tag) { g <- gc(); say("[GC %s] Vcells used %.0f MB | max used %.0f MB",
  tag, g[2, 2] * 8 / 1024^2, g[2, 3] * 8 / 1024^2) }
gctab("after compute")
param <- chat.cc@misc$.param$communication
layers <- unclass(chat.cc@net$cell$prob)
nnz <- vapply(layers, function(m) length(m@x), numeric(1))
say("[A] computeCommunProb: %.1fs | nLR1=%d | 总 nnz = %d (%.2f GB) | 每层中位 %.0f / 均值 %.0f / 最大 %d",
    t.comp, param$nLR1, sum(nnz), sum(nnz) * 12 / 1024^3, median(nnz), mean(nnz), max(nnz))

saveRDS(layers, file.path(bench_out, "mbm05_new_pre.rds"), compress = FALSE)
saveRDS(list(param = param, run.time = t.comp, layer.names = meta.new.names <- names(layers)),
        file.path(bench_out, "mbm05_new_meta.rds"))
say("[A] pre-layers 已保存（%.2f GB）",
    file.info(file.path(bench_out, "mbm05_new_pre.rds"))$size / 1024^3)

t0 <- Sys.time()
chat.filt <- filterProbability(chat.cc, nboot = 100, seed.use = 666L, thresh = 0.05)
t.filt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
gctab("after filter")
filt.layers <- unclass(chat.filt@net$cell$prob)
nnz.post <- vapply(filt.layers, function(m) length(m@x), numeric(1))
unchanged <- vapply(seq_along(layers), function(k)
  identical(filt.layers[[k]]@x, layers[[k]]@x) && identical(filt.layers[[k]]@i, layers[[k]]@i) &&
    identical(filt.layers[[k]]@p, layers[[k]]@p), logical(1))
say("[A] filterProbability: %.1fs | post 总 nnz = %d (%.2f GB, 保留率 %.1f%%) | cutoff=0 层 %d",
    t.filt, sum(nnz.post), sum(nnz.post) * 12 / 1024^3,
    100 * sum(nnz.post) / sum(nnz), sum(unchanged))
saveRDS(filt.layers, file.path(bench_out, "mbm05_new_post.rds"), compress = FALSE)
gctab("end")
writeLines(LOG, file.path(bench_out, "mbm05_A_log.txt"))
cat("SCENARIO A DONE\n")
