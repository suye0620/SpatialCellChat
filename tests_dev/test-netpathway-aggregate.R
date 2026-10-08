# Wave 2b: computeCommunProbPathway / aggregateNet / relabelSpatialCellChat 迁移
# 等价性契约：离散结果（pathways.sig 排序、count 整数、LR.sig 集合）逐位一致；
# 跨层浮点和按 1-ulp 契约（Plan A 先例），测试断言 ≤1e-12。
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
same <- function(a, b) identical(as.numeric(a), as.numeric(b))

## ---- fixture helpers ----
nC <- 30L
cells <- paste0("c", seq_len(nC))
set.seed(2026)
group <- factor(rep(c("A", "B", "C"), c(12L, 10L, 8L)), levels = c("A", "B", "C"))
names(group) <- cells
expr <- matrix(runif(5 * nC) * (runif(5 * nC) > 0.5), nrow = 5, ncol = nC,
               dimnames = list(paste0("gene", 1:5), cells))
expr <- as(expr, "dgCMatrix")
meta <- data.frame(label = group, row.names = cells)

mk_chat <- function() {
  chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")
  chat@assay$signaling <- expr
  chat@LR$LRsig <- data.frame(
    ligand = c("gene1", "gene3", "gene1", "gene2"),
    receptor = c("gene2", "gene4", "gene4", "gene5"),
    pathway_name = c("P1", "P1", "P2", "P3"),
    row.names = c("LR1", "LR2", "LR3", "LR4")
  )
  chat
}

mk_cell_layer <- function(nnz, seed, n.explicit.zeros = 0L) {
  set.seed(seed)
  idx <- sample.int(nC * nC, nnz)
  j <- (idx - 1L) %/% nC + 1L
  i <- (idx - 1L) %% nC + 1L
  m <- sparseMatrix(i = i, j = j, x = runif(nnz, 0.1, 2), dims = c(nC, nC), giveCsparse = TRUE)
  if (n.explicit.zeros > 0L) m@x[sample.int(nnz, n.explicit.zeros)] <- 0
  m
}

## 手工组级 fixture：K=3、4 层（P1: LR1+LR2, P2: LR3, P3: LR4 全零候选）
gd <- function(vals) as(matrix(vals, 3L, 3L), "dgCMatrix")
g1 <- gd(c(0.5, 0.2, 0.1, 0.3, 0.4, 0.2, 0.1, 0.1, 0.9))
g2 <- gd(c(0.7, 0.1, 0.0, 0.2, 0.3, 0.1, 0.4, 0.2, 0.1))
g3 <- gd(c(0.9, 0.0, 0.2, 0.1, 0.5, 0.3, 0.0, 0.6, 0.2))
g4 <- gd(rep(0, 9))   # P3：全零 pathway 候选
pv1 <- gd(rep(0.01, 9)); pv1[2, 3] <- 0.5          # g1[2,3] 被门控清零
pv2 <- gd(rep(0.01, 9)); pv2[1, 2] <- 0.05         # == thresh 被 >= 门控清零（锚点）
pv3 <- gd(rep(0.01, 9)); pv3[1, 1] <- 0.05         # == thresh
pv4 <- gd(rep(0.01, 9))

gate <- function(pm, pv) { pm <- as.matrix(pm); pv <- as.matrix(pv); pm[pv >= 0.05] <- 0; pm }
P1.ref <- gate(g1, pv1) + gate(g2, pv2)
P2.ref <- gate(g3, pv3)
# 锚点：P1 (1,2) = g1[1,2] 恰好（g2[1,2] 被 pval==thresh 清零，且 g1[1,2] 自身 pval=0.01 保留）
stopifnot(P1.ref[1, 2] == as.matrix(g1)[1, 2], as.matrix(g2)[1, 2] > 0,
          sum(P1.ref) > sum(P2.ref), sum(P2.ref) > 0)

chatA <- mk_chat()
prob.group <- SparseChatArray(list(LR1 = g1, LR2 = g2, LR3 = g3, LR4 = g4))
dimnames(prob.group) <- list(levels(group), levels(group), c("LR1", "LR2", "LR3", "LR4"))
pval.group <- SparseChatArray(list(LR1 = pv1, LR2 = pv2, LR3 = pv3, LR4 = pv4))
dimnames(pval.group) <- list(levels(group), levels(group), c("LR1", "LR2", "LR3", "LR4"))
chatA@net$group$prob <- prob.group
chatA@net$group$pval <- pval.group

cl <- list(LR1 = mk_cell_layer(120L, 11L, n.explicit.zeros = 3L),
           LR2 = mk_cell_layer(60L, 12L),
           LR3 = mk_cell_layer(40L, 13L),
           LR4 = mk_cell_layer(0L, 14L))   # 空层 -> P3 剔除
cell.prob <- SparseChatArray(cl)
dimnames(cell.prob) <- list(cells, cells, names(cl))
chatA@net$cell$prob <- cell.prob

## ---- 1. computeCommunProbPathway do.group：门控/求和/排序/剔除 ----
chatP <- computeCommunProbPathway(chatA)
gp <- chatP@netP$group$prob
check("1a netP$group$prob P1 bitwise == manual gated sum", same(gp[["P1"]], P1.ref))
check("1b netP$group$prob P2 bitwise == manual gated sum", same(gp[["P2"]], P2.ref))
check("1c gating anchor P1(1,2) == g1[1,2] exactly (LR2 zeroed at pval==thresh)",
      identical(as.numeric(gp[1, 2, "P1"]), as.numeric(as.matrix(g1)[1, 2])))
check("1d pathways order P1 before P2 (descending total)",
      identical(chatP@netP$pathways, c("P1", "P2")))
check("1e all-zero pathway P3 excluded", !("P3" %in% chatP@netP$pathways))
check("1f netP$group$prob dimnames (groups, groups, pathways)",
      identical(dimnames(gp)[[1]], levels(group)) &&
      identical(dimnames(gp)[[2]], levels(group)) &&
      identical(dimnames(gp)[[3]], c("P1", "P2")))
check("1g netP$tmp absent", is.null(chatP@netP$tmp))
check("1h validObject after computeCommunProbPathway", validObject(chatP))

## ---- 2. pval NULL -> 不门控 ----
chatB <- chatA; chatB@net$group$pval <- NULL
chatP2 <- computeCommunProbPathway(chatB)
gp2 <- chatP2@netP$group$prob
check("2a pval NULL: P1 == raw g1+g2 (no gating)", same(gp2[["P1"]], as.matrix(g1) + as.matrix(g2)))
check("2b pval NULL: P1(1,2) includes LR2 contribution",
      identical(as.numeric(gp2[1, 2, "P1"]),
                as.numeric(as.matrix(g1)[1, 2]) + as.numeric(as.matrix(g2)[1, 2])))

## ---- 3. do.cell：pathway 聚合与排序 ----
cp <- chatP@netP$cell$prob
check("3a netP$cell$prob P1 == LR1+LR2 cell layers bitwise",
      same(cp[["P1"]], as.matrix(cl$LR1) + as.matrix(cl$LR2)))
check("3b netP$cell$prob P2 == LR3 layer bitwise", same(cp[["P2"]], as.matrix(cl$LR3)))
check("3c pathways.cell == c('P1','P2') descending strength, empty P3 dropped",
      identical(chatP@netP$pathways.cell, c("P1", "P2")))
check("3d cell prob dims nC x nC x 2",
      identical(dim(cp), c(nC, nC, 2L)) &&
      identical(dimnames(cp)[[1]], cells) && identical(dimnames(cp)[[2]], cells))

## ---- 4. do.group=FALSE：cell-only 写回且 netP$group 保持 list ----
chatC <- mk_chat()
chatC@net$cell$prob <- cell.prob
chatP3 <- computeCommunProbPathway(chatC, do.group = FALSE)
check("4a do.group=FALSE: netP$group$prob NULL", is.null(chatP3@netP$group$prob))
check("4b do.group=FALSE: pathways NULL, cell fields present",
      is.null(chatP3@netP[["pathways"]]) && identical(chatP3@netP$pathways.cell, c("P1", "P2")))
check("4c do.group=FALSE: validObject (netP$group remains a list)",
      is.list(chatP3@netP$group) && validObject(chatP3))

## ---- 5. object=NULL 返回 list 形状 ----
np <- computeCommunProbPathway(object = NULL, net = chatA@net, pairLR.use = chatA@LR$LRsig)
check("5a object=NULL returns isomorphic list",
      is.list(np) && identical(np$pathways, c("P1", "P2")) &&
      inherits(np$group$prob, "SparseChatArray") &&
      inherits(np$cell$prob, "SparseChatArray") &&
      identical(np$pathways.cell, c("P1", "P2")))
check("5b object=NULL P1 matches object path bitwise", same(np$group$prob[["P1"]], P1.ref))

## ---- 6. aggregateNet 默认分支：两步门控 count/weight + LR.sig + 细胞级 ----
chatN <- aggregateNet(chatA)
gcnt <- chatN@net$group$count; gwgt <- chatN@net$group$weight
cnt.ref <- matrix(0, 3, 3); wgt.ref <- matrix(0, 3, 3)
for (pm0 in list(list(g1, pv1), list(g2, pv2), list(g3, pv3), list(g4, pv4))) {
  vm <- as.matrix(pm0[[2]]); vm[as.matrix(pm0[[1]]) == 0] <- 1
  pm <- as.matrix(pm0[[1]]); pm[vm >= 0.05] <- 0
  wgt.ref <- wgt.ref + pm
  cnt.ref <- cnt.ref + (pm != 0)
}
check("6a net$group$count bitwise == manual two-step gated counts", same(gcnt, cnt.ref))
check("6b net$group$weight bitwise == manual two-step gated weights", same(gwgt, wgt.ref))
check("6c net$group$count/weight dims + dimnames",
      identical(dimnames(gcnt), list(levels(group), levels(group))) &&
      identical(dimnames(gwgt), list(levels(group), levels(group))))
check("6d net$group$LR.sig == layers with gated sum > 0",
      identical(chatN@net$group$LR.sig, c("LR1", "LR2", "LR3")))
check("6e net$cell$count == stored-entry counts (explicit zeros count as links)", {
  cnt.c <- matrix(0, nC, nC)
  for (m in cl) {
    idx <- cbind(m@i + 1L, rep.int(seq_len(ncol(m)), diff(m@p)))
    if (nrow(idx)) for (t in seq_len(nrow(idx))) cnt.c[idx[t, 1], idx[t, 2]] <- cnt.c[idx[t, 1], idx[t, 2]] + 1
  }
  same(chatN@net$cell$count, cnt.c) &&
    chatN@net$cell$count@x[which(chatN@net$cell$count@x != 0)][1] > 0
})
check("6f net$cell$weight bitwise == plain layer sums", {
  w.c <- matrix(0, nC, nC)
  for (m in cl) w.c <- w.c + as.matrix(m)
  same(chatN@net$cell$weight, w.c)
})
check("6g net$cell$LR.sig excludes empty LR4",
      identical(chatN@net$cell$LR.sig, c("LR1", "LR2", "LR3")))
check("6h legacy net$count/net$weight/net$LR.sig no longer written",
      is.null(chatN@net$count) && is.null(chatN@net$weight) && is.null(chatN@net$LR.sig))
check("6i validObject after aggregateNet", validObject(chatN))

## ---- 7. aggregateNet subset 分支：signaling 过滤（count=行数, weight=prob 和） ----
sub <- aggregateNet(chatA, signaling = "P1", return.object = FALSE)
df.ref <- data.frame(
  source = rep(levels(group)[rep(1:3, each = 3)], times = 2),
  target = rep(levels(group)[rep(1:3, each = 1, times = 3)], times = 2)
)
# 手工复刻：逐层门控（单步 >=） -> pm>0 行 -> P1 层 -> 按 (source,target) 汇总
rows <- list()
for (k in list(list("LR1", g1, pv1), list("LR2", g2, pv2))) {
  pm <- gate(k[[2]], k[[3]])
  idx <- which(pm > 0, arr.ind = TRUE)
  rows[[k[[1]]]] <- data.frame(source = levels(group)[idx[, 1]], target = levels(group)[idx[, 2]],
                               prob = pm[idx], stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)
key <- paste(df$source, df$target, sep = "_")
cnt.ref7 <- tapply(rep(1, nrow(df)), key, sum)
wgt.ref7 <- tapply(df$prob, key, sum)
key2 <- paste(sub$source, sub$target, sep = "_")
check("7a subset count == significant LR rows per pair", {
  ok <- TRUE
  for (kk in names(cnt.ref7)) {
    st <- strsplit(kk, "_")[[1]]
    ok <- ok && identical(as.numeric(sub$group$count[st[1], st[2]]), as.numeric(cnt.ref7[kk]))
  }
  ok && identical(sum(as.numeric(sub$group$count)), sum(as.numeric(cnt.ref7)))
})
check("7b subset weight == sum of significant probs per pair", {
  ok <- TRUE
  for (kk in names(cnt.ref7)) {
    st <- strsplit(kk, "_")[[1]]
    ok <- ok && identical(as.numeric(sub$group$weight[st[1], st[2]]), as.numeric(wgt.ref7[kk]))
  }
  ok
})
check("7c subset weight sums == manual P1 gated total",
      identical(sum(as.numeric(sub$group$weight)), sum(P1.ref)))

## ---- 8. subset sources.use/targets.use（数字与字符） ----
sub2 <- aggregateNet(chatA, signaling = "P1", sources.use = 1, return.object = FALSE)
check("8a numeric sources.use=1 keeps only group A rows",
      identical(rownames(sub2$group$count), "A") && identical(colnames(sub2$group$count), levels(group)))
sub3 <- aggregateNet(chatA, signaling = "P1", sources.use = "A", targets.use = c("B", "C"), return.object = FALSE)
check("8b character sources/targets subset 1x2 block", identical(dim(sub3$group$count), c(1L, 2L)) ||
      (all(rownames(sub3$group$count) %in% "A") && all(colnames(sub3$group$count) %in% c("B", "C"))))
check("8c numeric/character subsets bitwise equal",
      identical(as.numeric(sub2$group$count["A", c("B", "C")]), as.numeric(sub3$group$count["A", c("B", "C")])))
check("8d remove.isolate=FALSE keeps full level dims",
      identical(dim(sub$group$count), c(3L, 3L)))

## ---- 9. subset 空结果 -> 逐字基线报错 ----
err <- tryCatch({
  aggregateNet(chatA, signaling = "P3", return.object = FALSE); NULL
}, error = function(e) conditionMessage(e))
check("9a empty subset result stops with baseline message",
      is.character(err) && grepl("No significant signaling interactions are inferred based on the input!", err, fixed = TRUE))

## ---- 10. relabelSpatialCellChat 端到端 ----
mk_e2e_chat <- function() {
  set.seed(2026)
  nC2 <- 40L
  grp <- factor(rep(c("A", "B", "C"), c(15L, 12L, 13L)))
  expr2 <- matrix(runif(6 * nC2) * (runif(6 * nC2) > 0.5), nrow = 6, ncol = nC2,
                  dimnames = list(paste0("gene", 1:6), paste0("cell", seq_len(nC2))))
  expr2 <- as(expr2, "dgCMatrix")
  meta2 <- data.frame(label = grp, row.names = colnames(expr2))
  chat <- createSpatialCellChat(expr2, meta = meta2, group.by = "label", datatype = "RNA")
  chat@assay$signaling <- expr2
  chat@LR$LRsig <- data.frame(
    ligand = c("gene1", "gene5", "C1"),
    receptor = c("gene2", "gene6", "gene3"),
    pathway_name = c("P1", "P1", "P2"),
    row.names = c("LR1", "LR2", "LR3")
  )
  chat@DB$complex <- data.frame(subunit_1 = "gene5", subunit_2 = "gene6", row.names = "C1")
  layers2 <- list(
    LR1 = mk_cell_layer2(200, 1L, 12L, nC2),
    LR2 = mk_cell_layer2(80, 2L, 5L, nC2),
    LR3 = mk_cell_layer2(0, 3L, 0L, nC2)
  )
  pa <- SparseChatArray(layers2)
  dimnames(pa) <- list(colnames(expr2), colnames(expr2), names(layers2))
  chat@net$cell$prob <- pa
  chat
}
mk_cell_layer2 <- function(nnz, seed, n.zeros, n.d) {
  set.seed(seed)
  idx <- sample.int(n.d * n.d, nnz)
  j <- (idx - 1L) %/% n.d + 1L
  i <- (idx - 1L) %% n.d + 1L
  m <- sparseMatrix(i = i, j = j, x = runif(nnz, 0.1, 2), dims = c(n.d, n.d), giveCsparse = TRUE)
  if (n.zeros > 0L) m@x[sample.int(nnz, n.zeros)] <- 0
  m
}
chatR1 <- mk_e2e_chat()
options(future.globals.maxSize = Inf)
invisible(capture.output(
  chatR1 <- relabelSpatialCellChat(chatR1, labelSet = c("B" = "AB"), nboot = 20)
))
check("10a relabel e2e validObject", validObject(chatR1))
check("10b merged levels in netP$group$prob dimnames",
      identical(dimnames(chatR1@netP$group$prob)[[1]], c("A", "AB", "C")) &&
      identical(dimnames(chatR1@netP$group$prob)[[2]], c("A", "AB", "C")))
check("10c net$group$count/weight present after relabel chain",
      !is.null(chatR1@net$group$count) && !is.null(chatR1@net$group$weight) &&
      identical(dimnames(chatR1@net$group$weight)[[1]], c("A", "AB", "C")))
check("10d no net$tmp / netP$tmp anywhere in relabel chain",
      is.null(chatR1@net$tmp) && is.null(chatR1@netP$tmp))
check("10e netP$pathways present as char vector",
      is.character(chatR1@netP$pathways) && length(chatR1@netP$pathways) >= 1L)

## ---- 11. 确定性：同 seed 重跑 relabel 链逐位一致 ----
chatR2 <- mk_e2e_chat()
invisible(capture.output(
  chatR2 <- relabelSpatialCellChat(chatR2, labelSet = c("B" = "AB"), nboot = 20)
))
check("11a deterministic net$group$weight (bitwise)",
      identical(chatR1@net$group$weight, chatR2@net$group$weight))
check("11b deterministic netP$group$prob (bitwise)",
      identical(chatR1@netP$group$prob, chatR2@netP$group$prob))
check("11c deterministic net$group$count / pathways.sig / LR.sig",
      identical(chatR1@net$group$count, chatR2@net$group$count) &&
      identical(chatR1@netP$pathways, chatR2@netP$pathways) &&
      identical(chatR1@net$group$LR.sig, chatR2@net$group$LR.sig))

cat("All checks passed.\n")
