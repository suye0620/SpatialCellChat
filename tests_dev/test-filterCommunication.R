# Wave 2a: filterCommunication 迁移（net$group / net$cell + SparseChatArray 逐层操作）
# 等价口径：过滤规则（min.cells / min.links / min.cells.sr）与旧代码语义一致，
# 写回路径从 net$prob/net$pval/net$prob.cell/net$tmp$prob.cell 切到
# net$group$prob/net$group$pval/net$cell$prob。
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

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

## ---- fixture: 30 细胞 3 组（A=15, B=12, C=3），2 个 LR 层 ----
set.seed(2026)
nC <- 30L
group <- factor(rep(c("A", "B", "C"), c(15L, 12L, 3L)))
cells <- paste0("c", seq_len(nC))
names(group) <- cells
expr <- Matrix::rsparsematrix(5, nC, 0.4)
expr <- as(expr, "dgCMatrix")
rownames(expr) <- paste0("gene", seq_len(5))
colnames(expr) <- cells
meta <- data.frame(label = group, row.names = cells)

chat <- createSpatialCellChat(expr, meta = meta, group.by = "label",
                              datatype = "RNA")

mk_layer <- function(nnz, seed) {
  set.seed(seed)
  idx <- sample.int(nC * nC, nnz)
  j <- (idx - 1L) %/% nC + 1L
  i <- (idx - 1L) %% nC + 1L
  sparseMatrix(i = i, j = j, x = runif(nnz, 0.1, 1), dims = c(nC, nC))
}
mk_group_layer <- function(nnz, seed) {
  set.seed(seed)
  idx <- sample.int(9L, 9L)[seq_len(nnz)]   # K=3 组级 3x3，确定性取 nnz 个不同坐标
  j <- (idx - 1L) %/% 3L + 1L
  i <- (idx - 1L) %% 3L + 1L
  sparseMatrix(i = i, j = j, x = runif(nnz, 0.1, 1), dims = c(3L, 3L))
}

group.layers <- list(LR1 = mk_group_layer(8L, 1L), LR2 = mk_group_layer(6L, 2L))
prob.group <- SparseChatArray(group.layers)
dimnames(prob.group) <- list(levels(group), levels(group), names(group.layers))
pval.layers <- list(LR1 = mk_group_layer(8L, 3L), LR2 = mk_group_layer(6L, 4L))
pval.group <- SparseChatArray(pval.layers)
dimnames(pval.group) <- list(levels(group), levels(group), names(pval.layers))

chat@net$group$prob <- prob.group
chat@net$group$pval <- pval.group

cell.layers <- list(
  LR1 = mk_layer(200L, 11L),   # 丰富层：不应被过滤
  LR2 = mk_layer(3L, 12L),     # nnz=3 < min.links=5：被 min.links 剔除
  LR3 = NULL                   # 50 条目全部来自 2 个 sender：仅被 min.cells.sr 剔除
)
set.seed(13L)
sj <- sample.int(nC, 50L, replace = TRUE)
cell.layers$LR3 <- sparseMatrix(i = rep(1:2, each = 25L), j = sj,
                                x = runif(50L, 0.1, 1), dims = c(nC, nC))
prob.cell <- SparseChatArray(cell.layers)
dimnames(prob.cell) <- list(cells, cells, names(cell.layers))
chat@net$cell$prob <- prob.cell

## ---- 1. 组级过滤：min.cells=10 剔除 C 组（3 个细胞）----
out <- filterCommunication(chat, min.cells = 10, min.links = NULL, min.cells.sr = NULL)

check("group filter returns SparseChatArray",
      inherits(out@net$group$prob, "SparseChatArray"))
check("group filter keeps dimnames",
      identical(dimnames(out@net$group$prob), dimnames(chat@net$group$prob)))
check("group filter C row zeroed",
      all(out@net$group$prob[[1]][3L, ] == 0) && all(out@net$group$prob[[2]][3L, ] == 0))
check("group filter C col zeroed",
      all(out@net$group$prob[[1]][, 3L] == 0) && all(out@net$group$prob[[2]][, 3L] == 0))
check("group filter A-B block preserved",
      identical(out@net$group$prob[[1]][1:2, 1:2], chat@net$group$prob[[1]][1:2, 1:2]))
check("group filter pval synced (prob==0 -> 1)",
      all(out@net$group$pval[[1]][3L, ] == 1) && all(out@net$group$pval[[1]][, 3L] == 1))
# 旧语义：先清零 prob 的 C 行列，再用清零后的 prob 做全矩阵掩码 pval[prob==0] <- 1
ref.prob1 <- as.matrix(chat@net$group$prob[[1]])
ref.prob1[3L, ] <- 0
ref.prob1[, 3L] <- 0
ref.pval1 <- as.matrix(chat@net$group$pval[[1]])
ref.pval1[ref.prob1 == 0] <- 1
check("group filter pval matches baseline full-matrix masking",
      identical(as.matrix(out@net$group$pval[[1]]), ref.pval1))
check("group filter does not mutate input",
      length(chat@net$group$prob[[1]]@x) == 8L)

## ---- 2. 无 pval（do.permutation=FALSE 场景）不报错 ----
chat.nopval <- chat
chat.nopval@net$group$pval <- NULL
out2 <- filterCommunication(chat.nopval, min.cells = 10, min.links = NULL, min.cells.sr = NULL)
check("group filter without pval works",
      inherits(out2@net$group$prob, "SparseChatArray") && is.null(out2@net$group$pval))

## ---- 3. 无 group$prob 拒绝 ----
chat.bad <- chat
chat.bad@net$group$prob <- NULL
err <- tryCatch({ filterCommunication(chat.bad, min.cells = 10); NULL },
                error = function(e) conditionMessage(e))
check("rejection: group filter without group$prob errors",
      !is.null(err) && grepl("computeAvgCommunProb", err))

## ---- 4. 细胞级过滤：min.links + min.cells.sr ----
out3 <- filterCommunication(chat, min.cells = NULL, min.links = 5, min.cells.sr = 5)
check("cell filter returns SparseChatArray",
      inherits(out3@net$cell$prob, "SparseChatArray"))
check("cell filter keeps dimnames",
      identical(dimnames(out3@net$cell$prob), dimnames(chat@net$cell$prob)))
check("cell filter min.links removes sparse layer",
      length(out3@net$cell$prob[[2]]@x) == 0L)
check("cell filter min.cells.sr removes low-sender layer",
      length(out3@net$cell$prob[[3]]@x) == 0L)
check("cell filter rich layer untouched",
      identical(out3@net$cell$prob[[1]]@x, chat@net$cell$prob[[1]]@x) &&
        identical(out3@net$cell$prob[[1]]@i, chat@net$cell$prob[[1]]@i))
check("cell filter does not mutate input",
      length(chat@net$cell$prob[[2]]@x) == 3L)

## ---- 5. 只传 min.links（无 min.cells.sr）----
out4 <- filterCommunication(chat, min.cells = NULL, min.links = 5, min.cells.sr = NULL)
check("min.links only works",
      length(out4@net$cell$prob[[2]]@x) == 0L &&
        length(out4@net$cell$prob[[3]]@x) > 0L)

## ---- 6. 只传 min.cells.sr（无 min.links）----
out5 <- filterCommunication(chat, min.cells = NULL, min.links = NULL, min.cells.sr = 5)
check("min.cells.sr only works",
      length(out5@net$cell$prob[[3]]@x) == 0L &&
        length(out5@net$cell$prob[[2]]@x) == 0L)  # LR2 仅 3 条目 -> sender<=3 < 5，基线语义同样剔除

## ---- 7. 无 cell$prob 拒绝 ----
chat.bad2 <- chat
chat.bad2@net$cell$prob <- NULL
err2 <- tryCatch({ filterCommunication(chat.bad2, min.cells = NULL, min.links = 5); NULL },
                 error = function(e) conditionMessage(e))
check("rejection: cell filter without cell$prob errors",
      !is.null(err2) && grepl("computeCommunProb", err2))

## ---- 8. 全部 NULL 参数 = no-op ----
out6 <- filterCommunication(chat, min.cells = NULL, min.links = NULL, min.cells.sr = NULL)
check("all-NULL params no-op",
      identical(out6@net$group$prob, chat@net$group$prob) &&
        identical(out6@net$cell$prob, chat@net$cell$prob))

## ---- 9. 对象整体合法 ----
out.full <- filterCommunication(chat, min.cells = 10, min.links = 5, min.cells.sr = 5)
check("full filter validates", isTRUE(validateSpatialCellChat(out.full)))
check("full filter both levels applied",
      all(out.full@net$group$prob[[1]][3L, ] == 0) &&
        length(out.full@net$cell$prob[[2]]@x) == 0L &&
        length(out.full@net$cell$prob[[3]]@x) == 0L)

cat("\nAll filterCommunication Wave 2a checks passed.\n")
