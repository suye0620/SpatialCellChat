# Wave 1b: cpp_prob_layer kernel —— R oracle / 基线链 逐位对拍
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

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
same <- function(a, b) identical(as.numeric(a), as.numeric(b))

## R 向量化语义规范（oracle）：与 kernel 相同的乘法分组
oracle_v <- function(x0, row0, col0, L, R, Kh, n, cmask, fAG = NULL, fAN = NULL) {
  lr <- L[row0] * R[col0]
  h <- if (n == 1) lr / (Kh + lr) else lr^n / (Kh^n + lr^n)
  v <- x0 * h * cmask
  if (!is.null(fAG)) v <- v * (fAG[row0] * fAG[col0])
  if (!is.null(fAN)) v <- v * (fAN[row0] * fAN[col0])
  v
}

## ---- fixture: 8 细胞、稀疏对称模式 + 自通讯对角 ----
set.seed(11)
nC <- 8L
d <- matrix(0, nC, nC)
d[sample(nC * nC, 30)] <- runif(30, 0.5, 200)
d <- as(forceSymmetric(as(d, "dgCMatrix")), "dgCMatrix")
diag(d) <- runif(nC, 0.05, 0.3)          # 自通讯对角
x0 <- d@x
row0 <- d@i                               # 0-based
col0 <- rep.int(0:(nC - 1L), diff(d@p))
L <- runif(nC, 0, 1); L[sample(nC, 3)] <- 0
R <- runif(nC, 0, 1); R[sample(nC, 2)] <- 0
fAG <- pmin(runif(nC, 0, 1) + 0.2, 1)
fAN <- runif(nC, 0, 0.8)
mask1 <- rep(1, length(x0))
maskc <- numeric(length(x0)); maskc[sample(length(x0), 12)] <- 1

layer_from_kernel <- function(n, cmask, fAG = NULL, fAN = NULL, nthreads = 1L) {
  out <- cpp_prob_layer(x0, row0, col0, L, R, Kh = 0.5, n = n,
                        contact_mask = cmask,
                        fAG = if (is.null(fAG)) NULL else fAG,
                        fAN = if (is.null(fAN)) NULL else fAN,
                        nC = nC, nthreads = nthreads)
  new("dgCMatrix", i = out$i, p = out$p, x = out$x, Dim = c(nC, nC))
}
## kernel 层在 P.spatial 原位置上的取值（被丢弃的位置 = 0）
values_at_pattern <- function(lay) as.numeric(lay[cbind(row0 + 1L, col0 + 1L)])

## ---- 1. kernel vs R oracle：分支矩阵（n × contact × 因子组合）----
combos <- 0L
for (n in c(1, 2)) for (cm in list(mask1, maskc)) for (fa in list(NULL, fAG)) for (fn in list(NULL, fAN)) {
  lay <- layer_from_kernel(n, cm, fa, fn)
  v.o <- oracle_v(x0, row0 + 1L, col0 + 1L, L, R, 0.5, n, cm, fa, fn)
  check(sprintf("kernel vs oracle (n=%d, contact=%s, fAG=%s, fAN=%s)",
                n, identical(cm, maskc), !is.null(fa), !is.null(fn)),
        same(values_at_pattern(lay), v.o))
  combos <- combos + 1L
}
cat(sprintf("  (%d 组合)\n", combos))

## ---- 2. kernel vs 基线子函数链逐位对拍（n=1 + AGAN + contact 场景）----
# 基线链：crossprod(1×nC 稀疏行) → HillFunctionFordataLR → *P.spatial
#         → [contact] * adj → myElementwiseProduct(AG) → myElementwiseProduct(AN)
## contact 场景基线链：由 maskc scatter 出 adj 矩阵（值 1），插入 P1_Pspatial * adj 步骤
keep.c <- maskc == 1
adj.sparse <- Matrix::sparseMatrix(i = row0[keep.c] + 1L, j = col0[keep.c] + 1L,
                                   x = rep(1, sum(keep.c)), dims = c(nC, nC))
i.lr <- 3L   # 选一个旁分泌 LR
nLR1 <- 2L   # 让 i.lr > nLR1 → 接触依赖路径
AGv <- fAG; ANv <- fAN
# 基线 L/R 行（与 L/R 向量同源）
x_ <- as(matrix(L, nrow = 1), "TsparseMatrix")
y_ <- as(matrix(R, nrow = 1), "TsparseMatrix")
dataLR <- Matrix::crossprod(x_, y_)
P1 <- HillFunctionFordataLR(dataLR, Kh = 0.5, n = 1)
P1_Pspatial <- P1 * d
bl_layer <- P1_Pspatial * adj.sparse
bl_layer <- myElementwiseProduct(bl_layer, AGv)
bl_layer <- myElementwiseProduct(bl_layer, ANv)
# kernel 同场景：mask = adj 模式（本 fixture 直接用全 1 简化——AGAN 不受 mask 影响；
# 为覆盖 contact×AGAN 组合，用 maskc）
k_layer <- layer_from_kernel(1, maskc, fAG, fAN)
check("kernel layer vs baseline chain (slots i/p/x)",
      identical(k_layer@i, bl_layer@i) && identical(k_layer@p, bl_layer@p) &&
        identical(k_layer@x, bl_layer@x))
# 旁分泌场景（mask=1）：contact 不截断
k2 <- layer_from_kernel(1, mask1, fAG, fAN)
bl2 <- myElementwiseProduct(P1_Pspatial, AGv)
bl2 <- myElementwiseProduct(bl2, ANv)
check("kernel vs baseline chain (旁分泌, slots)", identical(k2@i, bl2@i) &&
      identical(k2@p, bl2@p) && identical(k2@x, bl2@x))

## ---- 3. 多线程确定性（nthreads 1 vs 4）----
k1 <- layer_from_kernel(1, maskc, fAG, fAN, nthreads = 1L)
k4 <- layer_from_kernel(1, maskc, fAG, fAN, nthreads = 4L)
check("nthreads 1 vs 4 bitwise", identical(k1@x, k4@x) && identical(k1@i, k4@i) &&
      identical(k1@p, k4@p))

## ---- 4. 全零层（L/R 全零）----
z.L <- numeric(nC)
out <- cpp_prob_layer(x0, row0, col0, z.L, R, Kh = 0.5, n = 1, contact_mask = mask1,
                      fAG = NULL, fAN = NULL, nC = nC)
check("all-zero L yields empty layer", length(out$x) == 0 && all(out$p == 0))

cat("\nAll Wave 1b kernel checks passed.\n")
