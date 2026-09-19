# Wave 1a: 子函数与基线逐位对拍（computeExpr_* 系列）
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/modeling.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
same <- function(a, b) identical(as.numeric(a), as.numeric(b))  # names 无关的逐位比较

## 注：基线 computeExpr_complex/_coreceptor 内部调用 my_future_sapply(hint.message=...),
## 与当前包装器签名漂移（先在缺陷，复合物/共受体路径在基线中本就不可运行，已记录）。
## 对拍 oracle 采用其内联语义原语：geometricMean(子矩阵) / ∏(1+e) / 1+Hill(e)。
## computeExpr_LR 与 _agonist/_antagonist 不依赖漂移包装器，可直接对拍。

## ---- fixture：12 细胞 × 9 基因（复合物亚基、cofactor 基因、零值模式）----
set.seed(7)
cells <- paste0("c", 1:12)
genes <- c("G1", "G2", "G3", "R1", "R2", "SU1", "SU2", "COF1", "COF2")
expr <- matrix(runif(length(genes) * length(cells)), nrow = length(genes),
               dimnames = list(genes, cells))
expr[expr < 0.45] <- 0
expr["SU2", 3:4] <- 0           # CPLX 在部分细胞表达为零
sig.dense <- expr
sig.sparse <- as(expr, "dgCMatrix")

complex_input <- data.frame(row.names = "CPLX",
                            subunit_1 = "SU1", subunit_2 = "SU2",
                            stringsAsFactors = FALSE)
complex_subunits <- .sc_complex_subunits(complex_input)
su.cplx <- unlist(complex_subunits[["CPLX"]])

pairLR <- data.frame(
  ligand = c("G1", "CPLX", "G2", "G3"),
  receptor = c("R1", "R2", "R2", "R1"),
  co_A_receptor = c("AG1", "AG2", NA, ""),
  co_I_receptor = c("AI1", NA, "AI1", ""),
  agonist = c("AG1", "AG1", NA, ""),
  antagonist = c(NA, "AN1", "AN1", ""),
  stringsAsFactors = FALSE
)
cofactor_input <- data.frame(row.names = c("AG1", "AG2", "AI1", "AN1"),
                             cofactor1 = c("COF1", "COF1", "COF2", "COF2"),
                             cofactor2 = c("", "COF1", "", ""),
                             stringsAsFactors = FALSE)

cache <- new.env(parent = emptyenv())
cfc <- new.env(parent = emptyenv())

## ---- 1. 表达行 vs 基线（单基因直接对拍；复合物用内联 oracle）----
bl_LR <- computeExpr_LR(c("G1", "G2", "G3"), sig.dense, complex_input)
for (j in seq_len(3L)) {
  g <- c("G1", "G2", "G3")[j]
  check(sprintf("expr row %s bitwise (vs baseline computeExpr_LR)", g),
        same(bl_LR[j, ], .sc_expr_row(g, sig.sparse, complex_subunits, cache)))
}
oracle.cplx <- geometricMean(as.matrix(sig.dense[su.cplx, , drop = FALSE]))
mine.cplx <- .sc_expr_row("CPLX", sig.sparse, complex_subunits, cache)
check("expr row CPLX bitwise (oracle = geometricMean over subunit submatrix)",
      same(oracle.cplx, mine.cplx))
check("complex zero-subunit semantics (c3/c4 -> 0)",
      as.numeric(mine.cplx)[3] == 0 && as.numeric(mine.cplx)[4] == 0)

## ---- 2. 基线 computeExpr_agonist/_antagonist 直接对拍（n=1 与 n=2）----
bl_ag1 <- computeExpr_agonist(sig.dense, pairLR, cofactor_input,
                              index.agonist = 1, Kh = 0.5, n = 1)
mine_ag1 <- .sc_cofactor_factor("AG1", cofactor_input, sig.sparse,
                                type = "agonist", Kh = 0.5, n = 1, cache = cfc)
check("agonist AG1 n=1 bitwise", same(bl_ag1, mine_ag1))
bl_ag2 <- computeExpr_agonist(sig.dense, pairLR, cofactor_input,
                              index.agonist = 2, Kh = 0.5, n = 1)
mine_ag2 <- .sc_cofactor_factor("AG1", cofactor_input, sig.sparse,
                                type = "agonist", Kh = 0.5, n = 1, cache = cfc)
check("agonist cache-hit (LR2 same AG1) bitwise", same(bl_ag2, mine_ag2))
bl_an2 <- computeExpr_antagonist(sig.dense, pairLR, cofactor_input,
                                 index.antagonist = 2, Kh = 0.5, n = 1)
mine_an2 <- .sc_cofactor_factor("AN1", cofactor_input, sig.sparse,
                                type = "antagonist", Kh = 0.5, n = 1, cache = cfc)
check("antagonist AN1 n=1 bitwise", same(bl_an2, mine_an2))
bl_agn2 <- computeExpr_agonist(sig.dense, pairLR, cofactor_input,
                               index.agonist = 1, Kh = 0.5, n = 2)
mine_agn2 <- .sc_cofactor_factor("AG1", cofactor_input, sig.sparse,
                                 type = "agonist", Kh = 0.5, n = 2, cache = cfc)
check("agonist AG1 n=2 (通用路径) bitwise", same(bl_agn2, mine_agn2))

## ---- 3. 共受体/激动剂/拮抗剂内联 oracle（多基因分支 + NA/空三种形态）----
oracle_coA <- function(name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(rep(1, 12L))
  g <- as.character(unlist(cofactor_input[name, grep("cofactor", colnames(cofactor_input))],
                           use.names = FALSE))
  g <- g[g != ""]
  if (!length(g)) return(rep(1, 12L))
  apply(1 + as.matrix(sig.dense[g, , drop = FALSE]), 2, prod)
}
check("co-A AG1 (single member) bitwise",
      same(oracle_coA("AG1"), .sc_cofactor_factor("AG1", cofactor_input, sig.sparse,
                                                  type = "A", cache = cfc)))
check("co-A AG2 (multi member) bitwise",
      same(oracle_coA("AG2"), .sc_cofactor_factor("AG2", cofactor_input, sig.sparse,
                                                  type = "A", cache = cfc)))
check("co-A NA -> neutral ones",
      same(oracle_coA(NA), .sc_cofactor_factor(NA, cofactor_input, sig.sparse,
                                               type = "A", cache = cfc) %||% rep(1, 12L)))
check("co-I AI1 与 co-A 分支语义一致",
      same(oracle_coA("AI1"), .sc_cofactor_factor("AI1", cofactor_input, sig.sparse,
                                                  type = "I", cache = cfc)))

## ---- 4. 缓存确定性 ----
check("expr row cache deterministic",
      identical(.sc_expr_row("G1", sig.sparse, complex_subunits, cache),
                .sc_expr_row("G1", sig.sparse, complex_subunits, cache)))

## ---- 5. 拒绝路径 ----
check("absent gene/complex rejected",
      inherits(tryCatch(.sc_expr_row("NOPE", sig.sparse, complex_subunits, cache),
                        error = identity), "error"))
check("complex with missing subunit rejected",
      {
        ci2 <- data.frame(row.names = "BAD", subunit_1 = "SU1", subunit_2 = "NOPE",
                          stringsAsFactors = FALSE)
        cs2 <- .sc_complex_subunits(ci2)
        inherits(tryCatch(.sc_expr_row("BAD", sig.sparse, cs2, cache),
                          error = identity), "error")
      })

cat("\nAll Wave 1a helper parity checks passed.\n")
