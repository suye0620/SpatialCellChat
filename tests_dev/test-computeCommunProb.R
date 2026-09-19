# Wave 1c: 新 computeCommunProb 主体 vs 冻结基线 oracle（fixture 逐位对拍）
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/reference-computeCommunProb-baseline.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

## ---- fixture: 20 细胞 5×4 网格（间距 10），9 基因，4 LR（分泌/接触/复合物/AGAN 全覆盖）----
set.seed(23)
cells <- paste0("c", 1:20)
genes <- c("G1", "G2", "G3", "R1", "R2", "SU1", "SU2", "COF1", "COF2")
expr <- matrix(round(runif(length(genes) * length(cells)), 3),
               nrow = length(genes), dimnames = list(genes, cells))
expr[expr < 0.4] <- 0
expr["SU2", c(3, 8, 13)] <- 0
coor <- expand.grid(x = seq(10, 50, 10), y = seq(10, 40, 10))
rownames(coor) <- cells
meta <- data.frame(label = factor(rep(c("A", "B", "C", "D"), each = 5)),
                   row.names = cells)

chat <- createSpatialCellChat(
  object = as(expr, "dgCMatrix"), meta = meta, group.by = "label",
  input.assay = "norm", datatype = "spatial",
  coordinates = coor, spatial.factors = list(ratio = 1, tol = 5)
)
chat@DB <- list(
  interaction = data.frame(
    ligand = c("G1", "CPLX", "G2", "G3"),
    receptor = c("R1", "R2", "R2", "R1"),
    agonist = c("", "", "", "AG1"),
    antagonist = c("", "", "", "AN1"),
    co_A_receptor = c("", "", "", "AG2"),
    co_I_receptor = c("", "", "", "AI1"),
    annotation = c("Secreted Signaling", "Secreted Signaling",
                   "Cell-Cell Contact", "Secreted Signaling"),
    interaction_name = c("G1-R1", "CPLX-R2", "G2-R2", "G3-R1"),
    pathway_name = c("P1", "P2", "P3", "P4"),
    stringsAsFactors = FALSE),
  complex = data.frame(row.names = "CPLX", subunit_1 = "SU1", subunit_2 = "SU2",
                       stringsAsFactors = FALSE),
  cofactor = data.frame(row.names = c("AG1", "AG2", "AI1", "AN1"),
                        cofactor1 = c("COF1", "COF1", "COF2", "COF2"),
                        cofactor2 = c("", "COF1", "", ""),
                        stringsAsFactors = FALSE)
)
LRsig <- chat@DB$interaction[, c("ligand", "receptor", "agonist", "antagonist",
                                 "co_A_receptor", "co_I_receptor", "annotation",
                                 "interaction_name", "pathway_name")]
# 真实流程中 LRsig 继承 subsetData/identifyOverExpressedInteractions 的 annotation 排序；
# fixture 必须复刻该约定，否则 nLR1 边界两边不一致
LRsig <- LRsig[order(factor(LRsig$annotation,
  levels = c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling",
             "Cell-Cell Contact"))), , drop = FALSE]
chat@LR$LRsig <- LRsig
assay(chat, "signaling") <- assay(chat, "norm")

## ---- 基线 oracle（plain 输入，冻结参考）----
oracle <- baseline_reference_commun_prob(
  data.signaling = assay(chat, "signaling"), LR.use = LRsig,
  coordinates = coor, spatial.factors = list(ratio = 1, tol = 5),
  complex_input = chat@DB$complex, cofactor_input = chat@DB$cofactor,
  group = chat@idents, Kh = 0.5, n = 1, distance.use = TRUE,
  interaction.range = 250, scale.distance = 1,
  use.AGAN = TRUE, contact.dependent = TRUE,
  contact.range = 10, contact.dependent.forced = FALSE)

## ---- 新主体（内存模式）----
chat.cc <- computeCommunProb(chat, scale.distance = 1, use.AGAN = TRUE,
                             contact.dependent = TRUE, contact.range = 10,
                             nthreads = 1L, verbose = FALSE)
check("new object validates", isTRUE(validObject(chat.cc)))
prob.new <- chat.cc@net$cell$prob
check("net$cell$prob is SparseChatArray with 4 layers",
      inherits(prob.new, "SparseChatArray") && length(prob.new) == nrow(LRsig))

## ---- 逐层逐位对拍 ----
for (k in seq_len(nrow(LRsig))) {
  nm <- rownames(LRsig)[k]
  oracle.layer <- as.matrix(oracle$net$prob.cell[, , k])
  new.layer <- as.matrix(unclass(prob.new)[[k]])
  check(sprintf("layer %s bitwise (values)", nm),
        identical(as.numeric(oracle.layer), as.numeric(new.layer)))
  check(sprintf("layer %s pattern", nm),
        identical(which(oracle.layer != 0), which(new.layer != 0)))
}
check("LR layer names match", identical(names(prob.new), rownames(LRsig)))

## ---- 多线程确定性 ----
chat.n4 <- computeCommunProb(chat, scale.distance = 1, use.AGAN = TRUE,
                             contact.dependent = TRUE, contact.range = 10,
                             nthreads = 4L, verbose = FALSE)
check("nthreads 1 vs 4 layers bitwise",
      identical(unclass(chat.cc@net$cell$prob), unclass(chat.n4@net$cell$prob)))

## ---- spill 模式 ----
sp <- file.path(tempdir(), "ccp_spill_test")
if (dir.exists(sp)) unlink(sp, recursive = TRUE)
chat.sp <- computeCommunProb(chat, scale.distance = 1, use.AGAN = TRUE,
                             contact.dependent = TRUE, contact.range = 10,
                             spill.dir = sp, verbose = FALSE)
check("spill mode leaves net$cell$prob NULL for stream-in",
      is.null(chat.sp@net$cell$prob))
check("spill mode object validates", isTRUE(validObject(chat.sp)))
sp.files <- list.files(sp, pattern = "\\.rds$")
check("spill files written (one per LR)", length(sp.files) == nrow(LRsig))
for (k in seq_len(nrow(LRsig))) {
  sp.layer <- readRDS(file.path(sp, paste0(names(prob.new)[k], ".rds")))
  check(sprintf("spill layer %d bitwise vs in-memory", k),
        identical(as.numeric(sp.layer@x), as.numeric(unclass(prob.new)[[k]]@x)))
}
unlink(sp, recursive = TRUE)

## ---- params/log 记录 ----
check("misc$.param$communication recorded",
      !is.null(chat.cc@misc$.param$communication) &&
        chat.cc@misc$.param$communication$nLR == nrow(LRsig) &&
        chat.cc@misc$.param$communication$Kh == 0.5)
check("operation logged",
      any(vapply(chat.cc@misc$.log,
                 function(e) identical(e[["function"]], "computeCommunProb"), logical(1L))))

cat("\nAll Wave 1c main-body checks passed.\n")

## =====================================================================
## Wave 1d: filterProbability 迁移（读 net$cell$prob / misc$.param$communication
## + spill 消费回填 + 高效分位捷径逐位等价 + ⑤ signif 门控）
## =====================================================================
nboot <- 15L
quantile.prob <- 0.95

## ---- (A) 拒绝路径：无 computeCommunProb 输出 ----
chat.bad <- chat.cc
chat.bad@net$cell$prob <- NULL
chat.bad@misc$.param$communication <- NULL
err.bad <- tryCatch({ filterProbability(chat.bad, nboot = 20); NULL },
                    error = function(e) conditionMessage(e))
check("rejection: filterProbability without layers errors",
      !is.null(err.bad) && grepl("computeCommunProb", err.bad))

## ---- (B) 内存模式：逐层逐位 vs 内联基线复刻 ----
prob.in <- unclass(chat.cc@net$cell$prob)
set.seed(666L)
d.spatial <- chat.cc@images$.distance$d.spatial
Matrix::diag(d.spatial) <- 1
adj.contact <- chat.cc@images$.distance$adj.contact
nC <- NROW(d.spatial)
nLR1 <- chat.cc@misc$.param$communication$nLR1
permutation <- replicate(nrow(LRsig), base::sample(x = 1:nC, size = nboot, replace = F))
ref.layers <- lapply(seq_along(prob.in), function(i) {
  d_spatial <- if (i <= nLR1) d.spatial else adj.contact
  Prob.cell.i <- prob.in[[i]]
  sample.cells <- permutation[, i, drop = TRUE]
  sample.prob.cell.i <- purrr::map(.x = sample.cells,
                                   .f = function(cell.index) {
                                     prob.index <- which(d_spatial[cell.index, , drop = TRUE] > 0)
                                     Prob.cell.i[cell.index, prob.index, drop = TRUE]
                                   }) %>% unlist()
  nboot.quantile <- quantile(sample.prob.cell.i, probs = quantile.prob)
  if (nboot.quantile == 0) return(Prob.cell.i)
  if (nboot.quantile > 0) {
    Prob.cell.i <- scMatrixTruncation(Prob.cell.i, cutoff = nboot.quantile,
                                      remain.cutoff.v = TRUE, repr = "C")
  }
  Prob.cell.i
})
chat.filt <- filterProbability(chat.cc, nboot = nboot, seed.use = 666L, thresh = 0.05)
filt.layers <- unclass(chat.filt@net$cell$prob)
check("filter output is SparseChatArray", inherits(chat.filt@net$cell$prob, "SparseChatArray"))
check("filter layer names preserved",
      identical(dimnames(chat.filt@net$cell$prob)[[3]], rownames(LRsig)))
for (k in seq_along(ref.layers)) {
  nm <- rownames(LRsig)[k]
  ref.m <- as.matrix(ref.layers[[k]])
  new.m <- as.matrix(filt.layers[[k]])
  check(sprintf("filter layer %s bitwise vs baseline", nm),
        identical(as.numeric(ref.m), as.numeric(new.m)) &&
          identical(which(ref.m != 0), which(new.m != 0)))
}

## ---- (C) spill 模式：读盘消费 -> 回填 net$cell$prob -> 清理 ----
sp2 <- file.path(tempdir(), "ccp_spill_fp")
if (dir.exists(sp2)) unlink(sp2, recursive = TRUE)
chat.sp2 <- computeCommunProb(chat, scale.distance = 1, use.AGAN = TRUE,
                              contact.dependent = TRUE, contact.range = 10,
                              spill.dir = sp2, verbose = FALSE)
check("spill2 files present pre-filter",
      length(list.files(sp2, pattern = "\\.rds$")) == nrow(LRsig))
chat.sp.filt <- filterProbability(chat.sp2, nboot = nboot, seed.use = 666L)
check("spill consumed: net$cell$prob materialized",
      inherits(chat.sp.filt@net$cell$prob, "SparseChatArray"))
check("spill consumed: files removed", !file.exists(sp2))
sp.param <- chat.sp.filt@misc$.param$communication
check("spill params cleared",
      is.null(sp.param$spill.dir) && is.null(sp.param$spill.backend) &&
        is.null(sp.param$spill.count))
check("spill filtered == in-memory filtered (bitwise)",
      identical(unclass(chat.sp.filt@net$cell$prob), filt.layers))
check("filterProbability logged",
      any(vapply(chat.sp.filt@misc$.log,
                 function(e) identical(e[["function"]], "filterProbability"), logical(1L))))

## ---- (D) ⑤ 门控：signif 数值语义 vs format 字符串语义 ----
dtemp <- data.frame(group = c("A", "A"), r1 = 1e-4, r2 = 0.5, r3 = 0.05)
gate.old <- 1 * (format(dtemp[, -1], digits = 1) >= 0.1)
gate.new <- 1 * (signif(dtemp[, -1], 1) >= 0.1)
check("⑤ format gate lets ghost tiny rate pass", gate.old[1, "r1"] == 1)
check("⑤ signif gate zeros ghost tiny rate", gate.new[1, "r1"] == 0)
check("⑤ normal rates unchanged by signif fix",
      gate.new[1, "r2"] == 1 && gate.new[1, "r3"] == 0 && gate.old[1, "r2"] == 1)

## ---- (E) 高效分位捷径：近零层 cutoff 必为 0 -> 与基线“quantile==0 返回原层”等价 ----
cells20 <- colnames(assay(chat.cc, "signaling"))
mk <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 4), x = c(0.5, 0.7, 0.9),
                           dims = c(20, 20))
tiny.list <- lapply(seq_along(filt.layers), function(k) mk)
names(tiny.list) <- rownames(LRsig)
ta <- SparseChatArray(tiny.list)
dimnames(ta) <- list(cells20, cells20, rownames(LRsig))
chat.tiny <- chat.cc
chat.tiny@net$cell$prob <- ta
chat.tiny.filt <- filterProbability(chat.tiny, nboot = nboot, seed.use = 666L)
tiny.out <- unclass(chat.tiny.filt@net$cell$prob)
for (k in seq_along(tiny.out)) {
  check(sprintf("zero-shortcut layer %d returns layer unchanged", k),
        identical(as.numeric(tiny.out[[k]]@x), as.numeric(mk@x)) &&
          identical(methods::as(tiny.out[[k]], "dgTMatrix")@i, methods::as(mk, "dgTMatrix")@i))
}

cat("\nAll Wave 1d filterProbability checks passed.\n")
