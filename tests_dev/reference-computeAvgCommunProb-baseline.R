# 冻结基线参考实现：pre-migration computeAvgCommunProb（modeling.R @ be88f30）的
# plain-input 快照，供 Plan A 重构对拍（与 test-computeCommunProb 的冻结基线模式一致）。
# 除以下两点外逐行保真：
#   1) object 读取替换为显式入参（prob.cell_ / dataLavg / dataRavg / group）；
#   2) my_future_sapply/lapply 替换为串行 sapply/lapply（并行只影响调度，不影响数值；
#      串行 = 单 worker 语义）。
# dataLavg/dataRavg 语义 = 旧 net$tmp$Lavg/Ravg 的第 i 行（门控内经 `1 * (dataLR > 0)`）。
# 传 0/1 矩阵时 `1 * (x > 0) == x`，与 v2 的 suppL/suppR 门控输入逐位一致。
computeAvgCommunProb_baseline_ref <- function(
    prob.cell_,
    dataLavg,
    dataRavg,
    group,
    avg.type = c("avg", "sum"),
    min.percent = 0.1,
    min.cells.sr = 5,
    do.permutation = TRUE,
    nboot = 100,
    seed.use = 1L,
    permutation = NULL
) {
  avg.type <- match.arg(avg.type)
  if (avg.type == "avg") {
    computeAvgCommunProb_LR <- computeAvgCommunProb_LR_Avg
  } else {
    computeAvgCommunProb_LR <- computeAvgCommunProb_LR_Sum
  }
  nC <- ncol(dataLavg)
  numCluster <- nlevels(group)
  pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster)

  Prob <- array(0, dim = c(numCluster, numCluster, length(prob.cell_)))
  Pval <- array(1, dim = c(numCluster, numCluster, length(prob.cell_)))
  dimnames(Prob) <- list(levels(group), levels(group), names(prob.cell_))
  dimnames(Pval) <- dimnames(Prob)

  set.seed(seed.use)

  # 串行版 map_dbl（数值同义：length(Mat@x)）
  prob.sum <- vapply(prob.cell_, function(Mat) length(Mat@x), numeric(1))
  LRsig.use.idx <- which(prob.sum > 0)
  if (length(LRsig.use.idx) < 1) {
    stop("Each LR pair does not have any cell-level links/interactions.")
  }

  Prob.avg_ <- sapply(seq_along(LRsig.use.idx), function(x) {
    i <- LRsig.use.idx[[x]]
    prob.cell.i <- prob.cell_[[i]]
    dataLR_temp <- cbind(dataLavg[i, ], dataRavg[i, ])
    Prob.avg <- computeAvgCommunProb_LR(
      prob.cell.i,
      group = group,
      dataLR = dataLR_temp,
      min.percent = min.percent,
      min.cells.sr = min.cells.sr
    )
    # 基线此处有 `Prob.avg[pval.colo > thresh.colo] <- 0`；pval.colo 恒为全零
    # 矩阵（colocalization.use = FALSE 路径），掩码恒为 FALSE，等价 no-op，略去。
    Prob.avg
  }, simplify = FALSE)
  for (x in seq_along(LRsig.use.idx)) {
    i <- LRsig.use.idx[[x]]
    Prob[, , i] <- Prob.avg_[[x]]
  }

  prob.sum <- apply(Prob > 0, 3, sum)
  LRsig.use.idx <- which(prob.sum > 0)
  if (do.permutation) {
    if (is.null(permutation)) permutation <- replicate(nboot, sample.int(nC, size = nC))
    Pval_ <- lapply(seq_along(LRsig.use.idx), function(x) {
      i <- LRsig.use.idx[[x]]
      prob.cell.i <- prob.cell_[[i]]
      dataLR_temp <- cbind(dataLavg[i, ], dataRavg[i, ])
      Pnull <- as.vector(Prob[, , i])
      Pboot <- sapply(
        X = 1:nboot,
        FUN = function(nE) {
          groupboot <- group[permutation[, nE]]
          Pboot.avg <- computeAvgCommunProb_LR(
            prob.cell.i,
            group = groupboot,
            dataLR = dataLR_temp,
            min.percent = min.percent,
            min.cells.sr = min.cells.sr
          )
          as.vector(Pboot.avg)
        }
      )
      Pboot <- matrix(unlist(Pboot), nrow = length(Pnull), ncol = nboot, byrow = FALSE)
      nReject <- rowSums(Pboot - Pnull > 0)
      p <- nReject / nboot
      matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    })
    for (x in seq_along(LRsig.use.idx)) {
      i <- LRsig.use.idx[[x]]
      Pval[, , i] <- Pval_[[x]]
    }
    Pval[Prob == 0] <- 1
  } else {
    Pval <- NULL
  }
  list(Prob = Prob, Pval = Pval)
}
