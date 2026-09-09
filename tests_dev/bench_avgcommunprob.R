## Benchmark: computeAvgCommunProb_LR_Avg current vs sparse-matrix version
suppressMessages({
  library(Matrix); library(stats)
})

set.seed(1)
nC <- 5000; k <- 10; nboot <- 100
prob <- Matrix::rsparsematrix(nC, nC, density = 200000 / nC^2)  # ~200k nnz
group <- factor(sample(paste0("G", 1:k), nC, replace = TRUE))
dataLR <- cbind(runif(nC), runif(nC))   # nC x 2 dense

## ---- current implementation (as in computeAvgCommunProb_LR_Avg) ----
cur <- function(prob, group, dataLR, min.percent = 0.1, min.cells.sr = 5) {
  cell.type.mat <- model.matrix(~ group - 1)
  dataLR_temp <- 1 * (dataLR > 0)
  dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
  dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1),
                                    matrix(dataLR_percent[, 2], nrow = 1))
  if (sum(Prob_percent) == 0) {
    Prob.avg <- Prob_percent
  } else {
    Prob.avg <- Matrix::crossprod(x = cell.type.mat, y = prob %*% cell.type.mat)
    prob@x <- rep.int(1, times = length(prob@x))
    Prob.scale.factor <- Matrix::crossprod(x = cell.type.mat, y = prob %*% cell.type.mat)
    Prob.avg <- Prob.avg / Prob.scale.factor
    Prob.avg[is.nan(Prob.avg)] <- 0
    sender.counts <- Matrix::rowSums(prob)
    receptor.counts <- Matrix::colSums(prob)
    cells.sr <- cbind(sender.counts, receptor.counts)
    cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
    cells.sr <- 1 * (cells.sr[, -1] >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1),
                                  matrix(cells.sr[, 2], nrow = 1))
    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }
  as.matrix(Prob.avg)
}

## ---- proposed: precomputed one-hot G + sparse products; permuted by row-index ----
G <- Matrix::sparse.model.matrix(~ group - 1)         # nC x k
cnt <- as.numeric(Matrix::colSums(G))
Lp <- as.numeric(dataLR[, 1] > 0); Rp <- as.numeric(dataLR[, 2] > 0)

prop <- function(prob, G, perm = NULL, min.percent = 0.1, min.cells.sr = 5) {
  Gp <- if (is.null(perm)) G else G[perm, , drop = FALSE]
  pctL <- as.numeric(Matrix::crossprod(Gp, Lp)) / cnt
  pctR <- as.numeric(Matrix::crossprod(Gp, Rp)) / cnt
  Prob_percent <- outer(pctL >= min.percent, pctR >= min.percent) * 1
  if (sum(Prob_percent) == 0) return(Prob_percent)
  Prob.avg <- Matrix::crossprod(Gp, prob %*% Gp)
  prob1 <- prob; prob1@x <- rep.int(1, length(prob1@x))
  Prob.scale <- Matrix::crossprod(Gp, prob1 %*% Gp)
  Prob.avg <- Prob.avg / Prob.scale
  Prob.avg[is.nan(Prob.avg)] <- 0
  sr <- as.numeric(Matrix::crossprod(Gp, Matrix::rowSums(prob)))
  cr <- as.numeric(Matrix::crossprod(Gp, Matrix::colSums(prob)))
  cells.sr <- outer(sr >= min.cells.sr, cr >= min.cells.sr) * 1
  as.matrix(Prob.avg * Prob_percent * cells.sr)
}

## numerical identity (no permutation)
a <- cur(prob, group, dataLR)
b <- prop(prob, G)
cat("identity (no perm):", isTRUE(all.equal(a, b, tolerance = 0)), "\n")

## single-call cost
t0 <- system.time(for (i in 1:20) cur(prob, group, dataLR))[["elapsed"]] / 20
t1 <- system.time(for (i in 1:20) prop(prob, G))[["elapsed"]] / 20
## permuted cost (permutation test inner loop)
perm <- replicate(nboot, sample.int(nC, size = nC))
t2 <- system.time(for (e in 1:20) for (nE in 1:nboot) prop(prob, G, perm = perm[, nE]))[["elapsed"]] / 20 / nboot
cat(sprintf("current  : %8.2f ms/call\n", t0 * 1000))
cat(sprintf("sparse   : %8.2f ms/call\n", t1 * 1000))
cat(sprintf("sparse+perm: %8.2f ms/call\n", t2 * 1000))
cat(sprintf("permutation total: current ~%.0f s vs sparse ~%.0f s (nLR=500, nboot=100)\n",
            500 * nboot * t0, 500 * nboot * t2))
