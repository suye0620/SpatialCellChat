## Benchmark v2: future.multisession globals-transfer vs compute;
## plus O(nnz) outer-product trick replacing the dense nC x nC dataLR.
suppressMessages({
  library(Matrix); library(future); library(future.apply)
})
options(future.globals.maxSize = 2 * 1024^3)

set.seed(1)
nC <- 5000
nLR <- 50
P.spatial <- Matrix::rsparsematrix(nC, nC, density = 200000 / nC^2)  # ~200k nnz
dataL <- matrix(runif(nLR * nC), nLR, nC)
dataR <- matrix(runif(nLR * nC), nLR, nC)

## current implementation: dense outer product + Hill + sparse multiply
FUN_current <- function(i) {
  x_ <- Matrix::Matrix(dataL[i, ], nrow = 1, sparse = TRUE)
  y_ <- Matrix::Matrix(dataR[i, ], nrow = 1, sparse = TRUE)
  dataLR <- Matrix::crossprod(x_, y_)                      # nC x nC dense outer product
  dataLR@x <- dataLR@x^2 / (0.5^2 + dataLR@x^2)            # Hill
  dataLR * P.spatial
}

## proposed: evaluate Hill only at P.spatial nonzeros (O(nnz)), numerically identical
FUN_sparse <- function(i) {
  ps <- P.spatial
  L <- dataL[i, ]; R <- dataR[i, ]
  ## CSC traversal: column j holds rows i[p[j]..p[j+1]-1]; lr = L[row] * R[col]
  col <- rep.int(seq_len(ncol(ps)), diff(ps@p))
  lr <- L[ps@i + 1L] * R[col]
  ps@x <- ps@x * (lr^2 / (0.5^2 + lr^2))
  ps
}

tm <- function(label, expr) {
  t <- system.time(expr)[["elapsed"]]
  cat(sprintf("%-48s %8.2f s\n", label, t))
  invisible(t)
}

## correctness: current vs sparse must match exactly
future::plan("sequential")
r1 <- lapply(1:5, FUN_current)
r2 <- lapply(1:5, FUN_sparse)
cat("numerical identity (current vs sparse):",
    all(mapply(function(a, b) identical(a@x, b@x), r1, r2)), "\n")

## sequential compute cost
cat("\n-- compute-only (sequential) --\n")
t_seq <- tm("current (dense outer product)", lapply(1:nLR, FUN_current))
t_sp  <- tm("sparse O(nnz) trick           ", lapply(1:nLR, FUN_sparse))

## multisession: globals-transfer overhead
cat("\n-- multisession 4 workers --\n")
future::plan("multisession", workers = 4)
t_m1 <- tm("current, default sched (globals x nLR)",
           future_lapply(1:nLR, FUN_current, future.seed = TRUE))
t_m2 <- tm("current, chunked (globals x ~workers)",
           future_lapply(1:nLR, FUN_current, future.seed = TRUE,
                         future.scheduling = ceiling(nLR / 4)))
t_m3 <- tm("sparse O(nnz), default sched",
           future_lapply(1:nLR, FUN_sparse, future.seed = TRUE))
t_m4 <- tm("sparse O(nnz), chunked",
           future_lapply(1:nLR, FUN_sparse, future.seed = TRUE,
                         future.scheduling = ceiling(nLR / 4)))

future::plan("sequential")
cat(sprintf("\nspeedups: sparse/seq %.1fx; chunked/default %.1fx; combined %.1fx\n",
            t_seq / t_sp, t_m1 / t_m2, t_seq / t_m4))
