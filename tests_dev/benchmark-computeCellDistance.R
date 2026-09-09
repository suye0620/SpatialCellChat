# ====== Benchmark: tree-backed spatial range query ======
# Run: Rscript tests_dev/benchmark-computeCellDistance.R
#
# The package implementation uses BiocNeighbors::VptreeParam for an exact
# radius query. The reference intentionally materializes a dense distance
# matrix to expose the memory/time cost of the pre-refactor approach.

source("renv/activate.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(BiocNeighbors)
})
source("R/utilities.R")
source("R/spatial.R")

dense_reference <- function(coordinates, interaction.range, contact.range,
                             ratio = NULL, tol = 0) {
  coordinates <- as.matrix(coordinates)
  threshold <- interaction.range + tol
  if (!is.null(ratio)) threshold <- threshold / ratio
  distances <- as.matrix(stats::dist(coordinates))
  keep <- distances > 0 & distances <= threshold
  distances[!keep] <- 0
  if (!is.null(ratio)) distances <- distances * ratio
  d.spatial <- Matrix::Matrix(distances, sparse = TRUE)
  adj.contact <- createCellCellContactMatrixFrom_dspatial(
    d.spatial, tol = tol, contact.threshold = contact.range
  )
  list(d.spatial = d.spatial, adj.contact = adj.contact)
}

set.seed(42)
coords_check <- matrix(runif(200L), ncol = 2L)
reference <- dense_reference(coords_check, interaction.range = 0.2,
                             contact.range = 0.1, ratio = NULL, tol = 0)
actual <- suppressMessages(computeCellDistance(
  coords_check, interaction.range = 0.2,
  contact.range = 0.1, ratio = NULL, tol = 0
))
if (!isTRUE(all.equal(unname(as.matrix(actual[["d.spatial"]])),
                     unname(as.matrix(reference[["d.spatial"]])))) ||
    !isTRUE(all.equal(unname(as.matrix(actual[["adj.contact"]])),
                     unname(as.matrix(reference[["adj.contact"]]))))) {
  stop("tree-backed and dense reference results differ", call. = FALSE)
}
cat("[PASS] tree-backed range query matches dense reference\n")

measure_seconds <- function(fun, repeats = 3L) {
  elapsed <- system.time({
    for (i in seq_len(repeats)) invisible(fun())
  })[["elapsed"]]
  elapsed / repeats
}

sizes <- c(250L, 500L, 1000L)
rows <- lapply(sizes, function(n) {
  coords <- matrix(runif(n * 2L, min = 0, max = 1000), ncol = 2L)
  dense_time <- measure_seconds(function() dense_reference(
    coords, interaction.range = 75, contact.range = 30,
    ratio = NULL, tol = 0
  ))
  tree_time <- measure_seconds(function() computeCellDistance(
    coords, interaction.range = 75, contact.range = 30,
    ratio = NULL, tol = 0
  ))
  data.frame(
    n = n,
    dense_seconds = dense_time,
    tree_seconds = tree_time,
    speedup = if (tree_time > 0) dense_time / tree_time else NA_real_
  )
})
benchmark <- do.call(rbind, rows)
print(benchmark, row.names = FALSE)
cat("Note: timings are machine-dependent; the key contract is exactness with bounded intermediate memory.\n")
