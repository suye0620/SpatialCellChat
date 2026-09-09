# Deterministic 100k-cell sparse grid smoke test.
# This records elapsed time and R-level allocation deltas; it does not claim
# process peak RSS, which is not portable on the Windows validation host.

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(sf)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

nx <- 316L
ny <- 317L
n <- nx * ny
cells <- paste0("cell_", seq_len(n))
coordinates <- cbind(
  x = rep(seq_len(nx), each = ny),
  y = rep.int(seq_len(ny), nx)
)
rownames(coordinates) <- cells
expr <- Matrix::Matrix(
  1,
  nrow = 1L,
  ncol = n,
  dimnames = list("gene_1", cells),
  sparse = TRUE
)
meta <- data.frame(
  label = factor(rep(c("A", "B"), length.out = n)),
  row.names = cells
)
chat <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = 0.5)
)

gc()
before <- gc()
timing <- system.time({
  grid <- makeGridSpatialCellChat(chat, cellsize = c(1, 1))
})
after <- gc()

n_grid <- nrow(grid@images$coordinates)
delta_mb <- (after[["Ncells", "used"]] - before[["Ncells", "used"]]) * 56 / 1024^2
check("100k grid object validates", isTRUE(validateSpatialCellChat(grid)))
check("100k grid retains sparse assay", inherits(grid@assay$norm, "dgCMatrix"))
check("100k grid has near-complete occupancy", n_grid >= 90000L)
check("100k grid retains one member per source point", all(grid@images$.grid$within.nGrid > 0L))
cat(sprintf(
  "100k sparse grid: input=%d output=%d elapsed=%.3f s Ncells_delta=%.1f MB\n",
  n, n_grid, unname(timing[["elapsed"]]), delta_mb
))
cat("All 100k sparse grid checks passed.\n")
