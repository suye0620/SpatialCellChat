# ====== Test: computeCellDistance ======
# Run: Rscript tests_dev/test-computeCellDistance.R

Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(Matrix)
  library(BiocNeighbors)
})
source("R/utilities.R")
source("R/spatial.R")

ok <- 0L
fail <- 0L
check <- function(desc, condition) {
  if (isTRUE(condition)) {
    cat("[PASS] ", desc, "\n", sep = "")
    ok <<- ok + 1L
  } else {
    cat("[FAIL] ", desc, "\n", sep = "")
    fail <<- fail + 1L
  }
}

coords <- matrix(
  c(0, 0,
    3, 0,
    0, 4,
    10, 0),
  ncol = 2L, byrow = TRUE,
  dimnames = list(paste0("cell", 1:4), c("x", "y"))
)
res <- computeCellDistance(
  coords, interaction.range = 5, contact.range = 3, ratio = NULL, tol = NULL
)
check("distance cache records normalized parameters", isTRUE(
  identical(res[[".parameters"]], list(
    interaction.range = 5,
    contact.range = 3,
    ratio = NULL,
    tol = 0
  ))
))
check("distance cache matcher rejects changed thresholds", !.spatial_distance_cache_matches(
  res, interaction.range = 4, contact.range = 3, ratio = NULL, tol = NULL
))
distance_expected <- matrix(
  c(0, 3, 4, 0,
    3, 0, 5, 0,
    4, 5, 0, 0,
    0, 0, 0, 0),
  nrow = 4L, byrow = TRUE
)
contact_expected <- matrix(
  c(1, 1, 0, 0,
    1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1),
  nrow = 4L, byrow = TRUE
)
check("ratio=NULL uses coordinate units directly", isTRUE(
  all.equal(unname(as.matrix(res[["d.spatial"]])), distance_expected)
))
check("contact adjacency includes diagonal and contact pairs", isTRUE(
  all.equal(unname(as.matrix(res[["adj.contact"]])), contact_expected)
))

scaled <- computeCellDistance(
  coords[1:2, , drop = FALSE],
  interaction.range = 5, contact.range = 5, ratio = 0.5, tol = 0
)
check("ratio scales returned distances", isTRUE(
  identical(as.numeric(scaled[["d.spatial"]][1, 2]), 1.5)
))

coords_3d <- cbind(coords, z = c(0, 0, 0, 1))
res_3d <- computeCellDistance(
  coords_3d, interaction.range = 5, contact.range = 3, ratio = NULL, tol = 0
)
check("three-dimensional coordinates are supported", identical(dim(res_3d[["d.spatial"]]), c(4L, 4L)))

check("zero ratio is rejected", inherits(
  try(computeCellDistance(coords, ratio = 0), silent = TRUE), "try-error"
))
check("negative tolerance is rejected", inherits(
  try(computeCellDistance(coords, tol = -1), silent = TRUE), "try-error"
))

cat(sprintf("All computeCellDistance checks passed: %d; failed: %d\n", ok, fail))
if (fail > 0L) quit(status = 1L)
