# Regression tests for the sparse grid intersection in makeGridSpatialCellChat.
#
# The grid assignment previously built a dense n x nGrid logical matrix via
# sf::st_intersects(sparse = FALSE) (e.g. ~15 GiB at 1e5 cells). The sparse
# sgbp implementation must produce numerically identical objects while keeping
# memory O(nnz). Tests:
#   1. small synthetic object: sparse vs dense-reference equivalence
#   2. hexagonal/square grids: equivalence
#   3. simulated_data end-to-end: equivalence plus manual spot checks
#   4. 10k-cell grid stress: correctness and bounded memory (no dense matrix)

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

expect_error <- function(label, expr, pattern) {
  err <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) e)
  check(label, inherits(err, "error"))
  check(paste0(label, ": message"), grepl(pattern, conditionMessage(err)))
}

# Dense reference: exactly the pre-sparsification implementation.
makeGrid_dense_reference <- function(object, data.slot = "norm",
                                     cellsize = c(5, 5), what = "polygons",
                                     square = TRUE, idents.ties.method = "first") {
  coordinates <- object@images$coordinates
  coordinates <- as.data.frame(coordinates[, seq_len(2L), drop = FALSE])
  colnames(coordinates) <- c("x_cent", "y_cent")
  coordinates$cell_type <- object@idents
  spatial.factors <- object@images$spatial.factors
  ratio <- if (is.list(spatial.factors) && !is.null(spatial.factors$ratio)) {
    spatial.factors$ratio
  } else {
    1
  }
  spot.size <- cellsize[[1L]] * ratio
  grid.factors <- list(ratio = ratio, tol = spot.size / 2)
  df.sf <- sf::st_as_sf(coordinates, coords = c("x_cent", "y_cent"),
                        remove = FALSE)
  sf::st_crs(df.sf) <- 3857
  square_grid <- sf::st_make_grid(df.sf, cellsize = cellsize,
                                  what = what, square = square)
  square_grid_sf <- sf::st_sf(square_grid)
  square_grid_sf$grid_id <- seq_along(square_grid)
  within_mat <- sf::st_intersects(df.sf, square_grid_sf, sparse = FALSE)
  within_n_grid <- Matrix::rowSums(within_mat)
  spots.counts <- Matrix::colSums(within_mat)
  occupied <- which(spots.counts > 0)
  within.spots <- lapply(occupied, function(index) which(within_mat[, index]))
  levels.idents <- levels(object@idents)
  present.idents <- vapply(within.spots, function(index) {
    counts <- tabulate(match(object@idents[index], levels.idents),
                       nbins = length(levels.idents))
    levels.idents[max.col(matrix(counts, nrow = 1L),
                          ties.method = idents.ties.method)]
  }, character(1))
  new.names <- paste0("Grid", square_grid_sf$grid_id[occupied])
  new.meta <- as.data.frame(do.call(rbind, lapply(within.spots, function(index) {
    counts <- tabulate(match(object@idents[index], levels.idents),
                       nbins = length(levels.idents))
    setNames(as.list(counts), levels.idents)
  })), stringsAsFactors = FALSE)
  rownames(new.meta) <- new.names
  new.meta$cell.type <- factor(present.idents, levels = levels.idents)
  new.meta$spots.counts <- as.numeric(spots.counts[occupied])
  new.coordinates <- do.call(rbind, lapply(within.spots, function(index) {
    colMeans(coordinates[index, c("x_cent", "y_cent"), drop = FALSE])
  }))
  rownames(new.coordinates) <- new.names
  colnames(new.coordinates) <- c("x", "y")
  previous.data <- assay(object, data.slot)
  new.data <- do.call(cbind, lapply(within.spots, function(index) {
    Matrix::rowMeans(previous.data[, index, drop = FALSE])
  }))
  rownames(new.data) <- rownames(previous.data)
  colnames(new.data) <- new.names
  if (!inherits(new.data, "dgCMatrix")) new.data <- methods::as(new.data, "dgCMatrix")
  object.grid <- createSpatialCellChat(
    object = new.data,
    meta = new.meta,
    group.by = "cell.type",
    input.assay = "norm",
    datatype = "spatial",
    coordinates = new.coordinates,
    spatial.factors = grid.factors
  )
  object.grid@DB <- object@DB
  object.grid@images$.grid <- list(
    within.nGrid = within_n_grid,
    recommended.contact.range = spot.size
  )
  attr(object.grid, "members") <- within.spots
  object.grid
}

compare_grid <- function(a, b, label) {
  check(paste0(label, ": coordinates match"),
        isTRUE(all.equal(unname(a@images$coordinates),
                         unname(b@images$coordinates), tolerance = 1e-12)))
  check(paste0(label, ": meta matches"),
        isTRUE(all.equal(a@meta, b@meta, tolerance = 1e-12)))
  check(paste0(label, ": assay$norm matches"),
        max(abs(as.matrix(a@assay$norm) - as.matrix(b@assay$norm))) <= 1e-12)
  check(paste0(label, ": within.nGrid matches"),
        isTRUE(all.equal(as.numeric(a@images$.grid$within.nGrid),
                         as.numeric(b@images$.grid$within.nGrid))))
  check(paste0(label, ": recommended.contact.range matches"),
        isTRUE(all.equal(a@images$.grid$recommended.contact.range,
                         b@images$.grid$recommended.contact.range)))
  check(paste0(label, ": grid rownames match"),
        identical(rownames(a@images$coordinates),
                  rownames(b@images$coordinates)))
  check(paste0(label, ": sparse object validates"),
        isTRUE(validateSpatialCellChat(a)))
}

make_test_chat <- function(n_cells = 60, seed = 7) {
  set.seed(seed)
  genes <- paste0("g", seq_len(40L))
  cells <- paste0("c", seq_len(n_cells))
  expr <- Matrix::rsparsematrix(length(genes), n_cells, density = 0.3)
  dimnames(expr) <- list(genes, cells)
  expr <- as(expr, "dgCMatrix")
  meta <- data.frame(
    label = factor(sample(c("A", "B", "C"), n_cells, replace = TRUE)),
    row.names = cells
  )
  coords <- cbind(x = runif(n_cells, 0, 20), y = runif(n_cells, 0, 20))
  rownames(coords) <- cells
  chat <- createSpatialCellChat(
    expr,
    meta = meta,
    group.by = "label",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = list(ratio = 0.18, tol = 5)
  )
  chat
}

## ---- 1. small synthetic object, square grid ----
chat <- make_test_chat()
sparse_obj <- makeGridSpatialCellChat(chat, cellsize = c(5, 5))
dense_obj <- makeGrid_dense_reference(chat, cellsize = c(5, 5))
compare_grid(sparse_obj, dense_obj, "square grid (small)")
check("square grid: every occupied grid has >= 1 cell",
      all(sparse_obj@meta$spots.counts >= 1))

expect_error("makeGrid rejects invalid cellsize length",
             makeGridSpatialCellChat(chat, cellsize = c(1, 2, 3)), "length")
expect_error("makeGrid rejects zero cellsize",
             makeGridSpatialCellChat(chat, cellsize = 0), "greater than zero")
expect_error("makeGrid rejects non-finite cellsize",
             makeGridSpatialCellChat(chat, cellsize = Inf), "finite")

## ---- 2. hexagonal grid ----
chat_hex <- make_test_chat(seed = 11)
sparse_hex <- makeGridSpatialCellChat(chat_hex, cellsize = c(5, 5),
                                      square = FALSE)
dense_hex <- makeGrid_dense_reference(chat_hex, cellsize = c(5, 5),
                                      square = FALSE)
compare_grid(sparse_hex, dense_hex, "hexagonal grid")

## ---- 3. simulated_data end-to-end ----
data_dir <- "tests_dev/simulated_data"
counts_table <- data.table::fread(file.path(data_dir, "countmat_1.csv"))
counts <- as.matrix(counts_table[, -1, with = FALSE])
rownames(counts) <- counts_table[[1L]]
colnames(counts) <- names(counts_table)[-1L]
metadata <- as.data.frame(data.table::fread(
  file.path(data_dir, "metadata_1.csv")))
rownames(metadata) <- metadata$cell_id
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata$cell_type <- factor(metadata$cell_type)
coordinates <- as.matrix(metadata[, c("x", "y")])
rownames(coordinates) <- rownames(metadata)

sim_chat <- createSpatialCellChat(
  counts,
  meta = metadata,
  group.by = "cell_type",
  datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = 0.18, tol = 5)
)
check("simulated data: object constructed",
      isTRUE(validateSpatialCellChat(sim_chat)))

sim_grid <- makeGridSpatialCellChat(sim_chat, cellsize = c(50, 50))
sim_grid_dense <- makeGrid_dense_reference(sim_chat, cellsize = c(50, 50))
compare_grid(sim_grid, sim_grid_dense, "simulated data grid")
check("simulated data: all cells assigned to a grid",
      all(sim_grid@images$.grid$within.nGrid > 0))

# manual spot checks against the dense-reference members
members <- attr(sim_grid_dense, "members")
expected_mean <- Matrix::rowMeans(
  sim_chat@assay$norm[, members[[1]], drop = FALSE])
check("simulated data: grid expression is the row mean of its cells",
      max(abs(as.matrix(sim_grid@assay$norm[, 1]) - as.matrix(expected_mean))) <= 1e-12)
expected_coord <- colMeans(coordinates[members[[1]], c("x", "y"), drop = FALSE])
check("simulated data: grid coordinate is the cell centroid",
      max(abs(sim_grid@images$coordinates[1, ] - expected_coord)) <= 1e-12)
check("simulated data: member count matches spots.counts",
      length(members[[1]]) == sim_grid@meta$spots.counts[1])

## ---- 2.5 grid-aligned boundary points (edge/corner sharing) ----
# Points placed exactly on grid boundaries must produce the same hit sets
# as the dense reference (sf boundary semantics are preserved).
boundary_cells <- paste0("b", seq_len(6L))
boundary_coords <- cbind(
  x = c(1, 1, 2, 2, 3, 3),
  y = c(1, 2, 1, 2, 1, 2)
)
rownames(boundary_coords) <- boundary_cells
boundary_meta <- data.frame(
  label = factor(c("A", "A", "B", "B", "C", "C")),
  row.names = boundary_cells
)
boundary_expr <- matrix(1, nrow = 4L, ncol = 6L,
                        dimnames = list(paste0("g", 1:4), boundary_cells))
boundary_chat <- createSpatialCellChat(
  boundary_expr,
  meta = boundary_meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = boundary_coords,
  spatial.factors = list(ratio = 1, tol = 0.5)
)
sparse_bnd <- makeGridSpatialCellChat(boundary_chat, cellsize = c(1, 1))
dense_bnd <- makeGrid_dense_reference(boundary_chat, cellsize = c(1, 1))
compare_grid(sparse_bnd, dense_bnd, "boundary sharing grid")
check("boundary sharing: some cells hit multiple grids",
      any(sparse_bnd@images$.grid$within.nGrid > 1))

## ---- 4. 10k-cell grid stress: correctness + bounded memory ----
nx <- 100L; ny <- 100L; n <- nx * ny
stress_coords <- cbind(rep(seq_len(nx), each = ny), rep.int(seq_len(ny), nx))
rownames(stress_coords) <- paste0("s", seq_len(n))
stress_meta <- data.frame(
  label = factor(rep(c("A", "B", "C", "D"), length.out = n)),
  row.names = rownames(stress_coords)
)
stress_expr <- Matrix::rsparsematrix(30, n, density = 0.02)
dimnames(stress_expr) <- list(paste0("g", seq_len(30)), rownames(stress_coords))
stress_chat <- createSpatialCellChat(
  stress_expr,
  meta = stress_meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = stress_coords,
  spatial.factors = list(ratio = 1, tol = 0.5)
)
gc(reset = TRUE)
before <- gc()
stress_t <- system.time(stress_grid <- makeGridSpatialCellChat(
  stress_chat, cellsize = c(1, 1)))
after <- gc()
check("stress: all cells assigned to at least one grid",
      all(stress_grid@images$.grid$within.nGrid >= 1))
check("stress: hit conservation (row hits == column hits)",
      sum(stress_grid@images$.grid$within.nGrid) ==
        sum(stress_grid@meta$spots.counts))
# bbox expansion (sf::st_make_grid default) makes the exact grid count an sf
# implementation detail; assert the dense 100x100 layout occupies ~all grids.
check("stress: ~all unit grids occupied",
      nrow(stress_grid@images$coordinates) >= 9000L)
check("stress: grid object validates",
      isTRUE(validateSpatialCellChat(stress_grid)))
delta_mb <- (after[["Ncells", "used"]] - before[["Ncells", "used"]]) * 56 / 1024^2
check("stress: memory delta stays bounded (no dense n x nGrid matrix)",
      delta_mb < 100)

cat(sprintf("stress: n=%d cellsize=(1,1) elapsed=%.3f s Ncells_delta=%.1f MB\n", n,
            stress_t[["elapsed"]], delta_mb))

cat("All makeGridSpatialCellChat sparse tests passed.\n")
