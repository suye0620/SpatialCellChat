# Regression and contract tests for scalable grid-size estimation.
# This script is intentionally implementation-independent for the nearest-distance
# reference: base R `dist()` is the correctness oracle.

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(sf)
  library(ggplot2)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/visualization.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
expect_error <- function(label, code, pattern = NULL) {
  error <- tryCatch({force(code); NULL}, error = identity)
  check(label, !is.null(error))
  if (!is.null(pattern)) {
    check(paste0(label, ": message"), grepl(pattern, conditionMessage(error)))
  }
  invisible(error)
}

make_grid_chat <- function(coordinates, ratio = 0.18) {
  coordinates <- as.matrix(coordinates)
  cells <- paste0("cell_", seq_len(nrow(coordinates)))
  rownames(coordinates) <- cells
  expr <- Matrix::Matrix(
    matrix(seq_len(3L * nrow(coordinates)), nrow = 3L,
           dimnames = list(paste0("g", 1:3), cells)), sparse = TRUE)
  meta <- data.frame(
    label = factor(rep(c("A", "B"), length.out = nrow(coordinates))),
    row.names = cells
  )
  createSpatialCellChat(
    expr, meta = meta, group.by = "label", datatype = "spatial",
    coordinates = coordinates,
    spatial.factors = list(ratio = ratio, tol = ratio)
  )
}

coords <- cbind(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1))
chat <- make_grid_chat(coords)
check("canonical coordinates are x/y", identical(colnames(chat@images$coordinates), c("x", "y")))

# The implementation-independent nearest-neighbor oracle.
reference <- min(as.vector(dist(coords)))
check("base-R reference is finite", is.finite(reference) && reference == 1)
check("grid resolver exists", exists(".sc_resolve_grid_size", mode = "function"))
if (exists(".sc_resolve_grid_size", mode = "function")) {
  resolved <- .sc_resolve_grid_size(coords, grid.resolution = 2,
                                    ratio = chat@images$spatial.factors$ratio)
  check("default size uses nearest-neighbor distance",
        isTRUE(all.equal(resolved$base.cellsize, reference, tolerance = 1e-12)))
  check("default size applies resolution",
        isTRUE(all.equal(resolved$effective.cellsize, c(2, 2), tolerance = 1e-12)))
  coords3d <- cbind(coords, z = c(0, 100, 100, 0))
  reference3d <- min(as.vector(dist(coords3d)))
  check("three-dimensional coordinates support nearest-neighbor size",
        isTRUE(all.equal(.sc_resolve_grid_size(coords3d, grid.resolution = 1)$base.cellsize,
                         reference3d, tolerance = 1e-12)))
  check("vector cellsize is preserved component-wise",
        isTRUE(all.equal(.sc_resolve_grid_size(coords, cellsize = c(3, 4),
                                               grid.resolution = 2)$effective.cellsize,
                          c(6, 8))))
  duplicate_coords <- rbind(coords, c(0, 0))
  duplicate <- .sc_resolve_grid_size(duplicate_coords, grid.resolution = 1)
  check("duplicate coordinates retain zero spacing", identical(duplicate$base.cellsize, 0))
}

one_point_explicit <- .sc_resolve_grid_size(matrix(c(0, 0), nrow = 1L), cellsize = 2)
check("explicit cellsize does not require nearest-neighbor estimate",
      identical(one_point_explicit$base.cellsize, 2))

expect_error("one point is rejected", .sc_resolve_grid_size(matrix(c(0, 0), nrow = 1L)),
             "at least two")
expect_error("invalid cellsize length is rejected",
             .sc_resolve_grid_size(coords, cellsize = c(1, 2, 3)), "length")
expect_error("zero cellsize is rejected",
             .sc_resolve_grid_size(coords, cellsize = 0), "greater than zero")
expect_error("non-finite cellsize is rejected",
             .sc_resolve_grid_size(coords, cellsize = NaN), "finite")
expect_error("invalid resolution is rejected",
             .sc_resolve_grid_size(coords, grid.resolution = 0), "greater than zero")

options(SpatialCellChat.verbose = 1L)
plot_output <- capture.output(plot <- computeGridSize(
  chat, cellsize = 1, grid.resolution = 1, do.plot = TRUE,
  what = "polygons", square = TRUE
), type = "message")
check("computeGridSize plot path returns ggplot", inherits(plot, "ggplot"))
check("normal output reports effective cellsize",
      any(grepl("Effective grid cellsize", plot_output, fixed = TRUE)))
check("normal output reports occupancy",
      any(grepl("occupied", plot_output, fixed = TRUE)))
check("normal output provides next-step hint",
      any(grepl("makeGridSpatialCellChat", plot_output, fixed = TRUE)))

invisible_output <- capture.output({
  invisible_result <- computeGridSize(chat, cellsize = 1,
                                      grid.resolution = 1, do.plot = FALSE)
}, type = "message")
check("computeGridSize diagnostic path returns NULL", is.null(invisible_result))
check("non-plot output reports resolution",
      any(grepl("resolution", invisible_output, ignore.case = TRUE)))

options(SpatialCellChat.verbose = 2L)
detailed_output <- capture.output(computeGridSize(
  chat, cellsize = 1, grid.resolution = 1, do.plot = TRUE,
  what = "polygons", square = FALSE
), type = "message")
check("detailed output reports point count",
      any(grepl("points", detailed_output, ignore.case = TRUE)))
check("detailed output reports generated grids",
      any(grepl("generated", detailed_output, ignore.case = TRUE)))
check("detailed output reports assignment status",
      any(grepl("assigned", detailed_output, ignore.case = TRUE)))

options(SpatialCellChat.verbose = 3L)
debug_output <- capture.output(computeGridSize(
  chat, cellsize = 1, grid.resolution = 1, do.plot = TRUE,
  what = "polygons", square = TRUE
), type = "message")
check("verbosity 3 reports sparse allocation path",
      any(grepl("Sparse membership", debug_output, fixed = TRUE)))

options(SpatialCellChat.verbose = 0L)
warning_output <- capture.output(computeGridSize(
  chat, cellsize = 1, grid.resolution = 1, do.plot = TRUE,
  what = "centers", square = TRUE
), type = "message")
check("verbosity 0 keeps assignment warnings visible",
      any(grepl("unassigned", warning_output, ignore.case = TRUE)))

 boundary_chat <- make_grid_chat(cbind(x = c(0, 1, 2), y = c(0, 1, 0)), ratio = 1)
check("boundary object validates", isTRUE(validateSpatialCellChat(boundary_chat)))
check("grid membership helper exists", exists(".sc_grid_membership", mode = "function"))

cat("All computeGridSize contract tests passed.\n")
