# Phase 0/1 tests for the final SpatialCellChat schema

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
})
source("R/SpatialCellChat_class.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

expr <- Matrix::rsparsematrix(5, 4, density = 0.5)
rownames(expr) <- paste0("gene", seq_len(nrow(expr)))
colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
meta <- data.frame(
  label = factor(c("A", "B", "A", "B")),
  row.names = colnames(expr)
)

chat_rna <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  datatype = "RNA"
)
check("RNA object uses the canonical 11 slots", identical(
  methods::slotNames(chat_rna),
  c("assay", "images", "meta", "idents", "features", "LR", "dr",
    "net", "netP", "DB", "misc")
))
check("RNA assay uses smooth instead of project", identical(
  names(chat_rna@assay), c("raw", "norm", "scale", "smooth", "signaling")
))
check("RNA image state uses canonical defaults", identical(
  chat_rna@images,
  list(
    coordinates = NULL,
    coordinate.system = NULL,
    spatial.factors = NULL,
    rasters = list(),
    fov = list(),
    .distance = list(),
    .grid = list()
  )
))
check("RNA object stores idents as a factor", is.factor(chat_rna@idents) &&
        identical(names(chat_rna@idents), colnames(expr)) &&
        identical(as.character(chat_rna@idents), as.character(meta$label)))
check("idents accessor returns the factor", is.factor(idents(chat_rna)) &&
        identical(idents(chat_rna), chat_rna@idents))
check("single mode has no datasets", identical(chat_rna@misc$.mode, "single") &&
        length(chat_rna@misc$.datasets) == 0L)
check("single mode records datatype and construction parameters", identical(
  chat_rna@misc$.datatype, "RNA"
) && identical(
  names(chat_rna@misc$.param),
  c("sample", "group.by", "input.assay", "normalize", "scale.factor", "do.log",
    "coordinate.system", "spatial.factors")
))
check("RNA object validates", isTRUE(validateSpatialCellChat(chat_rna)))

coords <- matrix(seq_len(8), ncol = 2L, byrow = TRUE)
rownames(coords) <- colnames(expr)
chat_spatial <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = coords,
  spatial.factors = list(ratio = 0.18, tol = 5)
)
check("spatial coordinates are stored with canonical axes", identical(
  rownames(chat_spatial@images$coordinates), colnames(expr)
) && identical(colnames(chat_spatial@images$coordinates), c("x", "y")))
check("spatial image state stores calibrated coordinate metadata", identical(
  chat_spatial@images$coordinate.system,
  list(unit = "pixel", calibrated = TRUE, y.direction = "up", origin = "lower-left")
) && identical(
  chat_spatial@images$spatial.factors, list(ratio = 0.18, tol = 5)
))
check("spatial distance cache is an images list", is.list(chat_spatial@images$.distance))
distance_matrix <- Matrix::Diagonal(ncol(expr))
dimnames(distance_matrix) <- list(colnames(expr), colnames(expr))
chat_with_distance <- chat_spatial
chat_with_distance@images$.distance$d.spatial <- distance_matrix
check("compatible spatial distance cache validates", isTRUE(
  validateSpatialCellChat(chat_with_distance)
))
check("spatial object validates", isTRUE(validateSpatialCellChat(chat_spatial)))

bad <- chat_spatial
bad@meta <- bad@meta[-1, , drop = FALSE]
check("invalid metadata is rejected", inherits(
  try(methods::validObject(bad), silent = TRUE), "try-error"
))

bad_coords <- chat_spatial
bad_coords@images$coordinates <- bad_coords@images$coordinates[-1, , drop = FALSE]
check("invalid spatial coordinates are rejected", length(
  validateSpatialCellChat(bad_coords, strict = FALSE)
) > 0L
)
bad_distance <- chat_spatial
bad_distance@images$.distance$d.spatial <- Matrix::Diagonal(ncol(expr) - 1L)
check("invalid spatial distance cache is rejected", length(
  validateSpatialCellChat(bad_distance, strict = FALSE)
) > 0L)


expr_duplicate <- expr
rownames(expr_duplicate)[2] <- rownames(expr_duplicate)[1]
check("duplicate gene names are rejected", inherits(
  try(createSpatialCellChat(expr_duplicate, meta = meta, group.by = "label"), silent = TRUE),
  "try-error"
))

chat_uncalibrated <- createSpatialCellChat(
  expr, meta = meta, group.by = "label", datatype = "spatial",
  coordinates = coords,
  coordinate.system = list(unit = "arbitrary", calibrated = FALSE),
  spatial.factors = list(ratio = NULL, tol = NULL)
)
check("uncalibrated spatial coordinates remain valid", isTRUE(
  validateSpatialCellChat(chat_uncalibrated)
) && identical(
  chat_uncalibrated@images$coordinate.system,
  list(unit = "arbitrary", calibrated = FALSE, y.direction = "up", origin = "lower-left")
) && identical(
  chat_uncalibrated@images$spatial.factors, list(ratio = NULL, tol = NULL)
))

check("calibrated coordinates require a ratio", inherits(
  try(createSpatialCellChat(
    expr, meta = meta, group.by = "label", datatype = "spatial",
    coordinates = coords,
    coordinate.system = list(unit = "pixel", calibrated = TRUE),
    spatial.factors = list(tol = 5)
  ), silent = TRUE), "try-error"
))
check("calibrated coordinates require a tolerance", inherits(
  try(createSpatialCellChat(
    expr, meta = meta, group.by = "label", datatype = "spatial",
    coordinates = coords,
    coordinate.system = list(unit = "pixel", calibrated = TRUE),
    spatial.factors = list(ratio = 0.18)
  ), silent = TRUE), "try-error"
))

bad_images <- chat_spatial
bad_images@images$histology <- list(image = array(1, dim = c(2L, 2L)))
check("obsolete histology records are rejected", length(
  validateSpatialCellChat(bad_images, strict = FALSE)
) > 0L)
coords_3d <- cbind(coords, z = seq_len(nrow(coords)))
chat_3d <- createSpatialCellChat(
  expr, meta = meta, group.by = "label", datatype = "spatial",
  coordinates = coords_3d,
  spatial.factors = list(ratio = 0.18, tol = 5)
)
check("three-dimensional coordinates use canonical axes", isTRUE(
  validateSpatialCellChat(chat_3d)
) && identical(colnames(chat_3d@images$coordinates), c("x", "y", "z")))

cat("All SpatialCellChat class checks passed.\n")
