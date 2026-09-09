# Visium histology image and spatialDimPlot development tests

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
  library(Seurat)
  library(ggplot2)
  library(png)
  library(jsonlite)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/visualization.R")

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
coords <- matrix(c(20, 10, 40, 20, 60, 30, 80, 40), ncol = 2L, byrow = TRUE)
rownames(coords) <- colnames(expr)
colnames(coords) <- c("x", "y")

chat <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 5)
)

plain_plot <- spatialDimPlot(chat, image = FALSE)
check("spatialDimPlot works without histology image", inherits(plain_plot, "ggplot"))

distance_plot <- spatialCCCDistPlot(
  chat, signaling.type = "Secreted",
  interaction.range = 100, contact.range = 30, tol = 0
)
check("spatialCCCDistPlot computes from canonical images", inherits(distance_plot, "ggplot"))
chat@images$.distance <- computeCellDistance(
  chat@images$coordinates, interaction.range = 100,
  contact.range = 30, ratio = 1, tol = 0
)
cached_distance_plot <- spatialCCCDistPlot(
  chat, signaling.type = "Contact",
  interaction.range = 100, contact.range = 30, tol = 0
)
check("spatialCCCDistPlot consumes images$.distance", inherits(cached_distance_plot, "ggplot"))
narrow_distance_plot <- spatialCCCDistPlot(
  chat, signaling.type = "Secreted",
  interaction.range = 5, contact.range = 3, tol = 0
)
check("spatialCCCDistPlot recomputes for changed thresholds", inherits(
  narrow_distance_plot, "ggplot"
) && all(narrow_distance_plot$data$x <= 5))

image_array <- array(0.8, dim = c(100L, 100L, 3L))
image_coordinates <- data.frame(
  barcodes = rownames(coords),
  tissue = 1L,
  row = seq_len(nrow(coords)),
  col = seq_len(nrow(coords)),
  imagerow = coords[, 1],
  imagecol = coords[, 2],
  row.names = rownames(coords),
  check.names = FALSE
)
raster <- .sc_normalize_raster(list(
  image = image_array,
  scale.factors = list(spot = 1, fiducial = 1, hires = 1, lowres = 1),
  spot.coordinates = image_coordinates,
  source = list(type = "test")
), rownames(coords))
images(chat, "rasters") <- list(main = raster)
check("raster image is stored in images", is.list(chat@images$rasters$main))
check("raster coordinates are aligned to cells", identical(
  rownames(chat@images$rasters$main$spot.coordinates), colnames(expr)
))
check("object remains valid with raster image", isTRUE(validateSpatialCellChat(chat)))

image_plot <- spatialDimPlot(chat, image = TRUE, image.alpha = 0.25)
check("spatialDimPlot adds raster background", inherits(image_plot, "ggplot") &&
        length(image_plot$layers) >= 2L)

if (requireNamespace("scatterpie", quietly = TRUE)) {
  proportion <- data.frame(A = c(1, 0, 0.5, 0), B = c(0, 1, 0.5, 1))
  proportion_plot <- spatialDimPlot(chat, proportion = proportion, image = TRUE)
  check("spatialDimPlot keeps proportion mode with image", inherits(proportion_plot, "ggplot"))
} else {
  cat("[SKIP] proportion mode requires scatterpie\n")
}

image_dir <- file.path(tempdir(), "spatialcellchat_visium_image")
dir.create(image_dir, showWarnings = FALSE)
png::writePNG(image_array, file.path(image_dir, "tissue_lowres_image.png"))
jsonlite::write_json(
  list(
    spot_diameter_fullres = 10,
    fiducial_diameter_fullres = 20,
    tissue_hires_scalef = 0.5,
    tissue_lowres_scalef = 0.25
  ),
  file.path(image_dir, "scalefactors_json.json"), auto_unbox = TRUE
)
utils::write.table(
  image_coordinates,
  file.path(image_dir, "tissue_positions.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)
read_image <- readSpatialImage(image_dir)
check("readSpatialImage reads Visium raster", is.array(read_image$image))
check("readSpatialImage reads scale factors", identical(
  names(read_image$scale.factors), c("spot", "fiducial", "hires", "lowres")
))
check("readSpatialImage reads spot coordinates", identical(
  rownames(read_image$spot.coordinates), rownames(image_coordinates)
))

chat_from_image <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  datatype = "spatial",
  coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 5),
  image = read_image
)
check("createSpatialCellChat accepts raster list", isTRUE(
  validateSpatialCellChat(chat_from_image)
))
check("createSpatialCellChat stores raster list", is.list(chat_from_image@images$rasters$main))

seurat_counts <- abs(as.matrix(expr))
seurat_counts[seurat_counts == 0] <- 1
seurat_object <- Seurat::CreateSeuratObject(counts = seurat_counts)
seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
seurat_object@meta.data$label <- meta$label
seurat_image <- Seurat::Read10X_Image(image.dir = image_dir, filter.matrix = FALSE)
seurat_image <- seurat_image[colnames(expr)]
SeuratObject::DefaultAssay(seurat_image) <- Seurat::DefaultAssay(seurat_object)
seurat_object[["slice1"]] <- seurat_image
chat_from_seurat <- createSpatialCellChat(
  seurat_object,
  group.by = "label",
  datatype = "spatial",
  coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 5)
)
check("createSpatialCellChat imports Seurat raster", is.list(
  chat_from_seurat@images$rasters$main
))
check("Seurat raster has image and coordinates", is.array(
  chat_from_seurat@images$rasters$main$image
) && is.data.frame(chat_from_seurat@images$rasters$main$spot.coordinates))

cat("All SpatialCellChat spatial image checks passed.\n")
