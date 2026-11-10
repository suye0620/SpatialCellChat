args <- commandArgs(TRUE)
if (length(args) != 3L) stop('usage: Rscript test-installed-visualization-smoke.R LIB CHAT OUTDIR')
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
suppressPackageStartupMessages(library(SpatialCellChat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(plotly))
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
chat <- readRDS(args[[2L]])
stopifnot(inherits(chat, 'SpatialCellChat'))
files <- character()
write_plot <- function(name, plot, width = 7, height = 5.5) {
  file <- file.path(args[[3L]], paste0(name, '.pdf'))
  if (inherits(plot, 'Heatmap') || inherits(plot, 'HeatmapList')) {
    grDevices::pdf(file, width = width, height = height)
    tryCatch(ComplexHeatmap::draw(plot), finally = grDevices::dev.off())
  } else {
    ggplot2::ggsave(file, plot, width = width, height = height, limitsize = FALSE)
  }
  stopifnot(file.exists(file), file.info(file)$size > 1000)
  files <<- c(files, file)
}
write_plot('spatialC2CDistPlot-geometry', spatialC2CDistPlot(chat, mode = 'geometry'))
write_plot('netVisual_heatmap', netVisual_heatmap(chat), width = 7, height = 7)
feature <- rownames(chat@assay$signaling)[[1L]]
write_plot('spatialFeaturePlot', spatialFeaturePlot(chat, features = feature), width = 7, height = 5.5)
coords <- chat@images$coordinates
layers <- t(as.matrix(chat@assay$norm[head(rownames(chat@assay$norm), 3L), , drop = FALSE]))
rownames(layers) <- rownames(coords)
colnames(layers) <- c('layer1', 'layer2', 'layer3')
stack <- plotly_spatialStackPlot(chat, layers = layers, title.name = 'installed smoke')
stack.file <- file.path(args[[3L]], 'plotly_spatialStackPlot.html')
htmlwidgets::saveWidget(stack, stack.file, selfcontained = FALSE, libdir = 'plotly_files')
license.paths <- c(
  'plotly-latest.min.js' = system.file('htmlwidgets', 'lib', 'plotlyjs', 'LICENSE', package = 'plotly'),
  'typedarray.min.js' = system.file('htmlwidgets', 'lib', 'typedarray', 'LICENSE', package = 'plotly')
)
for (asset in list.files(file.path(args[[3L]], 'plotly_files'), recursive = TRUE, full.names = TRUE)) {
  base <- basename(asset)
  if (base %in% names(license.paths)) {
    license.file <- file.path(dirname(asset), 'LICENSE')
    stopifnot(file.copy(license.paths[[base]], license.file, overwrite = TRUE))
  }
}
stopifnot(file.exists(stack.file), file.info(stack.file)$size > 1000)
files <- c(files, stack.file)
summary <- data.frame(file = basename(files), status = 'PASS', stringsAsFactors = FALSE)
write.csv(summary, file.path(args[[3L]], 'summary.csv'), row.names = FALSE)
cat('Installed visualization smoke passed:', paste(basename(files), collapse = ', '), '\n')
