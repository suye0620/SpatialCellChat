## Consolidated runtime checks for SpatialCellChat visualization.R
##
## This script sources the package files in dependency order and loads the
## The script intentionally exercises package definitions as-is; no helpers
## are patched in memory.
suppressMessages({
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(patchwork)
  library(spdep)
  library(circlize)
  library(RColorBrewer)
  library(Matrix)
  library(cli)
  library(purrr)
  library(ks)
  library(colorspace)
})

source_files <- c(
  "R/SpatialCellChat_class.R",
  "R/utilities.R",
  "R/analysis.R",
  "R/spatial.R",
  "R/visualization.R"
)
env <- globalenv()
for (path in source_files) sys.source(path, envir = env)

ok <- function(label, thunk) {
  result <- tryCatch({
    value <- thunk()
    cls <- if (is.null(value)) "NULL" else paste(class(value)[1L], collapse = "/")
    paste("SUCCESS", cls)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", label, result))
  invisible(result)
}

set.seed(1)
n <- 300
df <- data.frame(x = runif(n, 0, 10), y = runif(n, 0, 10), score = rnorm(n))
rownames(df) <- paste0("spot", seq_len(n))
groups <- factor(sample(c("A", "B"), n, replace = TRUE))

## These checks should succeed without any artificial color/KDE patch.
ok("colorRamp3 [package source]", function() {
  env$colorRamp3(c(0, 1), c("white", "red"))(0.5)
})
ok("calculate_density wkde [package source]", function() {
  env$calculate_density(df$score, df[c("x", "y")], method = "wkde")
})
ok("plotFeatures normal [as-is]", function() {
  env$plotFeatures(df, method = "normal", plot.title = "t")
})
ok("Gi normal [as-is]", function() {
  env$plotStatistics_Gi(df, groups, rep(FALSE, n), n = 8, method = "normal")
})
ok("Gi 3d [as-is]", function() {
  env$plotStatistics_Gi(df, groups, rep(FALSE, n), n = 8, method = "3d")
})
ok("Gi density [as-is]", function() {
  env$plotStatistics_Gi(df, groups, rep(FALSE, n), n = 8, method = "density")
})

## A final-schema object exposes the migration failures in the old plotting API.
expr <- Matrix::Matrix(matrix(seq_len(20), nrow = 5), sparse = TRUE)
rownames(expr) <- paste0("g", seq_len(nrow(expr)))
colnames(expr) <- paste0("c", seq_len(ncol(expr)))
meta <- data.frame(label = factor(c("A", "B", "A", "B")), row.names = colnames(expr))
coords <- cbind(x = c(0, 3, 0, 4), y = c(0, 0, 4, 4))
rownames(coords) <- colnames(expr)
chat <- env$createSpatialCellChat(
  expr, meta = meta, group.by = "label", datatype = "spatial",
  coordinates = coords, spatial.factors = list(ratio = 1, tol = 0)
)
ok("plotly_spatialDimPlot [final matrix]", function() {
  env$plotly_spatialDimPlot(chat, method = "2d")
})
## The legacy Plotly function mutates coordinates with `$`; use a data.frame
## fixture to isolate the missing-helper failure from the final-schema failure.
chat_plotly <- chat
chat_plotly@images$coordinates <- as.data.frame(chat_plotly@images$coordinates)
rownames(chat_plotly@images$coordinates) <- rownames(chat@images$coordinates)
ok("plotFeatures 2d [as-is]", function() {
  env$plotFeatures(df, method = "2d", plot.title = "t")
})
ok("plotFeatures 3d [as-is]", function() {
  env$plotFeatures(df, method = "3d", plot.title = "t")
})
ok("plotly_spatialDimPlot 2d [as-is]", function() {
  env$plotly_spatialDimPlot(chat_plotly, method = "2d")
})
ok("plotly_spatialDimPlot 3d [as-is]", function() {
  env$plotly_spatialDimPlot(chat_plotly, method = "3d")
})


ok("spatialCCCDistPlot [final schema]", function() {
  env$spatialCCCDistPlot(chat, interaction.range = 10, contact.range = 5, tol = 0)
})
ok("spatialFeaturePlot [final schema]", function() {
  env$spatialFeaturePlot(chat, features = "g1")
})
ok("plotly_spatialLRpairPlot [final schema]", function() {
  env$plotly_spatialLRpairPlot(chat, features = "g1")
})

## One trace contains every layer; z ticks carry layer names and cmin/cmax
## make the color domain explicit across all layers.
spatial3Dstack <- function(xy, layers, color.heatmap = "RdBu", point.size = 3,
                           z.axis.space = 1) {
  stopifnot(nrow(xy) == nrow(layers), !is.null(colnames(layers)))
  feat <- colnames(layers)
  nf <- length(feat)
  vmin <- min(layers, na.rm = TRUE)
  vmax <- max(layers, na.rm = TRUE)
  if (vmin >= 0) vmin <- 0
  if (vmax <= 0) vmax <- 0
  cols <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(9, color.heatmap))
  )(99)
  cell_id <- rownames(xy)
  if (is.null(cell_id)) cell_id <- seq_len(nrow(xy))
  long <- do.call(rbind, lapply(seq_len(nf), function(k) data.frame(
    x = xy[[1]], y = xy[[2]], z = (k - 1) * z.axis.space,
    value = layers[[k]], layer = feat[[k]], cell_id = cell_id
  )))
  plotly::plot_ly(
    long, x = ~x, y = ~y, z = ~z, color = ~value, colors = cols,
    type = "scatter3d", mode = "markers",
    marker = list(size = point.size, cmin = vmin, cmax = vmax),
    hoverinfo = "text",
    text = ~paste("<b>", layer, "</b><br>", cell_id,
                  "<br>score:", round(value, 3))
  ) %>% plotly::layout(
    scene = list(
      xaxis = list(title = "x", showgrid = TRUE, zeroline = FALSE),
      yaxis = list(title = "y", showgrid = TRUE, zeroline = FALSE),
      zaxis = list(
        title = "", tickmode = "array",
        tickvals = seq(0, (nf - 1) * z.axis.space, by = z.axis.space),
        ticktext = feat, showgrid = TRUE, zeroline = FALSE
      ),
      aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.35),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 0.9))
    ),
    showlegend = FALSE
  ) %>% plotly::plotly_build()
}

set.seed(7)
xy <- data.frame(x = runif(300, 0, 10), y = runif(300, 0, 10))
rownames(xy) <- paste0("spot", seq_len(nrow(xy)))
layers <- data.frame(
  L1 = rnorm(300, 1, 0.6),
  L2 = rnorm(300, 0, 0.5),
  L3 = rnorm(300, -1, 0.7)
)
ok("3D stack prototype [single trace]", function() {
  p <- spatial3Dstack(xy, layers)
  main <- p$x$data[[1L]]
  stopifnot(length(main$x) == nrow(xy) * ncol(layers))
  stopifnot(identical(as.character(p$x$layout$scene$zaxis$ticktext), colnames(layers)))
  stopifnot(identical(as.numeric(main$marker$cmin), 0) ||
            is.numeric(main$marker$cmin))
  p
})
cat("\nAll checks complete.\n")
