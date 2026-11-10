% Generated from the roxygen comments in R/visualization.R
\name{netVisual_CommunFieldGrid}
\alias{netVisual_CommunFieldGrid}
\title{netVisual_CommunFieldGrid}
\usage{
netVisual_CommunFieldGrid(object, signaling, slot.name = "netP",
  pattern = c("incoming", "outgoing"), cellsize = NULL, what = "polygons",
  square = T, grid.resolution = NULL, grid.field.mean = c("median", "sum"),
  color.use = NULL, point.size = 2, image.alpha = 0.3, min.mag = 0,
  arrow.skip = NULL, arrow.line.color = "grey20", arrow.line.width = c(0.35, 1.4),
  arrow.line.alpha = 0.9, arrow.angle = 22.5, arrow.length = 0.5,
  arrow.type = "closed", lineend = "round", ..., arrow.scale = 1, legend.size = 2,
  legend.text.size = 8, title.name = NULL)
}
\arguments{
\item{object}{SpatialCellChat object}
\item{signaling}{A signaling pathway or ligand-receptor pair to visualize.}
\item{slot.name}{Use \code{"netP"} for pathway names or \code{"net"} for ligand-receptor pair names.}
\item{pattern}{\code{"outgoing"} or \code{"incoming"}.}
\item{cellsize}{Positive finite numeric vector of length 1 or 2. If \code{NULL}, the minimum non-self nearest-neighbor distance is estimated without a dense pairwise distance matrix.}
\item{what}{Grid geometry returned by \code{sf::st_make_grid}: \code{"polygons"}, \code{"corners"}, or \code{"centers"}.}
\item{square}{Logical; if \code{FALSE}, use hexagonal cells.}
\item{grid.resolution}{Positive finite multiplier applied component-wise to \code{cellsize}; defaults to 2.}
\item{grid.field.mean}{How vector fields within each grid cell are aggregated: \code{"median"} or \code{"sum"}.}
\item{color.use}{Colors for each cell group.}
\item{point.size}{Size of spots.}
\item{image.alpha}{Transparency of individual spots.}
\item{min.mag}{Minimum vector magnitude required for arrows.}
\item{arrow.skip}{Number of interpolated grid locations skipped between displayed arrows. \code{NULL} adaptively targets at most 12 arrows per axis.}
\item{arrow.line.color}{Color of arrows or streamlines.}
\item{arrow.line.width}{One or two positive values controlling width; a two-value range maps width to magnitude.}
\item{arrow.line.alpha}{Transparency of arrows or streamlines.}
\item{arrow.angle}{Arrowhead angle.}
\item{arrow.length}{Arrowhead length.}
\item{arrow.type}{Grid arrowhead type.}
\item{lineend}{Line ending style.}
\item{...}{Other parameters passed to the underlying metR geom.}
\item{arrow.scale}{Positive finite scalar multiplying displayed arrow shaft lengths after filtering; defaults to 1. Does not change stored field values.}
\item{legend.size}{Legend point size.}
\item{legend.text.size}{Legend text size.}
\item{title.name}{Plot title.}
}
\value{A ggplot2 object.}
\description{Aggregate and plot a communication vector field on a spatial grid.}
\details{Grid geometry uses the first two canonical \code{images$coordinates} axes after the existing visualization coordinate preparation. If \code{images$spatial.factors$ratio} is unavailable, coordinates remain in raw units and no physical scale is reported. Zero-length vectors are transparent; nonzero arrows retain the requested opacity. \code{arrow.scale} changes only displayed arrow shaft lengths and is shown in the caption when not 1; use the same value when comparing lengths across plots. Large values may make arrows overlap.}
