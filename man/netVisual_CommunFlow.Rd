% Generated from the roxygen comments in R/visualization.R
\name{netVisual_CommunFlow}
\alias{netVisual_CommunFlow}
\title{netVisual_CommunFlow}
\usage{
netVisual_CommunFlow(object, signaling, slot.name = "netP",
  pattern = c("incoming", "outgoing"), cellsize = NULL, what = "polygons",
  square = T, grid.resolution = NULL, grid.field.mean = c("median", "sum"),
  color.use = NULL, point.size = 2, image.alpha = 0.55, L = NULL, min.L = NULL,
  res = 0.5, n = NULL, jitter = 1, arrow.line.color = "grey25",
  arrow.line.width = c(0.25, 1.15), arrow.line.alpha = 0.72, arrow.angle = 18,
  arrow.length = 0.45, arrow.type = "closed", lineend = "round", ...,
  legend.size = 2, legend.text.size = 8, title.name = NULL)
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
\item{L}{Typical streamline length in coordinate units. \code{NULL} scales it to four effective grid-cell widths.}
\item{min.L}{Minimum visible streamline length. \code{NULL} uses one quarter of an effective grid-cell width.}
\item{res}{Resolution parameter for streamline integration.}
\item{n}{Streamline seeds per axis. \code{NULL} derives a bounded value from the occupied-grid count.}
\item{jitter}{Jitter of starting points.}
\item{arrow.line.color}{Color of arrows or streamlines.}
\item{arrow.line.width}{One or two positive values controlling width; a two-value range maps width to magnitude.}
\item{arrow.line.alpha}{Transparency of arrows or streamlines.}
\item{arrow.angle}{Arrowhead angle.}
\item{arrow.length}{Arrowhead length.}
\item{arrow.type}{Grid arrowhead type.}
\item{lineend}{Line ending style.}
\item{...}{Other parameters passed to the underlying metR geom.}
\item{legend.size}{Legend point size.}
\item{legend.text.size}{Legend text size.}
\item{title.name}{Plot title.}
}
\value{A ggplot2 object.}
\description{Aggregate and plot communication streamlines on a spatial grid.}
\details{Grid geometry uses the first two canonical \code{images$coordinates} axes after the existing visualization coordinate preparation. If \code{images$spatial.factors$ratio} is unavailable, coordinates remain in raw units and no physical scale is reported.

The default streamline layer checks the installed metR implementation for a zero-vector filtering defect. The affected metR 0.19.0 implementation receives a locally guarded correction; an affected, unverified implementation produces an error. The stored field and the metR namespace are not modified. A finite all-zero field renders the background without streamlines. An explicitly supplied custom stat is not replaced or bypassed.}
