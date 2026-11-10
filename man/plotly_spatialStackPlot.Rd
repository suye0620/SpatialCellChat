% Generated manually for the public Plotly stack entry point.
\name{plotly_spatialStackPlot}
\alias{plotly_spatialStackPlot}
\title{Stack spatial feature layers in a Plotly 3D view}
\usage{
plotly_spatialStackPlot(
  object,
  layers,
  color.heatmap = "RdBu",
  point.size = 3,
  z.axis.space = 1,
  title.name = NULL,
  colorbar.title = "value"
)
}
\arguments{
  \item{object}{A single-dataset \code{SpatialCellChat} object with spatial coordinates.}
  \item{layers}{A numeric matrix or data frame with cells in rows and named layers in columns. Row names must contain exactly the coordinate cell names; rows are aligned by name.}
  \item{color.heatmap}{A valid RColorBrewer palette name or a vector of colors.}
  \item{point.size}{A positive finite marker size.}
  \item{z.axis.space}{A positive finite distance between adjacent layer planes.}
  \item{title.name}{Optional plot title.}
  \item{colorbar.title}{Colorbar title.}
}
\value{
A built Plotly widget with one \code{scatter3d} trace and one shared colorbar.
}
\description{
Render numeric cell-by-layer values as horizontal spatial slices. The z-axis encodes layer order only; it is not a physical tissue-depth measurement.
}
\details{
The input rows are aligned by cell name. All layers share one global color domain and one colorbar. Missing values are omitted from the rendered trace; each layer must still contain at least one finite value.
}
\keyword{visualization}
