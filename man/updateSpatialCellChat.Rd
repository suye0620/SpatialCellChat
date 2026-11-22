% Generated from the updateSpatialCellChat documentation in R/merging.R.
\name{updateSpatialCellChat}
\alias{updateSpatialCellChat}
\title{Migrate an upstream SpatialCellChat object to the current schema}
\usage{
updateSpatialCellChat(object)
}
\arguments{
\item{object}{An upstream 15-slot \code{SpatialCellChat}, historical 15-slot
\code{SpatialChat}, or current 11-slot \code{SpatialCellChat} object.}
}
\value{
A validated 11-slot SpatialCellChat object from SpatialCellChatDev.
}
\description{
Converts an upstream 15-slot object or a historical \code{SpatialChat} object
with the same stored layout to the current 11-slot layout. A current object is
validated and returned unchanged. Unknown layouts are rejected.
}
\details{
For a single dataset, normalized, raw, scaled, projected, and signaling
expression layers are mapped to the current assay layers. Cell metadata,
identities, database, ligand-receptor metadata, dimensional reductions,
communication probabilities, p-values, counts, weights, pathway results,
centrality, and auxiliary network records are retained when their axes can be
aligned. Three-dimensional probability and p-value arrays are converted to
sparse matrix layers. Spatial coordinates and calibration are migrated.

Legacy distance caches lack the identity metadata required by the current
schema and are not reused; the original image records and options are retained
in \code{misc$.legacy}. The communication distance plot requires an inference
cache with provenance, so rerun \code{computeCommunProb()} to recreate it.

A legacy merged object can be converted only when its \code{@data} contains
normalized expression, as with the upstream \code{merge.data = TRUE} option.
If it does not, migrate the original single-dataset objects and merge them
again. Named axes with inconsistent cell, gene, or group identifiers are
rejected. Where the upstream constructor left an axis unnamed, migration
assumes its original order matches the expression data or group levels.

The upstream and development packages use the same S4 class name. When both are
loaded, use \code{SpatialCellChatDev::updateSpatialCellChat(object)} explicitly.
}
