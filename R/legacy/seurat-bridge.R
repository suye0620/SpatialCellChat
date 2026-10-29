# Archived: the two pre-consolidation Seurat bridges
#
# Archived on 2026-10-29, one day after asSeurat() replaced them. They are no
# longer part of the package: R only compiles top-level .R files under R/, so
# nothing in this directory is sourced by R CMD INSTALL. They are kept here as
# the record of the old argument vocabulary.
#
# Why they were retired:
#   - both are thin wrappers over asSeurat(); keeping them exported meant two
#     names and three axes (slot.name x pattern x merge.chat) for one operation;
#   - the package itself no longer calls them (identifyCellTopics and
#     identifyTopicCommunities call asSeurat() directly);
#   - mergeChatIntoSeurat(merge.chat = FALSE) returned "Spatial = incoming" plus
#     incoming/outgoing aliases, which asSeurat() deliberately narrowed to the
#     single requested direction.
#
# Argument mapping:
#   CellChat2Seurat(object, counts, group.by, assay, check.object)
#     -> asSeurat(object, slot.name = "assay", counts = counts,
#                 assay.name = assay, group.by = group.by,
#                 check.object = check.object)
#   mergeChatIntoSeurat(object.list, slot.name = "net", pattern = "incoming",
#                       normalize.chat = TRUE, normalization.method, merge.chat)
#     -> asSeurat(object.list, slot.name = slot.name,
#                 what = "merge" if merge.chat else pattern,
#                 normalize.chat = ..., normalization.method = ...)
#
# Original locations at archive time:
#   R/utilities.R  383  mergeChatIntoSeurat
#   R/utilities.R  667  CellChat2Seurat
#
# To restore, source this file after R/utilities.R (both call asSeurat(), which
# still exists) and re-add the two export() lines to NAMESPACE.

#' Merge SpatialCellChat objects (deprecated)
#'
#' Deprecated: use [asSeurat()] with `slot.name = "net"` or `"netP"` and
#' `what = "incoming"`, `"outgoing"` or `"merge"`.
#'
#' @inheritParams asSeurat
#' @param object.list A SpatialCellChat object or a named list of them.
#' @param pattern Which centrality direction to write.
#' @param normalize.chat Whether to run [Seurat::NormalizeData()].
#' @param normalization.method Passed to [Seurat::NormalizeData()].
#' @param merge.chat Combine both directions into one assay (`what = "merge"`).
#'
#' @return A Seurat object.
mergeChatIntoSeurat <- function(
    object.list,
    slot.name = "net",
    pattern = c("incoming","outgoing"),
    normalize.chat = TRUE,
    normalization.method = "LogNormalize",
    merge.chat = FALSE
){
  .Deprecated("asSeurat")
  slot.name <- match.arg(as.character(slot.name), c("net", "netP"))
  pattern <- match.arg(pattern)
  asSeurat(object.list,
           slot.name = slot.name,
           what = if (isTRUE(merge.chat)) "merge" else pattern,
           normalize.chat = normalize.chat,
           normalization.method = normalization.method,
           check.object = FALSE)
}

#' Convert a CellChat object into a Seurat object (deprecated)
#'
#' Deprecated: use [asSeurat()] with `slot.name = "assay"`, which also carries
#' `assay$raw`, `assay$scale`, `assay$signaling` and `assay$smooth` over.
#'
#' @inheritParams asSeurat
#' @param assay Name of the assay in the returned Seurat object.
#'
#' @return a Seurat object
CellChat2Seurat <- function(object,counts = NULL, group.by = NULL, assay = "Spatial", check.object = TRUE) {
  .Deprecated("asSeurat")
  asSeurat(object, slot.name = "assay", counts = counts, assay.name = assay,
           group.by = group.by, check.object = check.object)
}
