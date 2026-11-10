# Merged SpatialCellChat construction for the final 11-slot schema.

.sc_merge_names <- function(object.list, add.names) {
  if (!is.list(object.list) || length(object.list) < 2L)
    stop("object.list must contain at least two SpatialCellChat objects.", call. = FALSE)
  if (is.null(add.names)) add.names <- names(object.list)
  if (is.null(add.names)) add.names <- paste0("dataset", seq_along(object.list))
  if (!is.character(add.names) || length(add.names) != length(object.list) ||
      anyNA(add.names) || any(!nzchar(add.names)) || anyDuplicated(add.names))
    stop("add.names must contain unique, non-empty names for every dataset.", call. = FALSE)
  as.character(add.names)
}

.sc_merge_assert_inputs <- function(object.list, dataset.names) {
  valid <- vapply(object.list, methods::is, logical(1), class2 = "SpatialCellChat")
  if (!all(valid))
    stop("every element of object.list must be a SpatialCellChat object.", call. = FALSE)
  for (i in seq_along(object.list)) {
    object <- object.list[[i]]
    if (!identical(object@misc$.mode, "single"))
      stop("mergeSpatialCellChat accepts only single SpatialCellChat objects; dataset '",
           dataset.names[[i]], "' is not single.", call. = FALSE)
    validateSpatialCellChat(object)
  }

  reference <- object.list[[1L]]
  datatype <- reference@misc$.datatype
  genes <- rownames(reference@assay$norm)
  groups <- levels(reference@idents)
  assay.names <- names(reference@assay)
  feature.names <- names(reference@features)
  for (i in seq_along(object.list)[-1L]) {
    object <- object.list[[i]]
    if (!identical(object@misc$.datatype, datatype))
      stop("all datasets must have the same misc$.datatype.", call. = FALSE)
    if (!identical(rownames(object@assay$norm), genes))
      stop("all datasets must have identical assay$norm gene names and order.", call. = FALSE)
    if (!identical(levels(object@idents), groups))
      stop("all datasets must have identical idents levels in the same order.", call. = FALSE)
    if (!identical(names(object@assay), assay.names))
      stop("all datasets must expose the same assay layer names.", call. = FALSE)
    if (!identical(names(object@features), feature.names) ||
        !identical(rownames(object@features), rownames(reference@features)))
      stop("all datasets must have compatible features metadata.", call. = FALSE)
    if (!identical(object@DB, reference@DB))
      stop("all datasets must use identical DB metadata for a merged comparison.", call. = FALSE)
    for (layer in assay.names) {
      left <- reference@assay[[layer]]
      right <- object@assay[[layer]]
      if (is.null(left) != is.null(right))
        stop("assay layer '", layer, "' is present in only some datasets.", call. = FALSE)
      if (!is.null(left) &&
          !identical(rownames(left), rownames(right)))
        stop("assay layer '", layer, "' must preserve identical gene order.", call. = FALSE)
    }
  }
  invisible(reference)
}

.sc_merge_cell_names <- function(object.list, dataset.names, cell.prefix) {
  result <- lapply(seq_along(object.list), function(i) {
    cells <- colnames(object.list[[i]]@assay$norm)
    if (isTRUE(cell.prefix)) paste(dataset.names[[i]], cells, sep = "_") else cells
  })
  merged <- unlist(result, use.names = FALSE)
  if (anyNA(merged) || any(!nzchar(merged)) || anyDuplicated(merged))
    stop("merged cell names must be unique and non-empty; use cell.prefix = TRUE to disambiguate datasets.",
         call. = FALSE)
  names(result) <- dataset.names
  result
}

.sc_merge_rename_dimnames <- function(value, old.cells, new.cells) {
  if (is.null(value)) return(value)
  if (inherits(value, "SparseChatArray")) {
    dn <- dimnames(value)
    if (!is.null(dn[[1L]]) && identical(dn[[1L]], old.cells)) dn[[1L]] <- new.cells
    if (!is.null(dn[[2L]]) && identical(dn[[2L]], old.cells)) dn[[2L]] <- new.cells
    return(`dimnames<-`(value, dn))
  }
  if (is.matrix(value) || inherits(value, "Matrix") || is.array(value)) {
    dn <- dimnames(value)
    if (!is.null(dn)) {
      changed <- FALSE
      for (i in seq_along(dn)) {
        if (!is.null(dn[[i]]) && identical(dn[[i]], old.cells)) {
          dn[[i]] <- new.cells
          changed <- TRUE
        }
      }
      if (changed) return(`dimnames<-`(value, dn))
    }
    return(value)
  }
  if (is.data.frame(value)) {
    if (identical(rownames(value), old.cells)) rownames(value) <- new.cells
    return(value)
  }
  if (is.list(value)) {
    out <- lapply(value, .sc_merge_rename_dimnames,
                  old.cells = old.cells, new.cells = new.cells)
    names(out) <- names(value)
    return(out)
  }
  value
}

.sc_merge_rename_images <- function(images, old.cells, new.cells) {
  images <- .sc_merge_rename_dimnames(images, old.cells, new.cells)
  if (is.list(images$rasters)) {
    for (i in seq_along(images$rasters)) {
      spot <- images$rasters[[i]]$spot.coordinates
      if (is.data.frame(spot) && identical(rownames(spot), old.cells)) {
        rownames(spot) <- new.cells
        images$rasters[[i]]$spot.coordinates <- spot
      }
    }
  }
  images
}

.sc_merge_rename_snapshot <- function(object, old.cells, new.cells) {
  net <- .sc_merge_rename_dimnames(object@net, old.cells, new.cells)
  netP <- .sc_merge_rename_dimnames(object@netP, old.cells, new.cells)
  dr <- .sc_merge_rename_dimnames(object@dr, old.cells, new.cells)
  meta <- object@meta
  rownames(meta) <- new.cells
  idents <- object@idents
  names(idents) <- new.cells
  misc <- object@misc
  misc$.mode <- "single"
  misc$.datasets <- list()
  list(
    net = net,
    netP = netP,
    dr = dr,
    meta = meta,
    idents = idents,
    images = .sc_merge_rename_images(object@images, old.cells, new.cells),
    LR = object@LR,
    DB = object@DB,
    misc = misc
  )
}

.sc_merge_assay_layer <- function(object.list, layer, cell.names, keep) {
  if (!isTRUE(keep)) return(NULL)
  values <- lapply(object.list, function(object) object@assay[[layer]])
  if (all(vapply(values, is.null, logical(1)))) return(NULL)
  do.call(cbind, Map(function(value, cells) {
    if (is.null(value)) return(NULL)
    colnames(value) <- cells
    value
  }, values, cell.names))
}

.sc_merge_metadata <- function(object.list, dataset.names, cell.names) {
  columns <- unique(unlist(lapply(object.list, function(object) names(object@meta)),
                          use.names = FALSE))
  metadata <- lapply(seq_along(object.list), function(i) {
    value <- object.list[[i]]@meta
    missing <- setdiff(columns, names(value))
    for (name in missing) value[[name]] <- NA
    value <- value[, columns, drop = FALSE]
    rownames(value) <- cell.names[[i]]
    value$datasets <- dataset.names[[i]]
    value
  })
  result <- do.call(rbind, metadata)
  result$datasets <- factor(result$datasets, levels = dataset.names)
  result
}

.sc_merge_images <- function(object.list, cell.names, datatype) {
  images <- list(coordinates = NULL, coordinate.system = NULL,
                 spatial.factors = NULL, rasters = list(), fov = list(),
                 .distance = list(), .grid = list())
  if (!identical(datatype, "spatial")) return(images)
  coordinate.systems <- lapply(object.list, function(object) object@images$coordinate.system)
  spatial.factors <- lapply(object.list, function(object) object@images$spatial.factors)
  if (!all(vapply(coordinate.systems, identical, logical(1), coordinate.systems[[1L]])) ||
      !all(vapply(spatial.factors, identical, logical(1), spatial.factors[[1L]])))
    stop("spatial datasets must have identical coordinate systems and spatial factors.",
         call. = FALSE)
  coordinates <- lapply(seq_along(object.list), function(i) {
    value <- object.list[[i]]@images$coordinates
    rownames(value) <- cell.names[[i]]
    value
  })
  images$coordinates <- do.call(rbind, coordinates)
  images$coordinate.system <- coordinate.systems[[1L]]
  images$spatial.factors <- spatial.factors[[1L]]
  images
}

#' Merge SpatialCellChat objects under the final 11-slot schema.
#'
#' @param object.list A named list of at least two single SpatialCellChat objects.
#' @param add.names Optional dataset names. When NULL, list names are used, then
#'   `dataset1`, `dataset2`, and so on.
#' @param merge.data Whether to retain raw and scale assay layers in the merged
#'   top-level assay. Norm, signaling, and smooth layers are always retained.
#' @param cell.prefix Whether to prefix cells with their dataset name.
#' @param show.plot Retained for API compatibility; merged construction never
#'   opens a plot device.
#' @return A merged SpatialCellChat object with named `misc$.datasets` snapshots.
#' @export
mergeSpatialCellChat <- function(object.list, add.names = NULL, merge.data = FALSE,
                                 cell.prefix = TRUE, show.plot = TRUE) {
  if (length(merge.data) != 1L || !is.logical(merge.data) || is.na(merge.data))
    stop("merge.data must be TRUE or FALSE.", call. = FALSE)
  if (length(cell.prefix) != 1L || !is.logical(cell.prefix) || is.na(cell.prefix))
    stop("cell.prefix must be TRUE or FALSE.", call. = FALSE)
  if (length(show.plot) != 1L || !is.logical(show.plot) || is.na(show.plot))
    stop("show.plot must be TRUE or FALSE.", call. = FALSE)

  dataset.names <- .sc_merge_names(object.list, add.names)
  .sc_merge_assert_inputs(object.list, dataset.names)
  cell.names <- .sc_merge_cell_names(object.list, dataset.names, cell.prefix)
  datatype <- object.list[[1L]]@misc$.datatype
  assay.names <- names(object.list[[1L]]@assay)
  keep <- setNames(rep(TRUE, length(assay.names)), assay.names)
  if ("raw" %in% names(keep)) keep[["raw"]] <- isTRUE(merge.data)
  if ("scale" %in% names(keep)) keep[["scale"]] <- isTRUE(merge.data)
  assay <- lapply(assay.names, function(layer)
    .sc_merge_assay_layer(object.list, layer, cell.names, keep[[layer]]))
  names(assay) <- assay.names

  snapshots <- lapply(seq_along(object.list), function(i) {
    old.cells <- colnames(object.list[[i]]@assay$norm)
    .sc_merge_rename_snapshot(object.list[[i]], old.cells, cell.names[[i]])
  })
  names(snapshots) <- dataset.names

  misc <- list(
    .mode = "merged",
    .datatype = datatype,
    .param = list(
      merge = TRUE,
      datasets = dataset.names,
      cell.prefix = isTRUE(cell.prefix),
      merge.data = isTRUE(merge.data)
    ),
    .log = list(),
    .var.features = list(),
    .datasets = snapshots
  )
  reference <- object.list[[1L]]
  lr.tables <- lapply(object.list, function(object) object@LR$LRsig)
  if (!all(vapply(lr.tables, is.data.frame, logical(1))) ||
      !all(vapply(lr.tables, function(x) identical(names(x), names(lr.tables[[1L]])), logical(1))))
    stop("all datasets must provide compatible LRsig tables with the same columns.", call. = FALSE)
  lr.union <- unique(do.call(rbind, lr.tables))
  if ("interaction_name" %in% names(lr.union)) {
    lr.union <- lr.union[!duplicated(lr.union$interaction_name), , drop = FALSE]
    rownames(lr.union) <- lr.union$interaction_name
  }
  merged.LR <- reference@LR
  merged.LR$LRsig <- lr.union
  merged.idents <- factor(
    unlist(lapply(seq_along(object.list), function(i)
      as.character(object.list[[i]]@idents)), use.names = FALSE),
    levels = levels(reference@idents), labels = levels(reference@idents)
  )
  names(merged.idents) <- unlist(cell.names, use.names = FALSE)
  merged <- methods::new(
    Class = "SpatialCellChat",
    assay = assay,
    images = .sc_merge_images(object.list, cell.names, datatype),
    meta = .sc_merge_metadata(object.list, dataset.names, cell.names),
    idents = merged.idents,
    features = reference@features,
    LR = merged.LR,
    dr = list(),
    net = list(cell = list(), group = list()),
    netP = list(cell = list(), group = list()),
    DB = reference@DB,
    misc = misc
  )
  names(merged@idents) <- colnames(merged@assay$norm)
  methods::validObject(merged)
  merged
}

# Historical analysis functions call this name when given a list. Keep it
# internal and route all construction through the new producer.
mergeCellChat <- function(object.list, add.names = NULL, merge.data = FALSE,
                          cell.prefix = TRUE, show.plot = FALSE) {
  mergeSpatialCellChat(object.list = object.list, add.names = add.names,
                       merge.data = merge.data, cell.prefix = cell.prefix,
                       show.plot = show.plot)
}

.sc_subset_indices <- function(object, cells.use, idents.use, group.by, invert) {
  cells <- colnames(object@assay$norm)
  if (!is.null(cells.use)) {
    if (is.numeric(cells.use)) {
      if (anyNA(cells.use) || any(cells.use != as.integer(cells.use)) ||
          any(cells.use < 1L | cells.use > length(cells)))
        stop("cells.use contains invalid cell indices.", call. = FALSE)
      keep <- rep(FALSE, length(cells))
      keep[unique(as.integer(cells.use))] <- TRUE
    } else if (is.character(cells.use)) {
      if (any(!cells.use %in% cells))
        stop("cells.use contains unknown cell names.", call. = FALSE)
      keep <- cells %in% unique(cells.use)
    } else {
      stop("cells.use must be NULL, numeric indices, or cell names.", call. = FALSE)
    }
  } else {
    if (is.null(idents.use))
      stop("provide either cells.use or idents.use.", call. = FALSE)
    labels <- if (is.null(group.by)) {
      as.character(object@idents)
    } else {
      if (length(group.by) != 1L || !is.character(group.by) ||
          !group.by %in% names(object@meta))
        stop("group.by must name one metadata column.", call. = FALSE)
      as.character(object@meta[[group.by]])
    }
    if (is.numeric(idents.use)) {
      levels.use <- unique(sort(as.integer(idents.use)))
      values <- unique(labels)
      if (anyNA(levels.use) || any(levels.use < 1L | levels.use > length(values)))
        stop("idents.use contains invalid group indices.", call. = FALSE)
      idents.use <- values[levels.use]
    }
    if (!is.character(idents.use) || any(!idents.use %in% unique(labels)))
      stop("idents.use contains unknown group labels.", call. = FALSE)
    keep <- labels %in% unique(idents.use)
  }
  if (isTRUE(invert)) keep <- !keep
  if (!any(keep)) stop("the requested subset contains no cells.", call. = FALSE)
  which(keep)
}

.sc_subset_value_by_cells <- function(value, old.cells, keep) {
  if (is.null(value)) return(value)
  new.cells <- old.cells[keep]
  if (inherits(value, "SparseChatArray")) {
    dn <- dimnames(value)
    row.keep <- !is.null(dn[[1L]]) && identical(dn[[1L]], old.cells)
    col.keep <- !is.null(dn[[2L]]) && identical(dn[[2L]], old.cells)
    if (row.keep && col.keep) return(value[keep, keep, , drop = FALSE])
    if (row.keep) return(value[keep, , , drop = FALSE])
    if (col.keep) return(value[, keep, , drop = FALSE])
    return(value)
  }
  if (is.matrix(value) || inherits(value, "Matrix") || is.array(value)) {
    dn <- dimnames(value)
    if (is.null(dn)) return(value)
    out <- value
    for (axis in seq_along(dn)) {
      if (!identical(dn[[axis]], old.cells)) next
      index <- rep(list(TRUE), length(dim(out)))
      index[[axis]] <- keep
      out <- do.call(`[`, c(list(out), index, list(drop = FALSE)))
    }
    return(out)
  }
  if (is.data.frame(value) && identical(rownames(value), old.cells))
    return(value[keep, , drop = FALSE])
  if (is.list(value)) {
    out <- lapply(value, .sc_subset_value_by_cells, old.cells = old.cells,
                  keep = keep)
    names(out) <- names(value)
    return(out)
  }
  value
}

.sc_subset_images <- function(images, old.cells, keep) {
  if (!is.null(images$coordinates))
    images$coordinates <- images$coordinates[keep, , drop = FALSE]
  images$.distance <- list()
  images$.grid <- list()
  if (is.list(images$rasters)) {
    for (i in seq_along(images$rasters)) {
      spot <- images$rasters[[i]]$spot.coordinates
      if (is.data.frame(spot) && !is.null(rownames(spot)) &&
          all(rownames(spot) %in% old.cells)) {
        images$rasters[[i]]$spot.coordinates <-
          spot[match(old.cells[keep], rownames(spot)), , drop = FALSE]
      }
    }
  }
  images
}

.sc_subset_snapshot <- function(snapshot, keep) {
  old.cells <- rownames(snapshot$meta)
  new.cells <- old.cells[keep]
  snapshot$net <- list(cell = list(), group = list())
  snapshot$netP <- list(cell = list(), group = list())
  snapshot$meta <- snapshot$meta[keep, , drop = FALSE]
  snapshot$idents <- droplevels(snapshot$idents[keep])
  names(snapshot$idents) <- new.cells
  snapshot$dr <- .sc_subset_value_by_cells(snapshot$dr, old.cells, keep)
  snapshot$images <- .sc_subset_images(snapshot$images, old.cells, keep)
  snapshot
}

#' Subset a SpatialCellChat object by cells or group labels.
#'
#' Communication results are invalidated because filtering changes group sizes,
#' group membership, and spatial distance relationships. Recompute communication
#' after subsetting instead of silently reusing stale aggregates.
#' @param object A SpatialCellChat object.
#' @param cells.use Cell names or numeric cell indices.
#' @param idents.use Group labels or numeric group indices.
#' @param group.by Metadata column used to resolve idents.use.
#' @param invert Whether to keep the complement of the requested selection.
#' @param thresh Retained for API compatibility; it does not gate subsetting.
#' @return A valid SpatialCellChat object with communication results cleared.
#' @export
subsetSpatialCellChat <- function(object, cells.use = NULL, idents.use = NULL,
                                  group.by = NULL, invert = FALSE, thresh = 0.05) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (length(invert) != 1L || !is.logical(invert) || is.na(invert))
    stop("invert must be TRUE or FALSE.", call. = FALSE)
  if (length(thresh) != 1L || !is.numeric(thresh) || !is.finite(thresh))
    stop("thresh must be one finite numeric value.", call. = FALSE)
  keep <- .sc_subset_indices(object, cells.use, idents.use, group.by, invert)
  old.cells <- colnames(object@assay$norm)
  object@assay <- lapply(object@assay, function(value) {
    if (is.null(value)) return(NULL)
    value[, keep, drop = FALSE]
  })
  object@meta <- object@meta[keep, , drop = FALSE]
  object@idents <- droplevels(object@idents[keep])
  names(object@idents) <- old.cells[keep]
  object@dr <- .sc_subset_value_by_cells(object@dr, old.cells, keep)
  object@images <- .sc_subset_images(object@images, old.cells, keep)
  object@net <- list(cell = list(), group = list())
  object@netP <- list(cell = list(), group = list())
  if (identical(object@misc$.mode, "merged")) {
    datasets <- object@misc$.datasets
    dataset.labels <- as.character(object@meta$datasets)
    kept.names <- names(datasets)[vapply(names(datasets), function(name) {
      any(dataset.labels == name)
    }, logical(1))]
    if (length(kept.names) < 2L)
      stop("a merged subset must retain at least two datasets.", call. = FALSE)
    object@misc$.datasets <- lapply(kept.names, function(name) {
      snapshot <- datasets[[name]]
      snapshot.keep <- rownames(snapshot$meta) %in% rownames(object@meta)
      .sc_subset_snapshot(snapshot, snapshot.keep)
    })
    names(object@misc$.datasets) <- kept.names
    object@meta$datasets <- factor(dataset.labels, levels = kept.names)
  } else {
    object@misc$.datasets <- list()
  }
  .sc_validate_after_update(object)
}

#' Update a legacy SpatialCellChat object to the final schema.
#'
#' The current implementation accepts only objects already using the final
#' 11-slot class and validates them. Legacy S4 layouts are rejected explicitly
#' because guessing old slot meanings can silently corrupt communication data.
#' @param object A SpatialCellChat object.
#' @return A validated SpatialCellChat object.
#' @export
updateSpatialCellChat <- function(object) {
  if (!methods::is(object, "SpatialCellChat"))
    stop("legacy object conversion is not supported without an explicit schema; provide a SpatialCellChat object created by the current constructor.", call. = FALSE)
  .sc_validate_after_update(object)
}

#' Lift group labels to a common set across datasets.
#'
#' Use a named character vector as an old-to-new label map. An unnamed vector is
#' accepted only to reorder the existing levels without changing any labels.
#' Communication results are invalidated because group axes have changed.
#' @param object A single or merged SpatialCellChat object.
#' @param group.new A named old-to-new label map, or the existing group labels
#'   in a new order when no labels need to change.
#' @return A valid SpatialCellChat object with communication results cleared.
#' @export
liftSpatialCellChat <- function(object, group.new = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(group.new))
    stop("group.new must be supplied explicitly for the final schema.", call. = FALSE)
  if (!is.character(group.new) || !length(group.new) || anyNA(group.new) ||
      any(!nzchar(group.new)))
    stop("group.new must contain non-empty character labels.", call. = FALSE)
  old.levels <- levels(object@idents)
  if (is.null(names(group.new))) {
    if (anyDuplicated(group.new) || !setequal(group.new, old.levels))
      stop("unnamed group.new may only reorder all existing group labels; use a named old-to-new map to relabel.", call. = FALSE)
    group.levels <- group.new
    lift <- function(value) as.character(value)
  } else {
    if (anyNA(names(group.new)) || any(!nzchar(names(group.new))) ||
        anyDuplicated(names(group.new)) || !setequal(names(group.new), old.levels))
      stop("names(group.new) must map every existing group label exactly once.", call. = FALSE)
    group.levels <- unique(unname(group.new))
    lift <- function(value) unname(group.new[as.character(value)])
  }
  object@idents <- factor(lift(object@idents), levels = group.levels)
  names(object@idents) <- colnames(object@assay$norm)
  object@net <- list(cell = list(), group = list())
  object@netP <- list(cell = list(), group = list())
  if (identical(object@misc$.mode, "merged")) {
    for (name in names(object@misc$.datasets)) {
      snapshot <- object@misc$.datasets[[name]]
      snapshot$idents <- factor(lift(snapshot$idents), levels = group.levels)
      names(snapshot$idents) <- rownames(snapshot$meta)
      snapshot$net <- list(cell = list(), group = list())
      snapshot$netP <- list(cell = list(), group = list())
      object@misc$.datasets[[name]] <- snapshot
    }
  }
  .sc_validate_after_update(object)
}
