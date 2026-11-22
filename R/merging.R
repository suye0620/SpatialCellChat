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
                 .distance = list(), .grid = list(),
                 coordinate.scope = if (identical(datatype, "spatial"))
                   "dataset-local" else "shared")
  images
}

#' Merge SpatialCellChat objects under the final 11-slot schema.
#'
#' @param object.list A list of at least two single SpatialCellChat objects.
#' @param add.names Optional dataset names. When NULL, list names are used, then
#'   `dataset1`, `dataset2`, and so on.
#' @param merge.data Whether to retain raw and scale assay layers in the merged
#'   top-level assay. Norm, signaling, and smooth layers are always retained.
#' @param cell.prefix Whether to prefix cells with their dataset name. Defaults
#'   to TRUE to avoid cell-name collisions.
#' @param show.plot Retained for API compatibility; merged construction never
#'   opens a plot device. Defaults to TRUE for compatibility.
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

.sc_update_schema <- function(object) {
  if (!isS4(object))
    stop("object must be a SpatialCellChat S4 object.", call. = FALSE)
  class.name <- unname(class(object))[1L]
  # slotNames() can resolve a different same-named class definition from the
  # S4 cache. The object's own attributes are the persisted slot values.
  slots <- setdiff(names(attributes(object)), "class")
  legacy <- c("data.raw", "data", "data.signaling", "data.scale", "data.project",
              "images", "net", "netP", "meta", "idents", "DB", "LR",
              "var.features", "dr", "options")
  current <- c("assay", "images", "meta", "idents", "features", "LR", "dr",
               "net", "netP", "DB", "misc")
  if (length(slots) == length(legacy) && setequal(slots, legacy) &&
      (identical(class.name, "SpatialCellChat") ||
       (identical(class.name, "SpatialChat") &&
          identical(attr(attr(object, "class"), "package"), "SpatialChat"))))
    return("legacy")
  if (length(slots) == length(current) && setequal(slots, current) &&
      identical(class.name, "SpatialCellChat"))
    return("current")
  stop("Unrecognized SpatialCellChat slot layout; cannot safely migrate this object.",
       call. = FALSE)
}

.sc_update_layer <- function(value, norm, label, subset = FALSE) {
  if (!.sc_matrix_like(value)) stop(label, " must be a matrix.", call. = FALSE)
  if (!length(value) && (all(dim(value) == 0L) || (subset && nrow(value) == 0L))) return(NULL)
  genes <- rownames(norm); cells <- colnames(norm)
  if (is.null(rownames(value)) || is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) || anyDuplicated(colnames(value)) ||
      !setequal(colnames(value), cells))
    stop(label, " must have unique gene names and the same named cells as data.", call. = FALSE)
  if (subset) {
    if (!all(rownames(value) %in% genes))
      stop(label, " contains genes absent from data.", call. = FALSE)
    gene.order <- genes[genes %in% rownames(value)]
  } else {
    if (!setequal(rownames(value), genes))
      stop(label, " genes do not match data.", call. = FALSE)
    gene.order <- genes
  }
  value[gene.order, cells, drop = FALSE]
}

.sc_update_sparse_array <- function(value, nodes, label, second.nodes = nodes) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "SparseChatArray")) {
    dims <- dim(value); dn <- dimnames(value)
    layers <- unclass(value)
  } else if (inherits(value, "sparse3Darray")) {
    dims <- dim(value); dn <- dimnames(value)
    if (!dims[3L]) return(NULL)
    # The upstream `[` method scans every stored triplet for each layer. Split
    # indices once so a large spatial network is not scanned once per LR.
    stored <- unclass(value)
    if (length(stored$i) != length(stored$k) || length(stored$j) != length(stored$k) ||
        length(stored$x) != length(stored$k) || !is.numeric(stored$k) ||
        anyNA(stored$k) || any(!is.finite(stored$k)) ||
        any(stored$k < 1L | stored$k > dims[3L] | stored$k != floor(stored$k)))
      stop(label, " has invalid sparse triplet indices.", call. = FALSE)
    idx <- split(seq_along(stored$k), factor(stored$k, levels = seq_len(dims[3L])))
    layers <- lapply(idx, function(k) Matrix::sparseMatrix(
      i = stored$i[k], j = stored$j[k], x = stored$x[k], dims = dims[1:2]))
  } else if (is.array(value) && length(dim(value)) == 3L) {
    dims <- dim(value); dn <- dimnames(value)
    if (!dims[3L]) return(NULL)
    layers <- lapply(seq_len(dims[3L]), function(i) value[, , i, drop = TRUE])
  } else if (is.list(value) && length(value) &&
             all(vapply(value, .sc_matrix_like, logical(1)))) {
    dims <- c(dim(value[[1L]]), length(value))
    first.dn <- dimnames(value[[1L]]) %||% list(NULL, NULL)
    dn <- c(first.dn, list(names(value)))
    for (i in seq_along(value)) {
      other <- dimnames(value[[i]])
      for (j in seq_len(2L)) {
        if (!identical(other[[j]], dn[[j]]))
          stop(label, " contains layers with inconsistent dimnames.", call. = FALSE)
      }
    }
    layers <- value
  } else {
    stop(label, " must be a three-dimensional array or list of matrices.", call. = FALSE)
  }
  if (length(dims) != 3L || !identical(as.integer(dims[1:2]),
                                         as.integer(c(length(nodes), length(second.nodes)))))
    stop(label, " dimensions do not match its node axes.", call. = FALSE)
  if (is.null(dn)) dn <- list(NULL, NULL, NULL)
  if (length(dn) != 3L) stop(label, " has invalid dimnames.", call. = FALSE)
  expected <- list(nodes, second.nodes)
  indices <- lapply(seq_len(2L), function(i) {
    if (is.null(dn[[i]])) return(seq_along(expected[[i]]))
    if (anyDuplicated(dn[[i]]) || !setequal(dn[[i]], expected[[i]]))
      stop(label, " dimnames do not match its node axes.", call. = FALSE)
    match(expected[[i]], dn[[i]])
  })
  layers <- lapply(layers, function(layer) {
    if (!.sc_matrix_like(layer)) layer <- matrix(layer, nrow = dims[1L], ncol = dims[2L])
    if (!identical(as.integer(dim(layer)), as.integer(dims[1:2])))
      stop(label, " contains a layer with incompatible dimensions.", call. = FALSE)
    methods::as(methods::as(layer[indices[[1L]], indices[[2L]], drop = FALSE], "generalMatrix"), "CsparseMatrix")
  })
  layers <- lapply(layers, function(x) methods::as(x, "dMatrix"))
  SparseChatArray(layers, dimnames = list(nodes, second.nodes, dn[[3L]]))
}

.sc_update_square_matrix <- function(value, nodes, label) {
  if (!.sc_matrix_like(value) || !identical(as.integer(dim(value)), rep.int(length(nodes), 2L)))
    stop(label, " must be a square matrix aligned with its node axes.", call. = FALSE)
  dn <- dimnames(value)
  idx <- lapply(seq_len(2L), function(i) {
    if (is.null(dn[[i]])) return(seq_along(nodes))
    if (anyDuplicated(dn[[i]]) || !setequal(dn[[i]], nodes))
      stop(label, " dimnames do not match its node axes.", call. = FALSE)
    match(nodes, dn[[i]])
  })
  value <- value[idx[[1L]], idx[[2L]], drop = FALSE]
  dimnames(value) <- list(nodes, nodes)
  value
}

.sc_update_network <- function(old, cells, groups, label) {
  if (!is.list(old)) stop(label, " must be a list.", call. = FALSE)
  out <- list(cell = list(), group = list())
  map <- c("prob.cell" = "prob", "pval.cell" = "pval", "count.cell" = "count",
           "weight.cell" = "weight", "sum.cell" = "weight", "LR.sig.cell" = "LR.sig",
           "prob" = "prob", "pval" = "pval", "count" = "count",
           "weight" = "weight", "sum" = "weight", "LR.sig" = "LR.sig")
  for (name in names(map)) {
    value <- old[[name]]
    if (is.null(value)) next
    cell <- endsWith(name, ".cell")
    resolution <- if (cell) "cell" else "group"
    measure <- unname(map[[name]])
    if (!is.null(out[[resolution]][[measure]])) {
      if (!identical(out[[resolution]][[measure]], value))
        warning(label, " has both weight and sum; weight takes precedence; original sum is retained in misc$.legacy.", call. = FALSE)
      next
    }
    nodes <- if (cell) cells else groups
    if (measure %in% c("prob", "pval")) value <- .sc_update_sparse_array(value, nodes, paste0(label, "$", name))
    if (measure %in% c("count", "weight")) value <- .sc_update_square_matrix(value, nodes, paste0(label, "$", name))
    out[[resolution]][[measure]] <- value
  }
  for (resolution in c("cell", "group")) {
    prob <- out[[resolution]]$prob
    pval <- out[[resolution]]$pval
    if (!is.null(prob) && !is.null(pval) &&
        (!identical(dim(prob), dim(pval)) || !identical(dimnames(prob), dimnames(pval))))
      stop(label, "$", resolution, " probability and p-value axes do not match.", call. = FALSE)
  }
  keep <- setdiff(names(old), c(names(map), "field"))
  out[keep] <- old[keep]
  if (!is.null(old$group.merged)) {
    mapping <- old$group.merged
    if (!is.factor(mapping)) mapping <- factor(mapping)
    if (length(mapping) != length(groups) || anyNA(mapping))
      stop(label, "$group.merged must map each group.", call. = FALSE)
    if (!is.null(names(mapping))) {
      if (anyDuplicated(names(mapping)) || !setequal(names(mapping), groups))
        stop(label, "$group.merged names must match group names.", call. = FALSE)
      mapping <- mapping[groups]
    }
    names(mapping) <- groups
    out$group$group.merged <- mapping
    out$group.merged <- NULL
    for (measure in c("count.merged", "weight.merged")) {
      if (!is.null(old[[measure]])) {
        out$group[[measure]] <- .sc_update_square_matrix(
          old[[measure]], levels(mapping), paste0(label, "$", measure))
        out[[measure]] <- NULL
      }
    }
  }
  if (!is.null(old$field)) {
    if (!is.list(old$field)) stop(label, "$field must be a list.", call. = FALSE)
    if (length(old$field) && (is.null(names(old$field)) || anyNA(names(old$field)) ||
                              any(!nzchar(names(old$field))) || anyDuplicated(names(old$field))))
      stop(label, "$field must have unique names for its components.", call. = FALSE)
    out$cell$field <- lapply(names(old$field), function(name) {
      value <- old$field[[name]]
      dims <- dim(value)
      if (is.null(dims) || length(dims) != 3L || !dims[2L] %in% c(2L, 3L))
        stop(label, "$field$", name, " has an unsupported layout.", call. = FALSE)
      dn <- dimnames(value)
      axes <- if (!is.null(dn) && !is.null(dn[[2L]])) dn[[2L]] else c("x", "y", "z")[seq_len(dims[2L])]
      .sc_update_sparse_array(value, cells, paste0(label, "$field$", name), axes)
    })
    names(out$cell$field) <- names(old$field)
  }
  out
}

.sc_update_images <- function(old, datatype, cells) {
  images <- .sc_default_images()
  if (!identical(datatype, "spatial")) return(images)
  coords <- old$coordinates
  if (is.null(coords)) stop("Legacy spatial object has no coordinates.", call. = FALSE)
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || nrow(coords) != length(cells) ||
      !ncol(coords) %in% c(2L, 3L) || any(!is.finite(coords)))
    stop("Legacy coordinates must be a finite cells-by-2 or cells-by-3 matrix.", call. = FALSE)
  rn <- rownames(coords)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(coords))))) {
    rownames(coords) <- cells
  } else {
    if (anyDuplicated(rn) || !setequal(rn, cells))
      stop("Legacy coordinate rownames do not match the expression cells.", call. = FALSE)
    coords <- coords[cells, , drop = FALSE]
  }
  colnames(coords) <- c("x", "y", "z")[seq_len(ncol(coords))]
  factors <- old$spatial.factors %||% old$scale.factors %||% list()
  if (is.data.frame(factors)) factors <- as.list(factors)
  if (!is.list(factors)) stop("Legacy spatial factors must be a list or data frame.", call. = FALSE)
  system <- .sc_normalize_coordinate_system(NULL, factors)
  factors <- .sc_normalize_spatial_factors(factors, system)
  if (!"ratio" %in% names(factors)) factors["ratio"] <- list(NULL)
  if (!"tol" %in% names(factors)) factors["tol"] <- list(NULL)
  images$coordinates <- coords
  images$coordinate.system <- system
  images$spatial.factors <- factors
  images
}

.sc_update_merged <- function(old, datatype) {
  norm <- old$data
  if (!.sc_matrix_like(norm) || !length(norm) || is.null(rownames(norm)) ||
      is.null(colnames(norm)) || anyDuplicated(rownames(norm)) || anyDuplicated(colnames(norm)))
    stop("Legacy merged object has no complete named normalized expression matrix in @data. Upstream merge.data = FALSE does not save it; migrate the original single-dataset objects and merge them again.", call. = FALSE)
  cells <- colnames(norm); genes <- rownames(norm)
  meta <- old$meta
  if (!is.data.frame(meta) || !identical(rownames(meta), cells) ||
      !is.factor(meta$datasets) || anyNA(meta$datasets))
    stop("Legacy merged meta must align to @data and have a factor datasets column.", call. = FALSE)
  datasets <- levels(meta$datasets)
  if (length(datasets) < 2L || any(!nzchar(datasets)) || anyDuplicated(datasets) ||
      any(!datasets %in% as.character(meta$datasets)))
    stop("Legacy merged object must contain at least two named nonempty datasets.", call. = FALSE)
  if (!is.list(old$idents) || !is.factor(old$idents$joint) ||
      length(old$idents$joint) != length(cells) || anyNA(old$idents$joint))
    stop("Legacy merged idents$joint must contain one factor label per cell.", call. = FALSE)
  joint <- old$idents$joint
  if (!is.null(names(joint))) {
    if (anyDuplicated(names(joint)) || !setequal(names(joint), cells))
      stop("Legacy merged idents$joint names do not match expression cells.", call. = FALSE)
    joint <- joint[cells]
  }
  names(joint) <- cells
  snapshots <- lapply(datasets, function(dataset) {
    selected <- cells[as.character(meta$datasets) == dataset]
    identity <- old$idents[[dataset]]
    if (!is.factor(identity) || length(identity) != length(selected) || anyNA(identity))
      stop("Legacy dataset '", dataset, "' has invalid identities.", call. = FALSE)
    image <- old$images[[dataset]] %||% list()
    net <- old$net[[dataset]] %||% list()
    netP <- old$netP[[dataset]] %||% list()
    candidates <- list(names(identity),
                       if (!is.null(net$prob.cell)) dimnames(net$prob.cell)[[1L]] else NULL,
                       if (!is.null(netP$prob.cell)) dimnames(netP$prob.cell)[[1L]] else NULL,
                       rownames(image$coordinates))
    candidates <- Filter(function(x) !is.null(x) &&
                           !identical(x, as.character(seq_along(x))), candidates)
    original <- if (length(candidates)) candidates[[1L]] else selected
    if (length(original) != length(selected) || anyNA(original) || anyDuplicated(original))
      stop("Cannot determine legacy cell names for dataset '", dataset, "'.", call. = FALSE)
    prefixed <- paste0(dataset, "_", original)
    if (setequal(original, selected)) original <- original[match(selected, original)]
    else if (setequal(prefixed, selected)) original <- original[match(selected, prefixed)]
    else stop("Cannot align legacy dataset '", dataset, "' cell names to the merged object.", call. = FALSE)
    if (!is.null(names(identity))) identity <- identity[original]
    names(identity) <- original
    if (anyNA(identity) || !identical(as.character(identity), as.character(joint[selected])))
      stop("Legacy identities disagree for dataset '", dataset, "'.", call. = FALSE)
    lr <- old$LR[[dataset]] %||% list()
    if (is.null(lr$LRsig)) lr$LRsig <- data.frame()
    if (!is.data.frame(lr$LRsig)) stop("Legacy LRsig must be a data.frame.", call. = FALSE)
    vf <- old$var.features[[dataset]] %||% list()
    if (is.character(vf)) vf <- list(features = vf)
    local.misc <- list(.mode = "single", .datatype = datatype,
                       .param = list(sample = dataset), .log = list(),
                       .var.features = vf, .datasets = list(),
                       .legacy = list(images = image, net.sum = net$sum,
                                      net.sum.cell = net$sum.cell,
                                      netP.sum = netP$sum,
                                      netP.sum.cell = netP$sum.cell))
    local.meta <- meta[selected, , drop = FALSE]
    rownames(local.meta) <- original
    snapshot <- list(net = .sc_update_network(net, original, levels(identity), paste0("net$", dataset)),
                     netP = .sc_update_network(netP, original, levels(identity), paste0("netP$", dataset)),
                     dr = list(), meta = local.meta, idents = identity,
                     images = .sc_update_images(image, datatype, original),
                     LR = lr, DB = old$DB, misc = local.misc)
    snapshot$net <- .sc_merge_rename_dimnames(snapshot$net, original, selected)
    snapshot$netP <- .sc_merge_rename_dimnames(snapshot$netP, original, selected)
    snapshot$images <- .sc_merge_rename_images(snapshot$images, original, selected)
    rownames(snapshot$meta) <- selected
    names(snapshot$idents) <- selected
    snapshot
  })
  names(snapshots) <- datasets
  lr.tables <- lapply(snapshots, function(x) x$LR$LRsig)
  nonempty <- Filter(function(x) nrow(x) > 0L, lr.tables)
  if (length(nonempty) && !all(vapply(nonempty, function(x) identical(names(x), names(nonempty[[1L]])), logical(1))))
    stop("Legacy dataset LRsig tables have incompatible columns; migrate the original datasets separately.", call. = FALSE)
  lr.union <- if (length(nonempty)) unique(do.call(rbind, nonempty)) else data.frame()
  if ("interaction_name" %in% names(lr.union))
    lr.union <- lr.union[!duplicated(lr.union$interaction_name), , drop = FALSE]
  images <- .sc_default_images()
  images$coordinate.scope <- if (identical(datatype, "spatial")) "dataset-local" else "shared"
  misc <- list(.mode = "merged", .datatype = datatype,
               .param = list(legacy = old$options), .log = list(),
               .var.features = old$var.features, .datasets = snapshots,
               .legacy = list(images = old$images, options = old$options))
  result <- SpatialCellChat(
    assay = list(raw = .sc_update_layer(old$data.raw, norm, "data.raw"), norm = norm,
                 scale = .sc_update_layer(old$data.scale, norm, "data.scale"),
                 smooth = .sc_update_layer(old$data.project, norm, "data.project", TRUE),
                 signaling = .sc_update_layer(old$data.signaling, norm, "data.signaling", TRUE)),
    images = images, meta = meta, idents = joint, features = data.frame(row.names = genes),
    LR = list(LRsig = lr.union), dr = old$dr,
    net = list(cell = list(), group = list()), netP = list(cell = list(), group = list()),
    DB = old$DB, misc = misc)
  validateSpatialCellChat(result)
  .sc_validate_after_update(result)
}

#' Update an upstream SpatialCellChat object to the current schema
#'
#' Converts single-dataset 15-slot objects from the upstream SpatialCellChat
#' package and historical `SpatialChat` objects with the same stored layout to
#' the development package's 11-slot schema. Expression, metadata,
#' identities, LR, database, dimensional reductions and communication results
#' are retained. Legacy probability arrays are converted layer by layer to
#' SparseChatArray. Unverified distance caches are not reused and the original
#' image records and options are retained in `misc$.legacy` for inspection.
#' Rerun `computeCommunProb()` if an inference distance cache is needed.
#' Already-current objects are validated and returned unchanged. Legacy merged
#' objects require the normalized expression matrix in `data`; if upstream
#' merge.data was FALSE, migrate and merge the original datasets instead.
#' Named axes that disagree with the expression or group identities cause an
#' error. When an upstream axis has no names, its original order is assumed to
#' match the expression data or group levels. Both packages use the same S4
#' class name, so call `SpatialCellChatDev::updateSpatialCellChat(object)`
#' explicitly when they are both loaded.
#' @param object An upstream 15-slot `SpatialCellChat`, historical 15-slot
#'   `SpatialChat`, or current 11-slot `SpatialCellChat` object.
#' @return A validated 11-slot SpatialCellChat object.
#' @export
updateSpatialCellChat <- function(object) {
  schema <- .sc_update_schema(object)
  if (identical(schema, "current")) {
    validateSpatialCellChat(object)
    return(.sc_validate_after_update(object))
  }
  slots <- setdiff(names(attributes(object)), "class")
  old <- attributes(object)[slots]
  options <- old$options
  mode <- options$mode %||% "single"
  datatype <- options$datatype %||% "RNA"
  if (!datatype %in% c("RNA", "spatial")) stop("Unsupported legacy datatype: ", datatype, call. = FALSE)
  if (identical(mode, "merged")) return(.sc_update_merged(old, datatype))
  if (!identical(mode, "single")) stop("Unsupported legacy mode: ", mode, call. = FALSE)
  norm <- old$data
  if (!.sc_matrix_like(norm) || !length(norm) || is.null(rownames(norm)) ||
      is.null(colnames(norm)) || anyDuplicated(rownames(norm)) || anyDuplicated(colnames(norm)))
    stop("Legacy data must contain a nonempty normalized expression matrix with unique gene and cell names.", call. = FALSE)
  cells <- colnames(norm); genes <- rownames(norm)
  meta <- old$meta
  if (!is.data.frame(meta) || nrow(meta) != length(cells) || anyDuplicated(rownames(meta)) ||
      !setequal(rownames(meta), cells))
    stop("Legacy meta rownames must match the expression cell names.", call. = FALSE)
  meta <- meta[cells, , drop = FALSE]
  ident <- old$idents
  if (!is.factor(ident) || length(ident) != length(cells) || anyNA(ident))
    stop("Legacy idents must contain one non-missing factor label per cell.", call. = FALSE)
  if (!is.null(names(ident))) {
    if (anyDuplicated(names(ident)) || !setequal(names(ident), cells))
      stop("Legacy idents names must match the expression cells.", call. = FALSE)
    ident <- ident[cells]
  }
  names(ident) <- cells
  assay <- list(raw = .sc_update_layer(old$data.raw, norm, "data.raw"), norm = norm,
                scale = .sc_update_layer(old$data.scale, norm, "data.scale"),
                smooth = .sc_update_layer(old$data.project, norm, "data.project", TRUE),
                signaling = .sc_update_layer(old$data.signaling, norm, "data.signaling", TRUE))
  var.features <- old$var.features
  if (is.character(var.features)) var.features <- list(features = var.features)
  if (!is.list(var.features)) stop("Legacy var.features must be a list or character vector.", call. = FALSE)
  param <- list(communication = options$parameter %||% list(), legacy = options)
  misc <- list(.mode = "single", .datatype = datatype, .param = param,
               .log = list(), .var.features = var.features, .datasets = list(),
               .legacy = list(images = old$images, options = options,
                              net.sum = old$net$sum, net.sum.cell = old$net$sum.cell,
                              netP.sum = old$netP$sum, netP.sum.cell = old$netP$sum.cell))
  result <- SpatialCellChat(
    assay = assay, images = .sc_update_images(old$images, datatype, cells),
    meta = meta, idents = ident, features = data.frame(row.names = genes),
    LR = old$LR, dr = old$dr,
    net = .sc_update_network(old$net, cells, levels(ident), "net"),
    netP = .sc_update_network(old$netP, cells, levels(ident), "netP"),
    DB = old$DB, misc = misc)
  validateSpatialCellChat(result)
  .sc_validate_after_update(result)
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
