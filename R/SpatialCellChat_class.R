# ====== SparseChatArray ======
#' SparseChatArray: A 3D sparse array backed by list_of_dgCMatrix
#'
#' @description
#' Stores a stack of identically-dimensioned sparse matrices (dgCMatrix) as a list,
#' exposing a 3D array-like interface via S3 methods. Each element along the third
#' dimension is one dgCMatrix, enabling GPU-friendly access per layer.
#'
#' @param x a non-empty named list of dgCMatrix with identical dimensions
#' @param dimnames optional object-level dimnames: row, column, and layer names
#' @return An object of class SparseChatArray.
#' @export
SparseChatArray <- function(x, dimnames = NULL) {
  if (!is.list(x)) stop("x must be a list")
  if (length(x) == 0) stop("x must be non-empty")
  if (!all(vapply(x, inherits, logical(1), "dgCMatrix")))
    stop("all elements must be dgCMatrix")
  dims <- unique(lapply(x, dim))
  if (length(dims) > 1) stop("all matrices must have the same dimensions")
  dims <- dims[[1]]
  dn <- .SparseChatArray_dimnames(x, dims, dimnames)
  .new_SparseChatArray(x, dn)
}

#' Subset a sparse communication array
#'
#' @param x A SparseChatArray.
#' @param i Sender-axis selector.
#' @param j Receiver-axis selector.
#' @param k Signaling layer selector.
#' @param drop Whether to simplify the result when an axis has length one.
#' @return A SparseChatArray, or a dropped matrix or vector when requested.
#' @export
`[.SparseChatArray` <- function(x, i, j, k, drop = TRUE) {
  d <- dim(x)
  dn <- dimnames(x)
  missing_i <- missing(i)
  missing_j <- missing(j)
  missing_k <- missing(k)

  idx_i <- if (missing_i) NULL else .SparseChatArray_resolve_index(i, dn[[1]], d[1], "i")
  idx_j <- if (missing_j) NULL else .SparseChatArray_resolve_index(j, dn[[2]], d[2], "j")
  idx_k <- if (missing_k) seq_len(d[3]) else .SparseChatArray_resolve_index(k, dn[[3]], d[3], "k")

  layers <- unclass(x)[idx_k]
  layer_names <- .SparseChatArray_subset_dimname(dn[[3]], idx_k, FALSE)

  if (!missing_k && length(layers) == 1L && drop) {
    mat <- layers[[1]]
    dimnames(mat) <- dn[1:2]
    return(.SparseChatArray_subset_layer(mat, idx_i, idx_j, missing_i, missing_j, drop = drop))
  } else {
    out <- lapply(layers, function(mat) {
      .SparseChatArray_subset_layer(mat, idx_i, idx_j, missing_i, missing_j, drop = FALSE)
    })
    out_dimnames <- list(
      .SparseChatArray_subset_dimname(dn[[1]], idx_i, missing_i),
      .SparseChatArray_subset_dimname(dn[[2]], idx_j, missing_j),
      layer_names
    )
    return(.new_SparseChatArray(out, out_dimnames))
  }
}

#' Extract a sparse communication layer
#'
#' @param x A SparseChatArray.
#' @param i Integer or name identifying one signaling layer.
#' @param ... Reserved for compatibility with the subsetting generic.
#' @param exact Whether character layer names must match exactly.
#' @return The selected sparse matrix.
#' @export
`[[.SparseChatArray` <- function(x, i, ..., exact = TRUE) {
  if (length(list(...)) > 0L)
    stop("[[.SparseChatArray does not support unused extra arguments", call. = FALSE)
  if (!isTRUE(exact))
    stop("[[.SparseChatArray requires exact layer matching", call. = FALSE)
  if (missing(i))
    stop("[[.SparseChatArray requires a layer index", call. = FALSE)

  d <- dim(x)
  dn <- dimnames(x)
  idx <- .SparseChatArray_resolve_index(i, dn[[3]], d[3], "layer")
  layers <- unclass(x)[idx]
  layer_names <- .SparseChatArray_subset_dimname(dn[[3]], idx, FALSE)

  if (length(layers) == 1L) {
    mat <- layers[[1]]
    dimnames(mat) <- dn[1:2]
    return(mat)
  }

  .new_SparseChatArray(
    layers,
    list(dn[[1]], dn[[2]], layer_names)
  )
}

#' Get sparse communication array dimensions
#'
#' @param x A SparseChatArray.
#' @return Integer vector of sender, receiver, and signaling layer counts.
#' @export
`dim.SparseChatArray` <- function(x) {
  c(nrow(unclass(x)[[1]]), ncol(unclass(x)[[1]]), length(x))
}

#' Get sparse communication axis labels
#'
#' @param x A SparseChatArray.
#' @return A list of sender, receiver, and signaling layer names.
#' @export
`dimnames.SparseChatArray` <- function(x) {
  attr(x, "Dimnames", exact = TRUE) %||%
    .SparseChatArray_dimnames(unclass(x), dim(x), NULL)
}

#' Set sparse communication axis labels
#'
#' @param x A SparseChatArray.
#' @param value A list of three axis name vectors, or NULL.
#' @return The array with validated axis names.
#' @export
`dimnames<-.SparseChatArray` <- function(x, value) {
  dn <- .SparseChatArray_dimnames(unclass(x), dim(x), value)
  .new_SparseChatArray(unclass(x), dn)
}

#' Count sparse communication layers
#'
#' @param x A SparseChatArray.
#' @return The number of layers, including empty layers.
#' @export
`length.SparseChatArray` <- function(x) length(unclass(x))

#' Get signaling layer names
#'
#' @param x A SparseChatArray.
#' @return A character vector of layer names, or NULL.
#' @export
`names.SparseChatArray` <- function(x) dimnames(x)[[3]]

#' Set signaling layer names
#'
#' @param x A SparseChatArray.
#' @param value A character vector of layer names, or NULL.
#' @return The array with updated layer names.
#' @export
`names<-.SparseChatArray` <- function(x, value) {
  dn <- dimnames(x)
  dn[[3]] <- .SparseChatArray_validate_dimname(value, dim(x)[3], "names")
  .new_SparseChatArray(unclass(x), dn)
}

#' Print a sparse communication array summary
#'
#' @param x A SparseChatArray.
#' @param ... Additional arguments accepted by the print method.
#' @return The original array, invisibly.
#' @usage \method{print}{SparseChatArray}(x, ...)
#' @export
print.SparseChatArray <- function(x, ...) {
  d <- dim(x)
  dn <- dimnames(x)
  cat("SparseChatArray:", d[1], "x", d[2], "x", d[3], "\n")
  if (!is.null(dn[[3]])) {
    cat("Layers:", paste(head(dn[[3]], 6), collapse = ", "))
    if (length(dn[[3]]) > 6) cat(", ...")
    cat("\n")
  }
  nnz <- sum(vapply(unclass(x), function(m) length(m@x), integer(1)))
  total <- d[1] * d[2] * d[3]
  cat(sprintf("Sparsity: %.2f%%\n", 100 * (1 - nnz / total)))
  invisible(x)
}

#' Sum sparse communication probabilities by margin
#'
#' @return A vector of retained axis totals or a matrix of sender by receiver totals, depending on the margin.
#' @noRd
#' @export
#' Sum margins of a SparseChatArray
#' @param x a SparseChatArray
#' @param margin 1=row sums per layer, 2=col sums per layer, 3=sum across all layers
#' @export
marginSums <- function(x, margin = 3L) {
  UseMethod("marginSums")
}
#' Sum a sparse communication array by margin
#' @param x A SparseChatArray.
#' @param margin The axis to sum over, using 1 for senders, 2 for receivers, and 3 for signaling layers.
#' @return A vector for margin 1 or 2, or a sparse matrix when summing over the signaling layers with margin 3.
#' @noRd
marginSums.SparseChatArray <- function(x, margin = 3L) {
  d <- dim(x)
  if (margin == 3L) {
    if (length(x) == 1L) {
      res <- unclass(x)[[1L]]
    } else {
      cpp_sum_layers <- get("cpp_sum_layers", inherits = TRUE)
      if (!is.function(cpp_sum_layers)) {
        stop("cpp_sum_layers is required to sum multiple SparseChatArray layers",
             call. = FALSE)
      }
      res <- cpp_sum_layers(unclass(x))
    }
    dimnames(res) <- dimnames(x)[1:2]
    return(res)
  } else if (margin == 1L) {
    res <- vapply(unclass(x), function(m) Matrix::rowSums(m), numeric(d[1]))
    dimnames(res) <- list(dimnames(x)[[1]], dimnames(x)[[3]])
    return(res)
  } else if (margin == 2L) {
    res <- vapply(unclass(x), function(m) Matrix::colSums(m), numeric(d[2]))
    dimnames(res) <- list(dimnames(x)[[2]], dimnames(x)[[3]])
    return(res)
  } else {
    stop("margin must be 1, 2, or 3", call. = FALSE)
  }
}

#' List stored sparse communication entries
#'
#' @param x A SparseChatArray.
#' @param row.names Row names for the returned data frame.
#' @param ... Additional arguments accepted by the data frame method.
#' @return A data frame of stored entries and their sender, receiver, and layer indices.
#' @usage \method{as.data.frame}{SparseChatArray}(x, row.names = NULL, ...)
#' @export
as.data.frame.SparseChatArray <- function(x, row.names = NULL, ...) {
  d <- dim(x)
  dn <- dimnames(x)
  layers <- names(x) %||% as.character(seq_len(d[3]))
  rows <- dn[[1]] %||% as.character(seq_len(d[1]))
  cols <- dn[[2]] %||% as.character(seq_len(d[2]))

  result <- do.call(rbind, lapply(seq_len(d[3]), function(k) {
    mat <- unclass(x)[[k]]
    if (length(mat@x) == 0) return(NULL)
    ii <- mat@i + 1L
    jj <- rep(seq_len(ncol(mat)), diff(mat@p))
    data.frame(
      source = rows[ii],
      target = cols[jj],
      layer  = layers[k],
      value  = mat@x,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(result)) {
    result <- data.frame(
      source = character(), target = character(),
      layer = character(), value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  result
}

#' Transpose all sparse communication layers
#'
#' @param x A SparseChatArray.
#' @return A SparseChatArray with the sender and receiver axes exchanged.
#' @usage \method{t}{SparseChatArray}(x)
#' @export
`t.SparseChatArray` <- function(x) {
  out <- pbapply::pblapply(unclass(x), t)
  dn <- dimnames(x)
  .new_SparseChatArray(out, list(dn[[2]], dn[[1]], dn[[3]]))
}

#' Construct an array from validated sparse layers
#' @param x The value to inspect or transform.
#' @param dimnames The sender, receiver, and signaling axis names.
#' @return A SparseChatArray with the supplied axis labels.
#' @noRd
.new_SparseChatArray <- function(x, dimnames) {
  x <- lapply(x, .SparseChatArray_drop_matrix_dimnames)
  names(x) <- NULL
  structure(x, class = "SparseChatArray", Dimnames = dimnames)
}

#' Apply a transformation independently to every SparseChatArray layer.
#'
#' The outer object stores only the class marker and the object-level
#' `Dimnames` attribute; individual matrix dimnames are normalized by
#' `.new_SparseChatArray()`.
#' @param x The value to inspect or transform.
#' @param FUN A function applied to each sparse communication layer.
#' @return A SparseChatArray with the original axis labels and transformed layers.
#' @noRd
.sc_map_sparse_layers <- function(x, FUN) {
  dn <- dimnames(x)
  layers <- lapply(unclass(x), FUN)
  .new_SparseChatArray(layers, dn)
}

#' Resolve sparse communication axis names
#' @param x The value to inspect or transform.
#' @param dims The expected sender, receiver, and signaling extents.
#' @param dimnames The sender, receiver, and signaling axis names.
#' @return A list of validated axis names, inferred from layers when absent.
#' @noRd
.SparseChatArray_dimnames <- function(x, dims, dimnames = NULL) {
  if (length(dims) == 2L) dims <- c(dims, length(x))

  if (is.null(dimnames)) {
    return(list(
      .SparseChatArray_common_matrix_dimname(x, 1L, "rownames"),
      .SparseChatArray_common_matrix_dimname(x, 2L, "colnames"),
      names(x)
    ))
  }

  if (!is.list(dimnames) || length(dimnames) != 3L)
    stop("dimnames must be NULL or a list of length 3")

  list(
    .SparseChatArray_validate_dimname(dimnames[[1]], dims[1], "dimnames[[1]]"),
    .SparseChatArray_validate_dimname(dimnames[[2]], dims[2], "dimnames[[2]]"),
    .SparseChatArray_validate_dimname(dimnames[[3]], dims[3], "dimnames[[3]]")
  )
}

#' Check agreement of axis names across layers
#' @param x The value to inspect or transform.
#' @param margin The axis to sum over, using 1 for senders, 2 for receivers, and 3 for signaling layers.
#' @param label The field label used in validation diagnostics.
#' @return The common labels or NULL when no layer supplies them.
#' @noRd
.SparseChatArray_common_matrix_dimname <- function(x, margin, label) {
  values <- lapply(x, function(m) dimnames(m)[[margin]])
  present <- !vapply(values, is.null, logical(1))
  if (!any(present)) return(NULL)

  ref <- values[[which(present)[1]]]
  same <- vapply(values[present], identical, logical(1), ref)
  if (!all(same))
    stop("all matrix ", label, " must be identical; pass dimnames explicitly")

  ref
}

#' Validate one sparse communication axis
#' @param value The input or replacement value.
#' @param expected The expected number of axis labels.
#' @param label The field label used in validation diagnostics.
#' @return Completion is invisible. Invalid labels raise an error.
#' @noRd
.SparseChatArray_validate_dimname <- function(value, expected, label) {
  if (is.null(value)) return(NULL)
  if (!is.character(value)) value <- as.character(value)
  if (length(value) != expected)
    stop(label, " must be NULL or have length ", expected)
  value
}

#' Clear labels on one sparse layer
#' @param mat A matrix for the current layer.
#' @return The matrix with its row and column names removed.
#' @noRd
.SparseChatArray_drop_matrix_dimnames <- function(mat) {
  dimnames(mat) <- list(NULL, NULL)
  mat
}

#' Resolve a sparse communication subscript
#' @param index The index or names selecting positions along the axis.
#' @param dimname The existing names for the axis.
#' @param extent The length of the axis.
#' @param label The field label used in validation diagnostics.
#' @return Integer positions selected by the subscript.
#' @noRd
.SparseChatArray_resolve_index <- function(index, dimname, extent, label) {
  if (!is.character(index)) return(index)
  if (is.null(dimname))
    stop("cannot use character subscript for ", label, " without dimnames")

  resolved <- match(index, dimname)
  missing_index <- is.na(resolved) & !is.na(index)
  if (any(missing_index))
    stop("subscript out of bounds for ", label, ": ",
         paste(index[missing_index], collapse = ", "))

  resolved
}

#' Resolve a single sparse communication layer
#' @param index The index or names selecting positions along the axis.
#' @param dimname The existing names for the axis.
#' @param extent The length of the axis.
#' @param label The field label used in validation diagnostics.
#' @return A single integer layer position.
#' @noRd
.SparseChatArray_resolve_single_index <- function(index, dimname, extent, label) {
  resolved <- .SparseChatArray_resolve_index(index, dimname, extent, label)
  if (length(resolved) != 1L)
    stop("[[.SparseChatArray subscript ", label, " must select exactly one element",
         call. = FALSE)
  if (is.numeric(resolved) && (is.na(resolved) || resolved < 1L || resolved > extent))
    stop("subscript out of bounds for ", label, call. = FALSE)
  if (is.logical(resolved)) {
    resolved <- which(resolved)
    if (length(resolved) != 1L)
      stop("[[.SparseChatArray subscript ", label, " must select exactly one element",
           call. = FALSE)
  }
  resolved
}

#' Select axis labels for a subset
#' @param dimname The existing names for the axis.
#' @param index The index or names selecting positions along the axis.
#' @param missing_index Whether the subscript was omitted.
#' @return The labels corresponding to the selected axis positions, or NULL.
#' @noRd
.SparseChatArray_subset_dimname <- function(dimname, index, missing_index) {
  if (is.null(dimname)) return(NULL)
  if (missing_index) return(dimname)
  dimname[index]
}

#' Subscript the sender and receiver axes of a layer
#' @param mat A matrix for the current layer.
#' @param i Sender-axis selector.
#' @param j Receiver-axis selector.
#' @param missing_i Whether the sender subscript was omitted.
#' @param missing_j Whether the receiver subscript was omitted.
#' @param drop Whether to simplify the result when an axis has length one.
#' @return A subset matrix or a dropped vector or scalar.
#' @noRd
.SparseChatArray_subset_layer <- function(mat, i, j, missing_i, missing_j, drop) {
  if (missing_i && missing_j) {
    mat[, , drop = drop]
  } else if (missing_i) {
    mat[, j, drop = drop]
  } else if (missing_j) {
    mat[i, , drop = drop]
  } else {
    mat[i, j, drop = drop]
  }
}


# ====== SpatialCellChat S4 class ======

#' The 11-slot SpatialCellChat object
#'
#' Stores expression, spatial geometry, observation metadata, ligand-receptor
#' results, pathway results, and derived analyses. Mode and datatype are
#' recorded in the misc slot.
#' @name SpatialCellChat-class
#' @docType class
#' @aliases SpatialCellChat
#' @slot assay Named assay matrices: raw, norm, scale, smooth, and signaling.
#' @slot images Spatial coordinates, coordinate calibration, rasters, FOV data, and distance/grid caches.
#' @slot meta A data frame of cell/spot metadata.
#' @slot idents A factor defining the cell/spot group for each observation.
#' @slot features A data frame of gene-level metadata.
#' @slot LR LR interaction annotations and selected LR tables.
#' @slot dr Dimensional reductions and other embeddings.
#' @slot net LR-level communication results. Group probabilities and p-values are stored under \code{net$group}; cell probabilities are under \code{net$cell}; temporary computation data are under \code{net$tmp}.
#' @slot netP Pathway-level communication results with the same group/cell organization as \code{net}.
#' @slot DB The communication database and gene/complex annotations used in analysis.
#' @slot misc Mode, datatype, parameters, operation log, variable-feature results, and merged dataset snapshots.
#' @exportClass SpatialCellChat
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass setValidity
SpatialCellChat <- methods::setClass("SpatialCellChat",
  slots = c(
    assay    = "list",       # raw, norm, scale, smooth, signaling
    images   = "list",       # coordinates, coordinate.system, spatial.factors, rasters, fov, caches
    meta     = "data.frame", # 细胞级元数据
    idents   = "factor",     # cell group labels for all cells
    features = "data.frame", # 基因级元数据
    LR       = "list",       # LR 通讯条目信息
    dr       = "list",       # pca, umap, spatial and other embeddings
    net      = "list",       # LR-level communication results
    netP     = "list",       # pathway-level communication results
    DB       = "list",       # CellChatDB
    misc     = "list"        # mode, datatype, params, log and datasets
  )
)

#' Test for an ordinary or Matrix matrix
#' @param x The value to inspect or transform.
#' @return A logical scalar.
#' @noRd
.sc_matrix_like <- function(x) {
  inherits(x, "matrix") || inherits(x, "Matrix")
}

#' Validate a matrix against its axis contract
#' @param x The value to inspect or transform.
#' @param expected_dim The required vector of matrix dimensions.
#' @param expected_dimnames The required list of row and column names.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_matrix <- function(x, expected_dim, expected_dimnames, label) {
  errors <- character()
  if (is.null(x)) return(errors)
  if (!.sc_matrix_like(x)) {
    return(paste0(label, " must be a matrix or Matrix object"))
  }
  if (!identical(as.integer(dim(x)), as.integer(expected_dim))) {
    errors <- c(errors, paste0(label, " has incompatible dimensions"))
  }
  if (!is.null(expected_dimnames)) {
    actual <- dimnames(x)
    if (!identical(actual[[1]], expected_dimnames[[1]]) ||
        !identical(actual[[2]], expected_dimnames[[2]])) {
      errors <- c(errors, paste0(label, " dimnames do not match the object"))
    }
  }
  errors
}

# Gene-subset layers (signaling, smooth) share one contract: the object cells as
# columns in order, rows drawn from assay$norm and keeping its gene order.
#' Check a subset assay against object axes
#' @param x The value to inspect or transform.
#' @param gene_names Gene names in canonical assay order.
#' @param cell_names Observation names in canonical object order.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_subset_layer <- function(x, gene_names, cell_names, label) {
  if (!.sc_matrix_like(x)) {
    return(paste0("assay$", label, " must be a matrix or Matrix object"))
  }
  if (ncol(x) != length(cell_names) || !identical(colnames(x), cell_names)) {
    return(paste0("assay$", label, " must carry the object cells as columns"))
  }
  genes <- rownames(x)
  if (is.null(genes)) {
    return(paste0("assay$", label, " must have gene rownames"))
  }
  if (!all(genes %in% gene_names)) {
    return(paste0("assay$", label, " genes must be a subset of assay$norm genes"))
  }
  if (is.unsorted(match(genes, gene_names), strictly = TRUE)) {
    return(paste0("assay$", label, " rows must preserve the assay$norm gene order"))
  }
  character()
}

#' Check the assay list and its aligned matrices
#' @param assay The list of named expression layers.
#' @param gene_names Gene names in canonical assay order.
#' @param cell_names Observation names in canonical object order.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_assay <- function(assay, gene_names, cell_names) {
  if (!is.list(assay)) return("assay must be a list")
  required <- c("raw", "norm", "scale", "smooth", "signaling")
  subset_layers <- c("signaling", "smooth")
  errors <- character()
  missing <- setdiff(required, names(assay))
  if (length(missing)) {
    errors <- c(errors, paste0("assay is missing: ", paste(missing, collapse = ", ")))
  }
  if (is.null(assay$norm)) {
    errors <- c(errors, "assay$norm must be present")
  }
  expected <- c(length(gene_names), length(cell_names))
  expected_names <- list(gene_names, cell_names)
  for (name in setdiff(intersect(required, names(assay)), subset_layers)) {
    errors <- c(errors, .sc_validate_matrix(
      assay[[name]], expected, expected_names, paste0("assay$", name)))
  }
  # assay$signaling is the subsetData gene subset; assay$smooth is the
  # projectData projection of that same space. Both are gene subsets of
  # assay$norm and neither may reorder its rows.
  for (name in intersect(subset_layers, names(assay))) {
    if (!is.null(assay[[name]])) {
      errors <- c(errors, .sc_validate_subset_layer(
        assay[[name]], gene_names, cell_names, name))
    }
  }
  errors
}

#' Check the observation identities
#' @param idents The factor or identity list to validate.
#' @param cell_names Observation names in canonical object order.
#' @param mode The object mode or analysis mode required by the caller.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_idents <- function(idents, cell_names, mode) {
  if (!inherits(idents, "factor"))
    return("idents must be a factor")

  errors <- character()
  if (length(idents) != length(cell_names))
    errors <- c(errors, "idents length must equal the number of cells")
  if (!identical(names(idents), cell_names))
    errors <- c(errors, "names(idents) must match cell names")
  # Missing group labels propagate into group-level communication matrices and
  # make their row/column semantics ambiguous.  Reject them at the class
  # boundary instead of allowing a later aggregation step to fail silently.
  if (anyNA(idents))
    errors <- c(errors, "idents must not contain NA group labels")
  errors
}

#' Check the declared coordinate convention
#' @param coordinate.system The declared system in which spatial coordinates are stored.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_coordinate_system <- function(coordinate.system,
                                            label = "images$coordinate.system") {
  if (is.null(coordinate.system)) return(character())
  if (!is.list(coordinate.system))
    return(paste0(label, " must be a list"))

  errors <- character()
  unit <- coordinate.system$unit
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit))
    errors <- c(errors, paste0(label, "$unit must be one non-empty string"))

  calibrated <- coordinate.system$calibrated
  if (!is.logical(calibrated) || length(calibrated) != 1L || is.na(calibrated))
    errors <- c(errors, paste0(label, "$calibrated must be TRUE or FALSE"))
  if (isTRUE(calibrated) && identical(unit, "unknown"))
    errors <- c(errors, paste0(label, " cannot be calibrated when unit is 'unknown'"))

  y.direction <- coordinate.system$y.direction
  if (!is.null(y.direction) &&
      (!is.character(y.direction) || length(y.direction) != 1L ||
       !y.direction %in% c("up", "down")))
    errors <- c(errors, paste0(label, "$y.direction must be 'up' or 'down'"))
  errors
}

#' Check distance calibration against the coordinate convention
#' @param factors Spatial calibration factors for the declared coordinate system.
#' @param coordinate.system The declared system in which spatial coordinates are stored.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_spatial_factors <- function(factors, coordinate.system,
                                         label = "images$spatial.factors") {
  # An uncalibrated coordinate system is allowed for plotting and relative
  # topology, but it must not be mistaken for micrometer coordinates.  Such an
  # object therefore stores ratio = NULL and is rejected by distance routines
  # until the caller supplies a physical conversion explicitly.
  if (!is.list(factors)) return(paste0(label, " must be a list"))
  errors <- character()
  if (!all(c("ratio", "tol") %in% names(factors)))
    errors <- c(errors, paste0(label, " must contain ratio and tol"))

  ratio <- factors$ratio
  if (!is.null(ratio) &&
      (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 0))
    errors <- c(errors, paste0(label, "$ratio must be NULL or one positive finite number"))

  tol <- factors$tol
  if (!is.null(tol) &&
      (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0))
    errors <- c(errors, paste0(label, "$tol must be NULL or one non-negative finite number"))

  calibrated <- if (is.list(coordinate.system)) isTRUE(coordinate.system$calibrated) else FALSE
  if (isTRUE(calibrated) && is.null(ratio))
    errors <- c(errors, paste0(label, "$ratio is required for calibrated coordinates"))
  if (isTRUE(calibrated) && is.null(tol))
    errors <- c(errors, paste0(label, "$tol is required for calibrated coordinates"))
  if (!isTRUE(calibrated) && !is.null(ratio))
    errors <- c(errors, paste0(label, "$ratio must be NULL when coordinates are uncalibrated"))
  errors
}
#' Check the raster coordinate transform
#' @param transform The affine mapping connecting raster and spatial coordinates.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_transform <- function(transform,
                                   label = "images$rasters transform") {
  if (!is.list(transform)) return(paste0(label, " must be a list"))
  errors <- character()
  matrix.value <- transform$matrix
  if (!is.matrix(matrix.value) || !identical(dim(matrix.value), c(3L, 3L)) ||
      !is.numeric(matrix.value) || any(!is.finite(matrix.value)))
    errors <- c(errors, paste0(label, "$matrix must be a finite numeric 3 x 3 matrix"))
  for (field in c("from", "to")) {
    value <- transform[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
      errors <- c(errors, paste0(label, "$", field, " must be one non-empty string"))
  }
  errors
}

#' Check raster records and their observation references
#' @param rasters The named list of raster records.
#' @param cell_names Observation names in canonical object order.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_rasters <- function(rasters, cell_names) {
  if (is.null(rasters)) return(character())
  if (!is.list(rasters)) return("images$rasters must be a list")
  if (!length(rasters)) return(character())

  errors <- character()
  raster_names <- names(rasters)
  if (is.null(raster_names)) raster_names <- as.character(seq_along(rasters))
  for (i in seq_along(rasters)) {
    label <- paste0("images$rasters$", raster_names[[i]])
    raster <- rasters[[i]]
    if (!is.list(raster)) {
      errors <- c(errors, paste0(label, " must be a list"))
      next
    }
    image <- raster$image
    image_dim <- if (is.null(image)) NULL else dim(image)
    if (is.null(image) || !is.array(image) || length(image_dim) < 2L ||
        length(image_dim) > 3L || any(image_dim[seq_len(2L)] < 1L))
      errors <- c(errors, paste0(label, "$image must be a non-empty 2D raster or 3D RGB array"))
    if (!is.null(raster$scale.factors) && !is.list(raster$scale.factors))
      errors <- c(errors, paste0(label, "$scale.factors must be a list"))
    if (!is.null(raster$coordinate.system))
      errors <- c(errors, .sc_validate_coordinate_system(
        raster$coordinate.system, paste0(label, "$coordinate.system")))
    if (is.null(raster$transform)) {
      errors <- c(errors, paste0(label, "$transform is required"))
    } else {
      errors <- c(errors, .sc_validate_transform(
        raster$transform, paste0(label, "$transform")))
    }
    if (!is.null(raster$spot.coordinates) &&
        !is.data.frame(raster$spot.coordinates))
      errors <- c(errors, paste0(label, "$spot.coordinates must be a data.frame"))
    if (!is.null(raster$spot.coordinates) &&
        !is.null(rownames(raster$spot.coordinates)) &&
        anyDuplicated(rownames(raster$spot.coordinates)))
      errors <- c(errors, paste0(label, "$spot.coordinates rownames must be unique"))
    if (!is.null(raster$spot.radius) &&
        (!is.numeric(raster$spot.radius) || length(raster$spot.radius) != 1L ||
         !is.finite(raster$spot.radius) || raster$spot.radius < 0))
      errors <- c(errors, paste0(label, "$spot.radius must be one non-negative finite number"))
  }
  errors
}

#' Check the field-of-view metadata records
#' @param fov Per-field-of-view metadata.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_fov <- function(fov) {
  if (is.null(fov)) return(character())
  if (!is.list(fov)) return("images$fov must be a list")
  if (!length(fov)) return(character())

  errors <- character()
  fov_names <- names(fov)
  if (is.null(fov_names)) fov_names <- as.character(seq_along(fov))
  for (i in seq_along(fov)) {
    label <- paste0("images$fov$", fov_names[[i]])
    layer <- fov[[i]]
    if (!is.list(layer)) {
      errors <- c(errors, paste0(label, " must be a list"))
      next
    }
    if (!is.null(layer$centroids) && !is.data.frame(layer$centroids))
      errors <- c(errors, paste0(label, "$centroids must be a data.frame"))
    if (!is.null(layer$boundaries) && !is.list(layer$boundaries))
      errors <- c(errors, paste0(label, "$boundaries must be a list"))
    if (!is.null(layer$boundaries) && is.list(layer$boundaries)) {
      for (boundary_name in names(layer$boundaries))
        if (!is.data.frame(layer$boundaries[[boundary_name]]))
          errors <- c(errors, paste0(label, "$boundaries$", boundary_name,
                                     " must be a data.frame"))
    }
    if (!is.null(layer$molecules) && !is.data.frame(layer$molecules))
      errors <- c(errors, paste0(label, "$molecules must be a data.frame"))
  }
  errors
}

#' Check the sparse distance and contact cache
#' @param distance The distance cache record, or distances to plot in a density panel.
#' @param cell_names Observation names in canonical object order.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_distance_cache <- function(distance, cell_names) {
  if (!is.list(distance)) return("images$.distance must be a list")
  errors <- character()
  validate_matrix <- function(value, label) {
    if (is.null(value)) return(character())
    .sc_validate_matrix(value, c(length(cell_names), length(cell_names)),
                        list(cell_names, cell_names), label)
  }
  for (name in intersect(c("d.spatial", "adj.contact"), names(distance)))
    errors <- c(errors, validate_matrix(distance[[name]], paste0("images$.distance$", name)))

  errors
}

#' Check spatial data and cached geometry
#' @param images The spatial images list.
#' @param cell_names Observation names in canonical object order.
#' @param datatype The data type recorded in the object.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_images <- function(images, cell_names, datatype) {
  if (!is.list(images)) return("images must be a list")
  errors <- character()
  coordinates <- images$coordinates
  coordinate.system <- images$coordinate.system
  coordinate.scope <- images$coordinate.scope
  if (!is.null(coordinate.scope) &&
      (!is.character(coordinate.scope) || length(coordinate.scope) != 1L ||
       is.na(coordinate.scope) ||
       !coordinate.scope %in% c("shared", "dataset-local")))
    errors <- c(errors, "images$coordinate.scope must be 'shared' or 'dataset-local'")

  if (identical(datatype, "spatial")) {
    if (identical(coordinate.scope, "dataset-local")) {
      if (!is.null(coordinates) || !is.null(coordinate.system) ||
          !is.null(images$spatial.factors))
        errors <- c(errors, paste0("dataset-local coordinates must remain in misc$.datasets snapshots; ",
                                   "top-level spatial coordinates and calibration must be NULL"))
    } else {
    # Spatial objects must always have analysis coordinates.  Raster images and
    # FOV geometries are optional because many platforms provide coordinates
    # without a histology image (for example Slide-seq and Stereo-seq).
    if (!.sc_matrix_like(coordinates)) {
      errors <- c(errors, "spatial objects require images$coordinates")
    } else {
      values <- tryCatch(as.matrix(coordinates), error = function(e) NULL)
      if (is.null(values) || nrow(values) != length(cell_names) ||
          !ncol(values) %in% c(2L, 3L))
        errors <- c(errors, "images$coordinates must be cells x 2 or cells x 3")
      if (!is.null(values) && (!is.numeric(values) || any(!is.finite(values))))
        errors <- c(errors, "images$coordinates must contain finite numeric values")
      if (!identical(rownames(coordinates), cell_names))
        errors <- c(errors, "images$coordinates rownames must match cell names")
    }
    if (is.null(coordinate.system))
      errors <- c(errors, "spatial objects require images$coordinate.system")
    errors <- c(errors, .sc_validate_coordinate_system(coordinate.system))
    errors <- c(errors, .sc_validate_spatial_factors(
      images$spatial.factors, coordinate.system))
    }
  }
  if ("histology" %in% names(images))
    errors <- c(errors, "images$histology is obsolete; store the record in images$rasters")

  errors <- c(errors, .sc_validate_rasters(images$rasters, cell_names))
  errors <- c(errors, .sc_validate_fov(images$fov))
  errors <- c(errors, .sc_validate_distance_cache(images$.distance, cell_names))
  if (!is.null(images$.grid) && !is.list(images$.grid))
    errors <- c(errors, "images$.grid must be a list")
  unique(errors[nzchar(errors)])
}


#' Check alignment of communication layers and axes
#' @param net compute the centrality measures on a specific signaling network given by a 2 or 3 dimemsional array net
#' @param cell_names Observation names in canonical object order.
#' @param group_names Expected group labels in identity level order.
#' @param label The field label used in validation diagnostics.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_net <- function(net, cell_names, group_names, label) {
  if (!is.list(net)) return(paste0(label, " must be a list"))
  errors <- character()
  for (resolution in c("cell", "group")) {
    node_names <- if (resolution == "cell") cell_names else group_names
    branch <- net[[resolution]]
    if (is.null(branch)) {
      errors <- c(errors, paste0(label, "$", resolution, " is missing"))
      next
    }
    if (!is.list(branch)) {
      errors <- c(errors, paste0(label, "$", resolution, " must be a list"))
      next
    }
    for (measure in intersect(c("prob", "pval"), names(branch))) {
      value <- branch[[measure]]
      if (!inherits(value, "SparseChatArray")) {
        errors <- c(errors, paste0(label, "$", resolution, "$", measure,
                                   " must be a SparseChatArray"))
      } else if (!identical(as.integer(dim(value)[1:2]),
                            c(length(node_names), length(node_names)))) {
        errors <- c(errors, paste0(label, "$", resolution, "$", measure,
                                   " has incompatible dimensions"))
      } else {
        dn <- dimnames(value)
        if (!identical(dn[[1]], node_names) || !identical(dn[[2]], node_names))
          errors <- c(errors, paste0(label, "$", resolution, "$", measure,
                                     " dimnames do not match node names"))
      }
    }
    for (measure in intersect(c("count", "weight"), names(branch))) {
      errors <- c(errors, .sc_validate_matrix(
        branch[[measure]], c(length(node_names), length(node_names)),
        list(node_names, node_names), paste0(label, "$", resolution, "$", measure)))
    }
    if (resolution == "group" && !is.null(branch$group.merged)) {
      group.merged <- branch$group.merged
      if (!is.factor(group.merged) || length(group.merged) != length(group_names) ||
          anyNA(group.merged) || !identical(names(group.merged), group_names)) {
        errors <- c(errors, paste0(label, "$group$group.merged must be a named factor aligned to group names"))
      } else {
        merged.names <- levels(group.merged)
        for (measure in intersect(c("count.merged", "weight.merged"), names(branch)))
          errors <- c(errors, .sc_validate_matrix(
            branch[[measure]], c(length(merged.names), length(merged.names)),
            list(merged.names, merged.names), paste0(label, "$group$", measure)))
      }
    } else if (resolution == "group" &&
               any(c("count.merged", "weight.merged") %in% names(branch))) {
      errors <- c(errors, paste0(label, "$group merged matrices require group.merged"))
    }
  }
  errors
}

#' Check the individual snapshots of a merged object
#' @param object A SpatialCellChat object.
#' @return A character vector of validation problems.
#' @noRd
.sc_validate_merged_datasets <- function(object) {
  datasets <- object@misc$.datasets
  if (!is.list(datasets) || length(datasets) < 2L || is.null(names(datasets)) ||
      anyNA(names(datasets)) || any(!nzchar(names(datasets))) || anyDuplicated(names(datasets)))
    return("merged misc$.datasets must contain at least two uniquely named snapshots")
  if (!"datasets" %in% names(object@meta))
    return("merged meta must contain a datasets column")
  if (!identical(levels(object@meta$datasets), names(datasets)))
    return("merged meta$datasets levels must match misc$.datasets names")

  errors <- character()
  merged.cells <- rownames(object@meta)
  for (dataset in names(datasets)) {
    record <- datasets[[dataset]]
    label <- paste0("misc$.datasets[['", dataset, "']]")
    if (!is.list(record)) {
      errors <- c(errors, paste0(label, " must be a list"))
      next
    }
    cells <- if (is.data.frame(record$meta)) rownames(record$meta) else NULL
    if (is.null(cells) || anyDuplicated(cells)) {
      errors <- c(errors, paste0(label, "$meta must have unique cell rownames"))
      next
    }
    expected.cells <- merged.cells[as.character(object@meta$datasets) == dataset]
    if (!identical(cells, expected.cells))
      errors <- c(errors, paste0(label, "$meta cells must match the merged dataset partition in order"))
    if (!inherits(record$idents, "factor") || length(record$idents) != length(cells) ||
        !identical(names(record$idents), cells) || anyNA(record$idents))
      errors <- c(errors, paste0(label, "$idents must be a named factor aligned to snapshot cells"))
    group.names <- if (inherits(record$idents, "factor")) levels(record$idents) else character()
    errors <- c(errors, .sc_validate_net(record$net, cells, group.names,
                                          paste0(label, "$net")))
    errors <- c(errors, .sc_validate_net(record$netP, cells, group.names,
                                          paste0(label, "$netP")))
    if (!is.list(record$LR) || !is.data.frame(record$LR$LRsig))
      errors <- c(errors, paste0(label, "$LR$LRsig must be a data.frame"))
    errors <- c(errors, .sc_validate_images(record$images, cells, object@misc$.datatype))
  }
  unique(errors[nzchar(errors)])
}

#' Collect 11-slot schema validation problems
#' @param object A SpatialCellChat object.
#' @return A character vector of validation problems, empty for a valid object.
#' @noRd
.sc_validate_object <- function(object) {
  errors <- character()
  cell_names <- colnames(object@assay$norm)
  gene_names <- rownames(object@assay$norm)
  if (is.null(cell_names) || is.null(gene_names))
    return("assay$norm must have row and column names")

  # Duplicate identifiers make dimnames-based subsetting ambiguous.  The
  # constructor rejects them too, while this check protects objects assembled
  # manually or modified through direct slot access.
  if (anyDuplicated(cell_names))
    errors <- c(errors, "assay$norm column names must be unique")
  if (anyDuplicated(gene_names))
    errors <- c(errors, "assay$norm row names must be unique")

  errors <- c(errors, .sc_validate_assay(object@assay, gene_names, cell_names))
  if (nrow(object@meta) != length(cell_names))
    errors <- c(errors, "meta must have one row per cell")
  if (!identical(rownames(object@meta), cell_names))
    errors <- c(errors, "meta rownames must match assay$norm colnames")
  if (nrow(object@features) != length(gene_names) ||
      !identical(rownames(object@features), gene_names))
    errors <- c(errors, "features rownames must match assay$norm rownames")
  mode <- object@misc$.mode
  datatype <- object@misc$.datatype
  if (!identical(mode, "single") && !identical(mode, "merged"))
    errors <- c(errors, "misc$.mode must be 'single' or 'merged'")
  if (!identical(datatype, "RNA") && !identical(datatype, "spatial"))
    errors <- c(errors, "misc$.datatype must be 'RNA' or 'spatial'")
  if (identical(object@images$coordinate.scope, "dataset-local") &&
      !identical(mode, "merged"))
    errors <- c(errors, "images$coordinate.scope = 'dataset-local' is reserved for merged objects")
  if (identical(mode, "merged") && identical(datatype, "spatial") &&
      !identical(object@images$coordinate.scope, "dataset-local"))
    errors <- c(errors, "merged spatial objects must declare images$coordinate.scope = 'dataset-local'")
  errors <- c(errors, .sc_validate_idents(object@idents, cell_names, mode))
  errors <- c(errors, .sc_validate_images(object@images, cell_names, datatype))
  group_names <- levels(object@idents)
  errors <- c(errors, .sc_validate_net(object@net, cell_names, group_names, "net"))
  errors <- c(errors, .sc_validate_net(object@netP, cell_names, group_names, "netP"))
  if (!is.list(object@misc)) errors <- c(errors, "misc must be a list")
  if (is.list(object@misc)) {
    for (name in c(".param", ".log", ".var.features", ".datasets")) {
      if (!name %in% names(object@misc))
        errors <- c(errors, paste0("misc is missing ", name))
    }
    if (identical(mode, "single") && length(object@misc$.datasets) > 0)
      errors <- c(errors, "single objects must not contain misc$.datasets")
    if (identical(mode, "merged") && !is.list(object@misc$.datasets))
      errors <- c(errors, "merged misc$.datasets must be a list")
    if (identical(mode, "merged"))
      errors <- c(errors, .sc_validate_merged_datasets(object))
  }
  unique(errors[nzchar(errors)])
}
# ====== spatial image helpers ======

#' Create empty spatial storage
#' @return The default list of image and spatial geometry entries.
#' @noRd
.sc_default_images <- function() {
  list(
    coordinates = NULL,
    coordinate.system = NULL,
    spatial.factors = NULL,
    rasters = list(),
    fov = list(),
    .distance = list(),
    .grid = list()
  )
}

#' Normalize the input coordinate convention
#' @param coordinate.system The declared system in which spatial coordinates are stored.
#' @param spatial.factors Spatial calibration factors including distance conversion where applicable.
#' @return The supported coordinate-system label.
#' @noRd
.sc_normalize_coordinate_system <- function(coordinate.system, spatial.factors = NULL) {
  if (is.null(coordinate.system)) {
    ratio <- if (is.list(spatial.factors)) spatial.factors$ratio else NULL
    coordinate.system <- list(
      unit = if (!is.null(ratio)) "pixel" else "arbitrary",
      calibrated = !is.null(ratio),
      y.direction = "up",
      origin = "lower-left"
    )
  } else {
    if (!is.list(coordinate.system))
      stop("coordinate.system must be a list", call. = FALSE)
    coordinate.system <- as.list(coordinate.system)
    if (is.null(coordinate.system$unit)) coordinate.system$unit <- "unknown"
    if (is.null(coordinate.system$calibrated))
      coordinate.system$calibrated <- !is.null(spatial.factors$ratio)
    if (is.null(coordinate.system$y.direction)) coordinate.system$y.direction <- "up"
    if (is.null(coordinate.system$origin)) coordinate.system$origin <- "lower-left"
  }
  coordinate.system
}

#' Normalize distance calibration for the coordinate system
#' @param spatial.factors Spatial calibration factors including distance conversion where applicable.
#' @param coordinate.system The declared system in which spatial coordinates are stored.
#' @return The normalized spatial factors list.
#' @noRd
.sc_normalize_spatial_factors <- function(spatial.factors, coordinate.system) {
  if (is.null(spatial.factors)) spatial.factors <- list()
  if (!is.list(spatial.factors))
    stop("spatial.factors must be a list", call. = FALSE)
  if (is.null(spatial.factors$ratio) && identical(coordinate.system$unit, "um"))
    spatial.factors$ratio <- 1
  if (!"ratio" %in% names(spatial.factors)) spatial.factors$ratio <- NULL
  if (!"tol" %in% names(spatial.factors)) spatial.factors$tol <- NULL
  spatial.factors
}

#' Normalize an affine raster transform
#' @param transform The affine mapping connecting raster and spatial coordinates.
#' @return A list with from and to coordinate-frame names and a finite 3 by 3 affine matrix.
#' @noRd
.sc_normalize_transform <- function(transform = NULL) {
  if (is.null(transform)) transform <- list()
  if (is.matrix(transform)) transform <- list(matrix = transform)
  if (!is.list(transform) || is.null(transform$matrix))
    stop("raster transform must contain a 3 x 3 matrix", call. = FALSE)
  matrix.value <- if (.sc_matrix_like(transform$matrix)) {
    as.matrix(transform$matrix)
  } else {
    transform$matrix
  }
  if (!is.matrix(matrix.value) || !identical(dim(matrix.value), c(3L, 3L)) ||
      !is.numeric(matrix.value) || any(!is.finite(matrix.value)))
    stop("raster transform matrix must be a finite numeric 3 x 3 matrix", call. = FALSE)
  list(
    from = transform$from %||% "analysis",
    to = transform$to %||% "raster_pixel",
    matrix = matrix.value
  )
}

#' Normalize one spatial image record
#' @param image A package raster record or supported Seurat spatial image.
#' @param cell_names Observation names in canonical object order.
#' @return A normalized raster record, or NULL when no image is supplied.
#' @noRd
.sc_normalize_raster <- function(image, cell_names = NULL) {
  if (is.null(image)) return(NULL)
  raster <- NULL

  if (is.list(image) && !methods::is(image, "SpatialImage")) {
    # A plain list is the package-native interchange format.  `coordinates`
    # is accepted only as an input spelling and is renamed to
    # `spot.coordinates`; canonical analysis coordinates live at
    # images$coordinates and are never duplicated inside a raster record.
    raster <- image
    if (is.null(raster$image) && !is.null(raster$raster))
      raster$image <- raster$raster
    if (is.null(raster$image))
      stop("raster list must contain image", call. = FALSE)
    if (is.null(raster$spot.coordinates) && !is.null(raster$coordinates)) {
      raster$spot.coordinates <- raster$coordinates
      raster$coordinates <- NULL
    }
  } else if (inherits(image, "VisiumV1") || methods::is(image, "SpatialImage")) {
    # Seurat 4 exposes VisiumV1 slots; newer SeuratObject versions expose
    # additional SpatialImage subclasses.  Prefer semantic accessors when
    # available, then fall back to known slots so one package version does not
    # hard-code another version's internal class layout.
    raster <- list(source = list(type = class(image)[[1L]]))
    slots <- methods::slotNames(image)
    if ("image" %in% slots) raster$image <- methods::slot(image, "image")
    if (is.null(raster$image) && "raster" %in% slots) {
      raster.value <- methods::slot(image, "raster")
      if (is.list(raster.value) && !is.null(raster.value$image))
        raster$image <- raster.value$image
      else if (is.array(raster.value) || is.matrix(raster.value))
        raster$image <- raster.value
    }
    if (is.null(raster$image) && requireNamespace("SeuratObject", quietly = TRUE)) {
      raster.value <- tryCatch(
        SeuratObject::GetImage(image, mode = "raster"),
        error = function(e) NULL
      )
      if (is.array(raster.value) || is.matrix(raster.value))
        raster$image <- raster.value
    }
    if ("scale.factors" %in% slots)
      raster$scale.factors <- as.list(methods::slot(image, "scale.factors"))
    if ("coordinates" %in% slots)
      raster$spot.coordinates <- methods::slot(image, "coordinates")
    if (is.null(raster$spot.coordinates) && requireNamespace("SeuratObject", quietly = TRUE)) {
      raster$spot.coordinates <- tryCatch(
        SeuratObject::GetTissueCoordinates(image), error = function(e) NULL)
    }
    if ("spot.radius" %in% slots)
      raster$spot.radius <- methods::slot(image, "spot.radius")
    if (is.null(raster$spot.radius) && requireNamespace("SeuratObject", quietly = TRUE))
      raster$spot.radius <- tryCatch(SeuratObject::Radius(image), error = function(e) NULL)
    if (is.null(raster$image))
      stop("the Seurat spatial image does not expose a raster image", call. = FALSE)
  } else {
    stop("image must be a Seurat SpatialImage object or a raster list", call. = FALSE)
  }

  if (!is.array(raster$image) && !is.matrix(raster$image))
    stop("raster$image must be a raster array or matrix", call. = FALSE)
  if (!is.null(raster$spot.coordinates)) {
    if (!is.data.frame(raster$spot.coordinates))
      stop("raster$spot.coordinates must be a data.frame", call. = FALSE)
    if (is.null(rownames(raster$spot.coordinates)) &&
        "barcodes" %in% names(raster$spot.coordinates))
      rownames(raster$spot.coordinates) <- as.character(raster$spot.coordinates$barcodes)
    if (!is.null(cell_names)) {
      if (is.null(rownames(raster$spot.coordinates)))
        stop("raster$spot.coordinates must have rownames when cell_names are supplied", call. = FALSE)
      missing_cells <- setdiff(cell_names, rownames(raster$spot.coordinates))
      if (length(missing_cells))
        stop("raster$spot.coordinates is missing cells: ", paste(missing_cells, collapse = ", "), call. = FALSE)
      raster$spot.coordinates <- raster$spot.coordinates[cell_names, , drop = FALSE]
    }
  }
  if (is.null(raster$scale.factors)) raster$scale.factors <- list()
  if (is.null(raster$source)) raster$source <- list(type = "custom")
  if (is.null(raster$coordinate.system)) {
    raster$coordinate.system <- list(
      unit = "pixel", calibrated = FALSE,
      y.direction = "down", origin = "top-left"
    )
  }
  if (is.null(raster$coordinate.system$unit)) raster$coordinate.system$unit <- "pixel"
  if (is.null(raster$coordinate.system$calibrated)) raster$coordinate.system$calibrated <- FALSE
  if (is.null(raster$coordinate.system$y.direction)) raster$coordinate.system$y.direction <- "down"
  raster$transform <- .sc_normalize_transform(raster$transform %||% diag(3L))
  raster
}

#' Read spatial image files for Seurat
#' @param image_dir Directory containing the spatial image and scale-factor files.
#' @param image_name Image file name or name of a stored image, as appropriate to the input.
#' @param filter_matrix Whether to filter positions when reading the image.
#' @return A Seurat spatial image object.
#' @noRd
readSpatialImage <- function(image_dir,
                             image_name = "tissue_lowres_image.png",
                             filter_matrix = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Seurat is required to read Visium images", call. = FALSE)
  visium_image <- Seurat::Read10X_Image(
    image.dir = image_dir,
    image.name = image_name,
    filter.matrix = filter_matrix
  )
  .sc_normalize_raster(visium_image)
}

#' Extract raster and calibration from Seurat
#' @param object A SpatialCellChat object.
#' @param image_name Image file name or name of a stored image, as appropriate to the input.
#' @param cell_names Observation names in canonical object order.
#' @return A list containing raster data and geometry, or NULL if unavailable.
#' @noRd
extractSeuratRasterImage <- function(object, image_name = NULL, cell_names = NULL) {
  if (!methods::is(object, "Seurat") || !requireNamespace("SeuratObject", quietly = TRUE))
    return(NULL)
  image_names <- SeuratObject::Images(object)
  if (!is.null(image_name)) {
    if (!is.character(image_name) || length(image_name) != 1L)
      stop("image_name must be a single image name", call. = FALSE)
    image_names <- image_name
  }
  if (!length(image_names)) return(NULL)
  if (!image_names[[1L]] %in% SeuratObject::Images(object))
    stop("image name was not found in the Seurat object", call. = FALSE)
  image <- object[[image_names[[1L]]]]
  tryCatch(
    .sc_normalize_raster(image, cell_names),
    error = function(e) {
      # A Seurat object may contain an FOV/segmentation image without a raster.
      # It is valid spatial data, so automatic raster import must not make the
      # whole constructor fail; an explicitly requested image still reports the
      # actionable error to the caller.
      if (!is.null(image_name)) stop(conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

#' Extract coordinates from Seurat spatial data
#' @param object A SpatialCellChat object.
#' @return A coordinate matrix or NULL when positions are unavailable.
#' @noRd
.sc_extract_seurat_coordinates <- function(object) {
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  coordinates <- tryCatch(
    Seurat::GetTissueCoordinates(object, scale = NULL),
    error = function(e) NULL
  )
  if (is.null(coordinates) && requireNamespace("SeuratObject", quietly = TRUE))
    coordinates <- tryCatch(
      SeuratObject::GetTissueCoordinates(object), error = function(e) NULL)
  if (!is.data.frame(coordinates) || !nrow(coordinates)) return(NULL)

  id_col <- intersect(c("cell", "barcode", "barcodes", "spot"), names(coordinates))
  if (length(id_col)) rownames(coordinates) <- as.character(coordinates[[id_col[[1L]]]])
  candidates <- list(
    c("imagecol", "imagerow"),
    c("x", "y"),
    c("col", "row")
  )
  selected <- candidates[vapply(candidates, function(x) all(x %in% names(coordinates)), logical(1))]
  if (!length(selected)) return(NULL)
  result <- as.matrix(coordinates[, selected[[1L]], drop = FALSE])
  if (!is.numeric(result) || any(!is.finite(result))) return(NULL)
  colnames(result) <- c("x", "y")
  result
}
#' Extract stored numeric matrix values
#' @param x The value to inspect or transform.
#' @return A vector of matrix values, using the explicit values of sparse matrices.
#' @noRd
.sc_matrix_values <- function(x) {
  if (isS4(x) && "x" %in% methods::slotNames(x)) methods::slot(x, "x") else as.vector(x)
}

#' Validate a SpatialCellChat object against the final 11-slot schema.
#' @param object A SpatialCellChat object.
#' @param strict If TRUE, stop on the first validation failure.
#' @noRd
#' @return TRUE when valid; otherwise a character vector when strict is FALSE.
#' @export
validateSpatialCellChat <- function(object, strict = TRUE) {
  if (!methods::is(object, "SpatialCellChat"))
    stop("object must be a SpatialCellChat", call. = FALSE)
  errors <- .sc_validate_object(object)
  if (length(errors) && isTRUE(strict))
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  if (length(errors)) errors else TRUE
}

methods::setValidity("SpatialCellChat", function(object) {
  errors <- .sc_validate_object(object)
  if (length(errors)) errors else TRUE
})


# ====== show method ======

#' Print a SpatialCellChat summary
#' @param object A SpatialCellChat object.
#' @return The result of writing the summary to the console, invisibly.
#' @docType methods
setMethod(f = "show", signature = "SpatialCellChat", definition = function(object) {
  cell_names <- colnames(object@assay$norm)
  n_cells <- length(cell_names)
  n_genes <- nrow(object@assay$norm)
  mode <- object@misc$.mode %||% "unknown"
  datatype <- object@misc$.datatype %||% "unknown"
  n_groups <- if (inherits(object@idents, "factor")) {
    nlevels(object@idents)
  } else {
    0L
  }
  n_lr <- if (is.data.frame(object@LR$LRsig)) nrow(object@LR$LRsig) else 0L
  has_net <- any(lengths(object@net) > 0L) || any(lengths(object@netP) > 0L)

  cat("SpatialCellChat\n")
  cat(sprintf("  mode: %s; datatype: %s\n", mode, datatype))
  cat(sprintf("  cells: %d; genes: %d; groups: %d\n", n_cells, n_genes, n_groups))
  cat(sprintf("  LR entries: %d\n", n_lr))
  cat(sprintf("  communication: %s\n", if (has_net) "computed" else "not computed"))
  cat(sprintf("  logged operations: %d\n", length(object@misc$.log)))
  invisible(NULL)
})


# ====== createSpatialCellChat ======

#' Create a new SpatialCellChat object
#'
#' @param object A raw count matrix or normalized expression matrix, Seurat
#'   object, or SingleCellExperiment object. The supplied assay is declared
#'   by \code{input.assay}.
#' @param meta A data frame of cell metadata. If input is a Seurat or SCE
#'   object, taken from the object by default.
#' @param group.by Column name in meta defining cell groups. If omitted for a
#'   Seurat object, its active identities are copied into metadata.
#' @param input.assay Either \code{"norm"} for already normalized expression
#'   values or \code{"raw"} for non-negative count values.
#' @param normalize Whether raw input should be normalized into \code{assay$norm}.
#'   Defaults to TRUE for \code{input.assay = "raw"} and FALSE otherwise.
#'   Contradictory assay/normalize combinations are rejected.
#' @param scale.factor Per-cell scaling factor used for raw input normalization.
#' @param do.log Whether raw input normalization applies log1p after scaling.
#' @param datatype One of \code{"RNA"} or \code{"spatial"}.
#' @param coordinates A two- or three-column matrix of analysis coordinates.
#' @param coordinate.system A list with \code{unit}, \code{calibrated}, and
#'   optional orientation metadata. Use \code{unit = "arbitrary"} and
#'   \code{calibrated = FALSE} when no physical scale is known.
#' @param spatial.factors A list with \code{ratio} and \code{tol}. Calibrated
#'   spatial data require both values; uncalibrated data store NULL values.
#' @param image Optional Seurat SpatialImage object or raster record. It is
#'   stored under \code{images$rasters}, never as a duplicate histology object.
#' @param image_name Optional name of a Seurat image to import.
#' @param assay_name Assay to use from a Seurat object.
#' @param do.sparse Whether to convert the expression layers to dgCMatrix.
#'
#' @param sample Sample label or labels used to populate spatial sample metadata.
#' @return A SpatialCellChat object.
#' @export
#' @importFrom methods new
#'
#' @examples
#' \dontrun{
#' chat <- createSpatialCellChat(data, meta = meta, group.by = "labels")
#' chat <- createSpatialCellChat(seu, group.by = "clusters",
#'   datatype = "spatial", coordinates = coords, spatial.factors = sf)
#' }
createSpatialCellChat <- function(object,
                                  meta = NULL,
                                  group.by = NULL,
                                  input.assay = c("norm", "raw"),
                                  normalize = NULL,
                                  scale.factor = 10000,
                                  do.log = TRUE,
                                  datatype = c("RNA", "spatial"),
                                  coordinates = NULL,
                                  coordinate.system = NULL,
                                  spatial.factors = NULL,
                                  image = NULL,
                                  image_name = NULL,
                                  assay_name = NULL,
                                  do.sparse = TRUE,
                                  sample = "sample1") {
  input.assay <- match.arg(input.assay)
  datatype <- match.arg(datatype)
  if (is.null(normalize)) normalize <- identical(input.assay, "raw")
  if (length(normalize) != 1L || !is.logical(normalize) || is.na(normalize))
    stop("normalize must be TRUE, FALSE, or NULL", call. = FALSE)
  if (identical(input.assay, "raw") && !isTRUE(normalize))
    stop("raw input requires normalize = TRUE", call. = FALSE)
  if (identical(input.assay, "norm") && isTRUE(normalize))
    stop("normalized input requires normalize = FALSE", call. = FALSE)
  if (length(sample) != 1L || is.na(sample) || !nzchar(sample))
    stop("sample must be a single non-empty string", call. = FALSE)
  group.by.explicit <- !is.null(group.by)
  is_seurat <- methods::is(object, "Seurat")
  image.was.auto <- FALSE
  .cli("createSpatialCellChat", .type = "subheader")
  .cli("Input: {.val {input.assay}} expression layer, {.val {datatype}} mode",
       .type = "info")

  # ---- Extract expression matrix from input ----
  data <- NULL
  if (inherits(object, c("matrix", "Matrix", "dgCMatrix"))) {
    data <- object
    # Matrix input historically accepted a `labels` column as a convenience.
    # Keep that path, but distinguish it from an explicitly misspelled group.by
    # so a typo is no longer silently converted into all-unknown labels.
    if (is.null(group.by) && !is.null(meta) &&
        is.data.frame(meta) && "labels" %in% colnames(meta))
      group.by <- "labels"
  } else if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat is required for Seurat objects", call. = FALSE)
    if (is.null(assay_name)) assay_name <- Seurat::DefaultAssay(object)
    if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name))
      stop("assay_name must be one Seurat assay name", call. = FALSE)
    if (identical(assay_name, "integrated"))
      warning("The 'integrated' assay is not suitable; use 'RNA' or 'SCT'")
    layer.name <- if (identical(input.assay, "raw")) "counts" else "data"
    data <- tryCatch(
      tryCatch(
        Seurat::GetAssayData(object, assay = assay_name, layer = layer.name),
        error = function(layer.error)
          Seurat::GetAssayData(object, assay = assay_name, slot = layer.name)
      ),
      error = function(e) stop(
        "could not extract ", input.assay, " expression layer '", layer.name,
        "' from assay '", assay_name, "': ", conditionMessage(e), call. = FALSE)
    )
    if (is.null(data))
      stop("assay '", assay_name, "' has no '", layer.name, "' expression layer",
           call. = FALSE)

    if (is.null(meta)) meta <- object@meta.data
    if (identical(datatype, "spatial") && is.null(coordinates))
      coordinates <- .sc_extract_seurat_coordinates(object)
    if (identical(datatype, "spatial") && is.null(image) &&
        requireNamespace("SeuratObject", quietly = TRUE)) {
      image_names <- SeuratObject::Images(object)
      if (!is.null(image_name)) {
        image <- image_name
      } else if (length(image_names)) {
        image <- image_names[[1L]]
        image.was.auto <- TRUE
      }
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
      stop("SingleCellExperiment is required for SCE objects", call. = FALSE)
    assay.name <- if (identical(input.assay, "raw")) "counts" else "logcounts"
    if (!(assay.name %in% SummarizedExperiment::assayNames(object)))
      stop("SCE object must contain an assay named '", assay.name,
           "' for input.assay = '", input.assay, "'", call. = FALSE)
    data <- if (identical(assay.name, "counts"))
      SingleCellExperiment::counts(object) else SingleCellExperiment::logcounts(object)
    if (is.null(meta)) meta <- as.data.frame(SingleCellExperiment::colData(object))
    if (is.null(group.by))
      stop("group.by must be defined for SCE input", call. = FALSE)
  } else {
    stop("object must be an expression matrix, Seurat, or SingleCellExperiment object",
         call. = FALSE)
  }

  if (!is.matrix(data) && !inherits(data, "Matrix"))
    stop("the extracted expression data must be matrix-like", call. = FALSE)
  if (length(dim(data)) != 2L || any(dim(data) < 1L))
    stop("expression data must have at least one gene and one cell", call. = FALSE)
  if (is.null(rownames(data)) || is.null(colnames(data)))
    stop("expression data must have gene rownames and cell colnames", call. = FALSE)
  if (anyDuplicated(rownames(data)) || anyDuplicated(colnames(data)))
    stop("expression gene and cell names must be unique", call. = FALSE)
  data_values <- .sc_matrix_values(data)
  if (!is.numeric(data_values) || any(!is.finite(data_values)))
    stop("expression data must contain only finite numeric values", call. = FALSE)
  if (identical(input.assay, "raw") && length(data_values) && any(data_values < 0))
    stop("raw expression data must contain non-negative counts", call. = FALSE)
  data.raw <- NULL
  if (isTRUE(normalize)) {
    data.raw <- data
    data <- normalizeData(data.raw, scale.factor = scale.factor, do.log = do.log,
                          verbose = FALSE)
  }
  if (do.sparse && !inherits(data, "dgCMatrix")) data <- methods::as(data, "dgCMatrix")
  if (!is.null(data.raw) && do.sparse && !inherits(data.raw, "dgCMatrix"))
    data.raw <- methods::as(data.raw, "dgCMatrix")

  cell_names <- colnames(data)
  gene_names <- rownames(data)
  if (is.null(meta)) meta <- data.frame(row.names = cell_names)
  if (inherits(meta, c("matrix", "Matrix"))) meta <- as.data.frame(meta)
  if (!is.data.frame(meta)) stop("meta must be a data frame", call. = FALSE)
  if (nrow(meta) != length(cell_names))
    stop("meta must have one row per cell", call. = FALSE)
  if (is.null(rownames(meta))) {
    # No rownames means the input supplies no alignment information; assigning
    # the expression order is safe only after the row count was checked above.
    rownames(meta) <- cell_names
  } else {
    if (anyDuplicated(rownames(meta)))
      stop("meta rownames must be unique", call. = FALSE)
    if (!setequal(rownames(meta), cell_names))
      stop("meta rownames must contain exactly the expression cell names", call. = FALSE)
    meta <- meta[cell_names, , drop = FALSE]
  }

  if (is_seurat && is.null(group.by)) {
    # Seurat's active identities are not guaranteed to be materialized as an
    # `ident` metadata column.  Copy them explicitly so construction remains
    # correct for objects whose metadata contain no cluster column.
    active.ident <- tryCatch(as.character(Seurat::Idents(object)), error = function(e) NULL)
    active.names <- tryCatch(names(Seurat::Idents(object)), error = function(e) NULL)
    if (!is.null(active.ident) && length(active.ident) == length(cell_names)) {
      if (!is.null(active.names) && setequal(active.names, cell_names))
        active.ident <- active.ident[match(cell_names, active.names)]
      meta$.ident <- active.ident
      group.by <- ".ident"
    }
  }
  if (is.null(group.by) && "labels" %in% colnames(meta)) group.by <- "labels"

  images <- .sc_default_images()
  if (identical(datatype, "spatial")) {
    if (is.null(coordinates))
      stop("coordinates required for spatial data; provide a two- or three-column matrix",
           call. = FALSE)
    coordinates <- as.matrix(coordinates)
    if (nrow(coordinates) != length(cell_names) || !ncol(coordinates) %in% c(2L, 3L))
      stop("coordinates must have one row per cell and exactly two or three columns", call. = FALSE)
    if (!is.numeric(coordinates) || any(!is.finite(coordinates)))
      stop("coordinates must contain finite numeric values", call. = FALSE)
    if (is.null(rownames(coordinates))) {
      rownames(coordinates) <- cell_names
    } else {
      if (anyDuplicated(rownames(coordinates)))
        stop("coordinates rownames must be unique", call. = FALSE)
      if (!setequal(rownames(coordinates), cell_names))
        stop("coordinates rownames must contain exactly the expression cell names", call. = FALSE)
      coordinates <- coordinates[cell_names, , drop = FALSE]
    }
    colnames(coordinates) <- c("x", "y", "z")[seq_len(ncol(coordinates))]

    if (is.null(coordinate.system) && is.null(spatial.factors))
      stop("provide spatial.factors for calibrated coordinates or coordinate.system with unit = 'arbitrary'",
           call. = FALSE)
    coordinate.system <- .sc_normalize_coordinate_system(coordinate.system, spatial.factors)
    spatial.factors <- .sc_normalize_spatial_factors(spatial.factors, coordinate.system)
    factor_errors <- .sc_validate_spatial_factors(spatial.factors, coordinate.system)
    if (length(factor_errors)) stop(paste(factor_errors, collapse = "\n"), call. = FALSE)
    images$coordinates <- coordinates
    images$coordinate.system <- coordinate.system
    images$spatial.factors <- spatial.factors

    if (is.character(image) && is_seurat) {
      raster <- extractSeuratRasterImage(
        object, if (isTRUE(image.was.auto)) NULL else image, cell_names)
    } else if (!is.null(image)) {
      raster <- .sc_normalize_raster(image, cell_names)
    } else {
      raster <- NULL
    }
    if (!is.null(raster)) {
      raster_name <- if (!is.null(image_name)) image_name else "main"
      images$rasters[[raster_name]] <- raster
    }
  }

  if (is.null(group.by)) group.by <- "ident"
  if (!group.by %in% colnames(meta)) {
    if (group.by.explicit)
      stop("group.by column '" , group.by, "' was not found in meta", call. = FALSE)
    # A matrix without labels is still useful for structural work, but it must
    # be explicit that all cells share a placeholder group rather than silently
    # inheriting an arbitrary metadata column.
    meta[[group.by]] <- factor(rep("unknown", length(cell_names)), levels = "unknown")
  }
  labels <- meta[[group.by]]
  if (length(labels) != length(cell_names) || anyNA(labels) ||
      any(!nzchar(trimws(as.character(labels)))))
    stop("group.by labels must have one non-empty, non-missing value per cell", call. = FALSE)
  joint <- factor(as.character(labels))
  names(joint) <- cell_names

  misc <- list(
    .mode = "single",
    .datatype = datatype,
    .param = list(
      sample = sample,
      group.by = group.by,
      input.assay = input.assay,
      normalize = isTRUE(normalize),
      scale.factor = if (isTRUE(normalize)) as.numeric(scale.factor) else NULL,
      do.log = if (isTRUE(normalize)) isTRUE(do.log) else NULL,
      coordinate.system = coordinate.system,
      spatial.factors = spatial.factors
    ),
    .log = list(),
    .var.features = list(),
    .datasets = list()
  )
  chat <- methods::new(Class = "SpatialCellChat",
    assay = list(raw = data.raw, norm = data, scale = NULL,
                 smooth = NULL, signaling = NULL),
    images = images,
    meta = meta,
    idents = joint,
    features = data.frame(row.names = gene_names),
    LR = list(),
    dr = list(),
    net = list(cell = list(), group = list()),
    netP = list(cell = list(), group = list()),
    DB = list(),
    misc = misc)

  # `new()` invokes the S4 validity contract.  Keeping this check at the
  # construction boundary catches malformed dimensions before downstream
  # analysis allocates communication arrays.
  methods::validObject(chat)
  .cli("Object created: {nrow(data)} genes x {ncol(data)} cells; {nlevels(joint)} groups",
       .type = "success")
  .log_operation(chat, "createSpatialCellChat", params = list(
    datatype = datatype,
    group.by = group.by,
    input.assay = input.assay,
    normalize = isTRUE(normalize),
    scale.factor = if (isTRUE(normalize)) as.numeric(scale.factor) else NULL,
    do.log = if (isTRUE(normalize)) isTRUE(do.log) else NULL,
    sample = sample
  ))
}



# ====== idents accessor ======

#' Get/set cell group labels
#'
#' @param chat A SpatialCellChat object
#' @noRd
#' @return A factor of cell group labels
#' @export
idents <- function(chat) {
  if (!methods::is(chat, "SpatialCellChat"))
    stop("chat must be a SpatialCellChat", call. = FALSE)
  chat@idents
}

#' @param chat A SpatialCellChat object.
#' @param value The input or replacement value.
#' @return The object with the supplied identity factor.
#' @noRd
#' @export
"idents<-" <- function(chat, value) {
  if (!inherits(value, "factor") || length(value) != nrow(chat@meta))
    stop("idents must be a factor with one value per cell", call. = FALSE)
  if (is.null(names(value))) names(value) <- rownames(chat@meta)
  if (!identical(names(value), rownames(chat@meta)))
    stop("idents names must match meta rownames", call. = FALSE)
  chat@idents <- value
  methods::validObject(chat)
  chat
}


#' List cell identity levels
#'
#' @param chat A SpatialCellChat object.
#' @return A character vector of group levels, or a list of vectors for list-valued identities.
#' @export
level.SpatialCellChat <- function(chat) {
  lev <- levels(idents(chat))
  cat("Cell groups:", cli::col_red(paste(lev, collapse = ", ")), "\n")
  invisible(lev)
}


# ====== Phase 2 accessors ======

#' Check the accessor input class
#' @param object A SpatialCellChat object.
#' @return The validated SpatialCellChat object.
#' @noRd
.sc_assert_spatial_cell_chat <- function(object) {
  if (!methods::is(object, "SpatialCellChat"))
    stop("object must be a SpatialCellChat", call. = FALSE)
  object
}

#' Validate a SpatialCellChat replacement
#' @param object A SpatialCellChat object.
#' @return The updated SpatialCellChat object, after schema validation.
#' @noRd
.sc_validate_after_update <- function(object) {
  methods::validObject(object)
  object
}

#' Access expression layers.
#' @param object A SpatialCellChat object.
#' @param layer NULL for all layers, or one of raw, norm, scale, smooth, signaling.
#' @return The assay list or a single assay matrix.
#' @noRd
#' @export
assay <- function(object, layer = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(layer)) return(object@assay)
  if (length(layer) != 1L || !is.character(layer) ||
      !layer %in% names(object@assay))
    stop("layer must name one assay layer", call. = FALSE)
  object@assay[[layer]]
}

#' Replace expression layers.
#' @param object A SpatialCellChat object.
#' @param layer An assay layer name or the layer value used by the caller.
#' @param value An assay matrix when layer is set, or a named assay list.
#' @return The object with updated assay data.
#' @noRd
#' @export
`assay<-` <- function(object, layer = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(layer)) {
    if (!is.list(value)) stop("assay must be a list", call. = FALSE)
    object@assay <- value
  } else {
    if (length(layer) != 1L || !is.character(layer) ||
        !layer %in% names(object@assay))
      stop("layer must name one assay layer", call. = FALSE)
    object@assay[layer] <- list(value)
  }
  .sc_validate_after_update(object)
}

#' Access cell metadata.
#' @param columns NULL for all metadata, or character column names.
#' @param object A SpatialCellChat object.
#' @return The metadata data frame, optionally restricted to selected columns.
#' @noRd
#' @export
meta <- function(object, columns = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(columns)) return(object@meta)
  if (!is.character(columns) || any(!columns %in% colnames(object@meta)))
    stop("columns must refer to existing metadata columns", call. = FALSE)
  object@meta[, columns, drop = FALSE]
}

#' Replace cell metadata.
#' @param object A SpatialCellChat object.
#' @param value A data frame with one aligned row per observation.
#' @return The object with updated metadata.
#' @noRd
#' @export
`meta<-` <- function(object, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (!is.data.frame(value)) stop("meta must be a data.frame", call. = FALSE)
  object@meta <- value
  .sc_validate_after_update(object)
}

#' Access spatial data and caches.
#' @param key NULL for all image data, or a name inside the images list.
#' @param object A SpatialCellChat object.
#' @return The spatial data list or the named entry.
#' @noRd
#' @export
images <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@images)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@images))
    stop("key must name an existing images entry", call. = FALSE)
  object@images[[key]]
}

#' Replace spatial data and caches.
#' @param object A SpatialCellChat object.
#' @param key Name of a stored entry, or NULL to use the entire list.
#' @param value The input or replacement value.
#' @return The object with updated spatial data.
#' @noRd
#' @export
`images<-` <- function(object, key = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("images must be a list", call. = FALSE)
    invalidate.distance <- !identical(object@images$coordinates, value$coordinates) ||
      !identical(object@images$spatial.factors, value$spatial.factors)
    object@images <- value
  } else {
    if (length(key) != 1L || !is.character(key))
      stop("key must be a single character name", call. = FALSE)
    previous.value <- object@images[[key]]
    object@images[[key]] <- value
    invalidate.distance <- key %in% c("coordinates", "spatial.factors") &&
      !identical(previous.value, value)
  }
  if (invalidate.distance) object@images$.distance <- list()
  .sc_validate_after_update(object)
}

#' Access ligand-receptor metadata.
#' @param key NULL for the complete LR list, or a name inside it.
#' @param object A SpatialCellChat object.
#' @return The LR annotation list or the named entry.
#' @noRd
#' @export
LR <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@LR)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@LR))
    stop("key must name an existing LR entry", call. = FALSE)
  object@LR[[key]]
}

#' Replace ligand-receptor metadata.
#' @param object A SpatialCellChat object.
#' @param key Name of a stored entry, or NULL to use the entire list.
#' @param value The input or replacement value.
#' @return The object with updated LR annotations.
#' @noRd
#' @export
`LR<-` <- function(object, key = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("LR must be a list", call. = FALSE)
    object@LR <- value
  } else {
    if (length(key) != 1L || !is.character(key))
      stop("key must be a single character name", call. = FALSE)
    object@LR[[key]] <- value
  }
  .sc_validate_after_update(object)
}

#' Access the ligand-receptor database.
#' @param key NULL for the complete DB list, or a name inside it.
#' @param object A SpatialCellChat object.
#' @return The communication database or the named entry.
#' @noRd
#' @export
DB <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@DB)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@DB))
    stop("key must name an existing DB entry", call. = FALSE)
  object@DB[[key]]
}

#' Replace the ligand-receptor database.
#' @param object A SpatialCellChat object.
#' @param key Name of a stored entry, or NULL to use the entire list.
#' @param value The input or replacement value.
#' @return The object with the updated database.
#' @noRd
#' @export
`DB<-` <- function(object, key = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("DB must be a list", call. = FALSE)
    object@DB <- value
  } else {
    if (length(key) != 1L || !is.character(key))
      stop("key must be a single character name", call. = FALSE)
    object@DB[[key]] <- value
  }
  .sc_validate_after_update(object)
}

#' Access a communication result.
#' @param slot.name Either net (LR-level) or netP (pathway-level).
#' @param resolution Either cell or group, or NULL for the whole slot.
#' @param measure Optional result name such as prob, pval, count, weight, centr, or field.
#' @param object A SpatialCellChat object.
#' @return The selected network branch or entry.
#' @noRd
#' @export
communication <- function(object, slot.name = c("net", "netP"),
                           resolution = NULL, measure = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  slot.name <- match.arg(slot.name)
  value <- methods::slot(object, slot.name)
  if (is.null(resolution)) {
    if (!is.null(measure)) stop("measure requires resolution", call. = FALSE)
    return(value)
  }
  resolution <- match.arg(resolution, c("cell", "group"))
  value <- value[[resolution]]
  if (is.null(measure)) return(value)
  if (length(measure) != 1L || !is.character(measure) ||
      !measure %in% names(value))
    stop("measure must name an existing communication result", call. = FALSE)
  value[[measure]]
}

#' Replace a communication result.
#' @param object A SpatialCellChat object.
#' @param slot.name Source slot: \code{"assay"}, \code{"net"} or \code{"netP"}.
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param value The input or replacement value.
#' @return The object with the updated network branch.
#' @noRd
#' @export
`communication<-` <- function(object, slot.name = c("net", "netP"),
                               resolution = NULL, measure = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  slot.name <- match.arg(slot.name)
  if (is.null(resolution)) {
    if (!is.null(measure) || !is.list(value))
      stop("a whole communication slot must be a list", call. = FALSE)
    methods::slot(object, slot.name) <- value
  } else {
    resolution <- match.arg(resolution, c("cell", "group"))
    if (is.null(measure)) {
      if (!is.list(value)) stop("a communication branch must be a list", call. = FALSE)
      branch <- methods::slot(object, slot.name)
      branch[[resolution]] <- value
      methods::slot(object, slot.name) <- branch
    } else {
      if (length(measure) != 1L || !is.character(measure) || !nzchar(measure))
        stop("measure must be a single non-empty name", call. = FALSE)
      branch <- methods::slot(object, slot.name)
      branch[[resolution]][[measure]] <- value
      methods::slot(object, slot.name) <- branch
    }
  }
  .sc_validate_after_update(object)
}

#' Access the misc list or analysis parameters.
#' @param key NULL for misc, or a key inside misc$.param.
#' @param object A SpatialCellChat object.
#' @return The parameter list or the named entry.
#' @noRd
#' @export
params <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  value <- object@misc$.param
  if (is.null(key)) return(value)
  if (length(key) != 1L || !is.character(key) || !key %in% names(value))
    stop("key must name an existing parameter", call. = FALSE)
  value[[key]]
}

#' Replace analysis parameters.
#' @param object A SpatialCellChat object.
#' @param key Name of a stored entry, or NULL to use the entire list.
#' @param value The input or replacement value.
#' @return The object with the updated parameter list.
#' @noRd
#' @export
`params<-` <- function(object, key = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("params must be a list", call. = FALSE)
    object@misc$.param <- value
  } else {
    if (length(key) != 1L || !is.character(key) || !nzchar(key))
      stop("key must be a single non-empty name", call. = FALSE)
    object@misc$.param[[key]] <- value
  }
  .sc_validate_after_update(object)
}

# ====== misc accessor (Seurat-compatible) ======

#' Get and set miscellaneous data
#'
#' @param object An object
#' @param ... Arguments passed to the getter or replacement method, including key.
#' @return The getter returns the miscellaneous data list or one requested entry. The setter returns the updated object.
#' @export misc
#' @concept data-access
misc <- function(object, ...) {
  UseMethod(generic = "misc", object = object)
}

#' @param value Data to add
#' @rdname misc
#' @export misc<-
"misc<-" <- function(object, ..., value) {
  UseMethod(generic = "misc<-", object = object)
}

#' @param key Optional name of the miscellaneous entry to retrieve or replace.
#' @rdname misc
#' @export
#' @method misc SpatialCellChat
misc.SpatialCellChat <- function(object, key = NULL, ...) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (!is.null(key)) object@misc[[key]] else object@misc
}

#' @rdname misc
#' @export
#' @method misc<- SpatialCellChat
"misc<-.SpatialCellChat" <- function(object, key = NULL, ..., value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("Use a named list")
    object@misc <- value
  } else {
    object@misc[[key]] <- value
  }
  .sc_validate_after_update(object)
}


#' Use a fallback for NULL
#' @param a The value to use if it is non-NULL.
#' @param b The fallback value.
#' @return The first argument when present, otherwise the fallback.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
