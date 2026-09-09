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
#' @return an object of class "SparseChatArray"
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

#' @export
`dim.SparseChatArray` <- function(x) {
  c(nrow(unclass(x)[[1]]), ncol(unclass(x)[[1]]), length(x))
}

#' @export
`dimnames.SparseChatArray` <- function(x) {
  attr(x, "Dimnames", exact = TRUE) %||%
    .SparseChatArray_dimnames(unclass(x), dim(x), NULL)
}

#' @export
`dimnames<-.SparseChatArray` <- function(x, value) {
  dn <- .SparseChatArray_dimnames(unclass(x), dim(x), value)
  .new_SparseChatArray(unclass(x), dn)
}

#' @export
`length.SparseChatArray` <- function(x) length(unclass(x))

#' @export
`names.SparseChatArray` <- function(x) dimnames(x)[[3]]

#' @export
`names<-.SparseChatArray` <- function(x, value) {
  dn <- dimnames(x)
  dn[[3]] <- .SparseChatArray_validate_dimname(value, dim(x)[3], "names")
  .new_SparseChatArray(unclass(x), dn)
}

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

#' @export
#' Sum margins of a SparseChatArray
#' @param x a SparseChatArray
#' @param margin 1=row sums per layer, 2=col sums per layer, 3=sum across all layers
#' @export
marginSums <- function(x, margin = 3L) {
  UseMethod("marginSums")
}
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

#' @export
`t.SparseChatArray` <- function(x) {
  out <- pbapply::pblapply(unclass(x), t)
  dn <- dimnames(x)
  .new_SparseChatArray(out, list(dn[[2]], dn[[1]], dn[[3]]))
}

.new_SparseChatArray <- function(x, dimnames) {
  x <- lapply(x, .SparseChatArray_drop_matrix_dimnames)
  names(x) <- NULL
  structure(x, class = "SparseChatArray", Dimnames = dimnames)
}

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

.SparseChatArray_validate_dimname <- function(value, expected, label) {
  if (is.null(value)) return(NULL)
  if (!is.character(value)) value <- as.character(value)
  if (length(value) != expected)
    stop(label, " must be NULL or have length ", expected)
  value
}

.SparseChatArray_drop_matrix_dimnames <- function(mat) {
  dimnames(mat) <- list(NULL, NULL)
  mat
}

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

.SparseChatArray_subset_dimname <- function(dimname, index, missing_index) {
  if (is.null(dimname)) return(NULL)
  if (missing_index) return(dimname)
  dimname[index]
}

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

.sc_matrix_like <- function(x) {
  inherits(x, "matrix") || inherits(x, "Matrix")
}

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

.sc_validate_assay <- function(assay, gene_names, cell_names) {
  if (!is.list(assay)) return("assay must be a list")
  required <- c("raw", "norm", "scale", "smooth", "signaling")
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
  for (name in intersect(required, names(assay))) {
    errors <- c(errors, .sc_validate_matrix(
      assay[[name]], expected, expected_names, paste0("assay$", name)))
  }
  errors
}

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

.sc_validate_images <- function(images, cell_names, datatype) {
  if (!is.list(images)) return("images must be a list")
  errors <- character()
  coordinates <- images$coordinates
  coordinate.system <- images$coordinate.system

  if (identical(datatype, "spatial")) {
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
  if ("histology" %in% names(images))
    errors <- c(errors, "images$histology is obsolete; store the record in images$rasters")

  errors <- c(errors, .sc_validate_rasters(images$rasters, cell_names))
  errors <- c(errors, .sc_validate_fov(images$fov))
  errors <- c(errors, .sc_validate_distance_cache(images$.distance, cell_names))
  if (!is.null(images$.grid) && !is.list(images$.grid))
    errors <- c(errors, "images$.grid must be a list")
  unique(errors[nzchar(errors)])
}


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
  }
  errors
}

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
  }
  unique(errors[nzchar(errors)])
}
# ====== spatial image helpers ======

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
.sc_matrix_values <- function(x) {
  if (isS4(x) && "x" %in% methods::slotNames(x)) methods::slot(x, "x") else as.vector(x)
}

#' Validate a SpatialCellChat object against the final 11-slot schema.
#' @param object A SpatialCellChat object.
#' @param strict If TRUE, stop on the first validation failure.
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

#' @param object A SpatialCellChat object.
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
#' @return A factor of cell group labels
#' @export
idents <- function(chat) {
  if (!methods::is(chat, "SpatialCellChat"))
    stop("chat must be a SpatialCellChat", call. = FALSE)
  chat@idents
}

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

#' @export
setIdent <- function(chat, ident.use = "ident") {
  if (!ident.use %in% colnames(chat@meta))
    stop("Column '", ident.use, "' not found in meta", call. = FALSE)
  value <- factor(chat@meta[[ident.use]])
  names(value) <- rownames(chat@meta)
  chat@idents <- value
  methods::validObject(chat)
  chat
}

#' @export
level.SpatialCellChat <- function(chat) {
  lev <- levels(idents(chat))
  cat("Cell groups:", cli::col_red(paste(lev, collapse = ", ")), "\n")
  invisible(lev)
}


# ====== Phase 2 accessors ======

.sc_assert_spatial_cell_chat <- function(object) {
  if (!methods::is(object, "SpatialCellChat"))
    stop("object must be a SpatialCellChat", call. = FALSE)
  object
}

.sc_validate_after_update <- function(object) {
  methods::validObject(object)
  object
}

#' Access expression layers.
#' @param object A SpatialCellChat object.
#' @param layer NULL for all layers, or one of raw, norm, scale, smooth, signaling.
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
#' @export
meta <- function(object, columns = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(columns)) return(object@meta)
  if (!is.character(columns) || any(!columns %in% colnames(object@meta)))
    stop("columns must refer to existing metadata columns", call. = FALSE)
  object@meta[, columns, drop = FALSE]
}

#' Replace cell metadata.
#' @export
`meta<-` <- function(object, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (!is.data.frame(value)) stop("meta must be a data.frame", call. = FALSE)
  object@meta <- value
  .sc_validate_after_update(object)
}

#' Access spatial data and caches.
#' @param key NULL for all image data, or a name inside the images list.
#' @export
images <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@images)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@images))
    stop("key must name an existing images entry", call. = FALSE)
  object@images[[key]]
}

#' Replace spatial data and caches.
#' @export
`images<-` <- function(object, key = NULL, value) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) {
    if (!is.list(value)) stop("images must be a list", call. = FALSE)
    object@images <- value
  } else {
    if (length(key) != 1L || !is.character(key))
      stop("key must be a single character name", call. = FALSE)
    object@images[[key]] <- value
  }
  .sc_validate_after_update(object)
}

#' Access ligand-receptor metadata.
#' @param key NULL for the complete LR list, or a name inside it.
#' @export
LR <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@LR)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@LR))
    stop("key must name an existing LR entry", call. = FALSE)
  object@LR[[key]]
}

#' Replace ligand-receptor metadata.
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
#' @export
DB <- function(object, key = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  if (is.null(key)) return(object@DB)
  if (length(key) != 1L || !is.character(key) || !key %in% names(object@DB))
    stop("key must name an existing DB entry", call. = FALSE)
  object@DB[[key]]
}

#' Replace the ligand-receptor database.
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

#' @rdname misc
#' @export misc
misc <- function(object, ...) {
  UseMethod(generic = "misc", object = object)
}

#' @rdname misc
#' @export misc<-
"misc<-" <- function(object, ..., value) {
  UseMethod(generic = "misc<-", object = object)
}

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


# ====== Internal helpers ======

.cli <- function(text, .type = "info", ..., .env = parent.frame()) {
  switch(.type,
    info      = cli::cli_alert_info(text, .envir = .env, ...),
    success   = cli::cli_alert_success(text, .envir = .env, ...),
    danger    = cli::cli_alert_danger(text, .envir = .env, ...),
    warning   = cli::cli_alert_warning(text, .envir = .env, ...),
    bullet    = cli::cli_bullets(c("*" = text), .envir = .env, ...),
    header    = cli::cli_h1(text, .envir = .env, ...),
    subheader = cli::cli_h2(text, .envir = .env, ...),
    text      = cli::cli_text(text, .envir = .env, ...),
    cli::cli_alert_info(text, .envir = .env, ...)
  )
}

.log_operation <- function(object, funcname, params = list()) {
  entry <- list(
    "function" = funcname,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    params = params
  )
  object@misc$.log <- c(object@misc$.log, list(entry))
  object
}

`%||%` <- function(a, b) if (is.null(a)) b else a
