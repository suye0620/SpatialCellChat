#' Plot cell-to-cell distances behind the inferred communication
#'
#' The single entry point for spatial distance distributions: `mode`
#' selects the records, `split.by` decides how they are cut into panels.
#'
#' @details
#' `mode` selects the records, `split.by` decides how they are cut into
#' panels. The full 2 x 2 is:
#'
#' \tabular{lllll}{
#'   # \tab split.by \tab mode \tab panel 1 \tab panel 2 \cr
#'   1 \tab "annotation" \tab "communication" \tab diffusible signaling \tab contact-dependent signaling \cr
#'   2 \tab "distance" \tab "communication" \tab all inferred cell pairs \tab those inside contact.range \cr
#'   3 \tab "annotation" \tab "geometry" \tab all cell pairs \tab contact-range pairs \cr
#'   4 \tab "distance" \tab "geometry" \tab all cell pairs \tab contact-range pairs \cr
#' }
#'
#' - Rows 1 and 2 are **disjoint** and **nested** respectively: the annotation
#'   split assigns every record to exactly one panel, the distance split makes
#'   panel 2 a subset of panel 1. The same cell pair may appear in more than one
#'   annotation panel when it has positive LR layers in more than one class.
#' - Rows 3 and 4 are **the same output**. Without the communication array there
#'   are no layers to classify, so `split.by` is ignored and a warning is raised
#'   when it was supplied explicitly.
#' - Rows 1 and 2 count **undirected cell-pair observations**: a pair with
#'   multiple positive LR layers contributes once per annotation panel. Rows 3
#'   and 4 count each geometric pair once.
#'
#' `mode = "communication"` uses positive values in the sparse
#' `net$cell$prob` layers, without testing statistical significance. A pair is
#' observed once per panel when either communication direction has at least one
#' positive LR layer. In particular, `computeCommunProb()` with `n = 1` produces no p-values at all, and
#' thresholding is a later step (`filterProbability()` or `filterCommunication()`).
#' Both modes use the strict upper triangle and exclude diagonal entries, so a
#' geometric pair is not counted twice and communication directions are merged
#' into one undirected pair before plotting.
#'
#' In the annotation split, layers are classified one by one from
#' `LR$LRsig$annotation` (not by the positional `nLR1` threshold of the old
#' implementation). Unrecognised annotations, and records whose layer cannot be
#' matched to `LR$LRsig`, fall into an "Other signaling" panel instead of being
#' dropped silently.
#'
#' Distances come from `images(object, ".distance")`. With `mode = "communication"`
#' that cache must already exist and match the coordinates and settings recorded
#' by `computeCommunProb()`; rerun inference when they change. With
#' `mode = "geometry"` a missing or stale cache is recomputed from
#' `images(object, "coordinates")` and `images(object, "spatial.factors")`, using
#' `interaction.range`, `contact.range` and `tol` when supplied — so the geometric
#' view works on an object that has never been through the inference. The x axis
#' always ends at the cached `interaction.range` plus 10.
#'
#' Self-pairs (distance 0) are dropped from every panel.
#'
#' @param object A SpatialCellChat object.
#' @param mode `"geometry"` (default) uses cell pairs from the distance cache;
#'   `"communication"` uses cell pairs with at least one positive
#'   ligand-receptor probability.
#' @param split.by How to cut the panels: `"annotation"` (diffusible vs
#'   contact-dependent) or `"distance"` (all records vs contact-range records).
#'   Ignored when `mode = "geometry"`.
#' @param density.alpha Transparency of the density fill.
#' @param interaction.range,contact.range,tol Distance settings used only when a
#'   cache has to be computed for `mode = "geometry"`. `NULL` falls back to
#'   250, 10 and the object's `spatial.factors$tol`.
#'
#' @return A `patchwork` object with one panel per group.
#' @export
spatialC2CDistPlot <- function(
    object,
    mode = c("geometry", "communication"),
    split.by = c("annotation", "distance"),
    density.alpha = 0.5,
    interaction.range = NULL,
    contact.range = NULL,
    tol = NULL
) {
  split.by.given <- !missing(split.by)
  mode <- match.arg(mode)
  split.by <- match.arg(split.by)
  if (split.by.given && identical(mode, "geometry")) {
    warning("`split.by` only applies when `mode = \"communication\"`; it is ignored here",
            call. = FALSE)
  }
  panels <- .sc_c2c_panels(object, split.by = split.by,
                            mode = mode,
                            density.alpha = density.alpha,
                            interaction.range = interaction.range,
                            contact.range = contact.range, tol = tol)
  patchwork::wrap_plots(panels, nrow = length(panels), ncol = 1)
}

#' Distance cache for the communication view: it must match the inference
#' @noRd
.sc_c2c_inference_cache <- function(object) {
  cache <- images(object, ".distance")
  if (is.null(cache) || is.null(cache$d.spatial))
    stop("`images(object, \".distance\")` is empty; run `computeCommunProb()` first, ",
         "or use `mode = \"geometry\"` for the geometric view", call. = FALSE)
  inference <- object@misc$.param$communication
  params <- inference[["distance.parameters"]]
  coordinates <- images(object, "coordinates")
  if (!is.list(inference) || !is.list(params) ||
      is.null(inference[["distance.coordinates"]]) ||
      !identical(inference[["distance.coordinates"]], cache[[".coordinates"]]) ||
      !identical(inference[["distance.spatial.factors"]],
                 images(object, "spatial.factors")) ||
      !.spatial_distance_cache_matches(
        cache,
        coordinates = coordinates,
        interaction.range = params[["interaction.range"]],
        contact.range = params[["contact.range"]],
        ratio = params[["ratio"]],
        tol = params[["tol"]]
      )) {
    stop("The distance cache or inferred probabilities are stale for the current spatial data; ",
         "rerun `computeCommunProb()` before using `mode = \"communication\"`",
         call. = FALSE)
  }
  cache
}

#' Distance cache for the geometric view: reuse when valid, compute otherwise
#' @noRd
.sc_c2c_geometry_cache <- function(object, interaction.range = NULL,
                                   contact.range = NULL, tol = NULL) {
  coordinates <- images(object, "coordinates")
  if (is.null(coordinates))
    stop("`images(object, \"coordinates\")` is empty; the geometric view needs coordinates",
         call. = FALSE)
  factors <- images(object, "spatial.factors")
  ratio <- if (is.list(factors)) factors$ratio else NULL
  if (is.null(tol) && is.list(factors)) tol <- factors$tol
  if (is.null(interaction.range)) interaction.range <- 250
  if (is.null(contact.range)) contact.range <- 10

  cache <- images(object, ".distance")
  if (.spatial_distance_cache_matches(cache, coordinates = coordinates,
                                      interaction.range = interaction.range,
                                      contact.range = contact.range, ratio = ratio,
                                      tol = tol)) {
    return(cache)
  }
  computeCellDistance(coordinates = coordinates, ratio = ratio,
                      interaction.range = interaction.range,
                      contact.range = contact.range, tol = tol)
}

#' X axis limit: the interaction range the distances were built with
#' @noRd
.sc_c2c_xmax <- function(cache) {
  params <- cache[[".parameters"]]
  range <- if (is.list(params)) params$interaction.range else NULL
  if (is.null(range) || !is.finite(range)) range <- 250
  range + 10
}

#' Build an upper-triangular distance frame from a pair mask
#' @noRd
.sc_pair_distance_frame <- function(pair.mask, cache, class = NULL) {
  d.spatial <- cache$d.spatial
  if (!identical(dim(pair.mask), dim(d.spatial)) ||
      !identical(dimnames(pair.mask), dimnames(d.spatial))) {
    stop("communication and distance matrices must have identical cell axes",
         call. = FALSE)
  }
  pair.mask <- Matrix::triu(Matrix::drop0(pair.mask), k = 1L)
  pair.distance <- Matrix::drop0(d.spatial * pair.mask)
  if (!length(pair.distance@x)) {
    out <- data.frame(distance = numeric(), contact = logical(),
                      stringsAsFactors = FALSE)
  } else {
    s <- summary(pair.distance)
    out <- data.frame(
      distance = s$x,
      contact = as.logical(cache$adj.contact[cbind(s$i, s$j)]),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(class)) out$class <- rep.int(class, nrow(out))
  out
}

#' Symmetric positive-probability mask for selected LR layers
#' @noRd
.sc_comm_pair_mask <- function(prob, layer.index = NULL) {
  if (is.null(layer.index)) {
    sums <- marginSums(prob, margin = 3L)
  } else {
    if (!length(layer.index)) return(NULL)
    sums <- marginSums(prob[, , layer.index, drop = FALSE], margin = 3L)
  }
  positive <- Matrix::drop0(sums > 0)
  Matrix::drop0(Matrix::triu(positive | Matrix::t(positive), k = 1L))
}

#' Cell-pair distances behind the plots
#'
#' Geometry uses one strict-upper-triangle distance per spatial pair.
#' Communication uses one row per undirected cell pair and annotation class.
#' @noRd
.sc_comm_distance_frame <- function(object, mode, cache, split.by = "distance") {
  if (identical(mode, "geometry")) {
    return(.sc_pair_distance_frame(
      Matrix::triu(cache$d.spatial, k = 1L), cache
    ))
  }
  if (!identical(mode, "communication"))
    stop("mode must be \"geometry\" or \"communication\"", call. = FALSE)

  prob <- object@net$cell$prob
  if (is.null(prob))
    stop("`net$cell$prob` is empty (spill mode keeps the layers on disk); ",
         "the communication view needs the in-memory probability array", call. = FALSE)

  layer.names <- dimnames(prob)[[3]]
  lr <- object@LR$LRsig
  ann.per.layer <- if (is.data.frame(lr) &&
                       all(c("interaction_name", "annotation") %in% colnames(lr))) {
    lr$annotation[match(layer.names, lr$interaction_name)]
  } else {
    rep(NA_character_, length(layer.names))
  }
  class.per.layer <- .sc_comm_annotation_class(ann.per.layer)

  if (identical(split.by, "distance")) {
    return(.sc_pair_distance_frame(.sc_comm_pair_mask(prob), cache))
  }

  class.order <- c("diffusible", "contact", "other")
  frames <- lapply(class.order, function(class) {
    idx <- which(class.per.layer == class)
    if (!length(idx)) return(NULL)
    .sc_pair_distance_frame(
      .sc_comm_pair_mask(prob, idx), cache, class = class
    )
  })
  frames <- Filter(Negate(is.null), frames)
  if (!length(frames)) {
    return(data.frame(distance = numeric(), contact = logical(),
                      class = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, frames)
}

#' Ligand-receptor annotation class used by the annotation split
#' @noRd
.sc_comm_annotation_class <- function(annotation) {
  diffusible <- c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")
  out <- rep("other", length(annotation))
  out[!is.na(annotation) & annotation %in% diffusible] <- "diffusible"
  out[!is.na(annotation) & annotation == "Cell-Cell Contact"] <- "contact"
  out
}

#' One distance-density panel
#' @noRd
.sc_comm_dist_panel <- function(distance, title, fill, xmax, density.alpha) {
  distance <- distance[!is.na(distance) & distance > 0]
  plot <- ggplot2::ggplot(data.frame(distance = distance)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0)),
                                breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, by = 20)) +
    ggplot2::labs(title = title, x = "Distance between cell pairs (um)", y = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
                   text = ggplot2::element_text(size = 10))
  if (length(distance) >= 2L && any(is.finite(distance))) {
    plot + ggplot2::geom_density(
      ggplot2::aes(x = distance, y = ggplot2::after_stat(density)),
      fill = fill, alpha = density.alpha, na.rm = TRUE)
  } else {
    plot + ggplot2::geom_blank()
  }
}

#' Build the panels shared by every distance plot
#' @noRd
.sc_comm_dist_panels <- function(object, split.by, mode, density.alpha, cache) {
  df <- .sc_comm_distance_frame(object, mode, cache, split.by = split.by)
  xmax <- .sc_c2c_xmax(cache)
  panels <- list()

  if (identical(mode, "geometry")) {
    # Rows 3 and 4 of the contract: geometry only, split.by is irrelevant.
    panels$all <- .sc_comm_dist_panel(df$distance, "All cell pairs",
                                      "#69b3a2", xmax, density.alpha)
    panels$contact.range <- .sc_comm_dist_panel(df$distance[df$contact],
                                                "Contact-range pairs",
                                                "#404080", xmax, density.alpha)
    return(panels)
  }

  if (identical(split.by, "distance")) {
    panels$all <- .sc_comm_dist_panel(df$distance, "All inferred pairs",
                                      "#69b3a2", xmax, density.alpha)
    if (any(df$contact)) {
      panels$contact.range <- .sc_comm_dist_panel(df$distance[df$contact],
                                                  "Contact-range pairs",
                                                  "#404080", xmax, density.alpha)
    }
    return(panels)
  }

  # Class-specific frames were deduplicated before this point; here we only split.
  has.diffusible <- any(df$class == "diffusible")
  has.contact <- any(df$class == "contact")
  if (has.diffusible) {
    panels$diffusible <- .sc_comm_dist_panel(df$distance[df$class == "diffusible"],
                                             "Diffusible signaling",
                                             "#69b3a2", xmax, density.alpha)
  }
  if (has.contact) {
    panels$contact <- .sc_comm_dist_panel(df$distance[df$class == "contact"],
                                          "Contact-dependent signaling",
                                          "#404080", xmax, density.alpha)
  }
  if ((has.diffusible || has.contact) && any(df$class == "other")) {
    panels$other <- .sc_comm_dist_panel(df$distance[df$class == "other"],
                                        "Other signaling",
                                        "#8c8c8c", xmax, density.alpha)
  }
  if (!length(panels)) {
    panels$all <- .sc_comm_dist_panel(df$distance, "All inferred pairs",
                                      "#69b3a2", xmax, density.alpha)
  }
  panels
}

#' Panels for a given mode
#' @noRd
.sc_c2c_panels <- function(object, split.by, mode, density.alpha,
                           interaction.range = NULL, contact.range = NULL,
                           tol = NULL) {
  object <- .sc_assert_spatial_cell_chat(object)
  mode <- match.arg(mode, c("geometry", "communication"))
  cache <- if (identical(mode, "communication")) {
    .sc_c2c_inference_cache(object)
  } else {
    .sc_c2c_geometry_cache(object, interaction.range, contact.range, tol)
  }
  .sc_comm_dist_panels(object, split.by, mode, density.alpha, cache)
}

#' @title Internal grid-size and membership helpers
#'
#' These helpers deliberately operate on the caller's coordinate frame and
#' already-created `sf` objects. They do not change coordinate orientation or
#' CRS metadata.
#' @noRd
.sc_resolve_grid_size <- function(
    coordinates,
    cellsize = NULL,
    grid.resolution = NULL,
    ratio = NULL
) {
  coordinates <- as.matrix(coordinates)
  if (!is.numeric(coordinates) || length(dim(coordinates)) != 2L ||
      ncol(coordinates) < 2L)
    stop("coordinates must be a numeric matrix with at least two columns", call. = FALSE)
  if (!nrow(coordinates))
    stop("coordinates must contain at least one row", call. = FALSE)
  if (any(!is.finite(coordinates)))
    stop("coordinates must contain finite numeric values", call. = FALSE)

  if (is.null(grid.resolution)) grid.resolution <- 2
  if (length(grid.resolution) != 1L || !is.numeric(grid.resolution) ||
      is.na(grid.resolution) || !is.finite(grid.resolution) ||
      grid.resolution <= 0)
    stop("grid.resolution must be a finite numeric scalar greater than zero",
         call. = FALSE)
  grid.resolution <- as.numeric(grid.resolution)

  if (!is.null(cellsize)) {
    if (!is.numeric(cellsize) || !length(cellsize) %in% c(1L, 2L))
      stop("cellsize must be a positive numeric vector of length 1 or 2",
           call. = FALSE)
    if (any(!is.finite(cellsize)) || any(cellsize <= 0))
      stop("cellsize must contain finite numeric values greater than zero",
           call. = FALSE)
    base.cellsize <- as.numeric(cellsize)
    source <- "user-supplied"
    nearest.distance.method <- NULL
  } else {
    if (nrow(coordinates) < 2L)
      stop("coordinates must contain at least two points", call. = FALSE)
    k <- min(2L, nrow(coordinates) - 1L)
    neighbors <- BiocNeighbors::findKNN(
      X = coordinates,
      k = k,
      get.index = TRUE,
      get.distance = TRUE,
      num.threads = 1L,
      BNPARAM = BiocNeighbors::VptreeParam()
    )
    neighbor.index <- as.matrix(neighbors$index)
    neighbor.distance <- as.matrix(neighbors$distance)
    row.index <- matrix(seq_len(nrow(coordinates)),
                        nrow = nrow(coordinates), ncol = ncol(neighbor.index))
    keep <- is.finite(neighbor.distance) &
      !is.na(neighbor.index) & neighbor.index != row.index
    if (!any(keep))
      stop("unable to find a non-self nearest neighbor for coordinates",
           call. = FALSE)
    base.cellsize <- min(neighbor.distance[keep])
    source <- "nearest-neighbor"
    nearest.distance.method <- "BiocNeighbors::findKNN"
  }

  effective.cellsize <- if (length(base.cellsize) == 1L) {
    rep(base.cellsize, 2L)
  } else {
    base.cellsize[seq_len(2L)]
  }
  effective.cellsize <- effective.cellsize * grid.resolution

  if (!is.null(ratio)) {
    if (length(ratio) != 1L || !is.numeric(ratio) || is.na(ratio) ||
        !is.finite(ratio) || ratio <= 0)
      stop("ratio must be a finite numeric scalar greater than zero",
           call. = FALSE)
    ratio <- as.numeric(ratio)
  }
  list(
    base.cellsize = base.cellsize,
    effective.cellsize = effective.cellsize,
    grid.resolution = grid.resolution,
    ratio = ratio,
    physical.cellsize = if (is.null(ratio)) NULL else effective.cellsize * ratio,
    calibrated = !is.null(ratio),
    source = source,
    nearest.distance.method = nearest.distance.method
  )
}
.sc_grid_membership <- function(points, grid) {
  hits <- sf::st_intersects(points, grid, sparse = TRUE)
  grid.index <- unlist(hits, use.names = FALSE)
  point.index <- rep.int(seq_along(hits), lengths(hits))
  n_grid <- length(sf::st_geometry(grid))
  grid.counts <- tabulate(grid.index, nbins = n_grid)
  list(
    hits = hits,
    point.counts = lengths(hits),
    grid.counts = grid.counts,
    point.index = point.index,
    grid.index = grid.index
  )
}

#' @title computeGridSize
#' @description
#' Estimate a grid cell size and optionally preview the resulting grid. This
#' function does not aggregate expression data or mutate the object.
#'
#' @details
#' Coordinates are read from the canonical `images$coordinates` matrix. Grid
#' geometry uses its first two analysis axes (`x`, `y`) without swapping or
#' transforming them. With `cellsize = NULL`, the minimum non-self
#' nearest-neighbor distance is estimated exactly without a dense pairwise
#' distance matrix. `grid.resolution` multiplies the resolved size
#' component-wise. When `images$spatial.factors$ratio` is unavailable, output
#' reports raw coordinate units and does not claim a physical scale.
#' Diagnostics use `getOption("SpatialCellChat.verbose", 1L)`: level 0 keeps
#' warnings and errors visible, level 1 emits the concise summary, level 2 adds
#' grid diagnostics, and level 3 reports sparse allocation details.
#'
#' @param object SpatialCellChat object with canonical spatial coordinates.
#' @param cellsize NULL or a positive finite numeric vector of length 1 or 2.
#'   NULL estimates the minimum non-self nearest-neighbor distance.
#' @param grid.resolution Positive finite numeric multiplier, defaulting to 2.
#' @param do.plot Boolean. If TRUE, return a `ggplot` preview; if FALSE,
#'   print diagnostics and return `NULL` invisibly.
#' @param what Character. One of `"polygons"`, `"corners"`, or `"centers"`.
#' @param square Boolean. If FALSE, create a hexagonal grid.
#' @return A `ggplot` when `do.plot = TRUE`; otherwise invisible `NULL`.
#' @export
computeGridSize <- function(
    object,
    cellsize = NULL,
    grid.resolution = NULL,
    do.plot = TRUE,
    what = "polygons",
    square = TRUE
) {
  object <- .sc_assert_spatial_cell_chat(object)
  coordinates <- object@images$coordinates
  if (is.null(coordinates) || ncol(coordinates) < 2L)
    stop("images$coordinates must contain at least x and y columns", call. = FALSE)
  coordinates <- as.matrix(coordinates[, seq_len(2L), drop = FALSE])
  if (!is.numeric(coordinates) || any(!is.finite(coordinates)))
    stop("images$coordinates must contain finite numeric values", call. = FALSE)
  colnames(coordinates) <- c("x_cent", "y_cent")

  spatial.factors <- object@images$spatial.factors
  ratio <- if (is.list(spatial.factors)) spatial.factors$ratio else NULL
  resolved <- .sc_resolve_grid_size(
    coordinates = coordinates,
    cellsize = cellsize,
    grid.resolution = grid.resolution,
    ratio = ratio
  )
  newcellsize <- resolved$effective.cellsize
  n_points <- nrow(coordinates)
  shape <- if (isTRUE(square)) "square" else "hexagonal"
  resolution <- resolved$grid.resolution
  base <- resolved$base.cellsize
  effective_x <- newcellsize[[1L]]
  effective_y <- newcellsize[[2L]]

  .cli("Grid-size preview", .type = "header")
  .cli("Input: {.val {n_points}} points; using {.val x/y} axes", .type = "info")
  .cli("Cellsize source: {.val {resolved$source}}; base={.val {base}}; resolution={.val {resolution}}",
       .type = "text", .verbose = 2L)
  .cli("Effective grid cellsize: x={.val {effective_x}}, y={.val {effective_y}}; resolution={.val {resolution}}",
       .type = "text", .verbose = 1L)
  if (resolved$calibrated) {
    physical_x <- resolved$physical.cellsize[[1L]]
    physical_y <- resolved$physical.cellsize[[2L]]
    .cli("Calibration: ratio={.val {ratio}}; physical cellsize={.val {physical_x}} x {.val {physical_y}}",
         .type = "text", .verbose = 1L)
  } else {
    .cli("Calibration: {.val uncalibrated}; physical scale is unavailable",
         .type = "warning", .verbose = 0L)
  }

  if (!isTRUE(do.plot))
    return(invisible(NULL))

  df.sf <- sf::st_as_sf(
    data.frame(coordinates, cell_type = object@idents),
    coords = c("x_cent", "y_cent"),
    remove = FALSE
  )
  sf::st_crs(df.sf) <- 3857
  square_grid <- sf::st_make_grid(
    df.sf,
    cellsize = newcellsize,
    what = what,
    square = square
  )
  square_grid_sf <- sf::st_sf(square_grid)
  square_grid_sf$grid_id <- seq_along(square_grid)

  membership <- .sc_grid_membership(df.sf, square_grid_sf)
  within_nGrid <- membership$point.counts
  grid_nspots <- membership$grid.counts
  n_grid <- length(grid_nspots)
  n_occupied <- sum(grid_nspots > 0L)
  n_empty <- n_grid - n_occupied
  n_unassigned <- sum(within_nGrid == 0L)
  n_multi_hit <- sum(within_nGrid > 1L)
  n_hits <- sum(within_nGrid)
  occupancy_pct <- if (n_grid) 100 * n_occupied / n_grid else 0
  .cli("Grid: {.val {shape}}/{.val {what}}; generated={.val {n_grid}}; occupied={.val {n_occupied}} ({sprintf('%.1f', occupancy_pct)}%); empty={.val {n_empty}}",
       .type = "text", .verbose = 1L)
  .cli("Assignment: {.val {n_points - n_unassigned}} assigned; {.val {n_unassigned}} unassigned; {.val {n_multi_hit}} multi-grid points; hits={.val {n_hits}}",
       .type = if (n_unassigned) "warning" else "success",
       .verbose = if (n_unassigned) 0L else 1L)
  .cli("Sparse membership: {.val {n_hits}} point-grid hits; no dense point-grid matrix allocated",
       .type = "debug", .verbose = 3L)
  .cli("Use `cellsize = c({effective_x}, {effective_y})` in `makeGridSpatialCellChat()`.",
       .type = "info", .verbose = 1L)

  if (is.null(object@idents)) {
    df.sf$cell_type <- factor("SpatialCellChatObj", levels = "SpatialCellChatObj")
  }
  color.use <- scPalette(nlevels(df.sf$cell_type))
  names(color.use) <- levels(df.sf$cell_type)
  ggplot() +
    ggplot2::geom_sf(data = square_grid) +
    ggplot2::geom_sf(data = df.sf, mapping = aes(color = cell_type)) +
    scale_color_manual(values = color.use, na.value = "grey40") +
    theme_minimal() + theme(
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    theme(legend.key = element_blank()) +
    labs(color = "Cell Type")
}


#' @title makeGridSpatialCellChat
#'
#' @description
#' Aggregate a SpatialCellChat object into spatial grid cells while retaining
#' cell identities, centroids, counts, and the selected assay layer.
#'
#' @details
#' Coordinates are read from canonical `images$coordinates`; grid geometry uses
#' the first two analysis axes (`x`, `y`) without swapping or transforming them.
#' Point-on-edge and point-on-corner membership follows
#' [sf::st_intersects()], so a point may contribute to multiple grids. The
#' `images$spatial.factors$ratio` value is retained for calibrated data;
#' uncalibrated data use raw coordinate units and do not imply physical units.
#' Invalid `cellsize` values (zero, negative, non-finite, or unsupported
#' lengths) are rejected.
#'
#' @param object SpatialCellChat object.
#' @param data.slot Assay layer used to make the grid. One of `norm`, `raw`,
#'   `scale`, `signaling`, or `smooth`; defaults to `norm`. The layer is read
#'   from `object@assay`.
#' @param cellsize Positive finite numeric vector of length 1 or 2. Grid
#'   geometry uses the first two canonical `images$coordinates` axes (`x`,
#'   `y`). Square cells use width and height; hexagonal cells use the distance
#'   between opposite edges. Use `computeGridSize()` for an estimate.
#' @param what Character. One of `"polygons"`, `"corners"`, or `"centers"`.
#' @param square Boolean. If FALSE, create a hexagonal grid.
#' @param idents.ties.method To determine the ident/cell group of each grid,
#'   use [base::max.col()] to choose the ident containing the most cells.
#'   See `ties.method` in [base::max.col()].
#'
#' @return SpatialCellChat object.
#' @export
makeGridSpatialCellChat <- function(
    object,
    data.slot = c("norm", "raw", "scale", "signaling", "smooth"),
    cellsize = c(5, 5),
    what = "polygons",
    square = TRUE,
    idents.ties.method = "first"
){
  object <- .sc_assert_spatial_cell_chat(object)
  data.slot <- match.arg(data.slot)
  coordinates <- object@images$coordinates
  if (is.null(coordinates) || ncol(coordinates) < 2L)
    stop("images$coordinates must contain at least x and y columns", call. = FALSE)
  coordinates <- as.data.frame(coordinates[, seq_len(2L), drop = FALSE])
  colnames(coordinates) <- c("x_cent", "y_cent")
  resolved <- .sc_resolve_grid_size(
    coordinates = as.matrix(coordinates[, c("x_cent", "y_cent"), drop = FALSE]),
    cellsize = cellsize,
    grid.resolution = 1
  )
  cellsize <- resolved$effective.cellsize
  coordinates$cell_type <- object@idents

  spatial.factors <- object@images$spatial.factors
  ratio <- if (is.list(spatial.factors) && !is.null(spatial.factors$ratio)) {
    spatial.factors$ratio
  } else {
    1
  }
  spot.size <- cellsize[[1L]] * ratio
  grid.factors <- list(ratio = ratio, tol = spot.size / 2)

  .cli("Making spatial grid", .type = "header")
  .cli("Input: {.val {nrow(coordinates)}} points; layer={.val {data.slot}}",
       .type = "info")
  if (is.null(spatial.factors$ratio)) {
    .cli("Calibration: {.val uncalibrated}; recommended contact range uses raw coordinate units",
         .type = "warning", .verbose = 0L)
  } else {
    .cli("Cellsize: {.val {cellsize}} units; physical spot size is {.val {spot.size}}",
         .type = "text", .verbose = 2L)
  }

  df.sf <- sf::st_as_sf(coordinates, coords = c("x_cent", "y_cent"),
                        remove = FALSE)
  sf::st_crs(df.sf) <- 3857
  square_grid <- sf::st_make_grid(df.sf, cellsize = cellsize,
                                  what = what, square = square)
  square_grid_sf <- sf::st_sf(square_grid)
  square_grid_sf$grid_id <- seq_along(square_grid)
  membership <- .sc_grid_membership(df.sf, square_grid_sf)
  within_n_grid <- membership$point.counts
  spots.counts <- membership$grid.counts
  occupied <- which(spots.counts > 0)
  n_grid <- length(spots.counts)
  n_unassigned <- sum(within_n_grid == 0L)
  n_multi_hit <- sum(within_n_grid > 1L)
  n_hits <- sum(within_n_grid)
  if (!length(occupied))
    stop("the requested grid does not contain any cells", call. = FALSE)
  within.spots <- lapply(occupied, function(grid_id) {
    membership$point.index[membership$grid.index == grid_id]
  })
  .cli("Grid: generated={.val {n_grid}}; occupied={.val {length(occupied)}}; empty={.val {n_grid - length(occupied)}}",
       .type = "text", .verbose = 1L)
  .cli("Assignment: {.val {nrow(coordinates) - n_unassigned}} assigned; {.val {n_unassigned}} unassigned; {.val {n_multi_hit}} multi-grid points; hits={.val {n_hits}}",
       .type = if (n_unassigned) "warning" else "success",
       .verbose = if (n_unassigned) 0L else 1L)
  .cli("Sparse membership retained {.val {n_hits}} point-grid hits for aggregation",
       .type = "debug", .verbose = 3L)
  levels.idents <- levels(object@idents)
  present.idents <- vapply(within.spots, function(index) {
    counts <- tabulate(match(object@idents[index], levels.idents),
                       nbins = length(levels.idents))
    levels.idents[max.col(matrix(counts, nrow = 1L),
                          ties.method = idents.ties.method)]
  }, character(1))

  new.names <- paste0("Grid", square_grid_sf$grid_id[occupied])
  new.meta <- as.data.frame(do.call(rbind, lapply(within.spots, function(index) {
    counts <- tabulate(match(object@idents[index], levels.idents),
                       nbins = length(levels.idents))
    setNames(as.list(counts), levels.idents)
  })), stringsAsFactors = FALSE)
  rownames(new.meta) <- new.names
  new.meta$cell.type <- factor(present.idents, levels = levels.idents)
  new.meta$spots.counts <- as.numeric(spots.counts[occupied])

  new.coordinates <- do.call(rbind, lapply(within.spots, function(index) {
    colMeans(coordinates[index, c("x_cent", "y_cent"), drop = FALSE])
  }))
  rownames(new.coordinates) <- new.names
  colnames(new.coordinates) <- c("x", "y")

  previous.data <- assay(object, data.slot)
  if (is.null(previous.data))
    stop("assay$", data.slot, " is empty", call. = FALSE)
  new.data <- do.call(cbind, lapply(within.spots, function(index) {
    Matrix::rowMeans(previous.data[, index, drop = FALSE])
  }))
  rownames(new.data) <- rownames(previous.data)
  colnames(new.data) <- new.names
  if (!inherits(new.data, "dgCMatrix")) new.data <- methods::as(new.data, "dgCMatrix")

  object.grid <- createSpatialCellChat(
    object = new.data,
    meta = new.meta,
    group.by = "cell.type",
    input.assay = "norm",
    datatype = "spatial",
    coordinates = new.coordinates,
    spatial.factors = grid.factors
  )
  object.grid@DB <- object@DB
  object.grid@images$.grid <- list(
    within.nGrid = within_n_grid,
    recommended.contact.range = spot.size
  )
  .cli("Making grid is done; recommended contact.range is stored in images$.grid",
       .type = "success")
  object.grid
}


.spatial_distance_cache_matches <- function(distance, coordinates,
                                            interaction.range, contact.range,
                                            ratio = NULL, tol = NULL) {
  if (!is.list(distance) ||
      !all(c("d.spatial", "adj.contact", ".parameters", ".coordinates") %in%
           names(distance)) ||
      !is.list(distance[[".parameters"]]))
    return(FALSE)
  coordinates <- tryCatch(as.matrix(coordinates), error = function(e) NULL)
  if (is.null(coordinates) || !is.matrix(coordinates) ||
      !is.numeric(coordinates) || any(!is.finite(coordinates)))
    return(FALSE)
  storage.mode(coordinates) <- "double"
  cached.coordinates <- distance[[".coordinates"]]
  if (!is.matrix(cached.coordinates) ||
      !identical(cached.coordinates, coordinates))
    return(FALSE)
  n_cells <- nrow(coordinates)
  for (name in c("d.spatial", "adj.contact")) {
    matrix.value <- distance[[name]]
    if (is.null(matrix.value) || !identical(dim(matrix.value), c(n_cells, n_cells)))
      return(FALSE)
    matrix.names <- dimnames(matrix.value)
    coordinate.names <- rownames(coordinates)
    if (!is.null(coordinate.names) &&
        (!identical(matrix.names[[1L]], coordinate.names) ||
         !identical(matrix.names[[2L]], coordinate.names)))
      return(FALSE)
  }
  params <- distance[[".parameters"]]
  expected <- list(
    interaction.range = interaction.range,
    contact.range = contact.range,
    ratio = ratio,
    tol = if (is.null(tol)) 0 else tol
  )
  same <- function(name) {
    actual <- params[[name]]
    wanted <- expected[[name]]
    if (is.null(actual) || is.null(wanted)) return(identical(actual, wanted))
    isTRUE(all.equal(as.numeric(actual), as.numeric(wanted), check.attributes = FALSE))
  }
  all(vapply(names(expected), same, logical(1)))
}

#' compute cell-cell distance and cell-cell contact adjacency matrix
#'
#' @description
#' Compute cell-cell distance based on the spatial coordinates and
#' generate cell-cell contact adjacency matrix with a contact.range restriction
#'
#' @details Distance builders emit no ordinary status messages. Inference callers
#' report their distance parameters and cache status under their own verbosity controls.
#'
#' @param coordinates Numeric matrix or data.frame; each row gives the spatial location of one cell/spot. Two- and three-dimensional coordinates are supported.
#' @param interaction.range Numeric(must positive). The maximum interaction/diffusion range of ligands. This hard threshold is used to filter out the connections between spatially distant cells
#' @param contact.range Numeric. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data).
#' Users can run the function `computeCellDistance` to get the center-to-center distance in the result's "d.spatial" key, which will help decide the value of `contact.range`.
#' @param tol Numeric. Add a non-negative tolerance when applying the interaction and contact thresholds, in the same physical unit as `interaction.range`. `NULL` is treated as zero.
#' Typically, `tol` should equal half the cell/spot diameter; for example, `65/2` for 10X Visium or `10/2` for Slide-seq.
#'
#' If the cell/spot size is not known, `tol` can be set to zero while `contact.range` is chosen from the returned center-to-center distances.
#'
#' @param ratio NULL or Numeric. Conversion factor from coordinate units to the physical unit used by `interaction.range` and `contact.range` (for example, micrometers). `NULL` means the coordinates are already expressed in that working unit.
#'
#' For example, setting `ratio = 0.18` indicates that 1 pixel equals 0.18um; distances are returned in micrometers after conversion.
#' For 10X Visium, this is the theoretical spot size (65um) divided by the full-resolution spot diameter in pixels (`spot.size.fullres`).
#'
#' @return List. The `d.spatial` and `adj.contact` keys store the cell-cell
#' distance and contact adjacency matrices. `.coordinates` records the exact
#' coordinate matrix used to build them, while `.parameters` records the
#' thresholds and scale used so downstream consumers do not silently reuse a
#' cache made from different coordinates or distance settings.
#'
#' @examples
computeCellDistance <- function (
    coordinates,
    interaction.range = 250,
    contact.range = 10,
    ratio = NULL,
    tol = 10/2
)
{
  if (is.data.frame(coordinates) || inherits(coordinates, "Matrix"))
    coordinates <- as.matrix(coordinates)
  if (!is.matrix(coordinates) || length(dim(coordinates)) != 2L ||
      !ncol(coordinates) %in% c(2L, 3L))
    stop("coordinates must be a numeric matrix with two or three columns", call. = FALSE)
  if (!is.numeric(coordinates) || any(!is.finite(coordinates)))
    stop("coordinates must contain finite numeric values", call. = FALSE)
  storage.mode(coordinates) <- "double"
  n_cells <- nrow(coordinates)
  if (n_cells < 1L)
    stop("coordinates must contain at least one row", call. = FALSE)

  validate_scalar <- function(value, name, lower = 0, strictly_positive = FALSE) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) ||
        (strictly_positive && value <= lower) ||
        (!strictly_positive && value < lower))
      stop(name, " must be a finite numeric scalar ",
           if (strictly_positive) "greater than " else "at least ", lower,
           call. = FALSE)
    as.numeric(value)
  }
  interaction.range <- validate_scalar(interaction.range, "interaction.range",
                                       strictly_positive = TRUE)
  contact.range <- validate_scalar(contact.range, "contact.range")
  tol <- if (is.null(tol)) 0 else validate_scalar(tol, "tol")
  if (!is.null(ratio))
    ratio <- validate_scalar(ratio, "ratio", strictly_positive = TRUE)

  # `interaction.range` and `contact.range` are physical distances.  When a
  # coordinate-to-physical conversion is absent, the caller has declared that
  # the coordinate unit is already the working unit.
  threshold <- interaction.range + tol
  if (!is.null(ratio)) threshold <- threshold / ratio

  neighbors <- BiocNeighbors::queryNeighbors(
    X = coordinates,
    query = coordinates,
    threshold = threshold,
    get.index = TRUE,
    get.distance = TRUE,
    num.threads = 1L,
    BNPARAM = BiocNeighbors::VptreeParam()
  )
  index <- neighbors$index
  distance <- neighbors$distance
  row_index <- rep.int(seq_len(n_cells), lengths(index))
  col_index <- unlist(index, use.names = FALSE)
  distance <- unlist(distance, use.names = FALSE)
  keep <- logical(length(distance))
  if (length(distance))
    keep <- row_index != col_index & distance > 0 & is.finite(distance)
  if (length(distance) && any(keep)) {
    i <- row_index[keep]
    j <- col_index[keep]
    x <- distance[keep]
    if (!is.null(ratio)) x <- x * ratio
  } else {
    i <- integer()
    j <- integer()
    x <- numeric()
  }
  coordinate_names <- rownames(coordinates)
  d.spatial <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(n_cells, n_cells),
    dimnames = list(coordinate_names, coordinate_names),
    repr = "C"
  )
  adj.contact <- createCellCellContactMatrixFrom_dspatial(
    d.spatial = d.spatial,
    tol = tol,
    contact.threshold = contact.range
  )
  list(
    d.spatial = d.spatial,
    adj.contact = adj.contact,
    .coordinates = coordinates,
    .parameters = list(
      interaction.range = interaction.range,
      contact.range = contact.range,
      ratio = ratio,
      tol = tol
    )
  )
}

#' create a cell-cell contact matrix from `d.spatial` object
#' @description
#' A function defined for `computeCellDistance`. We use the function `createSparseMatrixFrom_distObj` to convert the dist object into a CsparseMatrix object when running `computeCellDistance`.
#' Here we design the function to generate cell-cell contact adjacency matrix with a `contact.range` restriction and a `tol` tolerance passed from `computeCellDistance`'s parameters respectively
#'
#' @param d.spatial a CsparseMatrix (dsC Matrix) object, see in `computeCellDistance`:
#' d.spatial <- stats::dist(coordinates);
#' d.spatial <- createSparseMatrixFrom_distObj(d.spatial)
#' @param tol Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix,
#' passed from `computeCellDistance`'s same parameter.
#' By default `tol = NULL` means `tol` equals to the half value of cell/spot size in the unit of um.
#' @param contact.threshold Numeric. The interaction range threshold(Unit: microns) to restrict the contact-dependent signaling,
#' passed from `computeCellDistance`'s parameter `contact.range`.
#' If the distance of a pair of cells is larger than contact.threshold, regard the two cells as not in contact. Otherwise, regard the two cells as in contact.
#'
#' @details This matrix helper emits no ordinary status messages; callers own
#' distance and contact diagnostics.
#'
#' @return A `CsparseMatrix` object in `Matrix` Package
#' @export
#'
#' @examples
createCellCellContactMatrixFrom_dspatial <- function(d.spatial,contact.threshold,tol){

  adj.contact <- as(d.spatial,Class = "TsparseMatrix")
  # `min(adj.contact@x)` is the smallest value in the d.spatial matrix.
  # Take the larger value in `contact.threshold` and `min(adj.contact@x)` as the
  # interaction range threshold (Unit: microns) to restrict the contact-dependent signaling.
  # contact.range <- max(contact.threshold,min(adj.contact@x))
  contact.range <- contact.threshold

  cellcell.contact.index <- which(adj.contact@x <= (contact.range+tol))

  # update
  adj.contact@i <- adj.contact@i[cellcell.contact.index]
  adj.contact@j <- adj.contact@j[cellcell.contact.index]
  adj.contact@x <- rep.int(1,times = length(cellcell.contact.index))

  adj.contact <- as(adj.contact,Class = "CsparseMatrix")

  # fill up the diagonal, because d.spatial's diagonal is all-zero
  Matrix::diag(adj.contact) <- 1
  return(adj.contact)

}



#' Assessment of the colocalization between any cell groups
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param group a factor vector defining the labels of each cell/spot
#' @param idents.use the identity to use
#' @param symmetric whether producing a symmetric colocalization matrix
#' @param nboot number of permutation when assessing the colocalization
#' @param seed.use set a random seed. By default, set the seed to 1.
#'
#' @return A square matrix giving the computed p-values after permuting cell labels
#'
#' @export

computeColocalization <- function(coordinates, group = NULL, idents.use = NULL, symmetric = TRUE, nboot = 100, seed.use = 1L) {
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  if (!is.factor(group)) {
    stop("Please input the `group` as a factor!")
  }
  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use0 <- level.use
  if (is.null(idents.use)) {
    level.use <- level.use[level.use %in% unique(group)]
  } else {
    numCluster <- length(idents.use)
    level.use <- level.use[level.use %in% idents.use]
    level.use0 <- level.use
  }

  fdr = matrix(NaN,numCluster,numCluster)
  for (i in c(1:numCluster)){
    for (j in c(1:numCluster)){
      label_1 = level.use[i]
      label_2 = level.use[j]
      data_label_1 = coordinates[group==label_1,]
      data_label_2 = coordinates[group==label_2,]

      shuffle_sel = !(group %in% c(label_1,label_2))
      data_nont = coordinates[shuffle_sel,]

      coord2_dist = matrix(0,nrow=nrow(data_label_1),ncol=nrow(data_label_2))
      coord1 = cbind(data_label_1$x_cent,data_label_1$y_cent) # every row is a cell in label1
      coord2_shuf_dist_array = array(dim=c(nrow(data_label_1),nrow(data_label_2), nboot))

      for (idx2 in c(1:nrow(data_label_2))){
        coord2 = c(data_label_2$x_cent[idx2], data_label_2$y_cent[idx2])
        coord2_dist[,idx2] = raster::pointDistance(coord1,coord2,lonlat=F)
      }

      set.seed(seed.use)
      permutation <- replicate(nboot, sample.int(nrow(data_nont), size = nrow(data_label_1)))
      for (idx3 in c(1:nboot)){
        #data_label_1_shuf = data_nont[sample(c(1:nrow(data_nont)),nrow(data_label_1)),]
        data_label_1_shuf = data_nont[permutation[, idx3], ]
        coord1_shuf = cbind(data_label_1_shuf$x_cent,data_label_1_shuf$y_cent)
        for (idx2 in c(1:nrow(data_label_2))){
          coord2 = c(data_label_2$x_cent[idx2], data_label_2$y_cent[idx2])
          coord2_shuf_dist_array[,idx2,idx3] = raster::pointDistance(coord1_shuf,coord2,lonlat=F)
        }
      }


      label2_mean_dists = rowMeans(coord2_dist)
      label2_mean_dists_median = median(label2_mean_dists) # observed median of mean dists
      label2_shuf_medians = rep(0,nboot)

      fdr_count= c(0,0)
      for (idx4 in c(1: nboot)){
        label2_mean_dists_shuf = rowMeans(coord2_shuf_dist_array[,,idx4])
        label2_shuf_medians[idx4] = median(label2_mean_dists_shuf)
        if (label2_shuf_medians[idx4]>label2_mean_dists_median){
          fdr_count[1] = fdr_count[1]+1
        }
        else{
          fdr_count[2] = fdr_count[2]+1
        }
      }
      fdr[i,j] = fdr_count[2]/nboot
    }
  }
  if (symmetric) {
    for (i in 1:(numCluster-1)){
      for (j in (i+1):numCluster){
        if (fdr[i,j] > fdr[j,i]) {
          fdr[i,j] <- fdr[j,i]
        }
        fdr[j,i] <- fdr[i,j]
      }
    }

  }
  rownames(fdr) <- level.use0; colnames(fdr) <- level.use0
  return(fdr)
}




