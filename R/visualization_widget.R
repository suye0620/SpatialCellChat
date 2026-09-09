#' Prepare a communication-flow payload for an interactive widget
#'
#' This file is additive: the existing static visualization functions are kept
#' unchanged. The preparation code intentionally follows their coordinate frame
#' and grid semantics so the browser renderer is a view, not a second field
#' definition.

.sc_widget_check_scalar <- function(value, name, lower = -Inf, upper = Inf,
                                    integer = FALSE, inclusive.lower = TRUE) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value)
  if (valid && integer) valid <- value == as.integer(value)
  if (valid && inclusive.lower) valid <- value >= lower else if (valid) valid <- value > lower
  if (valid) valid <- value <= upper
  if (!valid) {
    bound <- if (is.finite(lower) && is.finite(upper)) {
      paste0(if (inclusive.lower) "[" else "(", lower, ", ", upper, "]")
    } else if (is.finite(lower)) {
      paste0(if (inclusive.lower) "at least " else "greater than ", lower)
    } else {
      paste0("at most ", upper)
    }
    stop(name, " must be one finite numeric value ", bound, call. = FALSE)
  }
  as.numeric(value)
}

.sc_widget_color_map <- function(labels, color.use) {
  groups <- levels(labels)
  if (is.null(color.use)) {
    color.use <- scPalette(length(groups))
    names(color.use) <- groups
  } else {
    if (!is.atomic(color.use) || !length(color.use))
      stop("color.use must be a non-empty named or group-ordered color vector", call. = FALSE)
    if (is.null(names(color.use))) {
      if (length(color.use) != length(groups))
        stop("unnamed color.use must have one color per identity level", call. = FALSE)
      names(color.use) <- groups
    }
    if (anyNA(names(color.use)) || any(!nzchar(names(color.use))))
      stop("color.use names must be non-empty identity levels", call. = FALSE)
    if (anyDuplicated(names(color.use)))
      stop("color.use names must be unique", call. = FALSE)
    missing <- setdiff(groups, names(color.use))
    if (length(missing))
      stop("color.use is missing identity levels: ", paste(missing, collapse = ", "),
           call. = FALSE)
    color.use <- color.use[groups]
  }
  color.use <- as.character(color.use)
  if (anyNA(color.use) || any(!nzchar(color.use)))
    stop("color.use must contain non-empty color values", call. = FALSE)
  stats::setNames(color.use, groups)
}

.sc_widget_assets <- function() {
  files <- c(
    "spatialcellchat-commun-flow.js",
    "spatialcellchat-commun-flow.css"
  )
  source.assets <- normalizePath(
    file.path(getwd(), "inst", "htmlwidgets"),
    winslash = "/", mustWork = FALSE
  )
  installed.assets <- system.file("htmlwidgets", package = "SpatialCellChat")
  asset.dir <- if (nzchar(installed.assets) &&
                   all(file.exists(file.path(installed.assets, files)))) {
    installed.assets
  } else if (all(file.exists(file.path(source.assets, files)))) {
    source.assets
  } else {
    stop(
      "cannot locate communication-flow widget assets; reinstall SpatialCellChat or run from the package root",
      call. = FALSE
    )
  }
  list(
    directory = asset.dir,
    from.source = identical(asset.dir, source.assets)
  )
}

.sc_widget_dependency <- function(assets) {
  htmltools::htmlDependency(
    name = "spatialcellchat-commun-flow",
    version = "0.1.0",
    src = assets$directory,
    script = "spatialcellchat-commun-flow.js",
    stylesheet = "spatialcellchat-commun-flow.css"
  )
}

.sc_widget_regular_grid <- function(x, y, u, v, valid, field.grid.max) {
  nx <- length(x)
  ny <- length(y)
  if (length(u) != nx * ny || length(v) != nx * ny || length(valid) != nx * ny)
    stop("interpolated field dimensions do not match the regular grid", call. = FALSE)

  ix <- unique(round(seq.int(1L, nx, length.out = min(nx, field.grid.max))))
  iy <- unique(round(seq.int(1L, ny, length.out = min(ny, field.grid.max))))
  u <- matrix(as.numeric(u), nrow = nx, ncol = ny)
  v <- matrix(as.numeric(v), nrow = nx, ncol = ny)
  valid <- matrix(as.logical(valid), nrow = nx, ncol = ny)
  list(
    x = as.numeric(x[ix]),
    y = as.numeric(y[iy]),
    u = as.numeric(u[ix, iy, drop = FALSE]),
    v = as.numeric(v[ix, iy, drop = FALSE]),
    valid = as.logical(valid[ix, iy, drop = FALSE])
  )
}

.sc_widget_zero_grid <- function(df.field, occupied.grid, field.grid.max) {
  x.range <- range(df.field$x_cent, finite = TRUE)
  y.range <- range(df.field$y_cent, finite = TRUE)
  if (diff(x.range) == 0) x.range <- x.range + c(-0.5, 0.5)
  if (diff(y.range) == 0) y.range <- y.range + c(-0.5, 0.5)
  x <- seq(x.range[1L], x.range[2L], length.out = min(2L, field.grid.max))
  y <- seq(y.range[1L], y.range[2L], length.out = min(2L, field.grid.max))
  points <- expand.grid(x = x, y = y)
  points.sf <- sf::st_as_sf(points, coords = c("x", "y"), remove = FALSE, crs = 3857)
  hits <- sf::st_intersects(points.sf, occupied.grid, sparse = TRUE)
  valid <- lengths(hits) > 0L
  list(x = x, y = y,
       u = numeric(length(points$x)), v = numeric(length(points$x)), valid = valid)
}

.sc_widget_interpolate <- function(current.gridded, occupied.grid, field.grid.max,
                                   df.field) {
  finite.u <- is.finite(current.gridded$u)
  finite.v <- is.finite(current.gridded$v)
  if (!all(finite.u & finite.v)) {
    current.gridded$u[!finite.u] <- 0
    current.gridded$v[!finite.v] <- 0
  }
  nonzero <- any(abs(current.gridded$u) > 0 | abs(current.gridded$v) > 0)
  if (!nonzero)
    return(.sc_widget_zero_grid(df.field, occupied.grid, field.grid.max))

  if (length(unique(current.gridded$x)) < 2L ||
      length(unique(current.gridded$y)) < 2L)
    stop("communication field interpolation requires at least two distinct x and y grid centers",
         call. = FALSE)

  if (!requireNamespace("oce", quietly = TRUE))
    stop("package 'oce' is required for communication flow interpolation; install it before rendering",
         call. = FALSE)
  u.se <- tryCatch(
    oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$u),
    error = function(e) stop("failed to interpolate the communication-flow u component: ",
                             conditionMessage(e), call. = FALSE)
  )
  v.se <- tryCatch(
    oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$v),
    error = function(e) stop("failed to interpolate the communication-flow v component: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (!identical(u.se$xg, v.se$xg) || !identical(u.se$yg, v.se$yg) ||
      !identical(dim(u.se$zg), dim(v.se$zg)))
    stop("u and v interpolation returned incompatible grids", call. = FALSE)

  nx <- length(u.se$xg)
  ny <- length(u.se$yg)
  if (!identical(dim(u.se$zg), c(nx, ny)))
    stop("interpolated communication field has inconsistent dimensions", call. = FALSE)

  points <- expand.grid(x = u.se$xg, y = u.se$yg)
  points.sf <- sf::st_as_sf(points, coords = c("x", "y"), remove = FALSE, crs = 3857)
  hits <- sf::st_intersects(points.sf, occupied.grid, sparse = TRUE)
  valid <- lengths(hits) > 0L
  u <- as.numeric(u.se$zg)
  v <- as.numeric(v.se$zg)
  finite <- is.finite(u) & is.finite(v)
  valid <- valid & finite
  u[!valid] <- 0
  v[!valid] <- 0

  .sc_widget_regular_grid(
    x = u.se$xg, y = u.se$yg, u = u, v = v, valid = valid,
    field.grid.max = field.grid.max
  )
}

.sc_prepare_commun_flow_widget <- function(
    object, signaling, slot.name, pattern, cellsize, square,
    grid.resolution, grid.field.mean, color.use, point.size, image.alpha,
    field.grid.max) {
  if (!methods::is(object, "SpatialCellChat"))
    stop("object must be a SpatialCellChat", call. = FALSE)
  if (!is.character(slot.name) || length(slot.name) != 1L ||
      !slot.name %in% c("net", "netP"))
    stop("slot.name must be 'net' or 'netP'", call. = FALSE)
  if (!is.character(pattern) || length(pattern) != 1L ||
      !pattern %in% c("incoming", "outgoing"))
    stop("pattern must be 'incoming' or 'outgoing'", call. = FALSE)
  if (!is.character(grid.field.mean) || length(grid.field.mean) != 1L ||
      !grid.field.mean %in% c("median", "sum"))
    stop("grid.field.mean must be 'median' or 'sum'", call. = FALSE)
  if (!is.character(signaling) || length(signaling) != 1L || is.na(signaling) ||
      !nzchar(signaling))
    stop("signaling must be one non-empty layer name", call. = FALSE)
  point.size <- .sc_widget_check_scalar(point.size, "point.size", lower = 0,
                                        inclusive.lower = FALSE)
  image.alpha <- .sc_widget_check_scalar(image.alpha, "image.alpha", lower = 0,
                                         upper = 1)
  field.grid.max <- .sc_widget_check_scalar(field.grid.max, "field.grid.max",
                                            lower = 2, upper = 256, integer = TRUE)
  field.grid.max <- as.integer(field.grid.max)

  result <- methods::slot(object, slot.name)
  cell <- result$cell
  field <- if (is.list(cell)) cell$field[[pattern]] else NULL
  if (!inherits(field, "SparseChatArray"))
    stop("communication fields must be stored as cell-level SparseChatArray results; run computeCommunField first",
         call. = FALSE)
  signaling.names <- dimnames(field)[[3L]]
  if (is.null(signaling.names) || !(signaling %in% signaling.names))
    stop("Please check the input 'signaling' and make sure it has been computed via `computeCommunField`.",
         call. = FALSE)

  field.matrix <- as.matrix(field[, , signaling, drop = TRUE])
  coordinates <- object@images$coordinates
  if (!is.matrix(coordinates) || ncol(coordinates) < 2L ||
      nrow(coordinates) != nrow(field.matrix))
    stop("images$coordinates must contain one row and at least two columns per cell",
         call. = FALSE)
  coordinates <- coordinates[, seq_len(2L), drop = FALSE]
  cell.names <- rownames(field.matrix)
  if (is.null(cell.names)) cell.names <- rownames(coordinates)
  if (is.null(cell.names) || is.null(rownames(coordinates)) ||
      !identical(rownames(coordinates), cell.names))
    stop("images$coordinates rownames must match communication cell names", call. = FALSE)
  if (!is.numeric(coordinates) || any(!is.finite(coordinates)))
    stop("images$coordinates must contain finite numeric values", call. = FALSE)
  if (anyDuplicated(data.frame(x = coordinates[, 1L], y = coordinates[, 2L])))
    stop("communication-flow widget requires unique cell coordinates", call. = FALSE)

  # Match netVisual_CommunFlow(): swap the canonical axes, then negate the
  # plotted y coordinate and vector component in the browser's spatial frame.
  temp.coordinates <- coordinates
  coordinates[, 1L] <- temp.coordinates[, 2L]
  coordinates[, 2L] <- temp.coordinates[, 1L]
  temp.field <- field.matrix
  field.matrix[, 1L] <- temp.field[, 2L]
  field.matrix[, 2L] <- temp.field[, 1L]
  colnames(coordinates) <- c("x_cent", "y_cent")
  colnames(field.matrix) <- c("dx", "dy")
  field.magnitude <- sqrt(field.matrix[, "dx"]^2 + field.matrix[, "dy"]^2)
  labels <- object@idents
  if (length(labels) != nrow(coordinates) || !identical(names(labels), cell.names))
    stop("object identities must be named and aligned to communication cell names", call. = FALSE)
  labels <- factor(as.character(labels), levels = levels(labels))
  df.field <- data.frame(
    cell_id = cell.names,
    x_cent = coordinates[, "x_cent"],
    y_cent = -coordinates[, "y_cent"],
    dx = field.matrix[, "dx"],
    dy = -field.matrix[, "dy"],
    mag = field.magnitude,
    label = labels,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  colors <- .sc_widget_color_map(labels, color.use)

  spatial.factors <- object@images$spatial.factors
  ratio <- if (is.list(spatial.factors)) spatial.factors$ratio else NULL
  resolved <- .sc_resolve_grid_size(
    coordinates = coordinates,
    cellsize = cellsize,
    grid.resolution = grid.resolution,
    ratio = ratio
  )
  grid <- sf::st_make_grid(
    sf::st_as_sf(df.field, coords = c("x_cent", "y_cent"), remove = FALSE, crs = 3857),
    cellsize = resolved$effective.cellsize,
    what = "polygons",
    square = square
  )
  grid.sf <- sf::st_sf(grid_id = seq_along(grid), geometry = grid)
  points.sf <- sf::st_as_sf(df.field, coords = c("x_cent", "y_cent"),
                            remove = FALSE, crs = 3857)
  membership <- .sc_grid_membership(points.sf, grid.sf)
  contained <- split(
    membership$point.index,
    factor(membership$grid.index, levels = seq_along(grid))
  )
  centers <- sf::st_coordinates(sf::st_centroid(grid.sf))
  occupied <- membership$grid.counts > 0L
  occupied.grid <- grid.sf[occupied, , drop = FALSE]
  if (!nrow(occupied.grid))
    stop("communication-flow grid contains no occupied cells", call. = FALSE)

  summarize <- if (identical(grid.field.mean, "median")) stats::median else base::sum
  u <- vapply(contained, function(index) {
    if (!length(index)) return(0)
    summarize(df.field$dx[index], na.rm = TRUE)
  }, numeric(1L))
  v <- vapply(contained, function(index) {
    if (!length(index)) return(0)
    summarize(df.field$dy[index], na.rm = TRUE)
  }, numeric(1L))
  current.gridded <- data.frame(
    x = centers[, 1L], y = centers[, 2L], u = u, v = v,
    obs = membership$grid.counts, stringsAsFactors = FALSE
  )
  current.gridded <- current.gridded[occupied, , drop = FALSE]
  grid.payload <- .sc_widget_interpolate(
    current.gridded, occupied.grid, field.grid.max, df.field
  )

  mass <- rep(NA_real_, nrow(df.field))
  probability <- if (is.list(cell)) cell$prob else NULL
  if (inherits(probability, "SparseChatArray") && signaling %in% names(probability)) {
    probability.layer <- probability[[signaling]]
    mass <- if (identical(pattern, "outgoing")) {
      Matrix::rowSums(probability.layer)
    } else {
      Matrix::colSums(probability.layer)
    }
    mass <- as.numeric(mass[match(cell.names, names(mass))])
    mass[!is.finite(mass)] <- 0
  }

  list(
    cell = list(
      id = as.character(df.field$cell_id),
      x = as.numeric(df.field$x_cent),
      y = as.numeric(df.field$y_cent),
      dx = as.numeric(df.field$dx),
      dy = as.numeric(df.field$dy),
      label = as.character(df.field$label),
      magnitude = as.numeric(df.field$mag),
      mass = mass,
      color = unname(colors[as.character(df.field$label)])
    ),
    field = grid.payload,
    bounds = list(
      xmin = min(c(df.field$x_cent, grid.payload$x)),
      xmax = max(c(df.field$x_cent, grid.payload$x)),
      ymin = min(c(df.field$y_cent, grid.payload$y)),
      ymax = max(c(df.field$y_cent, grid.payload$y))
    ),
    metadata = list(
      signaling = signaling,
      slot.name = slot.name,
      pattern = pattern,
      legendTitle = if (identical(pattern, "outgoing")) "Sources" else "Targets",
      current = list(
        color = "grey25",
        width = c(0.65, 1.55),
        alpha = 0.84
      ),
      field.mean = grid.field.mean,
      coordinate.unit = if (is.list(object@images$coordinate.system)) {
        object@images$coordinate.system$unit %||% "unknown"
      } else "unknown",
      semantics = "fixed communication vector field particle animation",
      direction = "moving particles show signal flow direction",
      current.meaning = "net directional current; opposing links may cancel"
    ),
    diagnostics = list(
      cell.count = nrow(df.field),
      grid.count = length(grid),
      occupied.grid.count = sum(occupied),
      grid.resolution = resolved$grid.resolution,
      effective.cellsize = resolved$effective.cellsize,
      calibrated = resolved$calibrated
    ),
    style = list(
      point.size = point.size,
      image.alpha = image.alpha,
      colors = as.list(colors)
    )
  )
}

#' Interactive particle animation for a communication flow field
#'
#' @description
#' Render one precomputed incoming or outgoing communication field as a local
#' Canvas 2D htmlwidget. The static layer shows identity-colored cell points
#' and a Sources or Targets identity legend. Moving particles are the only
#' flow-direction overlay and move over the same fixed field; this is not a
#' biological time-series or a physical transport-speed estimate.
#'
#' @param object A SpatialCellChat object with a field computed by
#'   \code{computeCommunField}.
#' @param signaling A single pathway or ligand-receptor field name.
#' @param slot.name Either \code{"netP"} or \code{"net"}.
#' @param pattern Either \code{"incoming"} or \code{"outgoing"}.
#' @param cellsize Positive numeric scalar or length-two vector. \code{NULL}
#'   uses the existing nearest-neighbor grid-size resolver.
#' @param square Logical; use square grid cells when \code{TRUE}, otherwise
#'   use the existing hexagonal grid geometry before regular interpolation.
#' @param grid.resolution Positive grid-size multiplier.
#' @param grid.field.mean Either \code{"median"} or \code{"sum"}.
#' @param color.use Optional named identity-color vector.
#' @param point.size Cell point radius in widget pixels before scale fitting.
#' @param image.alpha Cell point alpha; also adjustable in the browser control.
#' @param particle.count Number of browser particles, bounded at 5000.
#' @param particle.speed Visual coordinate-units-per-animation-second scale;
#'   defaults to 3 and is not a physical transport speed.
#' @param trail.length Number of trail fade steps; larger values leave longer
#'   particle trails.
#' @param field.grid.max Maximum grid extent per axis, from 2 through 256.
#' @param seed Non-negative integer random seed for deterministic particles.
#' @param zoom Initial browser zoom factor from 1 through 4.
#' @param width,height Optional htmlwidget dimensions.
#' @param title.name Optional widget title.
#' @return An htmlwidget object.
#' @export
#' @importFrom htmlwidgets createWidget
netVisual_CommunFlowWidget <- function(
    object,
    signaling,
    slot.name = "netP",
    pattern = c("incoming", "outgoing"),
    cellsize = NULL,
    square = TRUE,
    grid.resolution = NULL,
    grid.field.mean = c("median", "sum"),
    color.use = NULL,
    point.size = 2,
    image.alpha = 0.32,
    particle.count = 1000L,
    particle.speed = 3,
    trail.length = 24L,
    field.grid.max = 128L,
    seed = 1L,
    zoom = 1,
    width = NULL,
    height = NULL,
    title.name = NULL) {
  slot.name <- match.arg(slot.name, c("net", "netP"))
  pattern <- match.arg(pattern)
  grid.field.mean <- match.arg(grid.field.mean)
  if (!is.logical(square) || length(square) != 1L || is.na(square))
    stop("square must be TRUE or FALSE", call. = FALSE)
  particle.count <- .sc_widget_check_scalar(particle.count, "particle.count",
                                             lower = 1, upper = 5000, integer = TRUE)
  particle.speed <- .sc_widget_check_scalar(particle.speed, "particle.speed",
                                            lower = 0, inclusive.lower = FALSE)
  trail.length <- .sc_widget_check_scalar(trail.length, "trail.length",
                                          lower = 0, upper = 80, integer = TRUE)
  seed <- .sc_widget_check_scalar(seed, "seed", lower = 0, upper = 2147483647,
                                  integer = TRUE)
  zoom <- .sc_widget_check_scalar(zoom, "zoom", lower = 1, upper = 4)
  particle.count <- as.integer(particle.count)
  trail.length <- as.integer(trail.length)
  seed <- as.integer(seed)

  payload <- .sc_prepare_commun_flow_widget(
    object = object, signaling = signaling, slot.name = slot.name,
    pattern = pattern, cellsize = cellsize, square = square,
    grid.resolution = grid.resolution, grid.field.mean = grid.field.mean,
    color.use = color.use, point.size = point.size, image.alpha = image.alpha,
    field.grid.max = field.grid.max
  )
  if (is.null(title.name))
    title.name <- paste0(tools::toTitleCase(pattern), " communication flow of ", signaling)

  payload$options <- list(
    particleCount = particle.count,
    particleSpeed = particle.speed,
    trailLength = trail.length,
    seed = seed,
    zoom = zoom,
    title = as.character(title.name)
  )
  assets <- .sc_widget_assets()
  htmlwidgets::createWidget(
    name = "spatialcellchat-commun-flow",
    x = payload,
    width = width,
    height = height,
    sizingPolicy = htmlwidgets::sizingPolicy(
      viewer.padding = 0,
      browser.padding = 0,
      browser.fill = TRUE,
      knitr.figure = TRUE
    ),
    package = if (assets$from.source) "htmlwidgets" else "SpatialCellChat",
    dependencies = if (assets$from.source) list(.sc_widget_dependency(assets)) else NULL
  )
}
