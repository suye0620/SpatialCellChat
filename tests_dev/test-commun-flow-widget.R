# Contract tests for the additive communication-flow htmlwidget API.

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
  library(sf)
  library(dplyr)
  library(spatstat.sparse)
  library(htmlwidgets)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/analysis.R")
source("R/visualization.R")
source("R/visualization_widget.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

expect_error <- function(label, expression) {
  result <- try(force(expression), silent = TRUE)
  check(label, inherits(result, "try-error"))
}

cells <- paste0("c", seq_len(9L))
expr <- Matrix::Matrix(
  matrix(1, nrow = 3L, ncol = length(cells),
         dimnames = list(paste0("g", seq_len(3L)), cells)),
  sparse = TRUE
)
meta <- data.frame(
  label = factor(rep(c("A", "B"), length.out = length(cells))),
  row.names = cells
)
coords <- cbind(x = rep(0:2, each = 3L), y = rep(0:2, 3L))
rownames(coords) <- cells
chat <- createSpatialCellChat(
  expr, meta = meta, group.by = "label", datatype = "spatial",
  coordinates = coords, spatial.factors = list(ratio = 1, tol = 0.5)
)
prob.matrix <- Matrix::sparseMatrix(
  i = c(1L, 2L, 3L, 4L, 5L),
  j = c(2L, 3L, 4L, 5L, 6L),
  x = c(3, 2, 1, 2, 1),
  dims = c(length(cells), length(cells)),
  dimnames = list(cells, cells)
)
chat@netP <- list(
  cell = list(prob = SparseChatArray(
    list(signal = prob.matrix), dimnames = list(cells, cells, "signal")
  )),
  group = list()
)
chat <- computeCommunField(chat, signaling.name = "signal")
before <- chat
widget <- netVisual_CommunFlowWidget(
  chat, signaling = "signal", pattern = "outgoing",
  cellsize = c(1, 1), grid.resolution = 1,
  particle.count = 24L, particle.speed = 1.2,
  trail.length = 10L, field.grid.max = 32L, zoom = 2.2
)

check("valid call returns a spatialcellchat htmlwidget",
      inherits(widget, "htmlwidget") &&
        "spatialcellchat-commun-flow" %in% class(widget))
check("widget call does not mutate input",
      identical(chat, before))

dependency <- widget$dependencies[[1L]]
check("source widget embeds its local frontend assets",
      identical(attr(widget, "package"), "htmlwidgets") &&
        length(widget$dependencies) == 1L &&
        inherits(dependency, "html_dependency") &&
        identical(dependency$name, "spatialcellchat-commun-flow") &&
        identical(dependency$script, "spatialcellchat-commun-flow.js") &&
        identical(dependency$stylesheet, "spatialcellchat-commun-flow.css"))

render.dir <- tempfile("spatialcellchat-commun-flow-")
dir.create(render.dir)
render.file <- file.path(render.dir, "widget.html")
htmlwidgets::saveWidget(widget, render.file, selfcontained = FALSE, libdir = "lib")
rendered.html <- paste(readLines(render.file, warn = FALSE), collapse = "\n")
dependency.dir <- file.path(render.dir, "lib", "spatialcellchat-commun-flow-0.1.0")
check("source widget saves its local frontend assets",
      all(file.exists(file.path(dependency.dir, c(
        "spatialcellchat-commun-flow.js",
        "spatialcellchat-commun-flow.css"
      )))) &&
        grepl("spatialcellchat-commun-flow.js", rendered.html, fixed = TRUE) &&
        grepl("spatialcellchat-commun-flow.css", rendered.html, fixed = TRUE))
unlink(render.dir, recursive = TRUE)
check("payload preserves cell order and Flow axis transform",
      identical(widget$x$cell$id, cells) &&
        identical(widget$x$cell$x, as.numeric(coords[, "y"])) &&
        identical(widget$x$cell$y, -as.numeric(coords[, "x"])))
check("payload includes cell vectors and current legend",
      identical(widget$x$cell$dx, rep(0, length(cells))) == FALSE &&
        identical(widget$x$cell$dy, rep(0, length(cells))) == FALSE &&
        identical(widget$x$metadata$legendTitle, "Sources") &&
        identical(widget$x$metadata$current$width, c(0.65, 1.55)) &&
        identical(widget$x$metadata$current$alpha, 0.84))
check("payload has one regular field value and mask per grid point",
      length(widget$x$field$x) >= 2L && length(widget$x$field$y) >= 2L &&
        length(widget$x$field$u) == length(widget$x$field$v) &&
        length(widget$x$field$u) == length(widget$x$field$valid) &&
        is.logical(widget$x$field$valid))
check("payload is bounded by the requested grid cap",
      length(widget$x$field$x) <= 32L && length(widget$x$field$y) <= 32L)
check("payload includes mass separate from net current magnitude",
      identical(widget$x$cell$mass, c(3, 2, 1, 2, 1, 0, 0, 0, 0)) &&
        length(widget$x$cell$magnitude) == length(widget$x$cell$mass) &&
        all(c("mass", "magnitude") %in% names(widget$x$cell)))
check("payload does not contain an edge table",
      is.null(widget$x$edges) && is.null(widget$x$cell$prob))
check("browser controls receive bounded options",
      identical(widget$x$options$particleCount, 24L) &&
        identical(widget$x$options$trailLength, 10L) &&
        identical(widget$x$options$zoom, 2.2) &&
        identical(widget$x$metadata$semantics,
                  "fixed communication vector field particle animation"))

expect_error("missing precomputed field is rejected", {
  no.field <- chat
  no.field@netP$cell$field <- NULL
  netVisual_CommunFlowWidget(no.field, signaling = "signal", cellsize = c(1, 1))
})
expect_error("unknown signaling layer is rejected", {
  netVisual_CommunFlowWidget(chat, signaling = "missing", cellsize = c(1, 1))
})
expect_error("invalid particle count is rejected", {
  netVisual_CommunFlowWidget(chat, signaling = "signal", particle.count = 5001L,
                             cellsize = c(1, 1))
})
expect_error("invalid pattern is rejected", {
  netVisual_CommunFlowWidget(chat, signaling = "signal", pattern = "sideways",
                             cellsize = c(1, 1))
})
expect_error("invalid grid cap is rejected", {
  netVisual_CommunFlowWidget(chat, signaling = "signal", field.grid.max = 257L,
                             cellsize = c(1, 1))
})
expect_error("invalid zoom is rejected", {
  netVisual_CommunFlowWidget(chat, signaling = "signal", zoom = 4.1,
                             cellsize = c(1, 1))
})

expect_error("duplicate coordinates are rejected explicitly", {
  duplicate <- chat
  duplicate@images$coordinates[2L, ] <- duplicate@images$coordinates[1L, ]
  netVisual_CommunFlowWidget(duplicate, signaling = "signal", cellsize = c(1, 1))
})

zero <- chat
zero.field <- Matrix::sparseMatrix(
  i = integer(), j = integer(), x = numeric(),
  dims = c(length(cells), length(cells)), dimnames = list(cells, cells)
)
zero@netP$cell$prob <- SparseChatArray(
  list(signal = zero.field), dimnames = list(cells, cells, "signal")
)
zero <- computeCommunField(zero, signaling.name = "signal")
zero.widget <- netVisual_CommunFlowWidget(
  zero, signaling = "signal", cellsize = c(1, 1), particle.count = 8L
)
check("all-zero field returns a bounded dormant payload",
      all(zero.widget$x$field$u == 0) &&
        all(zero.widget$x$field$v == 0) &&
        length(zero.widget$x$field$valid) == length(zero.widget$x$field$u))

cat("All communication-flow widget contract checks passed.\n")
