# Smoke tests for grid sizing in communication-field visualizations.

setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/analysis.R")
source("R/visualization.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}


# Check the user-visible help signatures against the actual function formals.
for (fn in c("netVisual_CommunFieldGrid", "netVisual_CommunFlow")) {
  rd <- tools::parse_Rd(file.path("man", paste0(fn, ".Rd")))
  usage <- rd[[which(vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\usage"), logical(1)))]]
  call <- parse(text = paste(as.character(usage), collapse = ""))[[1L]]
  args <- as.list(call)[-1L]
  formal.args <- formals(get(fn))
  arg.names <- names(args)
  for (i in which(!nzchar(arg.names))) arg.names[i] <- as.character(args[[i]])
  check(paste0(fn, " help documents argument names and order"),
        identical(arg.names, names(formal.args)))
  check(paste0(fn, " help documents actual defaults"),
        all(vapply(which(nzchar(names(args))), function(i)
          identical(args[[i]], formal.args[[arg.names[i]]]), logical(1))))
  tools::checkRd(rd)
}


# Compare the same membership indices against the previous sf-subsetting
# implementation, including empty and multiply assigned boundary positions.
field.fixture <- sf::st_as_sf(
  data.frame(x = 1:5, y = 1:5, dx = c(0, 2, NA, Inf, -1),
             dy = c(NA, -2, 0, -Inf, 4)), coords = c("x", "y"))
groups.fixture <- list(integer(), 1L, c(1L, 2L), c(2L, 3L, 5L), 3L, c(4L, 5L))
for (stat in list(median = stats::median, sum = base::sum)) {
  expected.u <- vapply(groups.fixture, function(index)
    stat(field.fixture[index, , drop = FALSE]$dx, na.rm = TRUE), numeric(1))
  expected.v <- vapply(groups.fixture, function(index)
    stat(field.fixture[index, , drop = FALSE]$dy, na.rm = TRUE), numeric(1))
  actual <- .sc_grid_field_components(field.fixture, groups.fixture, stat)
  check("grid field aggregation preserves sf-subsetting results including empty and missing values",
        identical(actual$u, expected.u) && identical(actual$v, expected.v))
}

required <- c("oce", "metR")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("[SKIP] visualization execution requires missing package(s): ",
      paste(missing, collapse = ", "), "\n", sep = "")
  check("shared resolver remains available", exists(".sc_resolve_grid_size", mode = "function"))
  check("shared membership remains available", exists(".sc_grid_membership", mode = "function"))
  quit(save = "no", status = 0L)
}

n <- 9L
cells <- paste0("c", seq_len(n))
coords <- cbind(x = rep(0:2, each = 3), y = rep(0:2, 3))
rownames(coords) <- cells
expr <- Matrix::Matrix(matrix(1, nrow = 3L, ncol = n,
                              dimnames = list(paste0("g", 1:3), cells)), sparse = TRUE)
meta <- data.frame(label = factor(rep(c("A", "B"), length.out = n)), row.names = cells)
chat <- createSpatialCellChat(
  expr, meta = meta, group.by = "label", datatype = "spatial",
  coordinates = coords, spatial.factors = list(ratio = 1, tol = 0.5)
)
prob_matrix <- Matrix::sparseMatrix(
  i = c(1L, 2L, 3L), j = c(2L, 3L, 4L), x = c(3, 2, 1),
  dims = c(n, n), dimnames = list(cells, cells)
)
prob_sparse <- SparseChatArray(
  list(signal = prob_matrix),
  dimnames = list(cells, cells, "signal")
)
chat@netP <- list(cell = list(prob = prob_sparse), group = list())
chat_with_field <- computeCommunField(chat, signaling.name = "signal")
check("computeCommunField stores cell-level sparse fields",
      inherits(chat_with_field@netP$cell$field$outgoing, "SparseChatArray") &&
      inherits(chat_with_field@netP$cell$field$incoming, "SparseChatArray"))
check("computeCommunField does not mutate its input",
      is.null(chat@netP$cell$field))
chat <- chat_with_field
field.matrix <- as(chat@netP$cell$field$outgoing[["signal"]], "dgCMatrix")
chat@netP$cell$field$outgoing <- SparseChatArray(
  list(signal = field.matrix), dimnames = list(cells, colnames(field.matrix), "signal")
)

for (fn in c("netVisual_CommunField", "netVisual_CommunFieldZoomIn")) {
  result <- get(fn)(chat, signaling = "signal", pattern = "outgoing")
  check(paste0(fn, " reads computed cell-level field"), inherits(result, "ggplot"))
  output <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(output, result, width = 7, height = 5.5)
  check(paste0(fn, " actually renders PDF"), file.exists(output) && file.info(output)$size > 1000L)
  unlink(output)
}

options(SpatialCellChat.verbose = 2L)
for (fn in c("netVisual_CommunFieldGrid", "netVisual_CommunFlow")) {
  fun <- get(fn)
  reordered <- chat
  reordered@images$coordinates <- reordered@images$coordinates[rev(seq_len(n)), , drop = FALSE]
  check(paste0(fn, " rejects mismatched field and coordinate order"),
        inherits(try(fun(reordered, signaling = "signal"), silent = TRUE), "try-error"))
  input.before.explicit <- chat
  explicit <- capture.output(
    result <- fun(chat, signaling = "signal", pattern = "outgoing",
                  cellsize = c(1, 1), grid.resolution = 1, arrow.line.alpha = 0.6),
    type = "message"
  )
  check(paste0(fn, " explicit cellsize returns ggplot"), inherits(result, "ggplot"))
  check(paste0(fn, " does not mutate explicit render input"),
        identical(chat, input.before.explicit))
  check(paste0(fn, " reports shared size source"),
        any(grepl("Cellsize source", explicit, fixed = TRUE)))
  built <- ggplot_build(result)$data[[2L]]
  width.column <- if (fn == "netVisual_CommunFieldGrid") "size" else "linewidth"
  check(paste0(fn, " maps magnitude to visible width"),
        width.column %in% names(built) &&
          length(unique(round(built[[width.column]], 6))) > 1L)
  if (fn == "netVisual_CommunFieldGrid") {
    check("netVisual_CommunFieldGrid adaptively decimates arrows",
          nrow(built) <= 144L && any(grepl("Arrow sampling", explicit, fixed = TRUE)))
    check("nonzero vectors retain requested opacity",
          any(built$mag > 0) && all(built$alpha[built$mag > 0] == 0.6))
    # Rebuild the same interpolated layer with the former fixed alpha to ensure
    # magnitude, decimation, positions and size scales are unchanged.
    result$layers[[2L]]$mapping$alpha <- NULL
    result$layers[[2L]]$aes_params$alpha <- 0.6
    before <- ggplot_build(result)$data[[2L]]
    shared <- c("x", "y", "dx", "dy", "mag", "norm_mag", "size")
    check("arrow sampling and nonzero arrow geometry match the former drawing",
          identical(built[, shared], before[, shared]))
  } else {
    check("netVisual_CommunFlow adaptively limits streamline seeds",
          length(unique(built$group)) <= 9L &&
            any(grepl("seeds per axis=3", explicit, fixed = TRUE)))
  }
  input.before.implicit <- chat
  implicit <- capture.output(
    result <- fun(chat, signaling = "signal", pattern = "outgoing",
                  cellsize = NULL, grid.resolution = 1),
    type = "message"
  )
  check(paste0(fn, " nearest-neighbor cellsize returns ggplot"), inherits(result, "ggplot"))
  check(paste0(fn, " does not mutate inferred render input"),
        identical(chat, input.before.implicit))
  check(paste0(fn, " reports occupancy diagnostics"),
        any(grepl("occupied", implicit, fixed = TRUE)))
}


zero.chat <- chat
zero.matrix <- field.matrix
zero.matrix@x[] <- 0
zero.chat@netP$cell$field$outgoing <- SparseChatArray(
  list(signal = zero.matrix), dimnames = list(cells, colnames(zero.matrix), "signal")
)
zero.plot <- netVisual_CommunFieldGrid(zero.chat, signaling = "signal",
                                      pattern = "outgoing", cellsize = c(1, 1),
                                      grid.resolution = 1, arrow.skip = 0)
zero.pdf <- tempfile(fileext = ".pdf")
ggplot2::ggsave(zero.pdf, zero.plot, width = 7, height = 5.5)
check("an all-zero communication field renders without visible arrow marks",
      file.exists(zero.pdf) && nrow(ggplot_build(zero.plot)$data[[2L]]) > 0 &&
        all(ggplot_build(zero.plot)$data[[2L]]$alpha == 0))
unlink(zero.pdf)

cat("All grid visualization smoke tests passed.\n")
