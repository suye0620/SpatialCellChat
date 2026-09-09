# Smoke tests for grid sizing in communication-field visualizations.

setwd("F:/Rworkspace/SpatialCellChat")
source("renv/activate.R")
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

options(SpatialCellChat.verbose = 2L)
for (fn in c("netVisual_CommunFieldGrid", "netVisual_CommunFlow")) {
  fun <- get(fn)
  input.before.explicit <- chat
  explicit <- capture.output(
    result <- fun(chat, signaling = "signal", pattern = "outgoing",
                  cellsize = c(1, 1), grid.resolution = 1),
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

cat("All grid visualization smoke tests passed.\n")
