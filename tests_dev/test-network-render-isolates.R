# Actual PDF rendering and alignment checks for remove.isolate across network plots.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(ggplot2)
  library(igraph); library(colorspace); library(ggnetwork)
})
source("R/utilities.R"); source("R/visualization.R")
groups <- LETTERS[1:6]
colors <- setNames(c("red", "blue", "green", "orange", "purple", "brown"), groups)
weights <- setNames(1:6, groups)
coords <- data.frame(x = seq_along(groups), y = (seq_along(groups))^2)
labels <- factor(groups, levels = groups)
net <- matrix(1, 6, 6, dimnames = list(groups, groups))
diag(net) <- 0
original.graph <- graph_from_adjacency_matrix
captured.net <- NULL
captured.graph <- NULL
graph_from_adjacency_matrix <- function(adjmatrix, ...) {
  captured.net <<- adjmatrix
  original.graph(adjmatrix, ...)
}
plot.igraph <- function(x, ...) {
  captured.graph <<- x
  getFromNamespace("plot.igraph", "igraph")(x, ...)
}
out.dir <- "tests_dev/figures/P2/network-isolate-regression"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
probe <- function(method, case) {
  input <- net
  if (case == "mixed_loops") { diag(input) <- 1; input[1, 2] <- 0 }
  if (case == "isolated") { input[2, ] <- 0; input[, 2] <- 0 }
  if (case == "one_receiver") { input[2:3, ] <- 0; input[, 2:3] <- 0 }
  if (case == "default_colors") { input[2, ] <- 0; input[, 2] <- 0 }
  active <- if (case %in% c("isolated", "default_colors")) groups[-2] else if (case == "one_receiver") groups[-c(2, 3)] else groups
  expected.colors <- if (case == "default_colors") setNames(scPalette(length(groups)), groups) else colors
  expected <- input[active, active, drop = FALSE]
  args <- list(net = input, color.use = colors, vertex.weight = weights,
               remove.isolate = TRUE)
  if (case == "default_colors") args$color.use <- NULL
  if (method == "spatial") args <- c(args, list(coordinates = coords, labels = labels))
  if (method %in% c("hierarchy1", "hierarchy2")) {
    args$vertex.receiver <- c(1L, 2L, 3L)
    args$vertex.label.color <- setNames(c("red", "blue", "green", "orange", "purple", "brown"), groups)
    args$vertex.label.cex <- setNames(seq(0.5, 1, length.out = length(groups)), groups)
  }
  if (case == "filtered") { args$sources.use <- seq_along(groups); args$targets.use <- seq_along(groups) }
  fun <- get(paste0("netVisual_", method))
  file <- file.path(out.dir, paste0(method, "-", case, ".pdf"))
  captured.net <<- NULL; captured.graph <<- NULL
  relevant.warnings <- character()
  message <- tryCatch({
    grDevices::pdf(file, width = 8, height = 6)
    value <- tryCatch(withCallingHandlers(do.call(fun, args), warning = function(w) {
      if (grepl("longer argument not a multiple", conditionMessage(w), fixed = TRUE))
        relevant.warnings <<- c(relevant.warnings, conditionMessage(w))
    }), finally = grDevices::dev.off())
    stopifnot(!length(relevant.warnings))
    if (method == "spatial") {
      grDevices::pdf(file, width = 8, height = 6)
      tryCatch(print(value), finally = grDevices::dev.off())
      stopifnot(identical(as.character(value$data$node_idents), active),
                identical(unname(value$layers[[3L]]$aes_params$colour), unname(expected.colors[groups])),
                isTRUE(all.equal(unname(value$layers[[2L]]$aes_params$size),
                                 unname(as.numeric(weights[active]) / max(weights) * 15 + 1))))
    }
    if (method %in% c("circle", "spatial")) stopifnot(identical(as.matrix(captured.net), expected))
    if (method %in% c("hierarchy1", "hierarchy2")) {
      receiver <- which(active %in% groups[1:3])
      ordering <- if (method == "hierarchy1") c(receiver, setdiff(seq_along(active), receiver)) else c(setdiff(seq_along(active), receiver), receiver)
      n <- length(active)
      stopifnot(identical(unname(as.matrix(captured.net)[seq_len(n), n + seq_along(receiver), drop = FALSE]),
                          unname(expected[ordering, receiver, drop = FALSE])))
    }
    if (method != "spatial") {
      stopifnot(!is.null(captured.graph))
      vertex.names <- igraph::V(captured.graph)$name
      actual.colors <- igraph::V(captured.graph)$color
      real <- nzchar(vertex.names)
      stopifnot(identical(unname(actual.colors[real]), unname(expected.colors[vertex.names[real]])))
      real.size <- igraph::V(captured.graph)$size[real]
      expected.size <- as.numeric(weights[vertex.names[real]]) / max(weights) * 15 + if (method == "circle") 5 else 6
      stopifnot(isTRUE(all.equal(unname(real.size), unname(expected.size))))
      if (method %in% c("hierarchy1", "hierarchy2")) {
        receiver.names <- active[which(active %in% groups[1:3])]
        expected.names <- c(vertex.names[real], receiver.names)
        stopifnot(identical(unname(igraph::V(captured.graph)$label.color), unname(args$vertex.label.color[expected.names])),
                  isTRUE(all.equal(unname(igraph::V(captured.graph)$label.cex), unname(args$vertex.label.cex[expected.names]))))
      }
    }
    if (!file.exists(file) || file.info(file)$size < 1000L) stop("missing PDF")
    ""
  }, error = conditionMessage)
  if (nzchar(message) && file.exists(file)) unlink(file)
  data.frame(method = method, case = case, status = if (nzchar(message)) "ERROR" else "PDF-WRITTEN", detail = message)
}
summary <- do.call(rbind, lapply(c("circle", "spatial", "hierarchy1", "hierarchy2"), function(method) {
  cases <- if (method %in% c("hierarchy1", "hierarchy2")) c("none", "isolated", "one_receiver", "filtered", "default_colors") else if (method == "circle") c("none", "isolated", "default_colors", "mixed_loops") else c("none", "isolated", "default_colors")
  do.call(rbind, lapply(cases, function(case) probe(method, case)))
}))
# When all groups disappear, or all requested receivers disappear, report why.
zero <- net * 0
for (method in c("circle", "spatial", "hierarchy1", "hierarchy2")) {
  args <- list(net = zero, remove.isolate = TRUE)
  if (method == "spatial") args <- c(args, list(coordinates = coords, labels = labels))
  if (method %in% c("hierarchy1", "hierarchy2")) args$vertex.receiver <- 1:3
  message <- tryCatch({ do.call(get(paste0("netVisual_", method)), args); "" }, error = conditionMessage)
  stopifnot(grepl("No nonzero interactions remain", message, fixed = TRUE))
}
no.receiver <- net
no.receiver[1:3, ] <- 0; no.receiver[, 1:3] <- 0
for (method in c("hierarchy1", "hierarchy2")) {
  message <- tryCatch({ get(paste0("netVisual_", method))(no.receiver, vertex.receiver = 1:3, remove.isolate = TRUE); "" }, error = conditionMessage)
  stopifnot(grepl("No selected receiver groups remain", message, fixed = TRUE))
}
cat("Empty-network and absent-receiver errors passed.\n")

write.csv(summary, file.path(out.dir, "summary.csv"), row.names = FALSE)
print(summary)
if (any(summary$status != "PDF-WRITTEN")) quit(status = 1L)
