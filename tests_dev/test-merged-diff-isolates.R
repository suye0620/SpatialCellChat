# Check signed isolate pruning in actual merged differential network rendering.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
source("tests_dev/test-merge-spatialcellchat.R", local = FALSE)
out.dir <- "tests_dev/figures/merged/diff-isolate-edgecases"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
original.graph <- graph_from_adjacency_matrix
captured.net <- NULL
captured.graph <- NULL
captured.plot.args <- NULL
plot.igraph <- function(x, ...) {
  captured.graph <<- x
  captured.plot.args <<- list(...)
  getFromNamespace("plot.igraph", "igraph")(x, ...)
}
graph_from_adjacency_matrix <- function(adjmatrix, ...) {
  captured.net <<- adjmatrix
  original.graph(adjmatrix, ...)
}
probe <- function(case) {
  obj <- merged
  first <- obj@misc$.datasets$sample1$net$group$count
  second <- obj@misc$.datasets$sample2$net$group$count
  first[,] <- 0; second[,] <- 0
  if (case == "balanced") {
    first[1, 2] <- 1; first[2, 1] <- 1
    second[1, 1] <- 1; second[2, 2] <- 1
    retained <- rownames(first)[1:2]
  } else if (case == "no_isolates") {
    diag(second) <- 1
    retained <- rownames(first)
  } else if (case == "one_node") {
    second[1, 1] <- 1
    retained <- rownames(first)[1]
  } else if (case %in% c("isolated_middle", "isolated_default")) {
    first[1, 3] <- 1; second[3, 4] <- 2; second[4, 4] <- 3
    second[1, 1] <- 5; second[3, 1] <- 1
    retained <- rownames(first)[c(1, 3, 4)]
  } else if (case == "empty") {
    retained <- character()
  } else stop("unknown case")
  obj@misc$.datasets$sample1$net$group$count <- first
  obj@misc$.datasets$sample2$net$group$count <- second
  expected <- as.matrix(second - first)[retained, retained, drop = FALSE]
  file <- file.path(out.dir, paste0(case, ".pdf"))
  captured.net <<- NULL; captured.graph <<- NULL; captured.plot.args <<- NULL
  groups <- rownames(first)
  check.attributes <- case %in% c("isolated_middle", "isolated_default")
  colors <- setNames(c("red", "blue", "green", "orange"), groups)
  weights <- setNames(c(1, 2, 3, 8), groups)
  sizes <- setNames(c(.6, .7, .8, .9), groups)
  label.colors <- setNames(c("red", "blue", "green", "orange"), groups)
  shapes <- setNames(c("circle", "square", "circle", "square"), groups)
  relevant.warnings <- character()
  message <- tryCatch({
    grDevices::pdf(file, width = 8, height = 6)
    args <- list(object = obj, comparison = c(1, 2), measure = "count", remove.isolate = TRUE)
    if (check.attributes) {
      args$vertex.weight <- weights; args$vertex.label.cex <- sizes
      args$vertex.label.color <- label.colors; args$shape <- shapes
      if (case == "isolated_middle") args$color.use <- colors
    }
    tryCatch(withCallingHandlers(do.call(netVisual_diffInteraction, args), warning = function(w) {
      if (grepl("longer argument not a multiple", conditionMessage(w), fixed = TRUE))
        relevant.warnings <<- c(relevant.warnings, conditionMessage(w))
    }), finally = grDevices::dev.off())
    stopifnot(identical(as.matrix(captured.net), expected), !length(relevant.warnings))
    if (check.attributes) {
      expected.colors <- if (case == "isolated_default") setNames(scPalette(length(groups)), groups) else colors
      stopifnot(identical(igraph::V(captured.graph)$name, retained),
                identical(unname(igraph::V(captured.graph)$color), unname(expected.colors[retained])),
                isTRUE(all.equal(unname(igraph::V(captured.graph)$size), unname(weights[retained] / max(weights) * 15 + 5))),
                identical(unname(igraph::V(captured.graph)$label.color), unname(label.colors[retained])),
                identical(unname(igraph::V(captured.graph)$label.cex), unname(sizes[retained])),
                identical(unname(captured.plot.args$vertex.shape), unname(shapes[retained])))
    }
    ""
  }, error = conditionMessage)
  if (nzchar(message) && file.exists(file)) unlink(file)
  data.frame(case = case, retained = length(retained), status = if (nzchar(message)) "ERROR" else "PDF-WRITTEN", detail = message)
}
summary <- do.call(rbind, lapply(c("balanced", "no_isolates", "one_node", "isolated_middle", "isolated_default", "empty"), probe))
write.csv(summary, file.path(out.dir, "summary.csv"), row.names = FALSE)
print(summary)
if (any(summary$status[1:5] != "PDF-WRITTEN") || summary$status[6] != "ERROR" ||
    !grepl("No nonzero differential interactions", summary$detail[6], fixed = TRUE)) quit(status = 1L)
