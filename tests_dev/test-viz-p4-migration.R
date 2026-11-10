# Regression checks for interactive expression visualization on the new schema.
setwd(local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr)
  library(ggplot2); library(colorspace); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
source("tests_dev/fixtures/make_simulated_chat.R")

check <- function(label, result) {
  if (!isTRUE(result)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
expect_error <- function(label, expr) check(label, inherits(try(force(expr), silent = TRUE), "try-error"))
chat <- make_simulated_chat()
before <- chat
gene <- rownames(assay(chat, "norm"))[1L]
plots <- plotly_spatialLRpairPlot(chat, features = gene)
plot <- plotly::plotly_build(plots[[1L]])
traces <- plot$x$data
ids <- unlist(lapply(traces, function(x) x$customdata), use.names = FALSE)
check("feature builds an interactive trace for every cell", setequal(ids, rownames(chat@meta)))
check("plot leaves the input unchanged", identical(chat, before))
check("do.plot=FALSE returns the feature count", identical(plotly_spatialLRpairPlot(chat, features = gene, do.plot = FALSE), 1L))
html <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(plot, html, selfcontained = FALSE)
check("feature HTML is written", file.exists(html) && file.info(html)$size > 1000L)

reordered <- chat
reordered@images$coordinates <- reordered@images$coordinates[rev(seq_len(nrow(chat@meta))), , drop = FALSE]
plot2 <- plotly::plotly_build(plotly_spatialLRpairPlot(reordered, features = gene)[[1L]])
points <- plot2$x$data[[1L]]
expected <- as.numeric(assay(chat, "norm")[gene, points$customdata])
check("expression follows reordered cell names", isTRUE(all.equal(as.numeric(points$marker$color), expected)))
metadata.chat <- reordered
metadata.chat@meta$test_score <- seq_len(nrow(metadata.chat@meta))
metadata.plot <- plotly::plotly_build(plotly_spatialLRpairPlot(metadata.chat, features = "test_score")[[1L]])
metadata.points <- metadata.plot$x$data[[1L]]
check("numeric metadata follows reordered cell names",
      identical(as.numeric(metadata.points$marker$color),
                as.numeric(metadata.chat@meta[metadata.points$customdata, "test_score"])))
dim.plot <- plotly::plotly_build(plotly_spatialDimPlot(reordered))
dim.traces <- dim.plot$x$data
for (trace in dim.traces) {
  if (!length(trace$customdata)) next
  check(paste0("dimension labels follow cell names: ", trace$name),
        all(as.character(chat@idents[trace$customdata]) == trace$name))
}
bad <- chat
rownames(bad@images$coordinates)[1L] <- "missing-cell"
expect_error("mismatched cell names are rejected", plotly_spatialLRpairPlot(bad, features = gene))
expect_error("unknown feature is rejected", plotly_spatialLRpairPlot(chat, features = "unknown-feature"))
expect_error("multiple LR pairs are rejected", plotly_spatialLRpairPlot(chat, pairLR.use = c("G1-R1", "G2-R2")))
pair <- rownames(chat@LR$LRsig)[1L]
check("LR pair resolves ligand and receptor", length(plotly_spatialLRpairPlot(chat, pairLR.use = pair)) == 2L)
gated <- chat
pair <- intersect(dimnames(gated@net$cell$prob)[[3L]], dimnames(gated@net$group$prob)[[3L]])[1L]
zero <- function(m) { m@x[] <- 0; m }
cell.layers <- unclass(gated@net$cell$prob)
cell.index <- match(pair, dimnames(gated@net$cell$prob)[[3L]])
cell.layers[[cell.index]] <- zero(cell.layers[[cell.index]])
gated@net$cell$prob <- SparseChatArray(cell.layers, dimnames = dimnames(gated@net$cell$prob))
expect_error("explicit sparse zero is not active", plotly_spatialLRpairPlot(gated, pairLR.use = pair, enriched.only = TRUE))
cell.layers[[cell.index]][1L, 2L] <- 0.5
gated@net$cell$prob <- SparseChatArray(cell.layers, dimnames = dimnames(gated@net$cell$prob))
check("positive cell communication is active", length(plotly_spatialLRpairPlot(gated, pairLR.use = pair, enriched.only = TRUE)) == 2L)
group.layers <- unclass(gated@net$group$prob)
p.layers <- unclass(gated@net$group$pval)
group.index <- match(pair, dimnames(gated@net$group$prob)[[3L]])
group.layers[[group.index]] <- zero(group.layers[[group.index]])
group.layers[[group.index]][1L, 2L] <- 0.5
p.layers[[group.index]][1L, 2L] <- 0.05
gated@net$group$prob <- SparseChatArray(group.layers, dimnames = dimnames(gated@net$group$prob))
gated@net$group$pval <- SparseChatArray(p.layers, dimnames = dimnames(gated@net$group$pval))
check("group p-value boundary remains inclusive", length(plotly_spatialLRpairPlot(gated, pairLR.use = pair, enriched.only = TRUE, do.group = TRUE, thresh = 0.05)) == 2L)
expect_error("group p-value above threshold is excluded", plotly_spatialLRpairPlot(gated, pairLR.use = pair, enriched.only = TRUE, do.group = TRUE, thresh = 0.049))
coords <- as.data.frame(chat@images$coordinates)
scores <- data.frame(x = coords[, 1], y = coords[, 2], score = seq_len(nrow(coords)),
                     row.names = rownames(coords))
feature.3d <- plotly::plotly_build(plotFeatures(scores, method = "3d",
  cells.highlight = rownames(scores)[1:3], three.dim.z = 7, plot.title = "score"))
feature.traces <- feature.3d$x$data[1:2]
check("3D feature traces preserve all cell ids and highlight membership",
      setequal(unlist(lapply(feature.traces, `[[`, "customdata")), rownames(scores)) &&
        setequal(feature.traces[[1]]$customdata, rownames(scores)[1:3]))
check("3D feature traces use the requested z coordinate",
      all(vapply(feature.traces, function(trace) all(trace$z == 7), logical(1))))
gi.3d <- plotly::plotly_build(plotStatistics_Gi(scores, group.by = chat@idents[rownames(scores)],
  group.highlight = seq_len(nrow(scores)) > 3L, n = 4L, method = "3d", three.dim.z = 7))
gi.traces <- gi.3d$x$data[1:2]
check("3D Gi traces preserve all points and highlight membership",
      identical(vapply(gi.traces, function(trace) length(trace$x), integer(1)),
                c(3L, nrow(scores) - 3L)))
check("3D Gi traces use the requested z coordinate",
      all(vapply(gi.traces, function(trace) all(trace$z == 7), logical(1))))
app <- plotly_spatialLRpairPlot_shiny(chat, features = gene)
shiny::testServer(app, {
  session$flushReact()
  check("Shiny starts with all cells", identical(vals$keepcells, rownames(chat@meta)))
  session$setInputs(`plotly_selecting-A` = '[{"customdata":"c1"},{"customdata":"c2"}]')
  session$setInputs(update = 1)
  check("Shiny update uses first selection", identical(vals$keepcells, c("c1", "c2")))
  session$setInputs(`plotly_selecting-A` = '[{"customdata":"c3"}]')
  session$setInputs(update = 2)
  check("Shiny update uses latest selection", identical(vals$keepcells, "c3"))
  session$setInputs(rev = 1)
  check("Shiny reverses colors", isTRUE(vals$rev))
  session$setInputs(reset = 1)
  check("Shiny reset restores all cells", identical(vals$keepcells, rownames(chat@meta)))
})
cat("All P4 migration checks passed.\n")
