# Render the remaining P5 analysis figures for all three migration fixtures.
setwd(local({
  args <- commandArgs(FALSE)
  file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file)) dirname(dirname(normalizePath(file))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(igraph)
  library(future.apply); library(ggplot2); library(colorspace); library(Rcpp)
  library(circlize); library(ComplexHeatmap); library(spatstat.sparse)
  library(ggalluvial); library(NMF)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/merging.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
for (f in list.files("tests_dev/fixtures", pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(TRUE)
dataset <- if (length(args)) args[[1L]] else "simulated"
if (!dataset %in% c("simulated", "visium", "slideseq"))
  stop("dataset must be simulated, visium, or slideseq", call. = FALSE)
chat <- switch(dataset,
  simulated = make_simulated_chat(),
  visium = make_visium_chat(),
  slideseq = make_slideseq_chat()
)
if (identical(dataset, "slideseq")) options(future.globals.maxSize = 20 * 1024^3)

output.root <- file.path("tests_dev", "figures", "P5", "three-fixture-2026-11-11", dataset)
dir.create(output.root, recursive = TRUE, showWarnings = FALSE)
pathways <- dimnames(chat@netP$cell$prob)[[3L]]
if (!length(pathways)) stop("fixture has no pathway probability layers", call. = FALSE)
pathway <- pathways[[1L]]
features <- head(rownames(chat@assay$signaling), 2L)
if (length(features) < 2L) stop("fixture has fewer than two signaling features", call. = FALSE)

results <- list()
record <- function(name, draw, annotation = "") {
  file <- file.path(output.root, paste0(name, ".pdf"))
  message <- tryCatch({
    if (name == "netAnalysis-signalingRole-network") {
      grDevices::pdf(file, width = 7, height = 4)
      tryCatch(draw(), finally = grDevices::dev.off())
    } else {
      plot <- draw()
      if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
        grDevices::pdf(file, width = 8, height = 6)
        tryCatch(ComplexHeatmap::draw(plot), finally = grDevices::dev.off())
      } else if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) {
        ggplot2::ggsave(file, plot, width = 8, height = 6, limitsize = FALSE)
      } else {
        stop("unsupported plot class: ", paste(class(plot), collapse = ", "))
      }
    }
    if (!file.exists(file) || file.info(file)$size <= 1000L)
      stop("PDF missing or empty")
    ""
  }, error = function(e) {
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    conditionMessage(e)
  })
  status <- if (nzchar(message)) "FAIL" else "PASS"
  cat("[", status, "] ", dataset, "/", name,
      if (nzchar(message)) paste0(" | ", message) else "", "\n", sep = "")
  results[[name]] <<- data.frame(
    figure = name, status = status,
    detail = if (nzchar(message)) message else annotation,
    file = file.path("tests_dev", "figures", "P5", "three-fixture-2026-11-11",
                     dataset, paste0(name, ".pdf")),
    stringsAsFactors = FALSE
  )
}

rank.scores <- .sc_normalize_pattern_scores(.sc_group_prob_array(chat, "net"), "outgoing")
rank.max <- min(4L, min(sum(rowSums(rank.scores) > 0), ncol(rank.scores)) - 1L)
if (rank.max < 2L) stop("fixture has insufficient active group/signaling dimensions for rank selection")
rank.range <- seq.int(2L, rank.max)
record("selectK-net", function() {
  selectK(chat, slot.name = "net", pattern = "outgoing", k.range = rank.range,
          nrun = 2L, seed.use = 11L)
}, annotation = paste0("k.range=", paste(rank.range, collapse = ":")))
record("plotGeneExpression-bar-single", function() {
  plotGeneExpression(chat, features = features[[1L]], type = "bar")
})
record("plotGeneExpression-bar-multi", function() {
  plotGeneExpression(chat, features = features, type = "bar")
})
record("netVisual-aggregate-spatial", function() {
  netVisual_aggregate(chat, signaling = pathway, layout = "spatial")
})
record("netVisual-individual-spatial", function() {
  lr <- chat@LR$LRsig$interaction_name[
    chat@LR$LRsig$pathway_name == pathway
  ][[1L]]
  patchwork::wrap_plots(netVisual_individual(
    chat, signaling = pathway, pairLR.use = lr, layout = "spatial"
  ))
})
record("netAnalysis-signalingRole-scatter", function() {
  netAnalysis_signalingRole_scatter(
    chat, signaling = pathway, do.label = FALSE, show.legend = FALSE
  )
})
record("netAnalysis-signalingRole-network", function() {
  netAnalysis_signalingRole_network(
    chat, signaling = pathway, measure = c("outdeg", "indeg"),
    measure.name = c("Sender", "Receiver")
  )
})
record("netAnalysis-signalingRole-heatmap", function() {
  netAnalysis_signalingRole_heatmap(chat, signaling = pathway, pattern = "outgoing")
})

summary <- do.call(rbind, results)
summary$detail <- ifelse(nzchar(summary$detail), summary$detail,
                         paste0("pathway=", pathway, "; features=", paste(features, collapse = ",")))
utils::write.csv(summary, file.path(output.root, "summary.csv"), row.names = FALSE)
if (any(summary$status != "PASS")) quit(status = 1L)
