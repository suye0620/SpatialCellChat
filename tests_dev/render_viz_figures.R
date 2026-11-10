# Figure runner for the visualization migration.
#
# Renders each registered visualization function on each dataset and saves one
# PDF per (dataset, function) under tests_dev/figures/<phase>/<dataset>/.
# A failing figure does not abort the run; the summary decides the exit code.
#
# Usage:
#   Rscript tests_dev/render_viz_figures.R                     # all datasets, all registered figures
#   Rscript tests_dev/render_viz_figures.R --dataset=slideseq  # one dataset
#   Rscript tests_dev/render_viz_figures.R --phase=P0          # one phase
#   Rscript tests_dev/render_viz_figures.R --only=spatialC2CDistPlot

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
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp); library(ggplot2)
  library(colorspace); library(igraph); library(ggnetwork); library(circlize)
  library(ComplexHeatmap); library(ggalluvial)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/visualization.R")
source("R/visualization_widget.R")
source("R/modeling.R")
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
DATASET.filter <- get_arg("dataset", "all")
PHASE.filter <- get_arg("phase", "all")
ONLY.filter <- get_arg("only", NULL)
FRESH.fixture <- identical(get_arg("fresh", "false"), "true")

for (f in list.files("tests_dev/fixtures", pattern = "[.]R$", full.names = TRUE))
  source(f)

FIG.root <- "tests_dev/figures"
dir.create(FIG.root, showWarnings = FALSE, recursive = TRUE)

# dataset chat builders (cached; built on first use)
load_chat <- function(dataset) {
  switch(dataset,
         simulated = make_simulated_chat(cache = !FRESH.fixture),
         slideseq = if (FRESH.fixture) stop("Fresh Slide-seq fixture requires an explicit full inference workflow") else make_slideseq_chat(),
         visium = make_visium_chat(cache = !FRESH.fixture,
                                   distance.only = FRESH.fixture && identical(PHASE.filter, "P0")),
         stop("unknown dataset: ", dataset, call. = FALSE))
}
DATASETS <- c("simulated", "slideseq", "visium")

# figure registry: phase, dataset, figure name, expression (evaluated with `chat` bound)
FIGURES <- list(
  list(phase = "P0", dataset = "simulated", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", mode = "communication"))),
  list(phase = "P0", dataset = "simulated", name = "spatialC2CDistPlot-geometry",
       expr = quote(spatialC2CDistPlot(chat, mode = "geometry"))),
  list(phase = "P0", dataset = "simulated", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat, image = FALSE))),
  list(phase = "P0", dataset = "slideseq", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", mode = "communication"))),
  list(phase = "P0", dataset = "slideseq", name = "spatialC2CDistPlot-geometry",
       expr = quote(spatialC2CDistPlot(chat, mode = "geometry"))),
  list(phase = "P0", dataset = "slideseq", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat, image = FALSE))),
  list(phase = "P0", dataset = "visium", name = "spatialC2CDistPlot",
       expr = quote(spatialC2CDistPlot(chat, split.by = "annotation", mode = "communication"))),
  list(phase = "P0", dataset = "visium", name = "spatialC2CDistPlot-geometry",
       expr = quote(spatialC2CDistPlot(chat, mode = "geometry"))),
  list(phase = "P0", dataset = "visium", name = "spatialDimPlot",
       expr = quote(spatialDimPlot(chat)))
)

for (dataset in DATASETS) {
  FIGURES <- c(FIGURES, list(
    list(phase = "P1", dataset = dataset, name = "spatialFeaturePlot",
         expr = quote(spatialFeaturePlot(state$chat,
                      features = rownames(state$chat@assay$signaling)[1]))),
    list(phase = "P1", dataset = dataset, name = "spatialVisual_scoring",
         expr = quote({ensure_centrality(); spatialVisual_scoring(
           state$chat, signaling = state$signalings[1], do.group = FALSE)})),
    list(phase = "P1", dataset = dataset, name = "spatialTopicPlot",
         expr = quote({ensure_topics(); spatialTopicPlot(
           state$chat, slot.name = "netP", pattern = "outgoing")})),
    list(phase = "P1", dataset = dataset, name = "spatialGiPlot",
         expr = quote({ensure_centrality(); spatialGiPlot(
           state$chat, signaling.name = state$signalings[1], do.binary = FALSE,
           return.object = FALSE, show.plot = FALSE)})),
    list(phase = "P1", dataset = dataset, name = "spatialLeePlot",
         expr = quote(spatialLeePlot(state$chat, group.by2 = state$chat@idents,
                                     group.by1 = state$chat@idents, type = "heatmap")))
  ))
}

FIGURES <- c(FIGURES, list(
  list(phase = "P3", dataset = "slideseq", name = "netVisual_CommunFieldGrid-adjusted",
       expr = quote({ensure_field(); state$note <- paste(state$note, "point.size=0.22; image.alpha=0.4", sep="; ");
                     netVisual_CommunFieldGrid(state$chat, signaling = state$field.signal,
                                               pattern = "outgoing", point.size = 0.22, image.alpha = 0.4)})),
  list(phase = "P3", dataset = "slideseq", name = "netVisual_CommunFlow-adjusted",
       expr = quote({ensure_field(); state$note <- paste(state$note, "point.size=0.22; image.alpha=0.4; res=8", sep="; ");
                     netVisual_CommunFlow(state$chat, signaling = state$field.signal,
                                          pattern = "outgoing", point.size = 0.22, image.alpha = 0.4, res = 8)}))
))

for (dataset in DATASETS) {
  FIGURES <- c(FIGURES, list(
    list(phase = "P3", dataset = dataset, name = "netVisual_CommunField",
         expr = quote({ensure_field(); netVisual_CommunField(state$chat, signaling = state$field.signal,
                                                              pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netVisual_CommunFieldZoomIn",
         expr = quote({ensure_field(); netVisual_CommunFieldZoomIn(state$chat, signaling = state$field.signal,
                                                                    pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netVisual_CommunFieldGrid",
         expr = quote({ensure_field(); netVisual_CommunFieldGrid(state$chat, signaling = state$field.signal,
                                                                  pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netVisual_CommunFlow",
         expr = quote({ensure_field(); netVisual_CommunFlow(state$chat, signaling = state$field.signal,
                                                             pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netAnalysis_river",
         expr = quote({ensure_pattern(); netAnalysis_river(state$chat, slot.name = "netP", pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netAnalysis_dot",
         expr = quote({ensure_pattern(); netAnalysis_dot(state$chat, slot.name = "netP", pattern = "outgoing")})),
    list(phase = "P3", dataset = dataset, name = "netVisual_embedding",
         expr = quote({ensure_embedding(); netVisual_embedding(state$chat, type = "functional", pathway.remove = character())})),
    list(phase = "P3", dataset = dataset, name = "netVisual_embeddingZoomIn",
         expr = quote({ensure_embedding(); netVisual_embeddingZoomIn(state$chat, type = "functional", pathway.remove = character())})),
    list(phase = "P3", dataset = dataset, name = "netVisual_TopicComposition",
         expr = quote({ensure_topics(); netVisual_TopicComposition(state$chat,
                       y.group = t(state$chat@netP$topic$outgoing$cell), type = "dot")})),
    list(phase = "P3", dataset = dataset, name = "netVisual_TopicSignaling",
         expr = quote({ensure_topics(); netVisual_TopicSignaling(state$chat, slot.name = "netP",
                       pattern = "outgoing", feature_cat = "signaling")})),
    list(phase = "P3", dataset = dataset, name = "plotGeneExpression",
         expr = quote(plotGeneExpression(state$chat, features = rownames(state$chat@assay$signaling)[1], type = "dot")))
  ))
}

for (dataset in DATASETS) {
  FIGURES <- c(FIGURES, list(
    list(phase = "P2", dataset = dataset, name = "netVisual_heatmap",
         base.device = TRUE, expr = quote(netVisual_heatmap(state$chat))),
    list(phase = "P2", dataset = dataset, name = "netVisual_barplot",
         expr = quote(netVisual_barplot(state$chat))),
    list(phase = "P2", dataset = dataset, name = "netVisual_bubble",
         expr = quote(netVisual_bubble(state$chat, signaling = state$chat@netP$pathways[1]))),
    list(phase = "P2", dataset = dataset, name = "netVisual_chord_cell",
         base.device = TRUE,
         expr = quote(netVisual_chord_cell(state$chat, signaling = state$chat@netP$pathways[1]))),
    list(phase = "P2", dataset = dataset, name = "netVisual_chord_gene",
         base.device = TRUE,
         expr = quote(netVisual_chord_gene(state$chat, signaling = state$chat@netP$pathways[1], show.legend = FALSE)))
  ))
  for (layout in c("circle", "hierarchy", "spatial", "chord")) {
    FIGURES <- c(FIGURES, list(list(
      phase = "P2", dataset = dataset, name = paste0("netVisual-export-", layout),
      export = TRUE,
      expr = substitute({
        counts <- table(state$chat@LR$LRsig$pathway_name)
        candidates <- intersect(names(sort(counts)), state$chat@netP$pathways)
        selected <- NULL
        for (candidate in candidates) {
          selected <- tryCatch(.sc_viz_group_lr(state$chat, candidate, .05), error = function(e) NULL)
          if (!is.null(selected)) break
        }
        if (is.null(selected)) stop("No pathway has plottable group LR layers")
        cat("   export pathway:", candidate, "| LR layers:", length(selected$names), "\n")
        state$note <- paste0("pathway=", candidate, "; LR layers=", length(selected$names))
        receiver <- seq_len(max(1L, floor(nlevels(state$chat@idents) / 2L)))
        netVisual(state$chat, signaling = candidate, signaling.name = EXPORT_NAME,
                  layout = LAYOUT, vertex.receiver = receiver, out.format = "pdf")
      }, list(EXPORT_NAME = paste0("netVisual-export-", layout), LAYOUT = layout))
    )))
    for (kind in c("aggregate", "individual")) {
      FIGURES <- c(FIGURES, list(list(
        phase = "P2", dataset = dataset, name = paste0("netVisual_", kind, "-", layout),
        base.device = TRUE,
        expr = substitute({
          signaling <- state$chat@netP$pathways[1]
          receiver <- seq_len(max(1L, floor(nlevels(state$chat@idents) / 2L)))
          FUN(state$chat, signaling = signaling, layout = LAYOUT, vertex.receiver = receiver)
        }, list(FUN = as.name(paste0("netVisual_", kind)), LAYOUT = layout))
      )))
    }
  }
}

state <- new.env(parent = globalenv())
state$dataset <- NULL
state$chat <- NULL
ensure_centrality <- function() {
  if (!is.null(state$chat@netP$centr.cell)) return(invisible(NULL))
  if (identical(state$dataset, "slideseq")) {
    options(future.globals.maxSize = 20 * 1024^3)
  }
  available <- dimnames(state$chat@netP$cell$prob)[[3]]
  if (length(available) < 2L) stop("at least two pathway layers are needed for topic rank 2")
  state$signalings <- head(available, 2L)
  state$chat <- netAnalysis_computeCentrality(state$chat, slot.name = "netP",
                                               do.group = FALSE,
                                               signaling.name = state$signalings)
  invisible(NULL)
}
ensure_topics <- function() {
  ensure_centrality()
  if (is.null(state$chat@netP$topic$outgoing$cell)) {
    state$chat <- identifyCellTopics(state$chat, slot.name = "netP",
                                    pattern = "outgoing", topic.k = 2L,
                                    verbose = FALSE)
  }
  state$note <- paste0("outgoing; topic rank=2; pathways=", paste(state$signalings, collapse = ","))
  invisible(NULL)
}
ensure_field <- function() {
  if (is.null(state$field.signal)) {
    state$field.signal <- dimnames(state$chat@netP$cell$prob)[[3L]][1L]
    state$chat <- computeCommunField(state$chat, signaling.name = state$field.signal)
  }
  state$note <- paste0("pathway=", state$field.signal, "; outgoing field")
  invisible(NULL)
}
ensure_pattern <- function() {
  if (is.null(state$chat@netP$pattern$outgoing$pattern)) {
    set.seed(1)
    state$chat <- identifyCommunicationPatterns(state$chat, slot.name = "netP", pattern = "outgoing",
                                                k = 2, heatmap.show = FALSE)
  }
  state$note <- "outgoing; NMF rank=2"
  invisible(NULL)
}
ensure_embedding <- function() {
  if (is.null(state$chat@netP$similarity$functional$group$single)) {
    n <- dim(state$chat@netP$group$prob)[3L]
    if (n < 2L) stop("Embedding needs at least two pathways")
    set.seed(1)
    state$chat <- computeNetSimilarity(state$chat, slot.name = "netP", type = "functional", k = 2)
    state$chat <- netEmbedding(state$chat, slot.name = "netP", type = "functional",
                               pathway.remove = character(), umap.method = "uwot",
                               n_neighbors = min(10L, max(2L, n - 1L)), n_threads = 1L)
    n <- nrow(state$chat@netP$similarity$functional$dr$single)
    state$embedding.k <- min(2L, n - 1L)
    if (state$embedding.k < 1L) stop("Not enough embedded pathways for clustering")
    state$chat <- netClustering(state$chat, slot.name = "netP", type = "functional",
                                k = state$embedding.k, do.plot = FALSE)
  }
  state$note <- paste0("functional; pathways=", nrow(state$chat@netP$similarity$functional$dr$single),
                       "; kmeans k=", state$embedding.k)
  invisible(NULL)
}

results <- list()
for (spec in FIGURES) {
  if (!identical(DATASET.filter, "all") && !identical(DATASET.filter, spec$dataset)) next
  if (!identical(PHASE.filter, "all") && !identical(PHASE.filter, spec$phase)) next
  if (!is.null(ONLY.filter) && !identical(ONLY.filter, spec$name)) next

  tag <- paste(spec$phase, spec$dataset, spec$name, sep = "/")
  out.dir <- file.path(FIG.root, spec$phase, spec$dataset)
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  out.pdf <- file.path(out.dir, paste0(spec$name, ".pdf"))
  cat("== ", tag, "\n", sep = "")

  status <- "PASS"
  detail <- ""
  state$note <- ""
  t0 <- Sys.time()
  if (!identical(state$dataset, spec$dataset)) {
    state$chat <- NULL
    gc()
    state$dataset <- NULL
    state$signalings <- NULL
    state$field.signal <- NULL
    state$embedding.k <- NULL
    state$chat <- tryCatch(load_chat(spec$dataset),
                           error = function(e) { status <<- "FIXTURE-FAIL"; detail <<- conditionMessage(e); NULL })
    if (!is.null(state$chat)) state$dataset <- spec$dataset
  }
  chat <- state$chat
  p <- NULL
  if (!is.null(chat)) {
    p <- tryCatch({
      if (isTRUE(spec$base.device)) {
        width <- 9
        height <- 7
        if (grepl("^netVisual_individual-(circle|hierarchy|chord)$", spec$name)) {
          selected <- .sc_viz_group_lr(chat, chat@netP$pathways[1], 0.05, NULL)
          n <- length(selected$names)
          rows <- if (grepl("-hierarchy$", spec$name)) n else ceiling(n / min(n, 2L))
          height <- max(height, 3.5 * rows)
          if (height > 200 && rows * 2.2 <= 200) height <- rows * 2.2
          if (height > 200) stop("complete LR panel exceeds PDF device height limit: ", height, " inches")
          cat("   complete LR panels:", n, "| PDF size:", width, "x", height, "inches\n")
        }
        grDevices::pdf(out.pdf, width = width, height = height)
      }
      env <- list(chat = chat, state = state, ensure_centrality = ensure_centrality,
                  ensure_topics = ensure_topics, ensure_field = ensure_field,
                  ensure_pattern = ensure_pattern, ensure_embedding = ensure_embedding)
      result <- if (isTRUE(spec$export)) {
        old.dir <- getwd()
        tryCatch({
          setwd(normalizePath(out.dir))
          eval(spec$expr, envir = env)
        }, finally = setwd(old.dir))
      } else eval(spec$expr, envir = env)
      if (isTRUE(spec$base.device)) {
        plots <- if (is.list(result) && !inherits(result, "ggplot")) result else list(result)
        for (plot in plots) if (inherits(plot, "ggplot")) print(plot)
      }
      result
    }, error = function(e) {
      status <<- "FAIL"
      detail <<- conditionMessage(e)
      NULL
    }, finally = {
      if (isTRUE(spec$base.device) && grDevices::dev.cur() > 1L) grDevices::dev.off()
    })
  }
  if (isTRUE(spec$export) && status == "PASS") {
    layout <- sub("^netVisual-export-", "", spec$name)
    outputs <- file.path(out.dir, paste0(spec$name, "_", layout, "_", c("individual", "aggregate"), ".pdf"))
    if (!all(file.exists(outputs)) || any(file.info(outputs)$size <= 1000L)) {
      status <- "SAVE-FAIL"; detail <- "One or both exported PDFs are missing or empty"
    } else cat("   saved:", paste(outputs, collapse = ", "), "\n")
  } else if (isTRUE(spec$base.device) && status == "PASS") {
    if (!file.exists(out.pdf) || file.info(out.pdf)$size <= 1000L) {
      status <- "SAVE-FAIL"; detail <- "PDF was not written"
    } else cat("   saved:", out.pdf, "\n")
  } else if (inherits(p, "ggplot")) {
    tryCatch({
      width <- 7
      height <- 5.5
      if (identical(spec$name, "netVisual_bubble")) {
        if (!any(!is.na(p$data$interaction_name_2))) stop("bubble has no visible LR labels")
        width <- max(width, 3 + .3 * length(unique(stats::na.omit(p$data$source.target))))
        height <- max(height, 2 + .25 * length(unique(stats::na.omit(p$data$interaction_name_2))))
        cat("   bubble PDF size:", width, "x", height, "inches\n")
      }
      withCallingHandlers(
        ggplot2::ggsave(out.pdf, p, width = width, height = height, dpi = 150, limitsize = FALSE),
        error = function(e) if (identical(Sys.getenv("SC_VIZ_TRACE"), "1")) {
          print(summary(p$data))
          print(sys.calls())
        }
      )
      cat("   saved:", out.pdf, "\n")
    }, error = function(e) {
      status <<- "SAVE-FAIL"; detail <<- conditionMessage(e)
    })
  } else if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
    tryCatch({
      grDevices::pdf(out.pdf, width = 7, height = 5.5)
      tryCatch(ComplexHeatmap::draw(p), finally = grDevices::dev.off())
      cat("   saved:", out.pdf, "\n")
    }, error = function(e) { status <<- "SAVE-FAIL"; detail <<- conditionMessage(e) })
  } else if (status == "PASS") {
    status <- "FAIL"; detail <- "expression did not return a supported plot object"
  }
  if (status == "PASS" && identical(spec$dataset, "visium") &&
      identical(spec$name, "netVisual_chord_gene")) {
    status <- "VISUAL-REVIEW"
    detail <- "complete 88-LR pathway; sector labels overlap"
  } else if (status == "PASS" && identical(spec$dataset, "slideseq") &&
      spec$name %in% c("netVisual_CommunFieldGrid", "netVisual_CommunFlow")) {
    status <- "VISUAL-REVIEW"
    detail <- paste(state$note, "default dense background obscures vectors; flow integration warning requires review", sep = "; ")
  } else if (status == "PASS" && nzchar(state$note)) detail <- state$note
  cat("   -> ", status, if (nzchar(detail)) paste0(" | ", detail) else "", " (",
      round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s)\n", sep = "")
  results[[tag]] <- list(status = status, detail = detail, pdf = out.pdf)
}

df <- do.call(rbind, lapply(names(results), function(k) {
  data.frame(figure = k, status = results[[k]]$status, detail = results[[k]]$detail)
}))
cat("\n==== summary ====\n")
print(df, right = FALSE)
summary.name <- paste(c("summary", DATASET.filter, PHASE.filter,
                        if (!is.null(ONLY.filter)) ONLY.filter,
                        if (FRESH.fixture) "fresh"), collapse = "-")
utils::write.csv(df, file.path(FIG.root, paste0(summary.name, ".csv")), row.names = FALSE)
fails <- sum(df$status != "PASS")
cat("total:", nrow(df), "| pass:", nrow(df) - fails, "| fail:", fails, "\n")
if (fails > 0) quit(status = 1)
