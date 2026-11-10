# End-to-end remove.isolate regression on the three canonical fixtures.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(ggplot2)
  library(igraph); library(colorspace); library(ggnetwork); library(spatstat.sparse)
})
source("R/SpatialCellChat_class.R"); source("R/database.R"); source("R/utilities.R")
source("R/spatial.R"); source("R/modeling.R"); source("R/analysis.R"); source("R/visualization.R")
source("tests_dev/fixtures/make_simulated_chat.R"); source("tests_dev/fixtures/make_visium_pp1.R")
source("tests_dev/fixtures/make_slideseq_window.R")
results <- list()
for (dataset in c("simulated", "visium", "slideseq")) {
  chat <- switch(dataset, simulated = make_simulated_chat(), visium = make_visium_chat(), slideseq = make_slideseq_chat())
  chosen <- NULL
  for (signal in unique(chat@LR$LRsig$pathway_name)) {
    selection <- tryCatch(.sc_viz_group_lr(chat, signal, 0.05), error = function(e) NULL)
    if (is.null(selection)) next
    for (i in seq_along(selection$names)) {
      layer <- selection$prob[, , i, drop = FALSE][, , 1L]
      offdiag <- layer; diag(offdiag) <- 0
      active <- which(rowSums(layer != 0) + colSums(layer != 0) > 0)
      if (length(active) >= 2 && any(offdiag != 0)) {
        chosen <- list(signal = signal, lr = selection$names[i], receiver = active[1L]); break
      }
    }
    if (!is.null(chosen)) break
  }
  if (is.null(chosen)) stop("No pathway with a multi-group, non-autocrine LR layer in ", dataset)
  cat("Dataset:", dataset, "cells:", nrow(chat@meta), "pathway:", chosen$signal, "LR:", chosen$lr, "\n")
  out.dir <- file.path("tests_dev/figures/P2", dataset, "network-isolate-regression")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  for (mode in c("aggregate", "individual")) for (layout in c("circle", "spatial", "hierarchy")) {
    file <- file.path(out.dir, paste0(mode, "-", layout, ".pdf"))
    warning.messages <- character()
    detail <- tryCatch(withCallingHandlers({
      args <- list(object = chat, signaling = chosen$signal, layout = layout,
                   remove.isolate = TRUE, vertex.receiver = chosen$receiver, vertex.weight = NULL)
      if (mode == "individual") args$pairLR.use <- chosen$lr
      grDevices::pdf(file, width = 10, height = 6)
      tryCatch({
        value <- do.call(get(paste0("netVisual_", mode)), args)
        plots <- if (inherits(value, "ggplot")) list(value) else if (is.list(value) && !inherits(value, "recordedplot")) value else list()
        for (plot in plots) if (inherits(plot, "ggplot")) print(plot)
      }, finally = grDevices::dev.off())
      if (!file.exists(file) || file.info(file)$size < 1000L) stop("missing or empty PDF")
      ""
    }, warning = function(w) {
      warning.messages <<- c(warning.messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = conditionMessage)
    if (nzchar(detail) && file.exists(file)) unlink(file)
    results[[length(results) + 1L]] <- data.frame(dataset = dataset, cells = nrow(chat@meta), signal = chosen$signal,
      lr = chosen$lr, mode = mode, layout = layout, status = if (nzchar(detail)) "ERROR" else "PDF-WRITTEN", detail = detail,
      warnings = paste(unique(warning.messages), collapse = " | "), file = file)
  }
  rm(chat); gc()
}
summary <- do.call(rbind, results)
write.csv(summary, "tests_dev/figures/P2/network-isolate-three-fixtures-2026-11-10.csv", row.names = FALSE)
print(summary)
if (any(summary$status != "PDF-WRITTEN")) quit(status = 1L)
