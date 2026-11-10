# Check pathway-centrality node sizes by group name on the canonical fixtures.
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
  signal <- if (dataset == "simulated") "P1" else "TGFb"
  groups <- dimnames(.sc_viz_group_lr(chat, signal, 0.05)$prob)[[1L]]
  original <- chat@netP$centr
  if (!all(groups %in% dimnames(original)[[2L]]) || !signal %in% dimnames(original)[[3L]]) stop("fixture lacks group centrality for ", dataset)
  permuted <- chat
  permuted@netP$centr <- original[, rev(seq_len(dim(original)[2L])), , drop = FALSE]
  out.dir <- file.path("tests_dev/figures/P2", dataset, "centrality-spatial-regression")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  for (direction in c("incoming", "outgoing")) {
    measure <- if (direction == "incoming") "indeg" else "outdeg"
    values <- original[measure, groups, signal]
    if (!all(is.finite(values)) || max(values) <= 0) stop("fixture has no usable centrality values")
    expected <- as.numeric(values) / max(values) * 5 + 1
    normal.plot <- netVisual_aggregate(chat, signaling = signal, layout = "spatial", vertex.weight = direction)
    reordered.plot <- netVisual_aggregate(permuted, signaling = signal, layout = "spatial", vertex.weight = direction)
    stopifnot(identical(as.character(reordered.plot$data$node_idents), groups),
              isTRUE(all.equal(unname(normal.plot$layers[[2L]]$aes_params$size), expected)),
              isTRUE(all.equal(unname(reordered.plot$layers[[2L]]$aes_params$size), expected)),
              identical(chat@netP$centr, original))
    file <- file.path(out.dir, paste0(direction, "-reordered.pdf"))
    grDevices::pdf(file, width = 10, height = 6)
    tryCatch(print(reordered.plot), finally = grDevices::dev.off())
    if (!file.exists(file) || file.info(file)$size < 1000L) stop("missing or empty PDF")
    results[[length(results) + 1L]] <- data.frame(dataset=dataset, cells=nrow(chat@meta), signal=signal, direction=direction, status="PDF-WRITTEN", file=file)
  }
  missing <- permuted
  dimnames(missing@netP$centr)[[2L]][1L] <- ".missing-group"
  msg <- tryCatch({netVisual_aggregate(missing, signaling=signal, layout="spatial", vertex.weight="incoming"); ""}, error=conditionMessage)
  stopifnot(grepl("centrality must contain", msg, fixed=TRUE))
  rm(chat, permuted, missing); gc()
}
summary <- do.call(rbind, results)
write.csv(summary, "tests_dev/figures/P2/centrality-spatial-three-fixtures-2026-11-10.csv", row.names=FALSE)
print(summary)
