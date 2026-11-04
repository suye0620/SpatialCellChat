# Slide-seq fixture: a 900-cell window of the MBM05 puck taken through the full
# inference pipeline (proven recipe: interaction.range = 250, scale.distance = 9,
# n = 1). Cached to tests_dev/test_data/slideseq_window_chat.rds.
#
# Usage (after the package sources are loaded):
#   source("tests_dev/fixtures/make_slideseq_window.R")
#   chat <- make_slideseq_chat()

make_slideseq_chat <- function(
    rds = "tests_dev/test_data/MBM05Chat_pre_inference.rds",
    half = 380, max.cells = 900, lr.per.annotation = 40, seed = 1,
    cache = TRUE, cache.path = "tests_dev/test_data/slideseq_window_chat.rds") {
  if (isTRUE(cache) && file.exists(cache.path)) {
    return(readRDS(cache.path))
  }

  obj <- readRDS(rds)
  co <- obj@images$coordinates
  sf <- obj@images$spatial.factors

  cx <- median(co[, 1]); cy <- median(co[, 2])
  keep <- which(abs(co[, 1] - cx) <= half & abs(co[, 2] - cy) <= half)
  set.seed(seed)
  if (length(keep) > max.cells) keep <- sort(sample(keep, max.cells))

  lr <- obj@LR$LRsig
  picks <- unlist(lapply(split(seq_len(nrow(lr)), lr$annotation), function(ix) {
    sample(ix, min(length(ix), lr.per.annotation))
  }), use.names = FALSE)
  LR.use <- lr[sort(picks), , drop = FALSE]

  chat <- createSpatialCellChat(
    object = obj@assay$norm[, keep, drop = FALSE],
    meta = obj@meta[keep, , drop = FALSE],
    group.by = "rctd_cell_type",
    input.assay = "norm", datatype = "spatial",
    coordinates = co[keep, , drop = FALSE],
    spatial.factors = sf
  )
  chat@DB <- obj@DB
  chat@LR <- list(LRsig = LR.use)
  assay(chat, "signaling") <- obj@assay$signaling[, keep, drop = FALSE]
  assay(chat, "smooth") <- assay(chat, "signaling")

  chat <- computeCommunProb(chat, LR.use = LR.use, n = 1, interaction.range = 250,
                            scale.distance = 9, nthreads = 4L, verbose = FALSE)
  chat <- computeAvgCommunProb(chat, min.cells.sr = 1, do.permutation = TRUE,
                               nboot = 20, seed.use = 1L, nthreads = 4L,
                               verbose = FALSE)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = TRUE)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = TRUE)

  if (isTRUE(cache)) {
    dir.create(dirname(cache.path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(chat, cache.path)
  }
  chat
}
