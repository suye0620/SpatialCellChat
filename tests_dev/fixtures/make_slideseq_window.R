# Slide-seq fixture: the FULL MBM05 puck (29,536 cells) taken through the
# inference pipeline. Preferred source is the post-filter benchmark object
# (computeCommunProb -> filterProbability already applied); if it is missing the
# builder re-runs inference from MBM05Chat_pre_inference.rds. Cached to
# tests_dev/test_data/slideseq_full_chat.rds.
#
# Usage (after the package sources are loaded):
#   source("tests_dev/fixtures/make_slideseq_window.R")
#   chat <- make_slideseq_chat()

make_slideseq_chat <- function(
    pre = "tests_dev/test_data/MBM05Chat_pre_inference.rds",
    bench = "tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/new_postfilter_object.rds",
    lr.per.annotation = 40, seed = 1,
    cache = TRUE, cache.path = "tests_dev/test_data/slideseq_full_chat.rds") {
  if (isTRUE(cache) && file.exists(cache.path)) {
    return(readRDS(cache.path))
  }

  if (file.exists(bench)) {
    # full-puck object at the filterProbability stage: group-level results still
    # need the downstream steps (avg -> pathway -> aggregate -> centrality)
    chat <- readRDS(bench)
  } else {
    obj <- readRDS(pre)
    co <- obj@images$coordinates
    chat <- createSpatialCellChat(
      object = obj@assay$norm, meta = obj@meta, group.by = "rctd_cell_type",
      input.assay = "norm", datatype = "spatial",
      coordinates = co, spatial.factors = obj@images$spatial.factors
    )
    chat@DB <- obj@DB
    chat@LR <- obj@LR
    assay(chat, "signaling") <- obj@assay$signaling
    chat <- computeCommunProb(chat, n = 1, interaction.range = 250,
                              scale.distance = 9, nthreads = 4L, verbose = FALSE)
    chat <- filterProbability(chat)
  }

  if (is.null(chat@LR$LRsig) || !nrow(chat@LR$LRsig)) {
    stop("the slide-seq fixture has no LRsig; cannot run group-level steps",
         call. = FALSE)
  }
  # the full puck's data closures exceed future's 500 MiB default; the benchmark
  # scripts raise this the same way
  options(future.globals.maxSize = 20 * 1024^3)
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
