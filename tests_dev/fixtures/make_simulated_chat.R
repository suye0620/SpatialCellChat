# Simulated-data fixture: a small grid object taken through the full inference
# pipeline. Deterministic (seeded), fast (~5 s), and cached to
# tests_dev/test_data/simulated_chat.rds so the runner does not rebuild it.
#
# Usage (after the package sources are loaded):
#   source("tests_dev/fixtures/make_simulated_chat.R")
#   chat <- make_simulated_chat()

make_simulated_chat <- function(seed = 1, cache = TRUE,
                                cache.path = "tests_dev/test_data/simulated_chat.rds") {
  if (isTRUE(cache) && file.exists(cache.path)) {
    return(readRDS(cache.path))
  }

  set.seed(seed)
  cells <- paste0("c", seq_len(20))
  genes <- c("G1", "G2", "G3", "R1", "R2", "L3", "R3")
  expr <- matrix(round(runif(length(genes) * length(cells)), 3),
                 nrow = length(genes), dimnames = list(genes, cells))
  coor <- expand.grid(x = seq(10, 50, 10), y = seq(10, 40, 10))
  rownames(coor) <- cells
  meta.df <- data.frame(label = factor(rep(c("A", "B", "C", "D"), each = 5)),
                        row.names = cells)

  chat <- createSpatialCellChat(
    object = as(expr, "dgCMatrix"), meta = meta.df, group.by = "label",
    input.assay = "norm", datatype = "spatial",
    coordinates = coor, spatial.factors = list(ratio = 1, tol = 5)
  )
  chat@DB <- list(
    interaction = data.frame(
      ligand = c("G1", "G2", "L3"), receptor = c("R1", "R2", "R3"),
      agonist = c("", "", ""), antagonist = c("", "", ""),
      co_A_receptor = c("", "", ""), co_I_receptor = c("", "", ""),
      annotation = c("Secreted Signaling", "Cell-Cell Contact", "Secreted Signaling"),
      interaction_name = c("G1-R1", "G2-R2", "L3-R3"),
      pathway_name = c("P1", "P2", "P3"), stringsAsFactors = FALSE
    ),
    complex = data.frame(), cofactor = data.frame()
  )
  chat@LR$LRsig <- chat@DB$interaction[, c("ligand", "receptor", "agonist", "antagonist",
                                           "co_A_receptor", "co_I_receptor", "annotation",
                                           "interaction_name", "pathway_name")]
  rownames(chat@LR$LRsig) <- chat@LR$LRsig$interaction_name
  assay(chat, "signaling") <- assay(chat, "norm")
  assay(chat, "smooth") <- assay(chat, "signaling")

  chat <- computeCommunProb(chat, scale.distance = 9, n = 1, nthreads = 1L,
                            verbose = FALSE)
  chat <- computeAvgCommunProb(chat, min.cells.sr = 1, do.permutation = TRUE,
                               nboot = 20, seed.use = 1L, nthreads = 1L,
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
