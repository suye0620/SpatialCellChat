# Visium fixture: the psoriasis PP1 sample (905 spots) imported the way
# tests_dev/test_data/visium_PP1.R does, then taken through inference.
# Cached to tests_dev/test_data/visium_pp1_chat.rds (heavy: ~1-2 min build).
#
# Usage (after the package sources are loaded):
#   source("tests_dev/fixtures/make_visium_pp1.R")
#   chat <- make_visium_chat()

make_visium_chat <- function(
    rdata = "tests_dev/test_data/visium_human_psoriasis_PP1_correctedLocations202401_annotated.RData",
    json = "tests_dev/test_data/scalefactors_json.json",
    cache = TRUE, cache.path = "tests_dev/test_data/visium_pp1_chat.rds",
    distance.only = FALSE) {
  if (isTRUE(cache) && file.exists(cache.path)) {
    return(readRDS(cache.path))
  }

  suppressPackageStartupMessages(library(Seurat))
  load(rdata)  # provides w10x (Seurat) and the scalefactors json path metadata
  json.path <- if (file.exists(json)) json else {
    cand <- list.files("tests_dev/test_data", pattern = "scalefactors_json[.]json$",
                       full.names = TRUE, recursive = TRUE)
    cand[[1]]
  }

  # cell-type proportions from the deconvolution assay
  cell.type.decomposition <- w10x@assays[["predictions"]]@data |> head(-1)
  CellTypeDict <- c("Keratinocyte", "FIB", "SM", "MYL", "TC", "ENDO", "EG")
  names(CellTypeDict) <- c("Keratinocyte", "Fibroblast", "Smooth muscle",
                           "Myeloid cell", "T cell", "Endothelial", "Eccrine gland")
  rownames(cell.type.decomposition) <- CellTypeDict[rownames(cell.type.decomposition)]
  predictions.immune <- colSums(cell.type.decomposition[c("TC", "MYL"), ])
  cell.type.decomposition <- cell.type.decomposition |>
    Matrix::t() |> as.data.frame()
  cell.type.decomposition$Immune <- predictions.immune
  cell.type.decomposition <- cell.type.decomposition[, !colnames(cell.type.decomposition) %in% c("TC", "MYL")]
  cell.type.decomposition[cell.type.decomposition < 0.1] <- 0

  # spot labels: merge immune cells, rename KC
  cell.labels <- w10x@meta.data[, "final.clusters", drop = FALSE]
  cell.labels$final.clusters <- as.character(cell.labels$final.clusters)
  cell.labels$final.clusters[cell.labels$final.clusters %in% c("TC", "MYL")] <- "Immune"
  cell.labels$final.clusters[cell.labels$final.clusters == "KC"] <- "Keratinocyte"
  w10x@meta.data$final.ident <- as.factor(cell.labels$final.clusters)

  # coordinates: flip axes into the plot frame used by the original script
  spatial.locs <- Seurat::GetTissueCoordinates(w10x, scale = NULL,
                                               cols = c("imagerow", "imagecol"))
  spatial.locs[, 2] <- max(spatial.locs[, 2]) - spatial.locs[, 2]
  temp <- spatial.locs
  spatial.locs[, 1] <- temp[, 2]
  spatial.locs[, 2] <- temp[, 1]
  spatial.locs[, 1] <- max(spatial.locs[, 1]) - spatial.locs[, 1]
  colnames(spatial.locs) <- c("x", "y")

  scalefactors <- jsonlite::fromJSON(json.path)
  spot.size <- 65  # theoretical Visium spot diameter (um)
  conversion.factor <- spot.size / scalefactors$spot_diameter_fullres
  spatial.factors <- list(ratio = conversion.factor, tol = spot.size / 2)

  DefaultAssay(w10x) <- "Spatial"
  w10x <- Seurat::NormalizeData(w10x, verbose = FALSE)
  # group on the merged labels: final.ident's levels {FIB, SM, ENDO, EG,
  # Immune, Keratinocyte} are exactly the columns of the deconvolution table,
  # which the strict dimnames check in computeAvgCommunProb_Visium requires.
  chat <- createSpatialCellChat(
    object = w10x, group.by = "final.ident", assay = "Spatial",
    datatype = "spatial", coordinates = spatial.locs,
    spatial.factors = spatial.factors
  )
  # package data is not attached when sources are loaded directly
  if (!exists("CellChatDB.human")) load("data/CellChatDB.human.rda")
  chat@DB <- subsetDB(CellChatDB.human,
                      search = c("Secreted Signaling", "ECM-Receptor",
                                 "Cell-Cell Contact"), non_protein = FALSE)
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  # The refactored guard requires d * scale.distance >= 1; this sample's minimum
# spot-to-spot distance is ~99.7 um, so the original default 0.01 lands at 0.997.
# 0.011 is the smallest round value above the required 0.01003.
  chat <- computeCommunProb(chat, n = 1, scale.distance = 0.011,
                            nthreads = 4L, verbose = FALSE)
  # deconvolution path: group-level communication from spot-level probabilities,
  # as in tests_dev/test_data/visium_PP1.R
  # the deconvolution table is already cells x types; align its cell axis with
  # idents' names and its type axis with the group levels
  decon <- cell.type.decomposition[rownames(meta(chat)),
                                   levels(idents(chat)), drop = FALSE]
  chat <- computeAvgCommunProb_Visium(chat, cell.type.decomposition = as.matrix(decon),
                                      nboot = 20, do.permutation = FALSE)
  chat <- filterCommunication(chat, min.cells = 10)
  if (isTRUE(distance.only)) return(chat)
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
