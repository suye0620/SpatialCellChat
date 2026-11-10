# Regression checks for the real merged-object producer and snapshot contract.
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
  library(methods); library(Matrix); library(dplyr); library(ggplot2)
  library(cowplot); library(igraph); library(circlize); library(ComplexHeatmap)
  library(ggalluvial); library(Rcpp); library(colorspace)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/merging.R")
source("R/visualization.R")
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")

check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

one <- make_simulated_chat()
two <- make_simulated_chat()
two@misc$.param$sample <- "sample2"
lr.drop <- which(!two@LR$LRsig$pathway_name %in% two@netP$pathways)[1L]
two@LR$LRsig <- two@LR$LRsig[-lr.drop, , drop = FALSE]
two@net$group$count[1, 1] <- two@net$group$count[1, 1] + 1
two@net$group$weight[1, 1] <- two@net$group$weight[1, 1] + 1
two@netP$centr[3, 1, 1] <- two@netP$centr[3, 1, 1] + 1

merged <- mergeSpatialCellChat(list(sample1 = one, sample2 = two))
check("producer returns a merged SpatialCellChat", methods::is(merged, "SpatialCellChat") &&
        identical(merged@misc$.mode, "merged"))
check("merged object validates", isTRUE(validateSpatialCellChat(merged)))
check("dataset metadata and prefixed cells align",
      identical(levels(merged@meta$datasets), c("sample1", "sample2")) &&
        identical(rownames(merged@meta), colnames(merged@assay$norm)) &&
        all(startsWith(rownames(merged@meta)[seq_len(nrow(one@meta))], "sample1_")) &&
        all(startsWith(rownames(merged@meta)[nrow(one@meta) + seq_len(nrow(two@meta))], "sample2_")))
check("snapshot names and cell-level sparse axes align",
      identical(names(merged@misc$.datasets), c("sample1", "sample2")) &&
        identical(dimnames(merged@misc$.datasets[[1]]$net$cell$prob)[[1]],
                  paste0("sample1_", colnames(one@assay$norm))) &&
        identical(dimnames(merged@misc$.datasets[[2]]$net$cell$prob)[[1]],
                  paste0("sample2_", colnames(two@assay$norm))))
check("dataset-specific LRsig tables are retained and unioned",
      nrow(merged@misc$.datasets[[1]]$LR$LRsig) == nrow(one@LR$LRsig) &&
        nrow(merged@misc$.datasets[[2]]$LR$LRsig) == nrow(two@LR$LRsig) &&
        nrow(merged@LR$LRsig) == nrow(one@LR$LRsig))
check("group-level snapshot axes retain group names",
      identical(dimnames(merged@misc$.datasets[[1]]$net$group$prob)[[1]],
                levels(one@idents)))
check("rankNet comparison consumes the real producer",
      inherits(rankNet(merged, slot.name = "netP", mode = "comparison",
                       comparison = c(1, 2), measure = "count", return.data = TRUE)$gg.obj,
               "ggplot"))

similarity.merged <- computeNetSimilarityPairwise(merged, slot.name = "netP", comparison = c(1, 2), k = 2)
check("pairwise similarity consumes snapshot probabilities",
      is.matrix(similarity.merged@netP$similarity$functional$matrix[["1-2"]]) &&
        nrow(similarity.merged@netP$similarity$functional$matrix[["1-2"]]) == 4L)

similarity.merged <- netEmbedding(similarity.merged, slot.name = "netP", type = "functional", comparison = c(1, 2), umap.method = "uwot", n_neighbors = 2)
similarity.merged <- netClustering(similarity.merged, slot.name = "netP", type = "functional", comparison = c(1, 2), k = 2, do.plot = FALSE)
check("pairwise embedding reads merged snapshots",
      inherits(netVisual_embeddingPairwise(similarity.merged, slot.name = "netP", comparison = c(1, 2), do.label = FALSE, show.legend = FALSE), "ggplot"))

check("signaling role comparison reads snapshot centrality",
      inherits(netAnalysis_diff_signalingRole_scatter(merged, comparison = c(1, 2), do.label = FALSE, show.legend = FALSE), "ggplot"))
check("signaling change comparison reads snapshot centrality",
      inherits(netAnalysis_signalingChanges_scatter(merged, idents.use = "A", comparison = c(1, 2), do.label = FALSE, show.legend = FALSE), "ggplot"))

communication <- subsetCommunication(merged, slot.name = "netP")
check("merged subsetCommunication returns named snapshots",
      is.list(communication) && identical(names(communication), c("sample1", "sample2")) &&
        all(vapply(communication, nrow, integer(1)) > 0L))
check("merged barplot reads snapshot matrices",
      inherits(netVisual_barplot(merged, comparison = c(1, 2), measure = "count"), "ggplot"))
check("merged diff interaction reads snapshot matrices",
      inherits(netVisual_diffInteraction(merged, comparison = c(1, 2), measure = "count"), "recordedplot"))
check("merged bubble reads snapshot communication",
      inherits(netVisual_bubble(merged, comparison = c(1, 2), signaling = "P1", return.data = TRUE)$gg.obj,
               "ggplot"))
netVisual_heatmap(merged, comparison = c(1, 2), measure = "count")
check("merged heatmap reads snapshot matrices", TRUE)
bad <- two
bad@idents <- factor(ifelse(seq_along(bad@idents) == 1L, "different", as.character(bad@idents)))
names(bad@idents) <- colnames(bad@assay$norm)
check("group mismatch is rejected", inherits(try(mergeSpatialCellChat(
  list(sample1 = one, sample2 = bad)), silent = TRUE), "try-error"))

subset.one <- subsetSpatialCellChat(one, cells.use = colnames(one@assay$norm)[1:5])
check("single subset keeps aligned assay and metadata",
      ncol(subset.one@assay$norm) == 5L &&
        identical(colnames(subset.one@assay$norm), rownames(subset.one@meta)) &&
        isTRUE(validateSpatialCellChat(subset.one)))
subset.merged <- subsetSpatialCellChat(merged,
  cells.use = c(rownames(merged@meta)[1:5], rownames(merged@meta)[21:25]))
check("merged subset retains both datasets and clears stale networks",
      identical(levels(subset.merged@meta$datasets), c("sample1", "sample2")) &&
        length(subset.merged@misc$.datasets) == 2L &&
        !length(subset.merged@net$group) &&
        isTRUE(validateSpatialCellChat(subset.merged)))
lifted <- liftSpatialCellChat(merged,
  group.new = c(A = "AB", B = "AB", C = "C", D = "D"))
check("lift relabels group axes and invalidates networks",
      identical(levels(lifted@idents), c("AB", "C", "D")) &&
        !length(lifted@net$group) && isTRUE(validateSpatialCellChat(lifted)))
check("update validates final-schema objects",
      isTRUE(validateSpatialCellChat(updateSpatialCellChat(one))))

cat("All merged producer checks passed.\n")
