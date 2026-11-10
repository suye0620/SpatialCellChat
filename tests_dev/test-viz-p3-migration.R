# P3 group probabilities: preserve axes and exercise actual producers/renderers.
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
  library(methods); library(Matrix); library(cli); library(dplyr); library(igraph)
  library(future.apply); library(ggplot2); library(colorspace); library(Rcpp)
  library(circlize); library(ComplexHeatmap); library(spatstat.sparse); library(ggalluvial)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")
check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
chat <- make_simulated_chat()
prob <- chat@netP$group$prob
dense <- .sc_group_prob_array(chat, "netP")
check("group probability array preserves dimension names and layer values",
      identical(dimnames(dense), dimnames(prob)) &&
        all(vapply(seq_len(dim(prob)[3L]), function(i)
          isTRUE(all.equal(unname(dense[, , i]), unname(as.matrix(prob[[i]])))), logical(1))))
chat@meta$viz_test_group <- chat@idents
composition <- netVisual_TopicComposition(chat, x.group = "viz_test_group", y.group = "viz_test_group",
                                           type = "dot", remove.isolate = FALSE)
check("topic composition reads SpatialCellChat metadata", inherits(composition, "ggplot"))
review.dir <- "tests_dev/figures/P3/simulated"
dir.create(review.dir, recursive = TRUE, showWarnings = FALSE)
composition.file <- file.path(review.dir, "netVisual_TopicComposition-metadata.pdf")
ggplot2::ggsave(composition.file, composition, width = 7, height = 5.5)
check("metadata topic composition renders a PDF", file.info(composition.file)$size > 1000L)
scores <- rbind(first = rep(1, length(chat@idents)), second = rep(2, length(chat@idents)),
                isolated = rep(0, length(chat@idents)))
composition <- netVisual_TopicComposition(chat, y.group = scores, type = "dot", remove.isolate = TRUE)
check("topic isolate removal removes rows without dropping cell groups",
      !"isolated" %in% levels(composition$data$x) &&
        all(levels(chat@idents) %in% levels(composition$data$y)))

seurat.counts <- Matrix::Matrix(matrix(1, nrow = 3L, ncol = length(chat@idents),
  dimnames = list(chat@LR$LRsig$interaction_name, names(chat@idents))), sparse = TRUE)
seurat <- SeuratObject::CreateSeuratObject(seurat.counts)
embedding <- matrix(seq_len(ncol(seurat.counts) * 2L), ncol = 2L,
  dimnames = list(colnames(seurat.counts), c("Topic_1", "Topic_2")))
loadings <- matrix(seq_len(nrow(seurat.counts) * 2L), ncol = 2L,
  dimnames = list(rownames(seurat.counts), c("Topic_1", "Topic_2")))
seurat[["Topic_.outgoing"]] <- SeuratObject::CreateDimReducObject(
  embeddings = embedding, loadings = loadings, key = "Topic_", assay = "RNA")
seurat.plot <- netVisual_TopicSignaling(seurat, pattern = "outgoing", feature_cat = "signaling",
                                       db = chat@DB$interaction, topics = 1L)
seurat.file <- file.path(review.dir, "netVisual_TopicSignaling-Seurat.pdf")
ggplot2::ggsave(seurat.file, seurat.plot, width = 7, height = 5.5)
check("Seurat topic signaling renders with an interaction data frame",
      file.exists(seurat.file) && file.info(seurat.file)$size > 1000L)

pattern <- identifyCommunicationPatterns(chat, slot.name = "netP", pattern = "outgoing", k = 2,
                                         heatmap.show = FALSE)
check("NMF pattern producer stores group and signaling results",
      nrow(pattern@netP$pattern$outgoing$pattern$cell) > 0L &&
        nrow(pattern@netP$pattern$outgoing$pattern$signaling) > 0L)
for (fn in c("netAnalysis_river", "netAnalysis_dot")) {
  p <- get(fn)(pattern, slot.name = "netP", pattern = "outgoing")
  output <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(output, p, width = 7, height = 5.5)
  check(paste(fn, "renders a PDF"), file.exists(output) && file.info(output)$size > 1000L)
  unlink(output)
}

similarity <- computeNetSimilarity(chat, slot.name = "netP", type = "functional", k = 2)
similarity <- netEmbedding(similarity, slot.name = "netP", type = "functional",
                           pathway.remove = character(), umap.method = "uwot", n_neighbors = 2)
n.embedded <- nrow(similarity@netP$similarity$functional$dr$single)
if (n.embedded < 2L) stop("Need at least two embedded pathways for kmeans", call. = FALSE)
k <- min(2L, n.embedded - 1L)
cat("Embedding test pathways:", n.embedded, "| kmeans clusters:", k, "\n")
similarity <- netClustering(similarity, slot.name = "netP", type = "functional", k = k, do.plot = FALSE)
for (fn in c("netVisual_embedding", "netVisual_embeddingZoomIn")) {
  p <- get(fn)(similarity, slot.name = "netP", type = "functional", pathway.remove = character())
  output <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(output, p, width = 7, height = 5.5)
  check(paste(fn, "renders a PDF from actual producer chain"),
        file.exists(output) && file.info(output)$size > 1000L)
  unlink(output)
}
cat("All P3 group migration checks passed.\n")
