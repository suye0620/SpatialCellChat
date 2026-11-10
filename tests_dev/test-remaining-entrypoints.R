# Regression checks for omitted comparison/analysis entrances after the 11-slot migration.
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
  library(methods); library(Matrix); library(ggplot2); library(dplyr)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/merging.R")
source("R/analysis.R")
source("R/visualization.R")

check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

cells <- paste0("c", 1:4)
groups <- c("A", "B")
genes <- c("L1", "R1", "L2", "R2")
expr <- Matrix::Matrix(matrix(1, nrow = length(genes), ncol = length(cells),
                              dimnames = list(genes, cells)), sparse = TRUE)
meta <- data.frame(label = factor(c("A", "A", "B", "B"), levels = groups),
                   row.names = cells)
lr <- data.frame(
  interaction_name = c("LR1", "LR2"),
  ligand = c("L1", "L2"), receptor = c("R1", "R2"),
  pathway_name = c("P1", "P2"), annotation = c("Secreted Signaling", "Cell-Cell Contact"),
  stringsAsFactors = FALSE, row.names = c("LR1", "LR2")
)
make_matrix <- function(values) {
  as(Matrix::Matrix(matrix(values, nrow = 2L, byrow = TRUE,
                            dimnames = list(groups, groups)), sparse = TRUE),
     "dgCMatrix")
}
make_array <- function(layers) {
  SparseChatArray(layers, dimnames = list(groups, groups, names(layers)))
}
make_centr <- function(outgoing, incoming) {
  result <- array(0, dim = c(2L, 2L, 2L),
                  dimnames = list(c("outdeg", "indeg"), groups, c("P1", "P2")))
  result["outdeg", , ] <- matrix(outgoing, nrow = 2L)
  result["indeg", , ] <- matrix(incoming, nrow = 2L)
  result
}
make_chat <- function(outgoing, incoming, sample) {
  chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")
  chat@LR <- list(LRsig = lr)
  chat@DB <- list(
    interaction = lr,
    complex = data.frame(row.names = character()),
    geneInfo = data.frame(Symbol = genes, stringsAsFactors = FALSE)
  )
  lr.prob <- make_array(list(
    LR1 = make_matrix(c(1, 2, 0, 1)),
    LR2 = make_matrix(c(0, 1, 3, 0))
  ))
  lr.pval <- make_array(list(
    LR1 = make_matrix(c(0.01, 0.20, 0.03, 0.04)),
    LR2 = make_matrix(c(0.02, 0.01, 0.40, 0.02))
  ))
  path.prob <- make_array(list(
    P1 = make_matrix(c(1, 2, 0, 1)),
    P2 = make_matrix(c(0, 1, 3, 0))
  ))
  chat@net$group <- list(
    prob = lr.prob, pval = lr.pval,
    count = make_matrix(c(2, 3, 1, 2)),
    weight = make_matrix(c(1, 2, 3, 4))
  )
  chat@netP$group <- list(prob = path.prob)
  chat@netP$pathways <- c("P1", "P2")
  chat@netP$centr <- make_centr(outgoing, incoming)
  chat@misc$.param$sample <- sample
  methods::validObject(chat)
  chat
}

one <- make_chat(c(1, 2, 2, 5), c(10, 20, 20, 10), "one")
two <- make_chat(c(3, 5, 4, 8), c(30, 40, 18, 18), "two")

pairwise <- rankNetPairwise(one)
check("rankNetPairwise writes current group-level pairwise rankings",
      is.list(pairwise@net$pairwiseRank) &&
        identical(names(pairwise@net$pairwiseRank), groups) &&
        identical(pairwise@net$pairwiseRank$A$B$interaction_name, c("LR2", "LR1")))
check("rankNetPairwise preserves p-value ordering and sparse probabilities",
      identical(pairwise@net$pairwiseRank$A$B$pval, c(0.01, 0.20)) &&
        identical(pairwise@net$pairwiseRank$A$B$prob, c(1, 2)))
missing.pval <- one
missing.pval@net$group$pval <- NULL
check("rankNetPairwise rejects missing p-values instead of inventing evidence",
      grepl("requires net\\$group\\$pval", tryCatch({
        rankNetPairwise(missing.pval); ""
      }, error = conditionMessage)))
enriched <- identifyEnrichedInteractions(pairwise, from = "A", to = "B", thresh = 0.05)
check("identifyEnrichedInteractions consumes the migrated pairwise result",
      is.data.frame(enriched) && all(enriched$interaction_name %in% c("LR1", "LR2")))
check("invalid group names cannot silently re-pair endpoints",
      grepl("Invalid cell group", tryCatch({
        identifyEnrichedInteractions(pairwise, c("A", "missing"), c("missing", "B")); ""
      }, error = conditionMessage), fixed = TRUE))
check("empty and out-of-range group selections are rejected",
      nzchar(tryCatch({ identifyEnrichedInteractions(pairwise, integer(), integer()); "" },
                      error = conditionMessage)) &&
        nzchar(tryCatch({ identifyEnrichedInteractions(pairwise, 1, 3); "" },
                        error = conditionMessage)))
none.enriched <- identifyEnrichedInteractions(pairwise, "A", "B", thresh = 0)
check("empty threshold selection returns an empty result",
      is.data.frame(none.enriched) && nrow(none.enriched) == 0L)

merged <- mergeSpatialCellChat(list(one = one, two = two), show.plot = FALSE)
check("compareInteractions reads snapshot group matrices",
      identical(compareInteractions(merged, measure = "count")$data$count,
                c(8, 8)))
check("compareInteractions supports weight summaries and groups",
      identical(compareInteractions(merged, measure = "weight", group = c("x", "y"))$data$count,
                c(10, 10)) &&
        identical(as.character(compareInteractions(merged, group = c("x", "y"))$data$group), c("x", "y")))

pairwise.merged <- rankNetPairwise(merged)
check("rankNetPairwise stores dataset-local results in merged snapshots",
      all(vapply(pairwise.merged@misc$.datasets,
                 function(x) is.list(x$net$pairwiseRank), logical(1))))
pairwise.merged@misc$.datasets$two$net$pairwiseRank$A$B["LR1", "pval"] <- 0.01
pairwise.merged@misc$.datasets$two$net$pairwiseRank$A$B["LR2", "pval"] <- 0.20
merged.enriched <- identifyEnrichedInteractions(pairwise.merged, "A", "B", pair.only = FALSE)
check("merged enrichment consumes the matching dataset ranking",
      identical(names(merged.enriched), c("one", "two")) &&
        identical(merged.enriched$one$interaction_name, "LR2") &&
        identical(merged.enriched$two$interaction_name, "LR1"))
unranked <- pairwise.merged
unranked@misc$.datasets$two$net$pairwiseRank <- NULL
check("merged enrichment identifies a missing dataset ranking",
      grepl("dataset 'two'", tryCatch({
        identifyEnrichedInteractions(unranked, "A", "B"); ""
      }, error = conditionMessage), fixed = TRUE))

Y <- matrix(c(0, 0, 1, 0, 0, 0, 2, 0), nrow = 4L, byrow = TRUE,
            dimnames = list(c("P1--one", "P2--one", "P1--two", "P2--two"),
                            c("UMAP1", "UMAP2")))
merged@netP$similarity <- list(functional = list(dr = list(`1-2` = Y)))
ranked <- rankSimilarity(merged, type = "functional", comparison1 = c(1, 2))
check("rankSimilarity reads the merged embedding and aligns dataset suffixes",
      inherits(ranked, "ggplot") && identical(as.character(ranked$data$name), c("P1", "P2")) &&
        identical(as.numeric(ranked$data$dist), c(0, 1)))
permuted.embedding <- merged
permuted.embedding@netP$similarity$functional$dr$`1-2` <- Y[c(4, 2, 3, 1), , drop = FALSE]
check("rankSimilarity aligns permuted embedding rows by name",
      identical(rankSimilarity(permuted.embedding, type = "functional", comparison1 = c(1, 2))$data,
                ranked$data))
missing.embedding <- merged
missing.embedding@netP$similarity$functional$dr$`1-2` <- Y[1:2, , drop = FALSE]
check("rankSimilarity rejects a dataset absent from the embedding",
      grepl("does not contain layers", tryCatch({
        rankSimilarity(missing.embedding, type = "functional", comparison1 = c(1, 2)); ""
      }, error = conditionMessage), fixed = TRUE))

role.diff <- netAnalysis_diff_signalingRole_scatter(
  merged, comparison = c(1, 2), do.label = FALSE, show.legend = FALSE
)
check("differential signaling-role scatter reads snapshot centrality",
      inherits(role.diff, "ggplot") && nrow(role.diff$data) == length(groups))

changes <- netAnalysis_signalingChanges_scatter(
  merged, idents.use = "A", comparison = c(1, 2), do.label = FALSE, show.legend = FALSE
)
check("signaling-change scatter retains the historical row-sum filter",
      inherits(changes, "ggplot") && identical(rownames(changes$data), "P1") &&
        identical(as.numeric(changes$data$outgoing), 2) &&
        identical(as.numeric(changes$data$incoming), 20))

one.signal <- netAnalysis_diff_signalingRole_scatter(
  merged, comparison = c(1, 2), signaling = "P1", do.label = FALSE
)
check("differential scatter handles a single signaling layer",
      inherits(one.signal, "ggplot") && nrow(one.signal$data) == length(groups))
permuted <- merged
permuted@misc$.datasets$two$netP$centr <-
  permuted@misc$.datasets$two$netP$centr[, rev(groups), , drop = FALSE]
reordered.role <- netAnalysis_diff_signalingRole_scatter(
  permuted, comparison = c(1, 2), do.label = FALSE, show.legend = FALSE
)
check("differential scatter aligns centrality by group name",
      identical(role.diff$data, reordered.role$data))

coordinates <- cbind(x = seq_along(cells), y = seq_along(cells))
rownames(coordinates) <- cells
spatial.one <- createSpatialCellChat(expr, meta = meta, group.by = "label",
  datatype = "spatial", coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = 1))
spatial.one@LR <- one@LR
spatial.one@DB <- one@DB
spatial.merged <- mergeSpatialCellChat(list(one = spatial.one, two = spatial.one),
                                       show.plot = FALSE)
check("spatialVisual_scoring rejects merged dataset-local coordinates explicitly",
      grepl("dataset-local coordinates", tryCatch({
        spatialVisual_scoring(spatial.merged, signaling = "P1")
        ""
      }, error = conditionMessage)))

cat("All remaining comparison/analysis entry-point checks passed.\n")

