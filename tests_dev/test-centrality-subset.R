# Wave 3a: centrality / subsetCommunication migration
# Contract: final 11-slot schema, native R centrality reshape, and explicit netP semantics.
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
  library(methods); library(Matrix); library(cli); library(dplyr)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/modeling.R")
source("R/database.R")
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))

## ---- fixture: K=3, four LR layers, P1/P1/P2/P3 ----
K <- 3L
cells <- paste0("c", seq_len(12L))
groups <- c("A", "B", "C")
group <- factor(rep(groups, c(5L, 4L, 3L)), levels = groups)
names(group) <- cells
expr <- matrix(seq_len(5L * length(cells)) / 10,
               nrow = 5L, ncol = length(cells),
               dimnames = list(paste0("gene", seq_len(5L)), cells))
meta <- data.frame(label = group, row.names = cells)
chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "RNA")
chat@LR$LRsig <- data.frame(
  interaction_name = c("LR1", "LR2", "LR3", "LR4"),
  interaction_name_2 = c("LR1", "LR2", "LR3", "LR4"),
  ligand = c("gene1", "gene3", "gene1", "gene2"),
  receptor = c("gene2", "gene4", "gene4", "gene5"),
  pathway_name = c("P1", "P1", "P2", "P3"),
  annotation = c("a1", "a2", "a3", "a4"),
  row.names = c("LR1", "LR2", "LR3", "LR4"),
  stringsAsFactors = FALSE
)
# Deliberately no `evidence` column: BUG-A regression.

gd <- function(x) as(matrix(x, K, K), "dgCMatrix")
g1 <- gd(c(0.5, 0.2, 0.1, 0.3, 0.4, 0.2, 0.1, 0.1, 0.9))
g2 <- gd(c(0.7, 0.1, 0.0, 0.2, 0.3, 0.1, 0.4, 0.2, 0.1))
g3 <- gd(c(0.9, 0.0, 0.2, 0.1, 0.5, 0.3, 0.0, 0.6, 0.2))
g4 <- gd(rep(0, K * K))
pv1 <- gd(rep(0.01, K * K)); pv1[2, 3] <- 0.5
pv2 <- gd(rep(0.01, K * K)); pv2[1, 2] <- 0.05
pv3 <- gd(rep(0.01, K * K)); pv3[1, 1] <- 0.05
pv4 <- gd(rep(0.01, K * K))
layers <- list(LR1 = g1, LR2 = g2, LR3 = g3, LR4 = g4)
pvals <- list(LR1 = pv1, LR2 = pv2, LR3 = pv3, LR4 = pv4)
make_array <- function(x, dn3 = names(x), dn1 = groups, dn2 = groups) {
  out <- SparseChatArray(x)
  dimnames(out) <- list(dn1, dn2, dn3)
  out
}
chat@net$group$prob <- make_array(layers)
chat@net$group$pval <- make_array(pvals)
set.seed(20261010)
cell.layers <- lapply(seq_along(layers), function(i) {
  m <- matrix(runif(length(cells)^2) * (runif(length(cells)^2) > 0.55),
              nrow = length(cells), ncol = length(cells))
  as(m, "dgCMatrix")
})
names(cell.layers) <- names(layers)
chat@net$cell$prob <- make_array(cell.layers, dn1 = cells, dn2 = cells)

raw.layers <- lapply(layers, as.matrix)
gated.layers <- Map(function(pm, pv) {
  pm <- as.matrix(pm); pv <- as.matrix(pv)
  pv[pm == 0] <- 1
  pm[pv >= 0.05] <- 0
  pm
}, layers, pvals)

## ---- 1. group centrality, gating, shape, and write-back ----
central <- netAnalysis_computeCentrality(chat, do.group = TRUE, thresh = 0.05)
centr <- central@net[["centr"]]
check("1a group centrality shape and dimnames",
      identical(dim(centr), c(5L, K, 4L)) &&
      identical(dimnames(centr), list(
        c("outdeg_unweighted", "indeg_unweighted", "outdeg", "indeg", "page_rank"),
        groups, names(layers))))
check("1b unweighted outdegree equals gated row support", all(vapply(seq_along(gated.layers), function(i)
  identical(as.numeric(centr["outdeg_unweighted", , i]),
            as.numeric(Matrix::rowSums(gated.layers[[i]] > 0))), logical(1))))
check("1c unweighted indegree equals gated column support", all(vapply(seq_along(gated.layers), function(i)
  identical(as.numeric(centr["indeg_unweighted", , i]),
            as.numeric(Matrix::colSums(gated.layers[[i]] > 0))), logical(1))))
check("1d weighted outdegree and indegree equal gated sums", all(vapply(seq_along(gated.layers), function(i) {
  actual.out <- as.numeric(centr["outdeg", , i])
  actual.in <- as.numeric(centr["indeg", , i])
  ref.out <- as.numeric(rowSums(gated.layers[[i]]))
  ref.in <- as.numeric(colSums(gated.layers[[i]]))
  if (!isTRUE(all.equal(actual.out, ref.out, tolerance = 1e-12)) ||
      !isTRUE(all.equal(actual.in, ref.in, tolerance = 1e-12)))
    cat("1d layer", i, "actual out", paste(actual.out, collapse = ","), "ref out", paste(ref.out, collapse = ","),
        "actual in", paste(actual.in, collapse = ","), "ref in", paste(ref.in, collapse = ","), "\n")
  isTRUE(all.equal(actual.out, ref.out, tolerance = 1e-12)) &&
    isTRUE(all.equal(actual.in, ref.in, tolerance = 1e-12))
}, logical(1))))
check("1e page_rank equals direct igraph result", all(vapply(seq_along(gated.layers), function(i) {
  ref <- igraph::page_rank(igraph::graph_from_adjacency_matrix(gated.layers[[i]], "directed", weighted = TRUE))$vector
  isTRUE(all.equal(as.numeric(centr["page_rank", , i]), as.numeric(ref), tolerance = 1e-12))
}, logical(1))))
check("1f whole-slot write-back and validObject", validObject(central))
check("1g operation log records centrality parameters", {
  entry <- central@misc$.log[[length(central@misc$.log)]]
  identical(entry[["function"]], "netAnalysis_computeCentrality") &&
    identical(entry$params$slot.name, "net") && isTRUE(entry$params$do.group) &&
    identical(entry$params$degree.only, TRUE) && identical(entry$params$nrun, 4L)
})

## ---- 2. native R reshape equivalence ----
manual.centrality.layers <- lapply(gated.layers, computeCentralityLocal, degree.only = TRUE)
manual.reshape <- array(unlist(lapply(manual.centrality.layers, as.vector)), c(5L, K, 4L))
dimnames(manual.reshape) <- dimnames(centr)
check("2a native R reshape equals manual layer assembly", identical(unname(centr), unname(manual.reshape)))
central.null <- netAnalysis_computeCentrality(object = NULL, net = layers[[1]], degree.only = TRUE)
check("2b object=NULL 2D net returns one named layer", identical(dim(central.null), c(5L, K, 1L)) && identical(dimnames(central.null)[[3]], "net_temp"))
net3 <- array(c(unlist(lapply(layers, as.matrix))), c(K, K, 4L), dimnames = list(groups, groups, names(layers)))
central.3d <- netAnalysis_computeCentrality(object = NULL, net = net3, degree.only = TRUE)
check("2c named 3D net preserves layers", identical(dim(central.3d), c(5L, K, 4L)) && identical(dimnames(central.3d)[[3]], names(layers)))
central.filtered <- netAnalysis_computeCentrality(object = NULL, net = net3, signaling.name = c("LR3", "LR1"))
check("2d signaling.name filters by layer name", identical(dimnames(central.filtered)[[3]], c("LR3", "LR1")))
check("2e signaling.name invalid error", grepl("Please check the input `signaling.name`", err(netAnalysis_computeCentrality(object = NULL, net = net3, signaling.name = "NO")), fixed = TRUE))

## ---- 3. degree.only=FALSE ----
central.full <- netAnalysis_computeCentrality(object = NULL, net = net3, degree.only = FALSE)
check("3a full centrality measure names", identical(dimnames(central.full)[[1]], c(
  "outdeg_unweighted", "indeg_unweighted", "outdeg", "indeg", "hub", "authority", "eigen",
  "page_rank", "betweenness", "flowbet", "info")))
check("3b flowbet and info are numeric", is.numeric(central.full["flowbet", , ]) && is.numeric(central.full["info", , ]))
central.na <- netAnalysis_computeCentrality(object = NULL, net = layers[[1]], degree.only = NA)
check("3c degree.only NA normalizes to full measures", nrow(central.na) == 11L)

## ---- 4. cell-level centrality ----
central.cell <- netAnalysis_computeCentrality(chat, do.group = FALSE)
cellcentr <- central.cell@net[["centr.cell"]]
check("4a cell centrality reads net$cell$prob", identical(dim(cellcentr), c(5L, length(cells), 4L)) &&
      identical(dimnames(cellcentr)[[2]], cells) && identical(dimnames(cellcentr)[[3]], names(layers)))
check("4b cell centrality validObject", validObject(central.cell))

## ---- 5. netP centrality is unconditionally ungated ----
path.layers <- list(P1 = as(gated.layers[[1]] + gated.layers[[2]], "dgCMatrix"),
                    P2 = g3, P3 = g4)
chatP <- chat
chatP@netP$group$prob <- make_array(path.layers)
chatP@netP$group$pval <- make_array(list(P1 = pv1, P2 = pv3, P3 = pv4))
centralP <- netAnalysis_computeCentrality(chatP, slot.name = "netP", do.group = TRUE, thresh = 0.05)
check("5a netP ignores an injected pval array", all(vapply(seq_along(path.layers), function(i) {
  actual <- as.numeric(centralP@netP$centr["outdeg", , i])
  ref <- as.numeric(rowSums(as.matrix(path.layers[[i]])))
  isTRUE(all.equal(actual, ref, tolerance = 1e-12))
}, logical(1))))
check("5b netP centrality layer names", identical(dimnames(centralP@netP$centr)[[3]], names(path.layers)))

## ---- 6. subsetCommunication net: manual row construction and filters ----
df <- subsetCommunication(chat, thresh = 0.05)
expected.cols <- c("source", "target", "ligand", "receptor", "prob", "pval", "interaction_name",
                   "interaction_name_2", "pathway_name", "annotation", "evidence")
check("6a subset net returns expected columns with missing evidence filled", identical(colnames(df), expected.cols) &&
      is.character(df$evidence) && all(is.na(df$evidence)))
manual.rows <- do.call(rbind, Filter(Negate(is.null), lapply(seq_along(gated.layers), function(i) {
  idx <- which(gated.layers[[i]] > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(NULL)
  data.frame(source = groups[idx[, 1]], target = groups[idx[, 2]], interaction_name = names(layers)[i],
             prob = gated.layers[[i]][idx], pval = as.matrix(pvals[[i]])[idx], stringsAsFactors = FALSE)
})))
actual.key <- paste(df$source, df$target, df$interaction_name, sep = "|")
manual.key <- paste(manual.rows$source, manual.rows$target, manual.rows$interaction_name, sep = "|")
if (!identical(actual.key, manual.key)) {
  cat("actual keys:", paste(head(actual.key, 12), collapse = ";"), "\nmanual keys:", paste(head(manual.key, 12), collapse = ";"), "\n")
}
check("6b subset net support and row ordering equal manual construction", identical(actual.key, manual.key))
check("6c subset net probability and pvalue equal manual construction",
      identical(as.numeric(df$prob), as.numeric(manual.rows$prob)) &&
      identical(as.numeric(df$pval), as.numeric(manual.rows$pval)))
check("6d pval >= thresh entries are absent", !any(df$interaction_name == "LR2" & df$source == "A" & df$target == "B"))
check("6e signaling filter uses searchPair", {
  d <- subsetCommunication(chat, signaling = "P1")
  nrow(d) > 0L && all(d$pathway_name == "P1")
})
check("6f pairLR interaction_name filter", {
  d <- subsetCommunication(chat, pairLR.use = data.frame(interaction_name = "LR1"))
  nrow(d) > 0L && all(d$interaction_name == "LR1")
})
check("6g pairLR pathway_name switches to netP", {
  d <- subsetCommunication(chatP, pairLR.use = data.frame(pathway_name = "P1"))
  nrow(d) > 0L && all(d$pathway_name == "P1") && all(is.na(d$pval))
})
check("6h sources/targets accept numeric and character selectors", {
  d1 <- subsetCommunication(chat, sources.use = 1, targets.use = "B")
  d2 <- subsetCommunication(chat, sources.use = "A", targets.use = 2)
  nrow(d1) > 0L && nrow(d2) > 0L && all(d1$source == "A") && all(d1$target == "B") &&
    all(d2$source == "A") && all(d2$target == "B")
})
check("6i empty support has the established stop", grepl("No significant signaling interactions are inferred based on the input!", err(
  subsetCommunication(chat, pairLR.use = data.frame(interaction_name = "NO_SUCH_INTERACTION"))), fixed = TRUE))
check("6j DEG guard has the established stop", grepl("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'", err(
  subsetCommunication(chat, ligand.pvalues = 0.05)), fixed = TRUE))

## ---- 7. subsetCommunication netP: no pval, pathway join, aggregation ----
dfP <- subsetCommunication(chatP, slot.name = "netP")
check("7a netP returns pathway schema and explicit NA pval", identical(colnames(dfP), c("source", "target", "pathway_name", "prob", "pval")) &&
      nrow(dfP) > 0L && all(is.na(dfP$pval)))
path.expected <- do.call(rbind, Filter(Negate(is.null), lapply(seq_along(path.layers), function(i) {
  idx <- which(as.matrix(path.layers[[i]]) > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) return(NULL)
  data.frame(source = groups[idx[, 1]], target = groups[idx[, 2]], pathway_name = names(path.layers)[i],
             prob = as.matrix(path.layers[[i]])[idx], stringsAsFactors = FALSE)
})))
path.expected <- path.expected[order(path.expected$source, path.expected$target, path.expected$pathway_name), ]
path.actual <- dfP[order(dfP$source, dfP$target, dfP$pathway_name), ]
check("7b netP pathway join and probability equal manual layer rows", {
  identical(paste(path.actual$source, path.actual$target, path.actual$pathway_name, sep = "|"),
            paste(path.expected$source, path.expected$target, path.expected$pathway_name, sep = "|")) &&
    isTRUE(all.equal(as.numeric(path.actual$prob), as.numeric(path.expected$prob), tolerance = 0))
})
check("7c netP multi-LR pathway uses pathway_name join", all(dfP$pathway_name %in% names(path.layers)))

## ---- 8. explicit regressions and determinism ----
check("8a new schema mode NULL path does not crash (BUG-B)", nrow(subsetCommunication(chat)) == nrow(df))
check("8b numeric source overflow has explicit error (BUG-E)", grepl("beyond the number of cell groups", err(
  subsetCommunication(chat, sources.use = K + 1L)), fixed = TRUE))
check("8c numeric target overflow has explicit error (BUG-E)", grepl("beyond the number of cell groups", err(
  subsetCommunication(chat, targets.use = K + 1L)), fixed = TRUE))
check("8d invalid centrality slot is rejected by match.arg", grepl("should be one of", err(
  netAnalysis_computeCentrality(chat, slot.name = "foo")), fixed = TRUE))
legacy <- structure(list(), class = "CellChat")
check("8e legacy/non-SpatialCellChat input fails with clear assert", identical(
  err(netAnalysis_computeCentrality(legacy)), "object must be a SpatialCellChat"))
central.repeat <- netAnalysis_computeCentrality(object = NULL, net = net3, degree.only = TRUE)
check("8f centrality is deterministic", identical(unname(central.3d), unname(central.repeat)))
check("8g subset result remains valid after centrality write-back", validObject(central))

cat("All checks passed.\n")
