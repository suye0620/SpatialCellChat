# Run in separate R processes with argument dev-first or upstream-first.
args <- commandArgs(trailingOnly = TRUE)
order <- if (length(args)) args[[1L]] else "dev-first"
suppressPackageStartupMessages({
  if (identical(order, "upstream-first")) {
    library(SpatialCellChat)
    library(SpatialCellChatDev)
  } else {
    library(SpatialCellChatDev)
    library(SpatialCellChat)
  }
})
check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); "" }, error = conditionMessage)
  grepl(pattern, message, fixed = TRUE)
}
update <- SpatialCellChatDev::updateSpatialCellChat
validate <- SpatialCellChatDev:::validateSpatialCellChat
cells <- paste0("c", seq_len(4L))
genes <- paste0("g", seq_len(4L))
groups <- c("A", "B")
lrs <- c("LR1", "LR2")
expr <- matrix(seq_len(16L), 4L, dimnames = list(genes, cells))
meta <- data.frame(group = factor(c("A", "A", "B", "B")), row.names = cells)
make.old <- function(datatype = "RNA") {
  SpatialCellChat::createSpatialCellChat(
    expr, meta = meta, group.by = "group", datatype = datatype,
    coordinates = if (datatype == "spatial") cbind(x_cent = 1:4, y_cent = c(2, 4, 6, 8)) else NULL,
    spatial.factors = if (datatype == "spatial") list(ratio = 0.5, tol = 2) else NULL)
}
rna <- make.old()
rna@var.features <- list(features = c("g1", "g3"), detail = data.frame(x = 1))
rna@data.raw <- rna@data
rna@data.scale <- matrix(seq_len(16L), 4L, dimnames = list(genes, cells))
rna@data.signaling <- rna@data[c("g3", "g1"), , drop = FALSE]
rna@data.project <- rna@data.signaling
rna@DB <- list(interaction = data.frame(ligand = "g1"))
rna@LR <- list(LRsig = data.frame(interaction_name = lrs, pathway_name = c("P1", "P2"),
                                  ligand = c("g1", "g2"), receptor = c("g3", "g4"),
                                  annotation = "Secreted Signaling", row.names = lrs))
rna@dr <- list(umap = matrix(seq_len(8L), 4L, dimnames = list(cells, c("UMAP1", "UMAP2"))))
rna@options$parameter <- list(n = 1, interaction.range = 100)
cell.prob <- spatstat.sparse::sparse3Darray(
  i = c(1L, 4L), j = c(2L, 3L), k = c(1L, 2L), x = c(0.2, 0.8),
  dims = c(4L, 4L, 2L), dimnames = list(cells, cells, lrs))
rna@net$prob.cell <- cell.prob
rna@net$pval.cell <- array(seq_len(32L) / 40, c(4, 4, 2), dimnames = list(cells, cells, lrs))
rna@net$prob <- array(seq_len(8L) / 10, c(2, 2, 2), dimnames = list(groups, groups, lrs))
rna@net$pval <- array(seq_len(8L) / 20, c(2, 2, 2), dimnames = list(groups, groups, lrs))
rna@net$count.cell <- matrix(seq_len(16L), 4L, dimnames = list(cells, cells))
rna@net$sum <- matrix(seq_len(4L), 2L, dimnames = list(groups, groups))
rna@net$group.merged <- factor(c("AB", "AB"))
rna@net$count.merged <- matrix(4, 1, 1, dimnames = list("AB", "AB"))
rna@net$weight.merged <- matrix(10, 1, 1, dimnames = list("AB", "AB"))
rna@net$centr <- array(seq_len(8L), c(2, 2, 2), dimnames = list(c("outdeg", "indeg"), groups, lrs))
rna@net$centr.cell <- array(seq_len(16L), c(2, 4, 2), dimnames = list(c("outdeg", "indeg"), cells, lrs))
rna@net$field <- list(outgoing = array(seq_len(16L), c(4, 2, 2), dimnames = list(cells, c("x", "y"), lrs)))
rna@net$tmp <- list(Lavg = matrix(1:8, nrow = 2L))
rna@netP$prob <- array(seq_len(4L), c(2, 2, 1), dimnames = list(groups, groups, "P1"))
rna@netP$prob.cell <- spatstat.sparse::sparse3Darray(
  i = 2L, j = 3L, k = 1L, x = 0.7, dims = c(4L, 4L, 1L), dimnames = list(cells, cells, "P1"))
rna@netP$pathways <- "P1"
check("upstream object remains a valid 15-slot object", length(methods::slotNames(rna)) == 15L &&
        isTRUE(methods::validObject(rna)))
out <- update(rna)
check("old 15 slots converted to 11-slot dev class", length(methods::slotNames(out)) == 11L &&
        "assay" %in% methods::slotNames(out) && !"data" %in% methods::slotNames(out))
check("S4 and explicit validation pass", isTRUE(methods::validObject(out)) && isTRUE(validate(out)))
check("norm and raw values preserved", identical(out@assay$norm, rna@data) && identical(out@assay$raw, rna@data.raw))
check("subset gene layers reordered with values", identical(rownames(out@assay$signaling), c("g1", "g3")) &&
        identical(as.matrix(out@assay$signaling), as.matrix(rna@data.signaling[c("g1", "g3"), ])) &&
        identical(out@assay$smooth, out@assay$signaling))
check("metadata identities DB LR dr options preserved", identical(out@meta, meta) &&
        identical(unname(out@idents), unname(rna@idents)) && identical(names(out@idents), cells) &&
        identical(out@DB, rna@DB) && identical(out@LR, rna@LR) && identical(out@dr, rna@dr) &&
        identical(out@misc$.var.features, rna@var.features) && identical(out@misc$.legacy$options, rna@options))
check("cell sparse probability values and names preserved", inherits(out@net$cell$prob, "SparseChatArray") &&
        identical(dimnames(out@net$cell$prob), dimnames(cell.prob)) &&
        isTRUE(all.equal(as.matrix(out@net$cell$prob[[1L]]), as.matrix(cell.prob[, , 1L]))) &&
        isTRUE(all.equal(as.matrix(out@net$cell$prob[[2L]]), as.matrix(cell.prob[, , 2L]))))
list.input <- rna
list.input@net$prob.cell <- list(LR1 = unname(as.matrix(cell.prob[, , 1L])),
                                 LR2 = unname(as.matrix(cell.prob[, , 2L])))
list.output <- update(list.input)@net$cell$prob
check("unnamed matrix layers receive known cell axes", identical(dimnames(list.output), list(cells, cells, lrs)) &&
        identical(unname(as.matrix(list.output[[2L]])), unname(as.matrix(cell.prob[, , 2L]))))
set.seed(23)
indices <- sample.int(90L * 90L * 40L, 10000L)
i <- ((indices - 1L) %% 90L) + 1L
j <- (((indices - 1L) %/% 90L) %% 90L) + 1L
k <- ((indices - 1L) %/% (90L * 90L)) + 1L
nodes <- paste0("n", seq_len(90L))
large <- spatstat.sparse::sparse3Darray(i = i, j = j, k = k, x = seq_along(i) / 10000,
  dims = c(90L, 90L, 40L), dimnames = list(nodes, nodes, paste0("L", seq_len(40L))))
converted <- SpatialCellChatDev:::.sc_update_sparse_array(large, nodes, "test sparse array")
check("many-layer sparse conversion retains all triplets", length(converted) == 40L &&
        sum(vapply(unclass(converted), function(layer) length(layer@x), integer(1))) == length(i) &&
        all(vapply(seq_len(20L), function(n) converted[[k[n]]][i[n], j[n]] == n / 10000, logical(1))))
corrupt <- large
corrupt$k[1L] <- 41L
check("invalid sparse layer indices are rejected",
      expect_error(SpatialCellChatDev:::.sc_update_sparse_array(corrupt, nodes, "corrupt"),
                   "invalid sparse triplet indices"))
check("dense pvalues and group probabilities preserve values and names",
        isTRUE(all.equal(as.matrix(out@net$cell$pval[[2L]]), rna@net$pval.cell[, , 2L])) &&
        identical(dimnames(out@net$group$prob), dimnames(rna@net$prob)) &&
        isTRUE(all.equal(as.matrix(out@net$group$prob[[1L]]), rna@net$prob[, , 1L])) &&
        isTRUE(all.equal(as.matrix(out@net$group$pval[[2L]]), rna@net$pval[, , 2L])))
check("counts weight centrality and field preserved", identical(out@net$cell$count, rna@net$count.cell) &&
        identical(out@net$group$weight, rna@net$sum) && identical(out@net$centr, rna@net$centr) &&
        identical(out@net$centr.cell, rna@net$centr.cell) &&
        isTRUE(all.equal(as.matrix(out@net$cell$field$outgoing[[1L]]), rna@net$field$outgoing[, , 1L])))
check("merged group summaries remain available", identical(out@net$group$count.merged, rna@net$count.merged) &&
        identical(out@net$group$weight.merged, rna@net$weight.merged) &&
        identical(names(out@net$group$group.merged), groups))
check("pathway results preserved", identical(out@netP$pathways, "P1") &&
        isTRUE(all.equal(as.matrix(out@netP$group$prob[[1L]]), rna@netP$prob[, , 1L])) &&
        identical(dimnames(out@netP$cell$prob), dimnames(rna@netP$prob.cell)))
communication <- SpatialCellChatDev::subsetCommunication(out, thresh = 0.5)
check("migrated results are usable by subsetCommunication", nrow(communication) == 8L &&
        setequal(communication$interaction_name, lrs) &&
        all(communication$ligand %in% c("g1", "g2")))
check("current object validates unchanged", identical(update(out), out))
historic <- rna
attr(historic, "var.features") <- c("g1", "g3")
check("historical character var.features migrated", identical(update(historic)@misc$.var.features,
                                                            list(features = c("g1", "g3"))))
wrong.class.package <- rna
attr(wrong.class.package, "class") <- structure("SpatialCellChat", package = "SpatialCellChatDev")
check("actual stored slots override a conflicting class definition",
      identical(update(wrong.class.package)@assay$norm, rna@data))
spatial <- make.old("spatial")
spatial@images$result.computeCellDistance <- list(d.spatial = matrix(99, 4, 4))
spatial@images$custom <- "retained"
sp <- update(spatial)
check("spatial coordinates and ratio/tol migrated", identical(rownames(sp@images$coordinates), cells) &&
        identical(colnames(sp@images$coordinates), c("x", "y")) &&
        identical(unname(sp@images$coordinates), unname(spatial@images$coordinates)) &&
        identical(sp@images$spatial.factors$ratio, 0.5) && identical(sp@images$spatial.factors$tol, 2))
spatial3d <- spatial
spatial3d@images$coordinates <- cbind(spatial3d@images$coordinates, z_cent = seq_along(cells))
check("three-dimensional legacy coordinates are preserved",
      identical(colnames(update(spatial3d)@images$coordinates), c("x", "y", "z")))
check("stale cache not reused and legacy data retained", !length(sp@images$.distance) &&
        identical(sp@misc$.legacy$images, spatial@images) && isTRUE(validate(sp)))
plot <- SpatialCellChatDev::spatialDimPlot(sp, image = FALSE)
plot.file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot.file)
print(plot)
grDevices::dev.off()
check("migrated spatial object renders a plot", file.exists(plot.file) && file.info(plot.file)$size > 1000)
spatial@images$scale.factors <- as.data.frame(spatial@images$spatial.factors)
spatial@images$spatial.factors <- NULL
check("legacy scale.factors data frame migrated", identical(update(spatial)@images$spatial.factors$tol, 2))
uncalibrated <- spatial
uncalibrated@images$scale.factors <- NULL
uncalibrated.result <- update(uncalibrated)
check("missing calibration remains explicitly uncalibrated",
      identical(uncalibrated.result@images$coordinate.system$calibrated, FALSE) &&
        all(c("ratio", "tol") %in% names(uncalibrated.result@images$spatial.factors)) &&
        isTRUE(validate(uncalibrated.result)))
spatial@net <- rna@net
spatial@netP <- rna@netP
spatial@LR <- rna@LR
spatial@data.signaling <- rna@data.signaling
legacy.merged <- suppressWarnings(SpatialCellChat::mergeSpatialCellChat(
  list(s1 = spatial, s2 = spatial), add.names = c("s1", "s2"),
  merge.data = TRUE, cell.prefix = TRUE, show.plot = FALSE))
merged <- update(legacy.merged)
snapshot <- merged@misc$.datasets$s1
check("merged object with expression migrates to validated snapshots",
      identical(merged@misc$.mode, "merged") && identical(names(merged@misc$.datasets), c("s1", "s2")) &&
        identical(merged@assay$norm, legacy.merged@data) && isTRUE(methods::validObject(merged)) &&
        isTRUE(validate(merged)))
check("merged dataset cell axes and probability values are preserved",
      identical(rownames(snapshot$meta), paste0("s1_", cells)) &&
        identical(rownames(snapshot$images$coordinates), paste0("s1_", cells)) &&
        identical(dimnames(snapshot$net$cell$prob)[[1L]], paste0("s1_", cells)) &&
        identical(dimnames(snapshot$net$cell$prob)[[3L]], lrs) &&
        isTRUE(all.equal(unname(as.matrix(snapshot$net$cell$prob[[1L]])),
                         unname(as.matrix(cell.prob[, , 1L])))))
legacy.merged@data <- matrix(numeric(), 0L, 0L)
check("merged object missing normalized data is rejected with recovery direction",
      expect_error(update(legacy.merged), "migrate the original single-dataset objects"))
check("non-SpatialCellChat input rejected", expect_error(update(list()), "must be a SpatialCellChat S4 object"))
legacy.root <- file.path("tests_dev", "test_data", "old_version_chat_objs",
                         "simulation data", "simulation data", "ChatObjects")
legacy.3cts.path <- file.path(legacy.root, "chat_3cts.rds")
legacy.4cts.path <- file.path(legacy.root, "chat_4cts.rds")
if (file.exists(legacy.3cts.path) && file.exists(legacy.4cts.path)) {
  legacy.3cts <- readRDS(legacy.3cts.path)
  migrated.3cts <- update(legacy.3cts)
  check("historical SpatialChat 3cts fixture migrates without its defining package",
        identical(unname(class(legacy.3cts))[1L], "SpatialChat") &&
          ncol(migrated.3cts@assay$norm) == 300L &&
          identical(levels(migrated.3cts@idents), c("A", "B", "C")) &&
          isTRUE(validate(migrated.3cts)))
  legacy.4cts <- readRDS(legacy.4cts.path)
  migrated.4cts <- update(legacy.4cts)
  check("historical SpatialChat 4cts fixture migrates stored inference and centrality",
        identical(unname(class(legacy.4cts))[1L], "SpatialChat") &&
          ncol(migrated.4cts@assay$norm) == 500L &&
          inherits(migrated.4cts@net$cell$prob, "SparseChatArray") &&
          identical(dimnames(migrated.4cts@net$cell$prob)[[3L]],
                    rownames(migrated.4cts@LR$LRsig)) &&
          isTRUE(validate(migrated.4cts)))
  historical.centrality <- attr(legacy.4cts, "net", exact = TRUE)$centr.cell
  check("migration preserves the stored 4cts centrality exactly",
        identical(migrated.4cts@net$centr.cell, historical.centrality))
  recomputed.4cts <- SpatialCellChatDev::netAnalysis_computeCentrality(
    migrated.4cts, slot.name = "net", do.group = FALSE
  )
  recomputed.centrality <- recomputed.4cts@net$centr.cell
  axes <- dimnames(historical.centrality)
  check("recomputed 4cts centrality has the same named measures, cells and LR layers",
        identical(dim(recomputed.centrality), dim(historical.centrality)) &&
          all(mapply(setequal, axes, dimnames(recomputed.centrality))))
  recomputed.centrality <- recomputed.centrality[axes[[1L]], axes[[2L]], axes[[3L]], drop = FALSE]
  check("centrality recomputed on the preserved 4cts network matches by name",
        max(abs(recomputed.centrality - historical.centrality)) <= 1e-10)
} else {
  cat("[SKIP] Historical 3cts/4cts fixtures are unavailable under ", legacy.root, "\n", sep = "")
}
unknown <- rna
attr(unknown, "future.slot") <- 1
check("unknown slot layout rejected", expect_error(update(unknown), "Unrecognized SpatialCellChat slot layout"))
bad <- rna
bad@net$prob <- array(1, c(2, 2, 1), dimnames = list(c("A", "WRONG"), groups, "LR1"))
check("wrong network names rejected", expect_error(update(bad), "dimnames do not match"))
bad <- rna
dimnames(bad@net$pval)[[3L]] <- rev(lrs)
check("misaligned probability and p-value layers rejected",
      expect_error(update(bad), "probability and p-value axes do not match"))
bad <- spatial
bad@images$coordinates <- bad@images$coordinates[-1L, ]
check("misaligned coordinates rejected", expect_error(update(bad), "coordinates must be"))
cat("UPDATE_TESTS_OK order=", order, "\n", sep = "")
