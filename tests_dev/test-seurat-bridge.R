# asSeurat(): the Seurat bridge on the final 11-slot schema.
# Covers the three slot.name modes (assay / net / netP), the what axis
# (incoming / outgoing / merge), list input and the deprecated wrappers.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(pbapply)
  library(Seurat); library(SeuratObject)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err_msg <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
warn_msg <- function(expr) {
  msg <- NULL
  withCallingHandlers(expr, warning = function(w) {
    # Keep the first warning: downstream Seurat calls emit their own noise.
    if (is.null(msg)) msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  msg
}
# GetAssayData(slot=) is soft-deprecated in SeuratObject 5 but is the only
# spelling Seurat 4.4 understands; silence the notice.
gs <- function(obj, what, assay = "Spatial")
  suppressWarnings(Seurat::GetAssayData(obj, assay = assay, slot = what))

## ---- fixture ----------------------------------------------------------------
cells <- paste0("cell", seq_len(8L))
genes <- c("L1", "R1", "S1", "S2")
expr <- matrix(rpois(length(genes) * length(cells), 3), nrow = length(genes),
               dimnames = list(genes, cells))
meta.df <- data.frame(label = factor(rep(c("A", "B"), each = 4L)), row.names = cells)
coords <- cbind(x = seq_len(8L), y = rep(c(1, 2), 4L))
rownames(coords) <- cells

lr <- c("L1_R1", "L2_R2")
chat <- createSpatialCellChat(
  as(expr, "dgCMatrix"), meta = meta.df, group.by = "label", input.assay = "norm",
  datatype = "spatial", coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 1)
)
chat@DB <- list(
  interaction = data.frame(
    interaction_name = lr,
    pathway_name = c("PA", "PB"),
    ligand = c("L1", "L2"),
    receptor = c("R1", "R2"),
    annotation = rep("Secreted Signaling", 2L),
    stringsAsFactors = FALSE
  ),
  complex = data.frame(),
  cofactor = data.frame(),
  geneInfo = data.frame(Symbol = genes, stringsAsFactors = FALSE)
)

## ---- 1. slot.name = "assay": basic conversion -------------------------------
seu <- asSeurat(chat)
check("cell count preserved", ncol(seu) == ncol(assay(chat, "norm")))
check("cell names preserved", identical(colnames(seu), cells))
check("idents carry the CellChat groups",
      identical(as.character(Seurat::Idents(seu)), as.character(idents(chat))))
check("idents levels preserved", identical(levels(Seurat::Idents(seu)), levels(idents(chat))))
check("group.cellchat column written to metadata",
      identical(seu$group.cellchat, idents(chat)))
check("cell-level metadata carried over",
      identical(as.character(seu$label), as.character(meta.df$label)))
check("assay.name defaults to Spatial", identical(Seurat::DefaultAssay(seu), "Spatial"))

## ---- 2. layers follow the Seurat mapping ------------------------------------
counts.layer <- gs(seu, "counts")
check("counts equals assay$norm when no raw layer exists",
      isTRUE(all.equal(as.matrix(counts.layer), as.matrix(assay(chat, "norm")))))
check("feature names preserved", identical(rownames(counts.layer), genes))

chat.scale <- chat
assay(chat.scale, "scale") <- as.matrix(expr)
seu.scale <- asSeurat(chat.scale)
check("assay$scale becomes scale.data",
      isTRUE(all.equal(as.matrix(gs(seu.scale, "scale.data")), as.matrix(expr))))

chat.raw <- createSpatialCellChat(
  as(expr, "dgCMatrix"), meta = meta.df, group.by = "label", input.assay = "raw",
  datatype = "spatial", coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 1)
)
seu.rawlayers <- asSeurat(chat.raw)
check("assay$raw becomes counts",
      isTRUE(all.equal(unname(as.matrix(gs(seu.rawlayers, "counts"))),
                       unname(as.matrix(chat.raw@assay$raw)))))
check("assay$norm becomes data",
      isTRUE(all.equal(unname(as.matrix(gs(seu.rawlayers, "data"))),
                       unname(as.matrix(chat.raw@assay$norm)))))

chat.sub <- subsetData(chat, features = c("L1", "R1"))
seu.sub <- asSeurat(chat.sub)
subsets <- seu.sub[["Spatial"]]@misc[["subsets"]]
check("gene-subset layers are parked in assay misc",
      identical(names(subsets), "signaling") &&
        identical(rownames(subsets$signaling), c("L1", "R1")))

## ---- 3. spatial image -------------------------------------------------------
check("Spatial assay carries a SlideSeq image", length(Seurat::Images(seu)) == 1L)
img <- seu@images[[1]]
check("image keeps the coordinates",
      isTRUE(all.equal(unname(as.matrix(img@coordinates[, c("x", "y")])),
                       unname(as.matrix(coords[, c("x", "y")])))))

seu.rna <- asSeurat(chat, assay.name = "RNA")
check("RNA assay attaches no image", length(Seurat::Images(seu.rna)) == 0L)

## ---- 4. counts override / group.by / argument gating ------------------------
override <- matrix(1, nrow = 3L, ncol = length(cells),
                   dimnames = list(c("X1", "X2", "X3"), cells))
seu.ov <- asSeurat(chat, counts = override)
check("counts override replaces assay$norm",
      identical(rownames(seu.ov), c("X1", "X2", "X3")))
check("counts override with other genes does not mix in the norm data layer",
      identical(rownames(gs(seu.ov, "data")), c("X1", "X2", "X3")))

bad <- override
colnames(bad) <- NULL
check("counts without dimnames rejected",
      grepl("dimnames", err_msg(asSeurat(chat, counts = bad))))

seu.grp <- asSeurat(chat, group.by = "label")
check("group.by selects another metadata column",
      identical(as.character(Seurat::Idents(seu.grp)), as.character(meta.df$label)))

w.what <- warn_msg(asSeurat(chat, what = "merge"))
check("`what` is rejected outside the centrality modes",
      !is.null(w.what) && grepl("only applies", w.what))

## ---- 5. slot.name = "net": centrality directions ----------------------------
indeg <- matrix(seq_len(16) / 10, nrow = length(cells), dimnames = list(cells, lr))
outdeg <- matrix(seq_len(16) / 20, nrow = length(cells), dimnames = list(cells, lr))
centr <- array(0, dim = c(2L, length(cells), length(lr)),
               dimnames = list(c("indeg", "outdeg"), cells, lr))
centr["indeg", , ] <- indeg
centr["outdeg", , ] <- outdeg

chat.net <- chat
chat.net@net$centr.cell <- centr

seu.in <- asSeurat(chat.net, slot.name = "net", what = "incoming")
# Seurat rewrites "_" to "-" in feature names; the original stays in features.info.
seu.lr <- gsub("_", "-", lr)
check("net/incoming: features are the LR pairs", identical(rownames(seu.in), seu.lr))
check("net/incoming: counts carry indeg",
      isTRUE(all.equal(unname(as.matrix(gs(seu.in, "counts"))), unname(t(indeg)))))
check("net/incoming: idents still come from the CellChat groups",
      identical(as.character(Seurat::Idents(seu.in)), as.character(idents(chat))))
check("net: features.info carries the database annotation",
      identical(seu.in[["Spatial"]]@misc[["features.info"]]$pathway_name, c("PA", "PB")))

seu.out <- asSeurat(chat.net, slot.name = "net", what = "outgoing")
check("net/outgoing: counts carry outdeg",
      isTRUE(all.equal(unname(as.matrix(gs(seu.out, "counts"))), unname(t(outdeg)))))

seu.merge <- asSeurat(chat.net, slot.name = "net", what = "merge")
check("net/merge: both directions in one assay",
      identical(rownames(seu.merge), gsub("_", "-", c(paste0(lr, "-I"), paste0(lr, "-O")))))
check("net/merge: rows keep their source layer",
      identical(seu.merge[["Spatial"]]@misc[["features.info"]]$orig.features, c(lr, lr)))

seu.raw <- asSeurat(chat.net, slot.name = "net", what = "incoming", normalize.chat = FALSE)
check("normalize.chat = FALSE leaves data equal to counts",
      isTRUE(all.equal(as.matrix(gs(seu.raw, "data")), as.matrix(gs(seu.raw, "counts")))))
check("normalize.chat = TRUE rewrites the data layer",
      !isTRUE(all.equal(as.matrix(gs(seu.in, "data")), as.matrix(gs(seu.in, "counts")))))

check("net without centrality is rejected",
      grepl("netAnalysis_computeCentrality",
            err_msg(asSeurat(chat, slot.name = "net", what = "incoming"))))

w.counts <- warn_msg(asSeurat(chat.net, slot.name = "net", counts = override))
check("`counts` is rejected outside slot.name = 'assay'",
      !is.null(w.counts) && grepl("only applies", w.counts))

## ---- 6. slot.name = "netP" --------------------------------------------------
centr.p <- array(0, dim = c(2L, length(cells), 2L),
                 dimnames = list(c("indeg", "outdeg"), cells, c("PA", "PB")))
centr.p["indeg", , ] <- indeg
centr.p["outdeg", , ] <- outdeg
chat.np <- chat
chat.np@netP$centr.cell <- centr.p

seu.np <- asSeurat(chat.np, slot.name = "netP", what = "incoming")
check("netP: features are the pathways", identical(rownames(seu.np), c("PA", "PB")))
check("netP: pathway_name mirrors the feature name",
      identical(seu.np[["Spatial"]]@misc[["features.info"]]$pathway_name, c("PA", "PB")))

## ---- 7. list input ----------------------------------------------------------
lst <- list(A = chat.net, B = chat.net, C = chat.net)
seu.list <- suppressWarnings(asSeurat(lst, slot.name = "net", what = "incoming"))
check("every element of the list reaches the merge", ncol(seu.list) == 3L * length(cells))
check("Condition records the list names",
      identical(levels(seu.list$Condition), names(lst)))
check("Condition is per cell, not per object",
      identical(as.integer(table(seu.list$Condition)), rep(length(cells), 3L)))

## ---- 8. legacy names are out of the package ---------------------------------
check("CellChat2Seurat is not defined in R/",
      !exists("CellChat2Seurat", mode = "function"))
check("mergeChatIntoSeurat is not defined in R/",
      !exists("mergeChatIntoSeurat", mode = "function"))

legacy_env <- new.env(parent = globalenv())
sys.source("R/legacy/seurat-bridge.R", envir = legacy_env)
check("legacy file still defines both bridges",
      is.function(legacy_env$CellChat2Seurat) &&
        is.function(legacy_env$mergeChatIntoSeurat))
check("legacy CellChat2Seurat still delegates to asSeurat",
      ncol(legacy_env$CellChat2Seurat(chat, check.object = FALSE)) == length(cells))

cat("All Seurat bridge checks passed.\n")
