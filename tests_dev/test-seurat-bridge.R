# CellChat2Seurat(): the Seurat bridge on the final 11-slot schema.
# Covers the slot mapping (assay$norm / meta() / idents() / images$coordinates),
# the SlideSeq image hand-off, the counts override and the matrix round trip.
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

## ---- fixture ----------------------------------------------------------------
cells <- paste0("cell", seq_len(8L))
genes <- c("L1", "R1", "S1", "S2")
expr <- matrix(rpois(length(genes) * length(cells), 3), nrow = length(genes),
               dimnames = list(genes, cells))
meta.df <- data.frame(label = factor(rep(c("A", "B"), each = 4L)), row.names = cells)
coords <- cbind(x = seq_len(8L), y = rep(c(1, 2), 4L))
rownames(coords) <- cells

chat <- createSpatialCellChat(
  as(expr, "dgCMatrix"), meta = meta.df, group.by = "label", input.assay = "norm",
  datatype = "spatial", coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 1)
)

## ---- 1. basic conversion ----------------------------------------------------
seu <- CellChat2Seurat(chat, check.object = FALSE)
check("cell count preserved", ncol(seu) == ncol(assay(chat, "norm")))
check("cell names preserved", identical(colnames(seu), cells))
check("idents carry the CellChat groups",
      identical(as.character(Seurat::Idents(seu)), as.character(idents(chat))))
check("idents levels preserved", identical(levels(Seurat::Idents(seu)), levels(idents(chat))))
check("group.cellchat column written to metadata",
      identical(seu$group.cellchat, idents(chat)))
check("cell-level metadata carried over",
      identical(as.character(seu$label), as.character(meta.df$label)))

## ---- 2. matrix round trip ---------------------------------------------------
roundtrip <- Seurat::GetAssayData(seu, assay = "Spatial", slot = "counts")
check("matrix equals assay$norm",
      isTRUE(all.equal(as.matrix(roundtrip), as.matrix(assay(chat, "norm")))))
check("feature names preserved", identical(rownames(roundtrip), genes))

## ---- 3. spatial image -------------------------------------------------------
check("Spatial assay carries a SlideSeq image", length(Seurat::Images(seu)) == 1L)
img <- seu@images[[1]]
check("image keeps the coordinates",
      isTRUE(all.equal(unname(as.matrix(img@coordinates[, c("x", "y")])),
                       unname(as.matrix(coords[, c("x", "y")])))))

seu.rna <- CellChat2Seurat(chat, assay = "RNA", check.object = FALSE)
check("RNA assay attaches no image", length(Seurat::Images(seu.rna)) == 0L)

## ---- 4. counts override -----------------------------------------------------
override <- matrix(1, nrow = 3L, ncol = length(cells),
                   dimnames = list(c("X1", "X2", "X3"), cells))
seu.ov <- CellChat2Seurat(chat, counts = override, check.object = FALSE)
check("counts override replaces assay$norm",
      identical(rownames(seu.ov), c("X1", "X2", "X3")))

bad <- override
colnames(bad) <- NULL
check("counts without dimnames rejected",
      inherits(tryCatch(CellChat2Seurat(chat, counts = bad, check.object = FALSE),
                        error = function(e) e), "error"))

## ---- 5. group.by ------------------------------------------------------------
seu.grp <- CellChat2Seurat(chat, group.by = "label", check.object = FALSE)
check("group.by selects another metadata column",
      identical(as.character(Seurat::Idents(seu.grp)), as.character(meta.df$label)))

cat("All Seurat bridge checks passed.\n")
