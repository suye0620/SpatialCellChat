# Regression tests for spatialFeaturePlot on the current schema.
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
  library(methods); library(Matrix); library(dplyr); library(future.apply)
  library(cli); library(ggplot2)
})
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/database.R")
source("R/analysis.R")
source("R/modeling.R")
source("R/visualization.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err <- function(expr) tryCatch({ force(expr); "" }, error = conditionMessage)
cells <- paste0("c", 1:4)
genes <- c("L1", "R1", "L2", "R2", "L3", "R3")
expr <- Matrix(matrix(c(0, 1, 2, 3, 4, 5, 1, 0, 3, 2, 5, 4,
                        2, 3, 0, 1, 4, 5, 3, 2, 1, 0, 5, 4),
                      nrow = 6, dimnames = list(genes, cells)), sparse = TRUE)
coords <- data.frame(x = 1:4, y = 11:14, row.names = cells)
meta <- data.frame(label = factor(c("A", "A", "B", "B")),
                   score = c(10, 20, 30, 40), row.names = cells)
chat <- createSpatialCellChat(expr, meta = meta, group.by = "label", datatype = "spatial",
                              coordinates = coords, spatial.factors = list(ratio = 1, tol = 1))
chat@assay$signaling <- expr[1:4, , drop = FALSE]
lr <- data.frame(interaction_name = c("LR1", "LR2", "LR3"),
                 ligand = c("L1", "L2", "L3"), receptor = c("R1", "R2", "R3"),
                 pathway_name = c("P1", "P1", "P2"), stringsAsFactors = FALSE,
                 row.names = c("arbitrary1", "arbitrary2", "arbitrary3"))
chat@LR <- list(LRsig = lr)
chat@DB <- list(interaction = lr,
                complex = data.frame(subunit_1 = character(), subunit_2 = character(), row.names = character()),
                geneInfo = data.frame(Symbol = genes))
mat <- function(i = integer(), j = integer(), x = numeric(), n = 2L)
  as(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n)), "dgCMatrix")
group.prob <- SparseChatArray(list(LR1 = mat(1, 2, 0.4), LR2 = mat(2, 1, 0.6),
                                   LR3 = mat(1, 1, 0)),
                              dimnames = list(c("A", "B"), c("A", "B"), lr$interaction_name))
group.pval <- SparseChatArray(list(LR1 = mat(1, 2, 0.05), LR2 = mat(2, 1, 0.06),
                                   LR3 = mat()),
                              dimnames = dimnames(group.prob))
cell.prob <- SparseChatArray(list(LR1 = mat(1, 2, 0, 4L), LR2 = mat(3, 4, 0.1, 4L),
                                  LR3 = mat(n = 4L)),
                             dimnames = list(cells, cells, lr$interaction_name))
chat@net$group <- list(prob = group.prob, pval = group.pval, LR.sig = "LR2")
chat@net$cell <- list(prob = cell.prob, LR.sig = "LR2")

p <- spatialFeaturePlot(chat, features = "L1")
check("default signaling layer draws expression", inherits(p, "ggplot") &&
        identical(p$data$feature.data, as.numeric(expr["L1", ])))
check("legacy and current signaling names agree",
      identical(spatialFeaturePlot(chat, features = "L1", slot.name = "data.signaling")$data, p$data))
check("legacy and current norm names expose all genes",
      identical(spatialFeaturePlot(chat, features = "L3", slot.name = "data")$data,
                spatialFeaturePlot(chat, features = "L3", slot.name = "norm")$data))

reordered <- chat
reordered@images$coordinates <- reordered@images$coordinates[rev(cells), , drop = FALSE]
reordered@meta <- reordered@meta[rev(cells), , drop = FALSE]
check("coordinates are aligned to expression cells", identical(spatialFeaturePlot(reordered, features = "L1")$data, p$data))
pm <- spatialFeaturePlot(reordered, features = "score")
check("metadata is aligned to expression cells", identical(pm$data$feature.data, c(10, 20, 30, 40)))
bad <- reordered
rownames(bad@images$coordinates)[1] <- "wrong"
check("mismatched cell names fail explicitly", grepl("same unique cell names", err(spatialFeaturePlot(bad, features = "L1"))))
missing <- chat
missing@assay$signaling <- NULL
check("missing signaling layer names prerequisite", grepl("subsetData", err(spatialFeaturePlot(missing, features = "L1"))))

check("group selection retains p equal to threshold and excludes explicit zero",
      identical(.sc_spatial_feature_active_lr(chat, lr$interaction_name, TRUE, 0.05), "LR1"))
check("cell selection ignores group p-values and explicit zeros",
      identical(.sc_spatial_feature_active_lr(chat, lr$interaction_name, FALSE, 0.05), "LR2"))
no.pval <- chat
no.pval@net$group$pval <- NULL
check("group selection without p-values uses positive probability",
      identical(.sc_spatial_feature_active_lr(no.pval, lr$interaction_name, TRUE, 0.05), c("LR1", "LR2")))
bad.pval <- chat
dimnames(bad.pval@net$group$pval) <- list(c("B", "A"), c("A", "B"), lr$interaction_name)
check("mismatched p-value axes are rejected", grepl("matching names", err(.sc_spatial_feature_active_lr(bad.pval, "LR1", TRUE, 0.05))))

explicit <- spatialFeaturePlot(chat, pairLR.use = c("LR1", "LR2"), ncol = 1)
check("explicit group selection uses current threshold and returns all gene panels",
      inherits(explicit, "patchwork") && length(explicit$patches$plots) == 1L &&
        explicit$patches$plots[[1]]$labels$title == "L1" && explicit$labels$title == "R1")
pathway <- spatialFeaturePlot(chat, signaling = "P1", ncol = 1, thresh = 0.001)
check("pathway selection continues using LR.sig cache", inherits(pathway, "patchwork") &&
        pathway$patches$plots[[1]]$labels$title == "L2" && pathway$labels$title == "R2")
binary <- spatialFeaturePlot(chat, pairLR.use = "LR2", do.group = FALSE,
                             do.binary = TRUE, cutoff = 1)
check("binary explicit cell plot uses LR expression", inherits(binary, "ggplot") &&
        identical(as.character(binary$data$feature.data), c("Both", "Both", "None", "None")))
binary.path <- spatialFeaturePlot(chat, signaling = "P1", do.group = FALSE,
                                  do.binary = TRUE, cutoff = 1)
check("binary pathway matches LR by interaction name, not row name",
      identical(binary.path$data$feature.data, binary$data$feature.data))
complex.chat <- chat
complex.chat@LR$LRsig$ligand[1] <- "C1"
complex.chat@DB$complex <- data.frame(subunit_1 = "L1", subunit_2 = "L2", row.names = "C1")
complex.plot <- spatialFeaturePlot(complex.chat, pairLR.use = "LR1", do.binary = TRUE, cutoff = 0.5)
check("binary complex requires expression of every subunit",
      identical(as.character(complex.plot$data$feature.data), c("R1", "C1", "R1", "Both")))
check("continuous and binary plots build", inherits(ggplotGrob(p), "gtable") &&
        inherits(ggplotGrob(binary), "gtable") && inherits(patchwork::patchworkGrob(explicit), "gtable"))
check("empty LR selection errors", grepl("after filtering", err(spatialFeaturePlot(chat, pairLR.use = "LR3"))))
check("unknown LR reports an error when filtering is disabled",
      grepl("Unknown interaction_name", err(spatialFeaturePlot(chat, pairLR.use = "missing", enriched.only = FALSE))))
