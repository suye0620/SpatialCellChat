# Regression checks for pre-inference functions on the final object schema

setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(cli)
  library(pbapply)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

cells <- paste0("cell", seq_len(6L))
genes <- c("L1", "R1", "L2", "R2", "S1", "S2")
expr <- matrix(1, nrow = length(genes), ncol = length(cells),
               dimnames = list(genes, cells))
expr["L1", 1:3] <- 10
expr["R1", 4:6] <- 10
expr["S1", 1:3] <- 8
expr["S2", 4:6] <- 8
meta <- data.frame(
  label = factor(c("A", "A", "A", "B", "B", "B")),
  row.names = cells
)
coordinates <- cbind(x = seq_len(6L), y = rep(c(1, 2), 3L))
rownames(coordinates) <- cells

chat <- createSpatialCellChat(
  expr,
  meta = meta,
  group.by = "label",
  input.assay = "norm",
  datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = 1)
)
chat@DB <- list(
  interaction = data.frame(
    ligand = c("L1", "CPLX", "L2", "NOPE"),
    receptor = c("R1", "R2", "R2", "R1"),
    annotation = c("Cell-Cell Contact", "Secreted Signaling",
                   "ECM-Receptor", "Secreted Signaling"),
    interaction_name = c("L1-R1", "CPLX-R2", "L2-R2", "NOPE-R1"),
    stringsAsFactors = FALSE
  ),
  complex = data.frame(
    subunit_1 = "S1",
    subunit_2 = "S2",
    row.names = "CPLX",
    stringsAsFactors = FALSE
  )
)

chat <- subsetData(chat)
check("subsetData writes assay signaling layer", identical(
  rownames(chat@assay$signaling), genes
))
check("subsetData excludes genes absent from the database", !"NOPE" %in% rownames(chat@assay$signaling))
check("subsetData reads datatype from misc", identical(
  chat@DB$interaction$interaction_name,
  c("CPLX-R2", "NOPE-R1", "L2-R2", "L1-R1")
))
check("subsetData preserves assay norm", identical(
  dimnames(chat@assay$norm), dimnames(expr)
))

chat <- identifyOverExpressedGenes(
  chat,
  selection.method = "wilcox",
  thresh.p = 1,
  thresh.fc = 0
)
check("overexpressed genes are stored in misc var.features", is.list(chat@misc$.var.features) &&
        "features" %in% names(chat@misc$.var.features) &&
        "features.info" %in% names(chat@misc$.var.features))
check("overexpressed gene metadata is a data frame", is.data.frame(
  chat@misc$.var.features$features.info
))
check("overexpressed genes use signaling layer", all(
  chat@misc$.var.features$features %in% rownames(chat@assay$signaling)
))
markers <- identifyOverExpressedGenes(
  chat,
  selection.method = "wilcox",
  return.object = FALSE,
  thresh.p = 1,
  thresh.fc = 0
)
check("overexpressed gene return.object FALSE returns a data frame", is.data.frame(markers))

chat@misc$.var.features$features <- c("L1", "R1", "S1", "S2", "R2")
chat <- identifyOverExpressedInteractions(chat, variable.both = TRUE)
check("LRsig is stored in the LR slot", is.data.frame(chat@LR$LRsig))
check("complex and direct LR pairs use new layers and metadata", identical(
  chat@LR$LRsig$interaction_name,
  c("CPLX-R2", "L1-R1")
))
selected_pairs <- identifyOverExpressedInteractions(
  chat,
  variable.both = FALSE,
  return.object = FALSE
)
check("interaction return.object FALSE returns a data frame", is.data.frame(selected_pairs))
check("pre-inference object remains valid", isTRUE(validateSpatialCellChat(chat)))
grid_chat <- identifyOverExpressedGenes(
  chat,
  do.grid = TRUE,
  cellsize = c(2, 2),
  selection.method = "wilcox",
  thresh.p = 1,
  thresh.fc = 0
)
check("grid feature selection stores grid in images", is.list(grid_chat@images$.grid) &&
        inherits(grid_chat@images$.grid$object, "SpatialCellChat"))
check("grid feature selection preserves grid summary", is.numeric(
  grid_chat@images$.grid$recommended.contact.range
) && length(grid_chat@images$.grid$recommended.contact.range) == 1L &&
        is.numeric(grid_chat@images$.grid$within.nGrid))
check("grid object remains valid", isTRUE(validateSpatialCellChat(
  grid_chat@images$.grid$object
)))
cat("All pre-inference new-schema checks passed.\n")
