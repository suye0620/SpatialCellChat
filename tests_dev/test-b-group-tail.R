# The last two B-group helpers: computeAveExpr() and
# identifyOverExpressedLigandReceptor(). Both used to read slots the 11-slot
# class deleted, and the marker function additionally had three real bugs.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err_msg <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))

## ---- fixture ----------------------------------------------------------------
set.seed(41)
cells <- paste0("c", 1:12)
genes <- c("G1", "R1", "SU1", "SU2", "R2", "ZZZ")
expr <- matrix(round(runif(length(genes) * length(cells)), 3), nrow = length(genes),
               dimnames = list(genes, cells))
coor <- expand.grid(x = seq(10, 40, 10), y = seq(10, 30, 10))
rownames(coor) <- cells
meta.df <- data.frame(label = factor(rep(c("A", "B", "C"), each = 4)),
                      batch = factor(rep(c("b1", "b2"), each = 6)),
                      row.names = cells)

chat <- createSpatialCellChat(
  object = as(expr, "dgCMatrix"), meta = meta.df, group.by = "label",
  input.assay = "norm", datatype = "spatial",
  coordinates = coor, spatial.factors = list(ratio = 1, tol = 5)
)
assay(chat, "signaling") <- assay(chat, "norm")[c("G1", "R1", "SU1", "R2"), , drop = FALSE]

## ---- 1. computeAveExpr ------------------------------------------------------
ave.sig <- computeAveExpr(chat, type = "median")
check("computeAveExpr defaults to the signaling layer",
      setequal(rownames(ave.sig), c("G1", "R1", "SU1", "R2")))
check("computeAveExpr columns are the cell groups",
      identical(colnames(ave.sig), levels(idents(chat))))

expected <- tapply(as.numeric(as.matrix(assay(chat, "signaling"))["G1", ]),
                   as.character(idents(chat)), median)
check("computeAveExpr matches a manual per-group median",
      isTRUE(all.equal(as.numeric(ave.sig["G1", ]), as.numeric(expected))))

ave.norm <- computeAveExpr(chat, slot.name = "norm", type = "median")
check("computeAveExpr reads the full norm layer",
      nrow(ave.norm) == nrow(assay(chat, "norm")))

ave.batch <- computeAveExpr(chat, group.by = "batch", type = "median")
check("computeAveExpr honours group.by",
      identical(colnames(ave.batch), levels(meta(chat)[["batch"]])))

check("computeAveExpr rejects the old slot.name values",
      grepl("should be one of", err_msg(computeAveExpr(chat, slot.name = "data.signaling"))))

empty <- chat
assay(empty, "signaling") <- NULL
check("empty signaling error names subsetData()",
      grepl("subsetData", err_msg(computeAveExpr(empty))))

check("empty smooth error names projectData()",
      grepl("projectData", err_msg(computeAveExpr(chat, slot.name = "smooth"))))

assay(chat, "smooth") <- assay(chat, "signaling")
check("computeAveExpr reads the smooth layer once it exists",
      identical(rownames(computeAveExpr(chat, slot.name = "smooth", type = "median")),
                rownames(assay(chat, "smooth"))))

## ---- 2. identifyOverExpressedLigandReceptor ---------------------------------
chat@DB <- list(
  interaction = data.frame(
    ligand = c("G1", "CPLX"), receptor = c("R1", "R2"),
    annotation = c("Secreted Signaling", "Secreted Signaling"),
    interaction_name = c("G1-R1", "CPLX-R2"),
    pathway_name = c("P1", "P2"), stringsAsFactors = FALSE
  ),
  complex = data.frame(row.names = "CPLX", subunit_1 = "SU1", subunit_2 = "SU2",
                       stringsAsFactors = FALSE),
  cofactor = data.frame()
)
chat@misc$.var.features[["features.info"]] <- data.frame(
  features = c("G1", "R1", "SU1", "ZZZ"),
  pval = c(0.01, 0.02, 0.03, 0.90),
  group = c("A", "B", "C", "A"),
  stringsAsFactors = FALSE
)

oe <- identifyOverExpressedLigandReceptor(chat, features.name = "features")
tbl <- oe@misc$.var.features[["features.LR"]]
check("the LR table is written to misc$.var.features",
      is.data.frame(tbl) && nrow(tbl) >= 2L)
check("direct ligand/receptor features are kept", all(c("G1", "R1") %in% tbl$features))
check("complex subunits are reported under the complex name",
      "CPLX" %in% tbl$features && !"SU1" %in% tbl$features)
check("features outside the DB are dropped", !"ZZZ" %in% tbl$features)
check("object still validates", isTRUE(validObject(oe)))

only <- identifyOverExpressedLigandReceptor(chat, features.name = "features",
                                            features = c("R1", "ZZZ"),
                                            return.object = FALSE)
check("the features argument filters the input table",
      identical(unique(only$features), "R1"))

no.markers <- chat
no.markers@misc$.var.features <- list()
check("missing marker table names identifyOverExpressedGenes()",
      grepl("identifyOverExpressedGenes",
            err_msg(identifyOverExpressedLigandReceptor(no.markers))))

bad.db <- chat
bad.db@DB$interaction <- data.frame(foo = "x")
check("DB without ligand/receptor columns is rejected clearly",
      grepl("ligand and receptor",
            err_msg(identifyOverExpressedLigandReceptor(bad.db))))

no.complex <- chat
no.complex@DB$complex <- NULL
check("a missing complex table is tolerated",
      is.data.frame(identifyOverExpressedLigandReceptor(no.complex, return.object = FALSE)))

cat("All B-group tail checks passed.\n")
