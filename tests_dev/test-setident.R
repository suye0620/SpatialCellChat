# setIdent(): the cell-group axis must stay consistent.
# Covers permutation / rename / mixed, the different-partition branch (group-level
# results dropped, cell-level kept), the level map helper and updateClusterLabels.
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
source("R/analysis.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err_msg <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
msg_capture <- function(expr) {
  out <- NULL
  withCallingHandlers(expr, message = function(m) {
    if (is.null(out)) out <<- conditionMessage(m)
    invokeRestart("muffleMessage")
  })
  out
}

## ---- fixture: 20 cells, 4 groups, 3 LR --------------------------------------
set.seed(31)
cells <- paste0("c", 1:20)
genes <- c("G1", "G2", "R1", "R2", "L3", "R3")
expr <- matrix(round(runif(length(genes) * length(cells)), 3), nrow = length(genes),
               dimnames = list(genes, cells))
coor <- expand.grid(x = seq(10, 50, 10), y = seq(10, 40, 10))
rownames(coor) <- cells
meta.df <- data.frame(label = factor(rep(c("A", "B", "C", "D"), each = 5)),
                      row.names = cells)

chat <- createSpatialCellChat(
  object = as(expr, "dgCMatrix"), meta = meta.df, group.by = "label",
  input.assay = "norm", datatype = "spatial",
  coordinates = coor, spatial.factors = list(ratio = 1, tol = 5)
)
chat@DB <- list(
  interaction = data.frame(
    ligand = c("G1", "G2", "L3"), receptor = c("R1", "R2", "R3"),
    agonist = c("", "", ""), antagonist = c("", "", ""),
    co_A_receptor = c("", "", ""), co_I_receptor = c("", "", ""),
    annotation = c("Secreted Signaling", "Cell-Cell Contact", "Secreted Signaling"),
    interaction_name = c("G1-R1", "G2-R2", "L3-R3"),
    pathway_name = c("P1", "P2", "P3"), stringsAsFactors = FALSE
  ),
  complex = data.frame(), cofactor = data.frame()
)
chat@LR$LRsig <- chat@DB$interaction[, c("ligand", "receptor", "agonist", "antagonist",
                                         "co_A_receptor", "co_I_receptor", "annotation",
                                         "interaction_name", "pathway_name")]
rownames(chat@LR$LRsig) <- chat@LR$LRsig$interaction_name
assay(chat, "signaling") <- assay(chat, "norm")
chat0 <- chat          # no network yet

chat <- computeCommunProb(chat, scale.distance = 9, n = 1, nthreads = 1L, verbose = FALSE)
chat <- computeAvgCommunProb(chat, min.cells.sr = 1, do.permutation = TRUE, nboot = 20,
                             seed.use = 1L, nthreads = 1L, verbose = FALSE)
chat <- computeCommunProbPathway(chat)
chat <- aggregateNet(chat)
# do.group = TRUE is what writes the group-level `centr`; the default writes
# `centr.cell`, and `net$centr` would then partial-match it.
chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = TRUE)
chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = TRUE)

val <- function(arr, gi, gj, layer = 1L) as.numeric(arr[gi, gj, layer])
check("fixture has every group-level structure",
      inherits(chat@net$group$prob, "SparseChatArray") &&
        inherits(chat@netP$group$prob, "SparseChatArray") &&
        !is.null(chat@net$group$count) && !is.null(chat@net[["centr"]]))

## ---- 1. permutation ---------------------------------------------------------
old.prob <- chat@net$group$prob
old.centr <- chat@net[["centr"]]
new.order <- c("D", "C", "B", "A")
perm <- setIdent(chat, levels = new.order)
check("permutation: idents follow the requested order",
      identical(levels(idents(perm)), new.order))
check("permutation: net$group$prob is reordered",
      identical(dimnames(perm@net$group$prob)[[1]], new.order) &&
        identical(dimnames(perm@net$group$prob)[[2]], new.order))
check("permutation: netP$group$prob follows net",
      identical(dimnames(perm@netP$group$prob)[[1]], new.order))
check("permutation: count matrix is reordered",
      identical(dimnames(perm@net$group$count)[[1]], new.order))
check("permutation: centr dim2 is reordered",
      identical(dimnames(perm@net[["centr"]])[[2]], new.order))
check("permutation: values stay attached to their group",
      isTRUE(all.equal(val(old.prob, "A", "B"), val(perm@net$group$prob, "A", "B"))))
check("permutation: centr values stay attached to their group",
      isTRUE(all.equal(unname(old.centr[1, "A", 1]), unname(perm@net[["centr"]][1, "A", 1]))))
check("permutation: object validates", isTRUE(validObject(perm)))

## ---- 2. rename --------------------------------------------------------------
renamed.src <- chat
renamed.src@meta$newid <- factor(paste0("g", as.integer(idents(chat))))
ren <- setIdent(renamed.src, ident.use = "newid")
check("rename: idents use the new labels",
      identical(levels(idents(ren)), paste0("g", 1:4)))
check("rename: net layers are relabelled",
      identical(dimnames(ren@net$group$prob)[[1]], paste0("g", 1:4)) &&
        identical(dimnames(ren@net$group$prob)[[3]], dimnames(old.prob)[[3]]))
check("rename: values are untouched",
      isTRUE(all.equal(unname(as.matrix(old.prob[, , 1])),
                       unname(as.matrix(ren@net$group$prob[, , 1])))))
check("rename: names(idents) match the metadata rows",
      identical(names(idents(ren)), rownames(meta(ren))))
check("rename: object validates", isTRUE(validObject(ren)))

## ---- 3. rename + reorder in one call ---------------------------------------
mixed <- setIdent(renamed.src, ident.use = "newid", levels = c("g4", "g3", "g2", "g1"))
step <- setIdent(setIdent(renamed.src, ident.use = "newid"), levels = c("g4", "g3", "g2", "g1"))
check("mixed rename+reorder equals doing them in sequence",
      identical(levels(idents(mixed)), levels(idents(step))) &&
        isTRUE(all.equal(unname(as.matrix(mixed@net$group$prob[, , 1])),
                         unname(as.matrix(step@net$group$prob[, , 1])))))

## ---- 4. different partition: merge -----------------------------------------
merged <- chat
merged@meta$coarse <- factor(ifelse(as.character(idents(chat)) %in% c("A", "B"), "X", "Y"))
m.msg <- msg_capture(merged.out <- setIdent(merged, ident.use = "coarse"))
check("merge: idents collapse to the new groups",
      identical(levels(idents(merged.out)), c("X", "Y")))
check("merge: group-level results are dropped",
      !length(merged.out@net$group) && !length(merged.out@netP$group) &&
        is.null(merged.out@net$centr) && is.null(merged.out@netP$centr))
check("merge: cell-level probabilities are kept",
      inherits(merged.out@net$cell$prob, "SparseChatArray"))
check("merge: message names what to re-run",
      !is.null(m.msg) && grepl("computeAvgCommunProb", m.msg))
check("merge: object still validates", isTRUE(validObject(merged.out)))

## ---- 5. different partition: split -----------------------------------------
split.src <- chat
split.src@meta$fine <- factor(paste0(as.character(idents(chat)),
                                     rep(c("1", "2"), length.out = length(cells))))
msg_capture(split.out <- setIdent(split.src, ident.use = "fine"))
check("split: idents get the finer groups", nlevels(idents(split.out)) == 8L)
check("split: group-level results are dropped", !length(split.out@net$group))

## ---- 6. level map helper (empty levels) ------------------------------------
old.id <- factor(c(rep("A", 3), rep("B", 2)), levels = c("A", "B", "gone"))
new.id <- factor(c(rep("a", 3), rep("b", 2)), levels = c("b", "a", "gone2"))
mp <- .sc_group_level_map(old.id, new.id, levels(old.id), levels(new.id))
check("level map pairs observed levels by content", identical(unname(mp[c("A", "B")]), c("a", "b")))
check("level map pairs the empty level with the leftover one",
      identical(unname(mp["gone"]), "gone2"))

## ---- 7. no network yet ------------------------------------------------------
plain <- setIdent(chat0, levels = c("D", "C", "B", "A"))
check("setIdent works on a network-free object",
      identical(levels(idents(plain)), c("D", "C", "B", "A")) &&
        !length(plain@net$group))

## ---- 8. updateClusterLabels -------------------------------------------------
ul <- updateClusterLabels(chat, new.cluster.name = c("w", "x", "y", "z"))
check("updateClusterLabels renames", identical(levels(idents(ul)), c("w", "x", "y", "z")))
check("updateClusterLabels keeps the results", inherits(ul@net$group$prob, "SparseChatArray") &&
        identical(dimnames(ul@net$group$prob)[[1]], c("w", "x", "y", "z")))
check("updateClusterLabels records the labels in meta",
      identical(as.character(meta(ul)[["new.labels"]]), as.character(idents(ul))))
check("updateClusterLabels validates", isTRUE(validObject(ul)))
check("updateClusterLabels rejects a many-to-one mapping",
      grepl("one-to-one",
            err_msg(updateClusterLabels(chat, old.cluster.name = c("A", "B", "C", "D"),
                                        new.cluster.name = c("m", "m", "n", "n")))))

## ---- 9. one definition only -------------------------------------------------
src <- list.files("R", pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
defs <- unlist(lapply(src, function(f) grep("^setIdent <- function", readLines(f, warn = FALSE), value = TRUE)))
check("exactly one setIdent definition in R/", length(defs) == 1L)

cat("All setIdent checks passed.\n")
