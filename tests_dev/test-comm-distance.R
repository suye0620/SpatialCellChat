# spatialC2CDistPlot(): the single distance-plot entry point.
# Covers the 2 x 2 contract of split.by x enriched.only, the geometric mode that
# needs no inference, the cache self-compute path and the three deprecated aliases.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(Rcpp); library(ggplot2)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/visualization.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err_msg <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
warn_msg <- function(expr) {
  msg <- NULL
  withCallingHandlers(expr, warning = function(w) {
    if (is.null(msg)) msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  msg
}
titles <- function(panels) vapply(panels, function(p) p$labels$title, character(1))

## ---- fixture: 5 x 3 grid, 3 groups, 3 LR (secreted / contact / unknown) -----
cells <- paste0("c", 1:15)
genes <- c("L1", "R1", "L2", "R2", "L3", "R3")
set.seed(7)
expr <- matrix(round(runif(length(genes) * length(cells)), 3), nrow = length(genes),
               dimnames = list(genes, cells))
coor <- expand.grid(x = seq(10, 50, 10), y = seq(10, 30, 10))
rownames(coor) <- cells
meta.df <- data.frame(label = factor(rep(c("A", "B", "C"), each = 5)), row.names = cells)

chat <- createSpatialCellChat(
  object = as(expr, "dgCMatrix"), meta = meta.df, group.by = "label",
  input.assay = "norm", datatype = "spatial",
  coordinates = coor, spatial.factors = list(ratio = 1, tol = 5)
)
chat@DB <- list(
  interaction = data.frame(
    ligand = c("L1", "L2", "L3"), receptor = c("R1", "R2", "R3"),
    agonist = c("", "", ""), antagonist = c("", "", ""),
    co_A_receptor = c("", "", ""), co_I_receptor = c("", "", ""),
    annotation = c("Secreted Signaling", "Cell-Cell Contact", "Custom Signaling"),
    interaction_name = c("L1-R1", "L2-R2", "L3-R3"),
    pathway_name = c("P1", "P2", "P3"), stringsAsFactors = FALSE
  ),
  complex = data.frame(), cofactor = data.frame()
)
chat@LR$LRsig <- chat@DB$interaction[, c("ligand", "receptor", "agonist", "antagonist",
                                         "co_A_receptor", "co_I_receptor", "annotation",
                                         "interaction_name", "pathway_name")]
rownames(chat@LR$LRsig) <- chat@LR$LRsig$interaction_name
assay(chat, "signaling") <- assay(chat, "norm")

## ---- 1. geometric mode needs no inference ----------------------------------
geo <- spatialC2CDistPlot(chat)                      # defaults: enriched.only = FALSE
check("default call works on a never-inferred object", inherits(geo, "patchwork"))
geo.panels <- .sc_c2c_panels(chat, "annotation", FALSE, 0.5)
check("geometric mode always yields the two geometry panels",
      identical(names(geo.panels), c("all", "contact.range")))
check("geometric panel titles", identical(unname(titles(geo.panels)),
                                          c("All cell pairs", "Contact-range pairs")))
check("x axis ends at the cached interaction.range + 10",
      identical(geo.panels[[1]]$scales$scales[[2]]$limits, c(0, 260)))

## ---- 2. rows 3 and 4 are the same output ------------------------------------
p3 <- .sc_c2c_panels(chat, "annotation", FALSE, 0.5)
p4 <- .sc_c2c_panels(chat, "distance", FALSE, 0.5)
check("annotation + FALSE equals distance + FALSE (names)", identical(names(p3), names(p4)))
check("annotation + FALSE equals distance + FALSE (titles)", identical(titles(p3), titles(p4)))

w.ignored <- warn_msg(spatialC2CDistPlot(chat, split.by = "distance"))
check("supplying split.by with enriched.only = FALSE warns",
      !is.null(w.ignored) && grepl("ignored", w.ignored))
check("no warning when split.by is left to its default",
      is.null(warn_msg(spatialC2CDistPlot(chat))))

## ---- 3. cache self-compute path --------------------------------------------
check("no usable distance cache before anything runs",
      !.spatial_distance_cache_matches(images(chat, ".distance"),
                                       interaction.range = 250, contact.range = 10,
                                       ratio = 1, tol = 5))
cache <- .sc_c2c_geometry_cache(chat)
check("geometry cache is built on demand",
      is.list(cache) && all(c("d.spatial", "adj.contact", ".parameters") %in% names(cache)))
check("built cache records the default interaction range",
      identical(cache$.parameters$interaction.range, 250))
cache.wide <- .sc_c2c_geometry_cache(chat, interaction.range = 400)
check("explicit interaction.range reaches computeCellDistance",
      identical(cache.wide$.parameters$interaction.range, 400))

## ---- 4. inference + enriched modes -----------------------------------------
chat.cc <- computeCommunProb(chat, scale.distance = 9, n = 1, nthreads = 1L, verbose = FALSE)
cache.cc <- images(chat.cc, ".distance")

df.en <- .sc_comm_distance_frame(chat.cc, TRUE, cache.cc)
nnz <- sum(vapply(unclass(chat.cc@net$cell$prob), function(m) length(m@x), integer(1)))
check("enriched frame has one row per non-zero probability", nrow(df.en) == nnz && nnz > 0)
check("enriched frame carries the documented columns",
      identical(names(df.en), c("i", "j", "k", "layer", "annotation", "distance", "contact")))

panels.ann <- .sc_c2c_panels(chat.cc, "annotation", TRUE, 0.5)
check("row 1: diffusible + contact + other panels",
      identical(names(panels.ann), c("diffusible", "contact", "other")))
check("row 1: unknown annotation lands in Other signaling",
      identical(unname(titles(panels.ann)["other"]), "Other signaling"))
valid <- sum(!is.na(df.en$distance) & df.en$distance > 0)
check("row 1: every plot-able record lands in exactly one panel",
      valid > 0 &&
        valid == sum(vapply(panels.ann, function(p) nrow(p$data), integer(1))))

panels.dist <- .sc_c2c_panels(chat.cc, "distance", TRUE, 0.5)
check("row 2: all inferred + contact-range panels",
      identical(names(panels.dist), c("all", "contact.range")))
check("row 2: the second panel is a subset of the first",
      nrow(panels.dist[[2]]$data) <= nrow(panels.dist[[1]]$data))

check("row 1 returns a patchwork through the public function",
      inherits(spatialC2CDistPlot(chat.cc, split.by = "annotation", enriched.only = TRUE),
               "patchwork"))
check("enriched mode without a cache is rejected with a hint",
      grepl("computeCommunProb",
            err_msg(spatialC2CDistPlot(chat, enriched.only = TRUE, split.by = "distance"))))

## ---- 5. only one name is left ------------------------------------------------
check("spatialCommDistPlot no longer exists",
      !exists("spatialCommDistPlot", mode = "function"))
check("spatialCCCDistPlot no longer exists",
      !exists("spatialCCCDistPlot", mode = "function"))
check("communicationDistPlot no longer exists",
      !exists("communicationDistPlot", mode = "function"))
check("communicationDistPlot2 no longer exists",
      !exists("communicationDistPlot2", mode = "function"))
check("spatialC2CDistPlot is the only distance-plot export",
      sum(grepl("^export\\(.*DistPlot\\)$", readLines("NAMESPACE"))) == 1L)

cat("All spatialC2CDistPlot checks passed.\n")
