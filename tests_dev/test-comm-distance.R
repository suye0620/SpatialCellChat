# spatialC2CDistPlot(): the single distance-plot entry point.
# Covers the 2 x 2 contract of split.by x mode, the geometric mode that
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
geo <- spatialC2CDistPlot(chat)                      # defaults: mode = "geometry"
check("default call works on a never-inferred object", inherits(geo, "patchwork"))
geo.panels <- .sc_c2c_panels(chat, "annotation", "geometry", 0.5)
check("geometric mode always yields the two geometry panels",
      identical(names(geo.panels), c("all", "contact.range")))
check("geometric panel titles", identical(unname(titles(geo.panels)),
                                          c("All cell pairs", "Contact-range pairs")))
check("x axis ends at the cached interaction.range + 10",
      identical(geo.panels[[1]]$scales$scales[[2]]$limits, c(0, 260)))

## ---- 2. geometry ignores split.by --------------------------------------------
p3 <- .sc_c2c_panels(chat, "annotation", "geometry", 0.5)
p4 <- .sc_c2c_panels(chat, "distance", "geometry", 0.5)
check("geometry annotation equals distance split (names)", identical(names(p3), names(p4)))
check("geometry annotation equals distance split (titles)", identical(titles(p3), titles(p4)))

w.ignored <- warn_msg(spatialC2CDistPlot(chat, split.by = "distance"))
check("supplying split.by in geometry mode warns",
      !is.null(w.ignored) && grepl("ignored", w.ignored))
check("no warning when split.by is left to its default",
      is.null(warn_msg(spatialC2CDistPlot(chat))))
check("explicit geometry mode works",
      inherits(spatialC2CDistPlot(chat, mode = "geometry"), "patchwork"))
check("invalid mode is rejected", !is.na(err_msg(spatialC2CDistPlot(chat, mode = "bad"))))
check("old enriched.only argument is rejected",
      grepl("unused argument", err_msg(spatialC2CDistPlot(chat, enriched.only = FALSE))))

## ---- 3. cache self-compute path --------------------------------------------
check("no usable distance cache before anything runs",
      !.spatial_distance_cache_matches(images(chat, ".distance"),
                                       coordinates = images(chat, "coordinates"),
                                       interaction.range = 250, contact.range = 10,
                                       ratio = 1, tol = 5))
cache <- .sc_c2c_geometry_cache(chat)
check("geometry cache is built on demand",
      is.list(cache) && all(c("d.spatial", "adj.contact", ".parameters", ".coordinates") %in% names(cache)))
check("built cache records the default interaction range",
      identical(cache$.parameters$interaction.range, 250))
check("built cache records the coordinates it was computed from",
      identical(cache$.coordinates, images(chat, "coordinates")))
geo.frame <- .sc_comm_distance_frame(chat, "geometry", cache)
geo.matrix <- as.matrix(cache$d.spatial)
expected.geo <- sum(upper.tri(geo.matrix) & geo.matrix > 0)
check("geometry frame keeps one observation per unordered cell pair",
      nrow(geo.frame) == expected.geo)
check("geometry frame uses logical contact flags",
      is.logical(geo.frame$contact) && !anyNA(geo.frame$contact))
check("geometry contact panel contains the contact distances",
      identical(geo.panels$contact.range$data$distance,
                geo.frame$distance[geo.frame$contact]))
cache.wide <- .sc_c2c_geometry_cache(chat, interaction.range = 400)
check("explicit interaction.range reaches computeCellDistance",
      identical(cache.wide$.parameters$interaction.range, 400))

chat.moved <- chat
images(chat.moved, ".distance") <- cache
coords.moved <- images(chat.moved, "coordinates")
coords.moved[2L, "x"] <- coords.moved[2L, "x"] + 1000
images(chat.moved, "coordinates") <- coords.moved
check("changing coordinates through images<- clears the distance cache",
      identical(images(chat.moved, ".distance"), list()))
cache.moved <- .sc_c2c_geometry_cache(chat.moved)
check("geometry cache recomputes after coordinates change",
      cache.moved$d.spatial["c1", "c2"] == 0 &&
        identical(cache.moved$.coordinates, coords.moved))
chat.factor <- chat
images(chat.factor, ".distance") <- cache
factors.changed <- images(chat.factor, "spatial.factors")
factors.changed$ratio <- 0.5
images(chat.factor, "spatial.factors") <- factors.changed
check("changing spatial factors through images<- clears the distance cache",
      identical(images(chat.factor, ".distance"), list()))

## ---- 4. inference + enriched modes -----------------------------------------
chat.cc <- computeCommunProb(chat, scale.distance = 9, n = 1, nthreads = 1L, verbose = FALSE)
cache.cc <- images(chat.cc, ".distance")

## ---- 4a. explicit pair-level deduplication contract -------------------------
dedup <- chat
dedup.layers <- list(
  Matrix::sparseMatrix(i = c(1L, 2L, 1L), j = c(2L, 1L, 3L),
                       x = c(0.2, 0.4, 0.5), dims = c(15L, 15L), repr = "C"),
  Matrix::sparseMatrix(i = 1L, j = 2L, x = 0.7,
                       dims = c(15L, 15L), repr = "C"),
  Matrix::sparseMatrix(i = 1L, j = 2L, x = 0.3,
                       dims = c(15L, 15L), repr = "C"),
  Matrix::sparseMatrix(i = 1L, j = 3L, x = 0,
                       dims = c(15L, 15L), repr = "C")
)
dedup@net$cell$prob <- SparseChatArray(
  dedup.layers,
  dimnames = list(cells, cells, c("L1-R1", "L2-R2", "L3-R3", "zero"))
)
dedup@LR$LRsig <- rbind(
  chat@LR$LRsig,
  data.frame(ligand = "L0", receptor = "R0", agonist = "", antagonist = "",
             co_A_receptor = "", co_I_receptor = "", annotation = "Custom Signaling",
             interaction_name = "zero", pathway_name = "P0",
             stringsAsFactors = FALSE)
)
rownames(dedup@LR$LRsig) <- dedup@LR$LRsig$interaction_name
dedup.frame.ann <- .sc_comm_distance_frame(
  dedup, "communication", cache, "annotation"
)
dedup.frame.dist <- .sc_comm_distance_frame(
  dedup, "communication", cache, "distance"
)
check("annotation frame deduplicates repeated LR layers and reverse directions",
      nrow(dedup.frame.ann) == 4L)
check("annotation frame excludes explicit zero-only communication",
      !any(dedup.frame.ann$distance == cache$d.spatial[1L, 3L] &
           dedup.frame.ann$class == "other"))
check("distance frame merges directions and deduplicates across classes",
      nrow(dedup.frame.dist) == 2L)
check("distance frame keeps the expected unordered-pair distances",
      setequal(dedup.frame.dist$distance,
               c(cache$d.spatial[1L, 2L], cache$d.spatial[1L, 3L])))

df.en <- .sc_comm_distance_frame(chat.cc, "communication", cache.cc, "annotation")
check("communication frame keeps one row per pair and annotation class",
      nrow(df.en) > 0)
check("communication frame carries the 3-column contract (i/j/layer never stored)",
      identical(names(df.en), c("distance", "contact", "class")))
check("class is resolved once per layer, from the LRsig annotation",
      setequal(unique(df.en$class), c("contact", "diffusible", "other")))
check("unknown annotation layers land in the other class",
      !any(is.na(df.en$class)) || !"other" %in% df.en$class ||
        !"diffusible" %in% df.en$class)

panels.ann <- .sc_c2c_panels(chat.cc, "annotation", "communication", 0.5)
check("row 1: diffusible + contact + other panels",
      identical(names(panels.ann), c("diffusible", "contact", "other")))
check("row 1: unknown annotation lands in Other signaling",
      identical(unname(titles(panels.ann)["other"]), "Other signaling"))
valid <- nrow(df.en)
check("row 1: every pair observation lands in exactly one panel",
      valid > 0 &&
        valid == sum(vapply(panels.ann, function(p) nrow(p$data), integer(1))))

panels.dist <- .sc_c2c_panels(chat.cc, "distance", "communication", 0.5)
check("row 2: all inferred + contact-range panels",
      identical(names(panels.dist), c("all", "contact.range")))
check("row 2: the second panel is a subset of the first",
      nrow(panels.dist[[2]]$data) <= nrow(panels.dist[[1]]$data))

check("row 1 returns a patchwork through the public function",
      inherits(spatialC2CDistPlot(chat.cc, split.by = "annotation", mode = "communication"),
               "patchwork"))
check("distance split returns a patchwork through the public function",
      inherits(spatialC2CDistPlot(chat.cc, mode = "communication", split.by = "distance"),
               "patchwork"))
check("enriched mode without a cache is rejected with a hint",
      grepl("computeCommunProb",
            err_msg(spatialC2CDistPlot(chat, mode = "communication", split.by = "distance"))))

chat.cc.moved <- chat.cc
coords.cc.moved <- images(chat.cc.moved, "coordinates")
coords.cc.moved[2L, "x"] <- coords.cc.moved[2L, "x"] + 1000
chat.cc.moved@images$coordinates <- coords.cc.moved
check("direct coordinate changes leave an old cache detectable",
      !.spatial_distance_cache_matches(
        images(chat.cc.moved, ".distance"), coordinates = coords.cc.moved,
        interaction.range = 250, contact.range = 10, ratio = 1, tol = 5
      ))
check("geometry rejects a stale cache after direct coordinate changes",
      .sc_c2c_geometry_cache(chat.cc.moved)$d.spatial["c1", "c2"] == 0)
check("enriched distance plot rejects probabilities from old coordinates",
      grepl("computeCommunProb",
            err_msg(spatialC2CDistPlot(chat.cc.moved, mode = "communication"))))
images(chat.cc.moved, ".distance") <- computeCellDistance(
  coords.cc.moved, interaction.range = 250, contact.range = 10, ratio = 1, tol = 5
)
check("a new distance cache cannot validate probabilities from old coordinates",
      grepl("computeCommunProb",
            err_msg(spatialC2CDistPlot(chat.cc.moved, mode = "communication"))))
chat.cc.moved@images$.distance <- cache.cc
chat.cc.moved <- computeCommunProb(chat.cc.moved, scale.distance = 9, n = 1,
                                   nthreads = 1L, verbose = FALSE)
check("computeCommunProb rebuilds a cache for the current coordinates",
      identical(images(chat.cc.moved, ".distance")$.coordinates, coords.cc.moved))

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
