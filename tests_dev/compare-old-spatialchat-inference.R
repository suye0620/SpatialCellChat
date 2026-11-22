args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:5 || !args[[1L]] %in% c("upstream", "dev") ||
    !args[[2L]] %in% c("3cts", "4cts") ||
    (length(args) >= 4L && !identical(args[[4L]], "--trace-permutations")) ||
    (length(args) == 5L && !identical(args[[1L]], "dev")))
  stop(paste("Usage: Rscript tests_dev/compare-old-spatialchat-inference.R",
             "<upstream|dev> <3cts|4cts> <output-dir> [--trace-permutations [upstream-trace.rds]]"))

engine <- args[[1L]]
dataset <- args[[2L]]
output.dir <- args[[3L]]
trace.permutations <- length(args) >= 4L
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
fixture.dir <- file.path("tests_dev", "test_data", "old_version_chat_objs",
                         "simulation data", "simulation data", "ChatObjects")
object <- readRDS(file.path(fixture.dir, paste0("chat_", dataset, ".rds")))

if (identical(engine, "upstream")) {
  suppressPackageStartupMessages(library(SpatialCellChat))
  attr(object, "class") <- structure("SpatialCellChat", package = "SpatialCellChat")
  package.version <- as.character(utils::packageVersion("SpatialCellChat"))
} else {
  suppressPackageStartupMessages(library(SpatialCellChatDev))
  object <- SpatialCellChatDev::updateSpatialCellChat(object)
  package.version <- as.character(utils::packageVersion("SpatialCellChatDev"))
}

database <- data.frame(
  ligand = c("L03", "L13", "L23", "L33", "L43", "L53", "L04", "L14", "L24", "L34", "L44", "L54"),
  receptor = c("R01", "R11", "R21", "R31", "R41", "R51", "R02", "R12", "R22", "R32", "R42", "R52"),
  interaction_name = c("L03-R01", "L13-R11", "L23-R21", "L33-R31", "L43-R41", "L53-R51", "L04-R02", "L14-R12", "L24-R22", "L34-R32", "L44-R42", "L54-R52"),
  pathway_name = c("L03-R01", "L13-R11", "L23-R21", "L33-R31", "L43-R41", "L53-R51", "L04-R02", "L14-R12", "L24-R22", "L34-R32", "L44-R42", "L54-R52"),
  annotation = if (identical(dataset, "3cts")) {
    rep("Cell-Cell Contact", 12L)
  } else {
    c("Secreted Signaling", "Secreted Signaling", "Cell-Cell Contact",
      "Secreted Signaling", "Secreted Signaling", "Secreted Signaling",
      rep("Secreted Signaling", 6L))
  },
  stringsAsFactors = FALSE
)
gene.info <- data.frame(
  Symbol = c("L03", "L04", "L13", "L23", "L33", "L43", "L53", "L14", "L24", "L34", "L44", "L54",
             "R01", "R02", "R11", "R21", "R31", "R41", "R51", "R12", "R22", "R32", "R42", "R52"),
  stringsAsFactors = FALSE
)
object@DB <- updateCellChatDB(db = database, gene_info = gene.info)
object <- subsetData(object)
if (identical(engine, "dev")) {
  setEnvironment(4, future.globals.maxSize = 10000 * 1024^2, verbose = 0L)
} else {
  setEnvironment(4, future.globals.maxSize = 10000 * 1024^2)
}
set.seed(1L)
object <- identifyOverExpressedGenes(object, selection.method = "moransi")
object <- identifyOverExpressedInteractions(object, variable.both = FALSE)
object <- computeCommunProb(
  object,
  distance.use = TRUE, use.AGAN = FALSE,
  contact.dependent.forced = FALSE, contact.dependent = TRUE,
  interaction.range = 250, tol = 1, scale.distance = 0.2,
  contact.range = 10
)

as.layers <- function(value) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "SparseChatArray")) return(unclass(value))
  if (length(dim(value)) == 3L)
    return(lapply(seq_len(dim(value)[3L]), function(i) value[, , i, drop = TRUE]))
  value
}
cell.prob <- function(x) {
  if (identical(engine, "dev")) x@net$cell$prob else x@net$prob.cell
}
group.prob <- function(x) {
  if (identical(engine, "dev")) x@net$group$prob else x@net$prob
}
group.pval <- function(x) {
  if (identical(engine, "dev")) x@net$group$pval else x@net$pval
}
features <- function(x) {
  if (identical(engine, "dev")) x@misc$.var.features else x@var.features
}
snapshot <- list(
  engine = engine,
  dataset = dataset,
  package.version = package.version,
  genes = features(object),
  LRsig = object@LR$LRsig,
  prob.compute = as.layers(cell.prob(object))
)

object <- filterProbability(object)
snapshot$prob.filterProbability <- as.layers(cell.prob(object))
if (identical(dataset, "3cts")) {
  object <- filterCommunication(object, min.cells = NULL, min.links = 5, min.cells.sr = 3)
} else {
  object <- filterCommunication(object, min.cells = NULL, min.links = 10, min.cells.sr = 10)
}
snapshot$prob.filterCommunication <- as.layers(cell.prob(object))

average <- computeAvgCommunProb
if (trace.permutations) {
  source(file.path("tests_dev", "legacy-permutation-trace.R"))
  trace <- new.env(parent = emptyenv())
  upstream.seed <- NULL
  if (length(args) == 5L) {
    reference <- readRDS(args[[5L]])
    stopifnot(identical(reference$result$engine, "upstream"),
              identical(reference$result$dataset, dataset))
    upstream.seed <- reference$permutations$seed.before.drawing
    stopifnot(is.integer(upstream.seed), length(upstream.seed) > 1L)
  }
  average <- trace_group_permutations(average, engine, trace, upstream.seed)
}
if (identical(dataset, "3cts")) {
  object <- average(object, do.permutation = TRUE, nboot = 100)
} else {
  object <- average(object, do.permutation = TRUE, min.percent = 0, nboot = 100)
}
object <- filterCommunication(object, min.cells = 10,
                               min.links = NULL, min.cells.sr = NULL)
snapshot$prob.final <- as.layers(cell.prob(object))
snapshot$group.prob <- as.layers(group.prob(object))
snapshot$group.pval <- as.layers(group.pval(object))

object <- netAnalysis_computeCentrality(object, slot.name = "net", do.group = FALSE)
snapshot$centrality.cell <- object@net$centr.cell
if (identical(dataset, "3cts")) {
  object <- netAnalysis_computeCentrality(object, slot.name = "net", do.group = TRUE)
  snapshot$centrality.group <- object@net$centr
} else {
  snapshot$centrality.group <- NULL
}

output.stem <- file.path(output.dir, paste0(engine, "-", dataset))
if (trace.permutations) {
  stopifnot(exists("value", envir = trace, inherits = FALSE))
  saveRDS(list(result = snapshot, permutations = trace$value),
          paste0(output.stem, "-permutation-trace.rds"))
} else {
  saveRDS(snapshot, paste0(output.stem, ".rds"))
}
cat("INFERENCE_OK", engine, dataset, "LR=", nrow(object@LR$LRsig), "\n")
