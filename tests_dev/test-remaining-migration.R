# Focused regression checks for remaining analysis and expression-chain migration.
setwd(local({
  args <- commandArgs(FALSE)
  file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file)) dirname(dirname(normalizePath(file))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(igraph)
  library(future.apply); library(ggplot2); library(colorspace); library(Rcpp)
  library(circlize); library(ComplexHeatmap); library(spatstat.sparse)
  library(ggalluvial); library(NMF)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/merging.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")

check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}

out.dir <- "tests_dev/figures/P5/simulated"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
results <- character()
pass <- function(label, value) {
  check(label, value)
  results <<- c(results, paste0("PASS: ", label))
}

chat <- make_simulated_chat()
pathway <- chat@netP$pathways[1]

netp <- subsetCommunication(chat, slot.name = "netP", signaling = pathway)
pass("netP explicit pathway selection returns pathway rows",
     nrow(netp) > 0L && all(netp$pathway_name == pathway) && all(is.na(netp$pval)))

merged <- mergeSpatialCellChat(
  list(one = make_simulated_chat(), two = make_simulated_chat()),
  show.plot = FALSE
)
merged.netp <- subsetCommunication(merged, slot.name = "netP", signaling = pathway)
pass("merged netP explicit pathway selection returns one table per snapshot",
     is.list(merged.netp) && length(merged.netp) == 2L &&
       all(vapply(merged.netp, nrow, integer(1L)) > 0L) &&
       all(vapply(merged.netp, function(x) all(x$pathway_name == pathway), logical(1L))))

k.plot <- selectK(chat, slot.name = "net", pattern = "outgoing", k.range = 2, nrun = 2)
pass("selectK reads sparse group layers and honors a small rank range",
     inherits(k.plot, "ggplot") && all(unique(k.plot$data$k) == 2))

chat@DB$geneInfo <- data.frame(
  Symbol = c("G1", "G2", "L3", "R1", "R2", "R3"),
  stringsAsFactors = FALSE
)
chat@DB$complex <- data.frame(row.names = character())
chat@misc$.var.features[["features.info"]] <- data.frame(
  clusters = c("A", "B", "A", "B"),
  features = c("G1", "R1", "G2", "R2"),
  pvalues = c(0.01, 0.02, 0.03, 0.04),
  logFC = c(1, 2, 3, 4),
  pct.1 = c(0.5, 0.6, 0.7, 0.8),
  pct.2 = c(0.1, 0.2, 0.3, 0.4),
  stringsAsFactors = FALSE
)
deg.net <- netMappingDEG(chat, "features")
pass("netMappingDEG reads misc var.features and maps DEG columns",
     nrow(deg.net) > 0L &&
       all(c("ligand.pvalues", "receptor.pvalues") %in% names(deg.net)) &&
       any(!is.na(deg.net$ligand.pvalues)))

marker.table <- function(effect, pvalue) {
  out <- expand.grid(
    clusters = c("A", "B", "C", "D"),
    features = c("G1", "G2", "L3", "R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )
  out$pvalues <- pvalue
  out$logFC <- effect
  out$pct.1 <- 0.8
  out$pct.2 <- 0.2
  out
}
source.one <- make_simulated_chat()
source.two <- make_simulated_chat()
for (source in c("source.one", "source.two")) {
  value <- get(source)
  value@DB$geneInfo <- chat@DB$geneInfo
  value@misc$.var.features[["features.info"]] <- marker.table(
    if (source == "source.one") 1 else 2,
    if (source == "source.one") 0.01 else 0.02
  )
  assign(source, value)
}
merged.deg <- mergeSpatialCellChat(
  list(one = source.one, two = source.two), show.plot = FALSE
)
snapshot.mapped <- netMappingDEG(merged.deg, "features")
snapshot.effects <- vapply(c("one", "two"), function(dataset) {
  unique(snapshot.mapped$ligand.logFC[snapshot.mapped$datasets == dataset])
}, numeric(1))
pass("merged DEG mapping uses the matching dataset snapshot",
     all(snapshot.effects == c(one = 1, two = 2)) &&
       all(snapshot.mapped$datasets %in% c("one", "two")))

tagged.one <- marker.table(11, 0.001)
tagged.one$datasets <- "one"
tagged.two <- marker.table(22, 0.002)
tagged.two$datasets <- "two"
merged.deg@misc$.var.features[["features.info"]] <- rbind(tagged.one, tagged.two)
tagged.mapped <- netMappingDEG(merged.deg, "features")
tagged.effects <- vapply(c("one", "two"), function(dataset) {
  unique(tagged.mapped$ligand.logFC[tagged.mapped$datasets == dataset])
}, numeric(1))
pass("merged tagged DEG table is filtered by dataset",
     all(tagged.effects == c(one = 11, two = 22)))

missing.deg <- merged.deg
missing.deg@misc$.var.features[["features.info"]] <- NULL
for (dataset in names(missing.deg@misc$.datasets))
  missing.deg@misc$.datasets[[dataset]]$misc$.var.features[["features.info"]] <- NULL
missing.error <- tryCatch({
  netMappingDEG(missing.deg, "features")
  ""
}, error = conditionMessage)
pass("merged DEG mapping rejects missing dataset-aligned evidence",
     grepl("Dataset 'one'.*dataset-tagged merged marker table", missing.error))

single.bar <- plotGeneExpression(chat, features = "G1", type = "bar")
multi.bar <- plotGeneExpression(chat, features = c("G1", "G2"), type = "bar")
ggplot2::ggsave(file.path(out.dir, "plotGeneExpression-bar-single.pdf"), single.bar,
                width = 7, height = 4.5)
ggplot2::ggsave(file.path(out.dir, "plotGeneExpression-bar-multi.pdf"), multi.bar,
                width = 8, height = 5)
pass("single-feature and multi-feature bar expression plots render",
     inherits(single.bar, "patchwork") && inherits(multi.bar, "patchwork") &&
       file.info(file.path(out.dir, "plotGeneExpression-bar-single.pdf"))$size > 1000 &&
       file.info(file.path(out.dir, "plotGeneExpression-bar-multi.pdf"))$size > 1000)

spatial.aggregate <- netVisual_aggregate(chat, signaling = pathway, layout = "spatial")
spatial.individual <- netVisual_individual(
  chat, signaling = pathway,
  pairLR.use = chat@LR$LRsig$interaction_name[chat@LR$LRsig$pathway_name == pathway][1],
  layout = "spatial"
)
ggplot2::ggsave(file.path(out.dir, "netVisual-aggregate-spatial.pdf"), spatial.aggregate,
                width = 7, height = 6)
ggplot2::ggsave(file.path(out.dir, "netVisual-individual-spatial.pdf"),
                patchwork::wrap_plots(spatial.individual), width = 7, height = 6)
pass("aggregate and individual spatial network layouts render",
     inherits(spatial.aggregate, "ggplot") && is.list(spatial.individual) &&
       length(spatial.individual) == 1L)

role.scatter <- netAnalysis_signalingRole_scatter(
  chat, signaling = pathway, do.label = FALSE, show.legend = FALSE
)
ggplot2::ggsave(file.path(out.dir, "netAnalysis-signalingRole-scatter.pdf"),
                role.scatter, width = 7, height = 5)
role.network.file <- file.path(out.dir, "netAnalysis-signalingRole-network.pdf")
grDevices::pdf(role.network.file, width = 7, height = 4)
tryCatch(
  netAnalysis_signalingRole_network(chat, signaling = pathway,
                                    measure = c("outdeg", "indeg"),
                                    measure.name = c("Sender", "Receiver")),
  finally = grDevices::dev.off()
)
role.heatmap <- netAnalysis_signalingRole_heatmap(chat, signaling = pathway,
                                                   pattern = "outgoing")
heatmap.file <- file.path(out.dir, "netAnalysis-signalingRole-heatmap.pdf")
grDevices::pdf(heatmap.file, width = 7, height = 5)
tryCatch(ComplexHeatmap::draw(role.heatmap), finally = grDevices::dev.off())
pass("signaling-role network, scatter, and heatmap consumers render",
     inherits(role.scatter, "ggplot") &&
       file.info(role.network.file)$size > 1000 &&
       inherits(role.heatmap, "Heatmap") && file.info(heatmap.file)$size > 1000)

summary.file <- "tests_dev/figures/summary-P5-remaining-migration.txt"
writeLines(c(results, "ALL P5 REMAINING MIGRATION CHECKS PASSED"), summary.file)
cat("All remaining migration checks passed.\n")
