# Group-network migration: threshold boundary, layer axes, and actual renderers.
setwd(local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(igraph)
  library(future.apply); library(ggplot2); library(colorspace); library(Rcpp); library(circlize)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
source("R/analysis.R")
source("R/visualization.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
source("tests_dev/fixtures/make_simulated_chat.R")
check <- function(label, value) {
  if (!isTRUE(value)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
chat <- make_simulated_chat()
labels <- levels(chat@idents)
lr <- data.frame(interaction_name = c("LR1", "LR2"), pathway_name = c("P", "P"),
                 ligand = c("L1", "L2"), receptor = c("R1", "R2"),
                 row.names = c("LR1", "LR2"))
mat <- function(i, j, x) as(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(4, 4)), "dgCMatrix")
chat@LR$LRsig <- lr
chat@net$group$prob <- SparseChatArray(list(mat(c(1, 2), c(2, 3), c(.4, .6)), mat(2, 1, .7)),
                                        dimnames = list(labels, labels, lr$interaction_name))
chat@net$group$pval <- SparseChatArray(list(mat(c(1, 2), c(2, 3), c(.05, .06)), mat(2, 1, .01)),
                                        dimnames = dimnames(chat@net$group$prob))
selected <- .sc_viz_group_lr(chat, "P", .05)
check("p equals threshold retained and greater removed",
      selected$prob[1, 2, 1] == .4 && selected$prob[2, 3, 1] == 0 && selected$prob[2, 1, 2] == .7)
single <- .sc_viz_group_lr(chat, "P", .05, "LR2")
check("single LR remains a three-dimensional array", identical(dim(single$prob), c(4L, 4L, 1L)))
without.p <- chat
without.p@net$group$pval <- NULL
check("missing p values leave probability unchanged", .sc_viz_group_lr(without.p, "P", .05)$prob[2, 3, 1] == .6)
check("group input stays unchanged", chat@net$group$prob[["LR1"]][2, 3] == .6)
check("heatmap LR threshold retains equality and removes larger p",
      .sc_viz_group_matrix(chat, "net", "LR1", thresh = .05)[1, 2] == .4 &&
        .sc_viz_group_matrix(chat, "net", "LR1", thresh = .05)[2, 3] == 0)
check("barplot LR layer keeps original probability",
      .sc_viz_group_matrix(chat, "net", "LR1")[2, 3] == .6)
chat@net$group$count <- Matrix::Matrix(matrix(c(0, 2, 1, 0, rep(0, 12)), 4, 4,
                                              dimnames = list(labels, labels)), sparse = TRUE)
check("aggregated count reads group branch",
      identical(.sc_viz_group_matrix(chat, "net", measure = "count"), as.matrix(chat@net$group$count)))
filtered <- subsetCommunication_internal(chat@net, lr[, c("ligand", "receptor", "pathway_name", "interaction_name")],
                                         labels, signaling = "P", thresh = .05)
check("pathway selection uses interaction names independent of metadata column order",
      nrow(filtered) == 1L && filtered$interaction_name[1] == "LR2" &&
        filtered$interaction_name_2[1] == "LR2")
check("default barplot constructs its data without source/target filters",
      inherits(netVisual_barplot(chat), "ggplot"))
edges <- .sc_viz_group_edges(chat, .05)
check("chord edges retain p equal to threshold and omit p above threshold",
      nrow(edges) == 2L && any(edges$interaction_name == "LR1" & edges$prob == .4) &&
        any(edges$interaction_name == "LR2" & edges$prob == .7))
pdf.file <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf.file)
aggregate <- netVisual_aggregate(chat, signaling = "P", layout = "circle")
individual <- netVisual_individual(chat, signaling = "P", pairLR.use = "LR1", layout = "circle")
chord <- netVisual_chord_cell(chat, signaling = "P", slot.name = "netP")
grDevices::dev.off()
check("circle aggregate and individual draw", inherits(aggregate, "recordedplot") &&
        length(individual) == 1L && file.info(pdf.file)$size > 0)
check("chord cell aggregated path draws", inherits(chord, "recordedplot"))
unlink(pdf.file)
old.dir <- getwd()
tryCatch({
  setwd(tempdir())
  netVisual(chat, signaling = "P", layout = "circle", out.format = "pdf")
  check("netVisual exports individual and aggregate PDFs",
        all(file.exists(c("P_circle_individual.pdf", "P_circle_aggregate.pdf"))) &&
          all(file.info(c("P_circle_individual.pdf", "P_circle_aggregate.pdf"))$size > 1000))
  netVisual(chat, signaling = "P", layout = "spatial", out.format = c("pdf", "svg", "png"), height = 3)
  spatial.files <- as.vector(outer(c("P_spatial_individual.", "P_spatial_aggregate."),
                                   c("pdf", "svg", "png"), paste0))
  check("netVisual spatial exports both plots in PDF, SVG and PNG",
        all(file.exists(spatial.files)) && all(file.info(spatial.files)$size > 1000))
}, finally = setwd(old.dir))
cat("All group-network migration checks passed.\n")
