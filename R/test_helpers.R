# ====== Testing infrastructure for SpatialCellChat ======
#
# Usage:
#   1. source("R/test_helpers.R")                  — load helpers
#   2. rs <- start_session()                        — start persistent R process
#   3. source_chat(rs)                              — source all R files into session
#   4. init_test_data(rs)                           — build mini object once
#   5. rs$run(function() {                          — test any function
#       chat <- readRDS(".test_cache/chat.rds")
#       # ... test something ...
#       "done"
#     })
#   6. source_chat(rs)                              — after editing R files, reload
#   7. rs$close()                                   — done, kill the process

#' Start a persistent background R session for iterative testing
start_session <- function() {
  rs <- callr::r_session$new(wait_timeout = 300000)
  rs$run(function() { "alive" })
  cat("[test session] started, pid:", rs$get_pid(), "\n")
  rs
}

#' Source all R/ files into the session (replaces devtools::load_all)
#'
#' Pre-loads key dependencies (Matrix, methods, etc.) then sources all R/ files.
#' Use this instead of devtools::load_all() to avoid dependency checks.
#' Call after any edit to R source files to reload.
source_chat <- function(rs) {
  loaded <- rs$run(function() {
    # preload essential packages needed by class definitions
    library(methods, quietly = TRUE)
    library(Matrix, quietly = TRUE)
    library(Seurat, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    library(igraph, quietly = TRUE)
    library(cli, quietly = TRUE)

    files <- sort(list.files("R", pattern = "\\.R$", full.names = TRUE))
    files <- setdiff(files, file.path("R", "test_helpers.R"))
    results <- character()
    for (f in files) {
      res <- tryCatch({source(f, local = FALSE); "OK"}, error = function(e) paste("FAIL:", conditionMessage(e)))
      results <- c(results, paste0(basename(f), ": ", res))
    }
    results
  })
  cat("[test session] sourcing results:\n")
  for (r in loaded) cat("  ", r, "\n")
}

#' Build a tiny test dataset and create a minimal SpatialCellChat object
#'
#' Creates: 50 cells, 100 genes, 5 cell groups, 10 LR pairs
#' Saves:   .test_cache/chat.rds
init_test_data <- function(rs) {
  rs$run(function() {
    set.seed(42)
    nC <- 50; nG <- 100; nGroups <- 5

    expr <- Matrix::rsparsematrix(nG, nC, density = 0.3)
    rownames(expr) <- paste0("Gene", 1:nG)
    colnames(expr) <- paste0("cell_", 1:nC)
    expr <- as(expr, "dgCMatrix")

    groups <- sample(paste0("Type", 1:nGroups), nC, replace = TRUE)
    meta <- data.frame(labels = groups, row.names = colnames(expr))
    coords <- expand.grid(x = 1:10, y = 1:5)[1:nC, ]

    DB <- list(
      interaction = data.frame(
        interaction_name = paste0("LR", 1:10),
        ligand     = c("Gene1", "Gene1", "Gene2", "Gene3", "Gene4",
                      "Gene5", "Gene6", "Gene7", "Gene2", "Gene8"),
        receptor   = c("Gene9", "Gene10","Gene10","Gene9", "Gene8",
                      "Gene7", "Gene6", "Gene5", "Gene4", "Gene3"),
        annotation = rep("Secreted Signaling", 10),
        pathway_name = c(rep("PathwayA", 3), rep("PathwayB", 3), rep("PathwayC", 4)),
        stringsAsFactors = FALSE
      ),
      complex   = data.frame(),
      cofactor  = data.frame(),
      geneInfo  = data.frame(Symbol = paste0("Gene", 1:10))
    )

    LR <- list(LRsig = DB$interaction)

    chat <- createSpatialCellChat(
      expr, meta = meta, group.by = "labels",
      datatype = "spatial",
      coordinates = as.matrix(coords),
      spatial.factors = list(ratio = 0.18, tol = 5)
    )
    chat@DB <- DB
    chat@LR <- LR

    dir.create(".test_cache", showWarnings = FALSE)
    saveRDS(chat, ".test_cache/chat.rds")
    cat("[test data] 50 cells, 100 genes, 5 groups, 10 LR pairs\n")
  })
}

#' Show available checkpoints
list_checkpoints <- function() {
  fs <- list.files(".test_cache", pattern = "\\.rds$")
  cat("Checkpoints:\n")
  for (f in fs) {
    info <- file.info(file.path(".test_cache", f))
    cat(sprintf("  %-40s %.1f KB\n", f, info$size / 1024))
  }
}
