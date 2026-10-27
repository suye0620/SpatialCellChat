# projectData(): PPI projection onto assay$smooth on the final 11-slot schema.
# Covers the slot mapping, the solve/iterative equivalence, the validation
# branches added with the migration, and the subsetData invalidation rule.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(pbapply)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
source("R/utilities.R")
source("R/spatial.R")

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
err_msg <- function(expr) {
  tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
}
warn_msg <- function(expr) {
  msgs <- character()
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  paste(msgs, collapse = " | ")
}

## ---- fixture: 6 cells, 2 groups, 5 genes -------------------------------------
cells <- paste0("cell", seq_len(6L))
genes <- c("L1", "R1", "S1", "S2", "S3")
expr <- matrix(1, nrow = length(genes), ncol = length(cells),
               dimnames = list(genes, cells))
expr["L1", 1:3] <- 10
expr["R1", 4:6] <- 10
meta <- data.frame(label = factor(c("A", "A", "A", "B", "B", "B")), row.names = cells)
coords <- cbind(x = seq_len(6L), y = rep(c(1, 2), 3L))
rownames(coords) <- cells

chat <- createSpatialCellChat(
  expr, meta = meta, group.by = "label", input.assay = "norm",
  datatype = "spatial", coordinates = coords,
  spatial.factors = list(ratio = 1, tol = 1)
)
chat@DB <- list(
  interaction = data.frame(ligand = "L1", receptor = "R1",
                           annotation = "Secreted Signaling",
                           interaction_name = "L1-R1", stringsAsFactors = FALSE),
  complex = data.frame(), cofactor = data.frame()
)
chat <- subsetData(chat, features = c("L1", "R1", "S1"))

## ---- adjacency: 3-gene chain -------------------------------------------------
ppi <- matrix(0, 3L, 3L, dimnames = list(c("L1", "R1", "S1"), c("L1", "R1", "S1")))
ppi["L1", "R1"] <- ppi["R1", "L1"] <- 1
ppi["R1", "S1"] <- ppi["S1", "R1"] <- 1

## ---- 1. slot mapping ---------------------------------------------------------
chat.solve <- projectData(chat, ppi, alpha = 0.5, solver = "solve")
smooth <- assay(chat.solve, "smooth")
check("smooth keeps the signaling gene space",
      identical(rownames(smooth), rownames(assay(chat, "signaling"))))
check("smooth carries the object cells as columns", identical(colnames(smooth), cells))
check("object validates after the write", isTRUE(validateSpatialCellChat(chat.solve)))
check("solver recorded in misc$.param$projection",
      identical(chat.solve@misc$.param$projection$solver, "solve"))
check("projection logged",
      any(vapply(chat.solve@misc$.log,
                 function(e) identical(e[["function"]], "projectData"), logical(1L))))
check("signaling left untouched",
      identical(assay(chat.solve, "signaling"), assay(chat, "signaling")))

## ---- 2. numerical contract ---------------------------------------------------
chat.iter <- projectData(chat, ppi, alpha = 0.5, solver = "iterative",
                         tol = 1e-10, max.iter = 2000L)
smooth.iter <- assay(chat.iter, "smooth")
check("iterative and solve agree to 1e-8", max(abs(smooth - smooth.iter)) < 1e-8)

anorm <- ppi / rowSums(ppi)
expected <- as.matrix(assay(chat, "signaling"))[rownames(ppi), , drop = FALSE]
manual <- solve(diag(nrow(ppi)) - 0.5 * anorm, 0.5 * expected)
check("solve matches the closed-form solution",
      isTRUE(all.equal(unname(smooth), unname(manual), tolerance = 1e-12)))

## ---- 3. invalidation ---------------------------------------------------------
chat.re <- subsetData(chat.solve, features = c("L1", "R1"))
check("subsetData clears a stale smooth layer", is.null(assay(chat.re, "smooth")))

## ---- 4. rejection branches ---------------------------------------------------
chat.nosig <- chat
assay(chat.nosig, "signaling") <- NULL
check("empty signaling rejected",
      grepl("subsetData", err_msg(projectData(chat.nosig, ppi))))
bad <- ppi[1:2, , drop = FALSE]
check("non-square adjacency rejected", grepl("square", err_msg(projectData(chat, bad))))
bad <- ppi; rownames(bad) <- c("L1", "S1", "R1")
check("permuted rownames/colnames rejected", grepl("identical", err_msg(projectData(chat, bad))))
bad <- ppi; bad["S1", ] <- 0; bad[, "S1"] <- 0
check("zero row/column rejected", grepl("zero", err_msg(projectData(chat, bad))))
bad <- ppi; dimnames(bad) <- list(c("X1", "X2", "X3"), c("X1", "X2", "X3"))
check("zero gene overlap rejected", grepl("overlap", err_msg(projectData(chat, bad))))
check("alpha = 0 rejected", grepl("\\(0, 1\\)", err_msg(projectData(chat, ppi, alpha = 0))))
check("alpha = 'auto' rejected",
      grepl("\\(0, 1\\)", err_msg(projectData(chat, ppi, alpha = "auto"))))
check("alpha >= 1 rejected", grepl("\\(0, 1\\)", err_msg(projectData(chat, ppi, alpha = 1))))

## ---- 5. low-overlap warning --------------------------------------------------
net <- c("L1", paste0("S", seq_len(39L)))
big <- matrix(0, 40L, 40L, dimnames = list(net, net))
for (i in seq_len(39L)) big[i, i + 1L] <- big[i + 1L, i] <- 1
check("low overlap warns",
      grepl("network genes", warn_msg(projectData(chat, big, alpha = 0.5))))

## ---- 6. non-convergence is an error, not a silent truncation -----------------
check("max.iter exhaustion raises an error",
      grepl("did not converge",
            err_msg(projectData(chat, ppi, alpha = 0.99, solver = "iterative",
                                tol = 1e-14, max.iter = 1L))))

cat("All projectData checks passed.\n")
