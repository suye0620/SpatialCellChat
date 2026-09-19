# ====== Test: cpp_sum_layers (Rcpp cross-layer sum) ======
# Run: Rscript tests_dev/test-cpp_sum_layers.R
#
# Tests: correctness vs R Reduce, edge cases (1 layer, empty?),
#        dimnames preservation

Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({ library(Matrix); library(Rcpp) })

cat("Compiling Rcpp...\n")
sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)
cat("OK\n\n")

ok <- 0; fail <- 0
check <- function(desc, expr, expected = NULL) {
  val <- tryCatch(expr, error = function(e) e)
  if (inherits(val, "error")) {
    cat(sprintf("  \u2716 %s\n    ERROR: %s\n", desc, conditionMessage(val)))
    fail <<- fail + 1
  } else if (!is.null(expected) && !isTRUE(all.equal(val, expected))) {
    cat(sprintf("  \u2716 %s\n    expected: %s\n    got:      %s\n", desc,
        deparse(expected)[1], deparse(val)[1]))
    fail <<- fail + 1
  } else {
    cat(sprintf("  \u2713 %s\n", desc))
    ok <<- ok + 1
  }
}

cat("=== cpp_sum_layers ===\n")

set.seed(42)
make_mat <- function(nr, nc, density = 0.3) {
  m <- rsparsematrix(nr, nc, density)
  rownames(m) <- paste0("R", seq_len(nr))
  colnames(m) <- paste0("C", seq_len(nc))
  m
}

# Basic: 2 layers
m1 <- make_mat(5, 5); m2 <- make_mat(5, 5)
res2 <- cpp_sum_layers(list(m1, m2))
ref2 <- Reduce(`+`, list(m1, m2))
check("2 layers: class dgCMatrix", inherits(res2, "dgCMatrix"), TRUE)
check("2 layers: dim", dim(res2), c(5L, 5L))
check("2 layers: matches R Reduce", res2, ref2)
check("2 layers: dimnames preserved",
      dimnames(res2), dimnames(m1))

# Single layer
res1 <- cpp_sum_layers(list(m1))
check("1 layer: dim", dim(res1), c(5L, 5L))
check("1 layer: matches source", res1, m1)
check("1 layer: dimnames", dimnames(res1), dimnames(m1))

# 3 layers
m3 <- make_mat(5, 5)
res3 <- cpp_sum_layers(list(m1, m2, m3))
ref3 <- Reduce(`+`, list(m1, m2, m3))
check("3 layers: matches R Reduce", res3, ref3)

# Larger: 100x100, 10 layers
set.seed(1)
layers <- lapply(1:10, function(...) make_mat(100, 100, 0.05))
res_large <- cpp_sum_layers(layers)
ref_large <- Reduce(`+`, layers)
check("100x100x10: matches R Reduce", res_large, ref_large)

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
