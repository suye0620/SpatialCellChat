# ====== Test: SparseChatArray methods ======
# Run: Rscript tests_dev/test-SparseChatArray.R
#
# Tests: constructor, dim, object-level dimnames, print, [, [[,
#        marginSums, as.data.frame, unsupported arithmetic, t

source("renv/activate.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(Rcpp)
})
source("R/SpatialCellChat_class.R")

set.seed(42)
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

cat("=== SparseChatArray ===\n")

# Build test data
m1 <- rsparsematrix(5, 5, 0.3); m2 <- rsparsematrix(5, 5, 0.3)
rownames(m1) <- colnames(m1) <- letters[1:5]
rownames(m2) <- colnames(m2) <- letters[1:5]
arr <- SparseChatArray(list(LR1 = m1, LR2 = m2))

# Constructor
check("constructor returns SparseChatArray",
      inherits(arr, "SparseChatArray"), TRUE)

# dim
check("dim = c(5,5,2)", dim(arr), c(5L, 5L, 2L))

# dimnames
check("dimnames has 3 elements", length(dimnames(arr)), 3L)
check("row dimnames correct", dimnames(arr)[[1]], letters[1:5])
check("column dimnames correct", dimnames(arr)[[2]], letters[1:5])
check("layer names correct", names(arr), c("LR1", "LR2"))
check("internal layer dimnames are stripped",
      all(vapply(unclass(arr), function(m) {
        is.null(rownames(m)) && is.null(colnames(m))
      }, logical(1))),
      TRUE)
check("internal list names are stripped", is.null(names(unclass(arr))), TRUE)

# length
check("length = 2", length(arr), 2L)

# print (no error)
check("print produces output", { capture.output(print(arr)); TRUE }, TRUE)

# indexing: single layer
check("[,,1] is dgCMatrix", inherits(arr[,,1], "dgCMatrix"), TRUE)
check("[,,1] dim = c(5,5)", dim(arr[,,1]), c(5L, 5L))
check("[,,1] restores row dimnames", rownames(arr[,,1]), letters[1:5])
check("[[ by layer name returns dgCMatrix",
      inherits(arr[["LR1"]], "dgCMatrix"),
      TRUE)
check("[[ by layer name restores matrix dimnames",
      list(rownames(arr[["LR1"]]), colnames(arr[["LR1"]])),
      dimnames(arr)[1:2])
check("[[ numeric vector returns SparseChatArray",
      inherits(arr[[1:2]], "SparseChatArray"),
      TRUE)
check("[[ character vector returns SparseChatArray",
      inherits(arr[[c("LR1", "LR2")]], "SparseChatArray"),
      TRUE)
check("[[ rejects extra positional indices",
      inherits(tryCatch(arr[[1, 2, 1]], error = function(e) e), "error"),
      TRUE)
check("[[ rejects incomplete multidimensional subscripts",
      inherits(tryCatch(arr[[1, 2]], error = function(e) e), "error"),
      TRUE)

# indexing: multi-layer
check("[1:2,,] has 2 layers", length(arr[1:2,,]), 2L)
check("[1:2,,1] is dgCMatrix", inherits(arr[1:2,,1], "dgCMatrix"), TRUE)
arr_subset <- arr[c("a", "c"), "b", "LR2", drop = FALSE]
check("subset keeps row dimnames", dimnames(arr_subset)[[1]], c("a", "c"))
check("subset keeps column dimnames", dimnames(arr_subset)[[2]], "b")
check("subset keeps layer dimnames", dimnames(arr_subset)[[3]], "LR2")

# explicit object-level dimnames
m3 <- rsparsematrix(3, 4, 0.3); m4 <- rsparsematrix(3, 4, 0.3)
dn_explicit <- list(paste0("r", 1:3), paste0("c", 1:4), c("A", "B"))
arr_explicit <- SparseChatArray(list(m3, m4), dimnames = dn_explicit)
check("explicit dimnames stored", dimnames(arr_explicit), dn_explicit)
check("explicit layer names exposed", names(arr_explicit), c("A", "B"))
check("explicit dimnames are restored on extraction",
      list(rownames(arr_explicit[, , "A"]), colnames(arr_explicit[, , "A"])),
      dn_explicit[1:2])

# marginSums (cross-layer sum)
s3 <- marginSums(arr, 3)
check("marginSums(3) is dgCMatrix", inherits(s3, "dgCMatrix"), TRUE)
check("marginSums(3) dim", dim(s3), c(5L, 5L))
check("marginSums(3) preserves dimnames", rownames(s3), letters[1:5])

s1 <- marginSums(arr, 1)
check("marginSums(1) nrow = 5", nrow(s1), 5L)
check("marginSums(1) ncol = 2", ncol(s1), 2L)

s2 <- marginSums(arr, 2)
check("marginSums(2) nrow = 5", nrow(s2), 5L)
check("marginSums(2) ncol = 2", ncol(s2), 2L)

# as.data.frame
df <- as.data.frame(arr)
check("as.data.frame has 4 cols", ncol(df), 4L)
check("as.data.frame col names",
      sort(colnames(df)), sort(c("source", "target", "layer", "value")))

# unsupported arithmetic
check("unsupported subtraction fails",
      inherits(tryCatch(arr - arr, error = function(e) e), "error"),
      TRUE)
check("unsupported plus fails",
      inherits(tryCatch(arr + arr, error = function(e) e), "error"),
      TRUE)

# t
arr_t <- t(arr)
check("t() is SparseChatArray", inherits(arr_t, "SparseChatArray"), TRUE)
check("t() dim reversed", dim(arr_t), c(5L, 5L, 2L))
check("t() preserves names", names(arr_t), names(arr))

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)



Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

benchmark_config <- list(
  run = TRUE,
  n = 10000L,
  n_layers = 200L,
  sparsity = 0.9999,
  seed = 42L,
  iterations = 3L
)

make_binary_sparse_matrix <- function(n, density) {
  nnz <- max(1L, as.integer(floor(n * n * density)))
  positions <- sample.int(n * n, nnz, replace = FALSE)
  Matrix::sparseMatrix(
    i = (positions - 1L) %% n + 1L,
    j = (positions - 1L) %/% n + 1L,
    x = rep.int(1, nnz),
    dims = c(n, n),
    repr = "C"
  )
}

run_sparse_chat_benchmark <- function(
    arr,
    iterations = 1L
) {
  stopifnot(
    inherits(arr, "SparseChatArray"),
    length(iterations) == 1L, iterations >= 1
  )
  if (!requireNamespace("microbenchmark", quietly = TRUE)) {
    stop("Package 'microbenchmark' is required to run this benchmark; install it first",
         call. = FALSE)
  }

  d <- dim(arr)
  input_nnz_by_layer <- vapply(unclass(arr), Matrix::nnzero, integer(1))
  input_nnz <- sum(input_nnz_by_layer)
  density <- input_nnz / prod(d)
  cpp_result <- marginSums(arr, margin = 3L)
  cpp_benchmark <- microbenchmark::microbenchmark(
    marginSums_cpp = marginSums(arr, margin = 3L),
    times = iterations
  )

  result <- list(
    n = d[[1L]],
    n_layers = d[[3L]],
    sparsity = 1 - density,
    density = density,
    nnz_per_layer = input_nnz_by_layer,
    input_nnz = input_nnz,
    output_nnz = Matrix::nnzero(cpp_result),
    cpp_benchmark = cpp_benchmark,
    input_size_mb = as.numeric(object.size(arr)) / 1024^2,
    output_size_mb = as.numeric(object.size(cpp_result)) / 1024^2
  )

  gc()
  reduce_result <- Reduce(`+`, unclass(arr))
  reduce_benchmark <- microbenchmark::microbenchmark(
    reduce_r = Reduce(`+`, unclass(arr)),
    times = iterations
  )
  result$reduce_benchmark <- reduce_benchmark
  result$results_equal <- isTRUE(all.equal(cpp_result, reduce_result))

  invisible(result)
}


if (isTRUE(benchmark_config$run)) {
  density <- 1 - benchmark_config$sparsity
  set.seed(benchmark_config$seed)
  generation_time <- system.time({
    layers <- lapply(
      seq_len(benchmark_config$n_layers),
      function(...) make_binary_sparse_matrix(benchmark_config$n, density)
    )
  })
  names(layers) <- paste0("layer", seq_len(benchmark_config$n_layers))
  benchmark_arr <- SparseChatArray(layers)

  benchmark_result <- run_sparse_chat_benchmark(
    arr = benchmark_arr,
    iterations = benchmark_config$iterations
  )
  benchmark_result$generation_time <- generation_time
  print(benchmark_result)
}
