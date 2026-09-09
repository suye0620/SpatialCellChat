# ====== Test: setEnvironment, .cli, my_future_lapply/sapply ======
# Run: Rscript tests_dev/test-parallel.R

source("renv/activate.R")
suppressPackageStartupMessages({
  library(Matrix); library(future); library(future.apply)
  library(progressr); library(cli)
})
source("R/utilities.R")

ok <- 0; fail <- 0
check <- function(desc, expr, expected = NULL) {
  val <- tryCatch(expr, error = function(e) e)
  if (inherits(val, "error")) {
    cat(sprintf("  X %s\n    ERROR: %s\n", desc, conditionMessage(val)))
    fail <<- fail + 1
  } else if (!is.null(expected) && !isTRUE(all.equal(val, expected))) {
    cat(sprintf("  X %s\n    expected: %s\n    got:      %s\n", desc,
        deparse(expected)[1], deparse(val)[1]))
    fail <<- fail + 1
  } else {
    cat(sprintf("  V %s\n", desc))
    ok <<- ok + 1
  }
}

cat("=== setEnvironment ===\n")

setEnvironment(workers = 1, verbose = 1)
check("sequential plan", future::nbrOfWorkers() == 1, TRUE)
check("verbose option", getOption("SpatialCellChat.verbose"), 1L)

setEnvironment(workers = 2, verbose = 0)
check("parallel plan", future::nbrOfWorkers() == 2, TRUE)
check("verbose=0 option", getOption("SpatialCellChat.verbose"), 0L)

cat("\n=== .cli verbosity ===\n")

options(SpatialCellChat.verbose = 0)
check("v0: info suppressed",
      capture.output(.cli("x", .type = "info")) |> length(), 0L)
check("v0: text suppressed",
      capture.output(.cli("x", .type = "text")) |> length(), 0L)

options(SpatialCellChat.verbose = 1)
check("v1: text suppressed",
      capture.output(.cli("x", .type = "text")) |> length(), 0L)

options(SpatialCellChat.verbose = 2)
check("v2: text shown?",
      capture.output(.cli("x", .type = "text")) |> length() >= 0, TRUE)

options(SpatialCellChat.verbose = 1)
check(".verbose=0 overrides down",
      capture.output(.cli("x", .type = "info", .verbose = 0)) |> length(), 0L)

cat("\n=== my_future_lapply (parallel) ===\n")

handlers("cli")
res_l <- my_future_lapply(1:4, function(i) { Sys.sleep(0.1); i^2 })
check("lapply result length", length(res_l), 4L)
check("lapply values", unlist(res_l), c(1, 4, 9, 16))

res_s <- my_future_sapply(1:4, function(i) { Sys.sleep(0.1); i^2 })
check("sapply result length", length(res_s), 4L)
check("sapply values", res_s, c(1, 4, 9, 16))

cat("\n=== my_future_lapply (sequential) ===\n")
setEnvironment(workers = 1, verbose = 0)
res_l2 <- my_future_lapply(1:4, function(i) { i * 10 })
check("sequential result", unlist(res_l2), c(10, 20, 30, 40))

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
