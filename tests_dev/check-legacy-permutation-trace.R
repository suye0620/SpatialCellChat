# Verify the 4cts trace and its controlled replay using the upstream RNG state.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tests_dev/check-legacy-permutation-trace.R <trace-dir>")
root <- args[[1L]]
upstream <- readRDS(file.path(root, "upstream-4cts-permutation-trace.rds"))
dev <- readRDS(file.path(root, "dev-4cts-permutation-trace.rds"))
matched <- readRDS(file.path(root, "matched-seed", "dev-4cts-permutation-trace.rds"))
u <- upstream$permutations
d <- dev$permutations
m <- matched$permutations
stopifnot(identical(u$axes, d$axes), identical(u$axes, m$axes),
          identical(as.character(u$group), as.character(d$group)),
          identical(d$group, m$group),
          identical(names(u$layers), names(d$layers)), identical(names(u$layers), names(m$layers)),
          u$min.percent == 0, u$min.cells.sr == 5, u$nboot == 100L,
          identical(u$min.percent, d$min.percent), identical(u$min.cells.sr, d$min.cells.sr),
          identical(d$seed.use, u$seed.use))

draw.differences <- colSums(u$permutation != d$permutation)
stopifnot(!identical(u$seed.before.drawing, d$seed.before.drawing),
          identical(which(draw.differences > 0L), 1L),
          identical(u$seed.before.drawing, m$seed.before.drawing),
          identical(u$permutation, m$permutation))
cat("Only draw 1 differs in the two original runs, with", draw.differences[[1L]],
    "of", nrow(u$permutation), "indices changed.\n")

observed.helper <- getFromNamespace("my_future_sapply", "SpatialCellChat")
future::plan(future::sequential)
set.seed(u$seed.use)
stopifnot(identical(.Random.seed, d$seed.before.drawing))
invisible(observed.helper(X = seq_along(u$axes[[3L]]), FUN = identity,
                          simplify = FALSE, strategy.message = FALSE))
stopifnot(identical(.Random.seed, u$seed.before.drawing))
set.seed(u$seed.use)
invisible(observed.helper(X = seq_along(u$axes[[3L]]), FUN = identity,
                          simplify = FALSE, strategy.message = FALSE, future.seed = FALSE))
stopifnot(identical(.Random.seed, d$seed.before.drawing))
cat("The upstream observed-aggregation helper reproduces the recorded RNG state.\n")

rows <- expand.grid(source = u$axes[[1L]], target = u$axes[[2L]], stringsAsFactors = FALSE)
flips <- lapply(names(u$layers), function(lr) {
  up.layer <- u$layers[[lr]]
  dev.layer <- d$layers[[lr]]
  matched.layer <- m$layers[[lr]]
  up.delta <- up.layer$Pboot - up.layer$Pnull
  dev.delta <- dev.layer$Pboot - dev.layer$Pnull
  matched.delta <- matched.layer$Pboot - matched.layer$Pnull
  stopifnot(identical(up.delta > 0, matched.delta > 0),
            identical(dev.layer$Pnull, matched.layer$Pnull))
  # Zero observed probabilities receive p=1 in either implementation.
  idx <- which((up.delta > 0) != (dev.delta > 0) & up.layer$Pnull > 0, arr.ind = TRUE)
  data.frame(LR = rep(lr, nrow(idx)), rows[idx[, "row"], , drop = FALSE], draw = idx[, "col"],
             upstream.delta = up.delta[idx], development.delta = dev.delta[idx])
})
flips <- do.call(rbind, flips)
stopifnot(nrow(flips) == 5L, all(flips$draw == 1L))
print(flips, row.names = FALSE)

pvalues <- function(run) unlist(lapply(run$result$group.pval, function(p) as.numeric(as.matrix(p))))
up.pval <- pvalues(upstream)
dev.pval <- pvalues(dev)
matched.pval <- pvalues(matched)
stopifnot(identical(up.pval, matched.pval), sum(up.pval != dev.pval) == nrow(flips),
          all((up.pval <= 0.05) == (dev.pval <= 0.05)))
cat("All stored group p-values match in the controlled replay.\n")
