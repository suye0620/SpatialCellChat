# Compare old sf row-subsetting against direct vector indexing under the same inputs.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages(library(sf))
source("R/visualization.R")
set.seed(4921)
n <- 2000L
field <- st_as_sf(data.frame(x = seq_len(n), y = seq_len(n),
                             dx = rnorm(n), dy = rnorm(n)), coords = c("x", "y"))
groups <- lapply(seq_len(6000L), function(i) {
  size <- sample.int(4L, 1L) - 1L
  if (size) sample.int(n, size) else integer()
})
results <- lapply(c("median", "sum"), function(name) {
  fn <- if (name == "median") stats::median else base::sum
  old.time <- system.time({
    old <- list(u = vapply(groups, function(i) fn(field[i, , drop = FALSE]$dx, na.rm = TRUE), numeric(1)),
                v = vapply(groups, function(i) fn(field[i, , drop = FALSE]$dy, na.rm = TRUE), numeric(1)))
  })[["elapsed"]]
  new.time <- system.time(actual <- .sc_grid_field_components(field, groups, fn))[["elapsed"]]
  if (!identical(old, actual)) stop("aggregation results differ for ", name)
  data.frame(statistic = name, points = n, grid_groups = length(groups),
             empty_groups = sum(lengths(groups) == 0L),
             old_seconds = old.time, new_seconds = new.time, equivalent = TRUE)
})
summary <- do.call(rbind, results)
write.csv(summary, "tests_dev/figures/P3/grid-aggregation-benchmark-2026-11-10.csv", row.names = FALSE)
print(summary)
