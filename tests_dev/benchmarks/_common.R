## tests_dev/benchmarks/_common.R —— 基准脚本公共引导
## 用法：各 bench 脚本头两行 = 可移植 setwd + source 本文件（约定从任意目录 Rscript 运行）
## 产物目录：SCC_BENCH_OUT 环境变量可覆盖，默认 tests_dev/benchmarks/out/（gitignored）

root <- getwd()
stopifnot(file.exists(file.path(root, "DESCRIPTION")),
          file.exists(file.path(root, "renv", "activate.R")))

if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
  if (requireNamespace("renv", quietly = TRUE)) renv::load(root) else source("renv/activate.R")
}
options(future.globals.maxSize = 3e10, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(purrr); library(Rcpp)
})
source("R/SpatialCellChat_class.R")
source("R/database.R")
suppressPackageStartupMessages(library(spatstat.sparse))
source("R/utilities.R")
source("R/spatial.R")
source("R/modeling.R")
Rcpp::sourceCpp("src/SpatialChat_Rcpp.cpp", rebuild = FALSE, showOutput = FALSE)

LOG <- character()
say <- function(...) { m <- sprintf(...); cat(m, "\n"); flush.console(); LOG <<- c(LOG, m) }
bench_out <- Sys.getenv("SCC_BENCH_OUT", unset = file.path(root, "tests_dev", "benchmarks", "out"))
if (!dir.exists(bench_out)) dir.create(bench_out, recursive = TRUE)
