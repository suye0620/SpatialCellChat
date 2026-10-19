# ====== Test: CLI conditions and sequential/parallel progress ======
# Run: Rscript tests_dev/test-parallel.R

Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(progressr)
  library(cli)
})
source("R/utilities.R")

local({
  option_names <- c("SpatialCellChat.verbose", "future.globals.maxSize",
                    "progressr.enable", "progressr.handlers")
  old_options <- setNames(lapply(option_names, getOption), option_names)
  old_plan <- future::plan("list")
  old_global <- progressr::handlers(global = NA)
  old_rng_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    progressr::handlers(global = FALSE)
    future::plan(old_plan)
    options(old_options)
    progressr::handlers(global = old_global)
    do.call(RNGkind, as.list(old_rng_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  progressr::handlers(global = FALSE)

  passed <- 0L
  failed <- character()
  check <- function(label, expr) {
    result <- tryCatch(force(expr), error = identity)
    if (isTRUE(result)) {
      passed <<- passed + 1L
    } else {
      detail <- if (inherits(result, "error")) conditionMessage(result) else "assertion failed"
      failed <<- c(failed, paste(label, detail, sep = ": "))
      cat("FAIL:", label, "-", detail, "\n")
    }
  }
  messages <- function(expr) capture.output(force(expr), type = "message")
  contains <- function(output, token) any(grepl(token, output, fixed = TRUE))

  # Observe rendered stderr/message output, not stdout or ANSI snapshots.
  options(SpatialCellChat.verbose = 0L)
  for (type in c("info", "success", "bullet", "header", "subheader", "text", "debug")) {
    check(paste("global quiet suppresses", type),
          identical(messages(.cli("quiet-token", .type = type)), character()))
  }
  for (type in c("warning", "danger")) {
    check(paste("quiet retains visual", type),
          contains(messages(.cli("visual-token", .type = type)), "visual-token"))
  }
  check("explicit zero is a threshold override, not a quiet flag",
        contains(messages(.cli("override-token", .verbose = 0L)), "override-token"))

  options(SpatialCellChat.verbose = 1L)
  value <- "interpolation-token"
  check("normal messages interpolate the caller environment",
        contains(messages(.cli("{value}")), value))
  check("normal mode excludes details",
        identical(messages(.cli("detail-token", .type = "text")), character()))
  check("an explicit lower threshold admits details",
        contains(messages(.cli("detail-token", .type = "text", .verbose = 1L)), "detail-token"))
  check("an explicit higher threshold suppresses ordinary messages",
        identical(messages(.cli("info-token", .verbose = 2L)), character()))
  options(SpatialCellChat.verbose = 2L)
  check("detailed mode emits details",
        contains(messages(.cli("detail-token", .type = "text")), "detail-token"))
  check("detailed mode excludes debug",
        identical(messages(.cli("debug-token", .type = "debug")), character()))
  options(SpatialCellChat.verbose = 3L)
  check("debug threshold emits its message",
        contains(messages(.cli("debug-token", .type = "debug")), "debug-token"))

  # Install a user-selected reporter, but no global handler and no outer
  # with_progress(). The wrappers must supply their own finite lifecycle.
  observe_progress <- function(expr, verbose = 1L, enable = TRUE) {
    saved <- options(SpatialCellChat.verbose = verbose, progressr.enable = enable)
    saved_handlers <- getOption("progressr.handlers")
    on.exit({
      options(saved)
      options(progressr.handlers = saved_handlers)
    }, add = TRUE)
    started <- updated <- finished <- numeric()
    recorder <- function(...) {
      progressr::make_progression_handler(
        "test-recorder",
        reporter = list(
          initiate = function(config, state, progression, ...) {
            started <<- c(started, config$max_steps)
          },
          update = function(config, state, progression, ...) {
            updated <<- c(updated, state$step)
          },
          finish = function(config, state, progression, ...) {
            finished <<- c(finished, state$step)
          }
        ),
        interval = 0, times = Inf, intrusiveness = 1
      )
    }
    progressr::handlers(recorder)
    result <- force(expr)
    list(value = result, started = started, updated = updated, finished = finished)
  }
  completed <- function(observed, n) {
    identical(as.integer(observed$started), as.integer(n)) &&
      identical(as.integer(observed$updated), seq_len(n)) &&
      identical(as.integer(observed$finished), as.integer(n))
  }
  silent <- function(observed) {
    identical(observed$started, numeric()) &&
      identical(observed$updated, numeric()) &&
      identical(observed$finished, numeric())
  }

  policies <- list(
    normal = list(verbose = 1L, progress = TRUE, enable = TRUE, visible = TRUE),
    detailed = list(verbose = 2L, progress = TRUE, enable = TRUE, visible = TRUE),
    local_quiet = list(verbose = 2L, progress = FALSE, enable = TRUE, visible = FALSE),
    global_quiet = list(verbose = 0L, progress = TRUE, enable = TRUE, visible = FALSE),
    batch_opt_out = list(verbose = 2L, progress = TRUE, enable = FALSE, visible = FALSE),
    session_policy = list(verbose = 1L, progress = TRUE, enable = NULL, visible = interactive())
  )
  inputs <- c(first = 2L, second = 4L, third = 6L, fourth = 8L)
  # A named argument called x must not collide with a progress wrapper's
  # iteration argument; unnamed arguments and ... must also reach FUN.
  calculate <- function(item, multiplier, x, ...) item * multiplier + x + sum(...)
  random_value <- function(item, shift) c(item = item, runif(1), rnorm(1)) + shift
  warning_value <- function(item) {
    if (item == 2L) {
      warning(structure(list(message = "warning-token", call = NULL, code = 41L),
                        class = c("test_apply_warning", "warning", "condition")))
    }
    item
  }
  error_value <- function(item) {
    stop(structure(list(message = "error-token", call = NULL, code = 42L),
                   class = c("test_apply_error", "error", "condition")))
  }

  for (workers in c(1L, 2L)) {
    strategy <- if (workers == 1L) "sequential" else "multisession"
    setup_messages <- messages(setEnvironment(workers = workers, verbose = 0L))
    check(paste(strategy, "environment setup honors quiet"),
          identical(setup_messages, character()))
    progressr::handlers(global = FALSE)

    for (kind in c("lapply", "sapply")) {
      apply_fn <- if (kind == "lapply") my_future_lapply else my_future_sapply
      expected <- lapply(inputs, calculate, 3L, x = 5L, extra = 7L)
      if (kind == "sapply") expected <- simplify2array(expected)

      for (policy_name in names(policies)) {
        policy <- policies[[policy_name]]
        label <- paste(strategy, kind, policy_name)
        check(paste(label, "values, names, arguments and completion"), {
          observed <- observe_progress(
            apply_fn(inputs, calculate, 3L, x = 5L, extra = 7L,
                     future.seed = TRUE, future.scheduling = Inf,
                     .progress = policy$progress),
            verbose = policy$verbose, enable = policy$enable
          )
          identical(observed$value, expected) &&
            if (policy$visible) completed(observed, length(inputs)) else silent(observed)
        })
      }

      for (small in list(integer(), 5L)) {
        check(paste(strategy, kind, "empty/single iteration does not report progress", length(small)), {
          observed <- observe_progress(apply_fn(small, identity))
          reference <- lapply(small, identity)
          if (kind == "sapply" && length(reference)) reference <- simplify2array(reference)
          identical(observed$value, reference) && silent(observed)
        })
      }

      for (seed in list(TRUE, 619L)) {
        set.seed(1234)
        expected_rng <- future.apply::future_lapply(inputs, random_value, shift = 9,
                                                   future.seed = seed, future.scheduling = 1)
        expected_state <- .Random.seed
        if (kind == "sapply") expected_rng <- simplify2array(expected_rng)
        for (policy_name in c("normal", "detailed", "local_quiet", "global_quiet", "batch_opt_out")) {
          policy <- policies[[policy_name]]
          check(paste(strategy, kind, "seed", seed, policy_name, "preserves future RNG"), {
            set.seed(1234)
            observed <- observe_progress(
              apply_fn(inputs, random_value, shift = 9,
                       future.seed = seed, future.scheduling = 1,
                       .progress = policy$progress),
              verbose = policy$verbose, enable = policy$enable
            )
            identical(observed$value, expected_rng) && identical(.Random.seed, expected_state)
          })
        }
      }

      for (policy_name in c("normal", "local_quiet", "global_quiet", "batch_opt_out")) {
        policy <- policies[[policy_name]]
        label <- paste(strategy, kind, policy_name)
        check(paste(label, "preserves actual warning conditions and values"), {
          warnings <- list()
          observed <- withCallingHandlers(
            observe_progress(
              apply_fn(inputs, warning_value, .progress = policy$progress),
              verbose = policy$verbose, enable = policy$enable
            ),
            warning = function(condition) {
              warnings[[length(warnings) + 1L]] <<- condition
              invokeRestart("muffleWarning")
            }
          )
          reference <- if (kind == "lapply") as.list(inputs) else inputs
          identical(observed$value, reference) && length(warnings) == 1L &&
            inherits(warnings[[1L]], "test_apply_warning") && identical(warnings[[1L]]$code, 41L)
        })
        check(paste(label, "preserves actual error conditions"), {
          condition <- tryCatch(
            observe_progress(
              apply_fn(inputs, error_value, .progress = policy$progress),
              verbose = policy$verbose, enable = policy$enable
            ),
            error = identity
          )
          inherits(condition, "test_apply_error") && identical(condition$code, 42L)
        })
      }
    }

    for (progress in c(FALSE, TRUE)) {
      for (use_names in c(FALSE, TRUE)) {
        check(paste(strategy, "sapply character names", use_names, "progress", progress), {
          observed <- observe_progress(my_future_sapply(c("a", "bbb"), nchar,
                                                       USE.NAMES = use_names, .progress = progress))
          reference <- future.apply::future_sapply(c("a", "bbb"), nchar,
                                                   USE.NAMES = use_names, future.seed = TRUE)
          identical(observed$value, reference)
        })
      }
      matrix_value <- function(item) matrix(item + 1:4, nrow = 2L,
                                           dimnames = list(c("a", "b"), c("c", "d")))
      for (simplify in list(FALSE, TRUE, "array")) {
        check(paste(strategy, "sapply simplify", simplify, "progress", progress), {
          observed <- observe_progress(my_future_sapply(inputs, matrix_value,
                                                       simplify = simplify, .progress = progress))
          reference <- future.apply::future_sapply(inputs, matrix_value,
                                                   simplify = simplify, future.seed = TRUE)
          identical(observed$value, reference)
        })
      }
      check(paste(strategy, "lapply unlist semantics with progress", progress), {
        observed <- observe_progress(my_future_lapply(inputs, matrix_value,
                                                     simplify = TRUE, .progress = progress))
        reference <- unlist(future.apply::future_lapply(inputs, matrix_value, future.seed = TRUE))
        identical(observed$value, reference)
      })
    }
  }

  cat(sprintf("\n%d passed, %d failed\n", passed, length(failed)))
  if (length(failed)) stop(paste(failed, collapse = "\n"), call. = FALSE)
})
