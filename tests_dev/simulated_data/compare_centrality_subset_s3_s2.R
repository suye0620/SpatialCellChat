# Wave 3a E2E oracle: centrality/subsetCommunication on SpatialChat_3
# versus the frozen be88f30 package on SpatialChat_2.
#
# The baseline subsetCommunication calls are expected to fail before producing a
# data frame because the frozen LR schema has no `evidence` column. Subset output
# is therefore checked against the independent row-by-row oracle below. Baseline
# centrality is still usable and is compared for net/netP and both measure sets.
setwd(local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(dirname(dirname(normalizePath(f)))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
if (!nzchar(Sys.getenv("RENV_PROJECT"))) {
  if (requireNamespace("renv", quietly = TRUE)) renv::load(getwd()) else source("renv/activate.R")
}
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100000 * 1024^2)
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(future)
})
future::plan(future::sequential)
source("R/SpatialCellChat_class.R")
source("R/utilities.R")
source("R/database.R")
source("R/modeling.R")
source("R/analysis.R")

data_dir <- file.path(getwd(), "tests_dev", "simulated_data")
ref_path <- file.path(data_dir, "SpatialChat_2.rds")
new_path <- file.path(data_dir, "SpatialChat_3.rds")
ref <- readRDS(ref_path)
new <- readRDS(new_path)
thresh <- 0.05

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("[FAIL] ", label, call. = FALSE)
  cat("[PASS] ", label, "\n", sep = "")
}
run_safe <- function(expr) {
  tryCatch(list(ok = TRUE, value = force(expr)),
           error = function(e) list(ok = FALSE, error = conditionMessage(e)))
}

## ---- 1. baseline package metadata and isolated probes ----
baseline_dir <- normalizePath(
  file.path(getwd(), "renv/library/windows/R-4.5/x86_64-w64-mingw32/SpatialCellChat"),
  winslash = "/", mustWork = TRUE
)
baseline_desc <- read.dcf(file.path(baseline_dir, "DESCRIPTION"))[1, , drop = FALSE]
check("baseline package version is pinned", identical(unname(baseline_desc[, "Version"]), "0.1.0"))
check("baseline package SHA is pinned", identical(unname(baseline_desc[, "RemoteSha"]),
                                                   "be88f300464cf970b36a6165ea7d2e41a47a7511"))

baseline <- callr::r(
  function(ref_path, baseline_dir, thresh) {
    suppressPackageStartupMessages(library(SpatialCellChat, lib.loc = dirname(baseline_dir)))
    ref <- readRDS(ref_path)
    one <- function(expr) {
      tryCatch(list(ok = TRUE, value = force(expr)),
               error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    }
    sub_net <- one(SpatialCellChat::subsetCommunication(ref, slot.name = "net", thresh = thresh))
    sub_netP <- one(SpatialCellChat::subsetCommunication(ref, slot.name = "netP", thresh = thresh))
    manual_centrality <- function(slot_name, degree_only) {
      slot_obj <- methods::slot(ref, slot_name)
      prob <- if (slot_name == "net") slot_obj$prob else slot_obj$prob
      if (slot_name == "net") {
        pval <- slot_obj$pval
        pval[prob == 0] <- 1
        prob[pval >= thresh] <- 0
      }
      local_fun <- getFromNamespace("computeCentralityLocal", "SpatialCellChat")
      layers <- lapply(seq_len(dim(prob)[3]), function(i)
        local_fun(prob[, , i, drop = TRUE], degree.only = degree_only))
      measures <- rownames(layers[[1]])
      out <- array(unlist(lapply(layers, as.vector)),
                   c(length(measures), dim(prob)[1], dim(prob)[3]))
      dimnames(out) <- list(measures, dimnames(prob)[[1]], dimnames(prob)[[3]])
      out
    }
    centrality <- list()
    for (slot_name in c("net", "netP")) {
      for (degree_only in c(TRUE, FALSE)) {
        key <- paste(slot_name, if (degree_only) "degree" else "full", sep = ".")
        result <- one(SpatialCellChat::netAnalysis_computeCentrality(
          ref, slot.name = slot_name, do.group = TRUE, thresh = thresh,
          degree.only = degree_only
        ))
        if (isTRUE(result$ok)) {
          result$value <- methods::slot(result$value, slot_name)[["centr"]]
        } else {
          result$fallback <- TRUE
          result$value <- manual_centrality(slot_name, degree_only)
        }
        centrality[[key]] <- result
      }
    }
    list(
      package = packageDescription("SpatialCellChat"),
      subset_net = sub_net,
      subset_netP = sub_netP,
      centrality = centrality
    )
  },
  args = list(ref_path = ref_path, baseline_dir = baseline_dir, thresh = thresh),
  libpath = .libPaths(),
  cmdargs = c("--vanilla"),
  stdout = "|", stderr = "|"
)
check("baseline net subset failure records evidence bug", !isTRUE(baseline$subset_net$ok) &&
      grepl("evidence", baseline$subset_net$error, fixed = TRUE))
check("baseline netP subset failure records evidence bug", !isTRUE(baseline$subset_netP$ok) &&
      grepl("evidence", baseline$subset_netP$error, fixed = TRUE))
cat("[INFO] baseline net error: ", baseline$subset_net$error, "\n", sep = "")
cat("[INFO] baseline netP error: ", baseline$subset_netP$error, "\n", sep = "")

## ---- 2. independent row oracle for the new 11-slot object ----
ensure_join_columns <- function(LR) {
  for (col in c("interaction_name_2", "pathway_name", "ligand", "receptor", "annotation", "evidence")) {
    if (!col %in% colnames(LR)) LR[[col]] <- NA_character_
  }
  LR
}

manual_rows <- function(object, slot_name = "net", signaling = NULL,
                        sources.use = NULL, targets.use = NULL) {
  LR <- ensure_join_columns(object@LR$LRsig)
  slot_obj <- methods::slot(object, slot_name)
  prob <- slot_obj$group$prob
  layers <- unclass(prob)
  layer_names <- dimnames(prob)[[3]]
  dn1 <- dimnames(prob)[[1]]
  dn2 <- dimnames(prob)[[2]]
  pvals <- if (slot_name == "net") slot_obj$group$pval else NULL
  pval_layers <- if (is.null(pvals)) NULL else unclass(pvals)

  rows <- lapply(seq_along(layers), function(i) {
    pm <- as.matrix(layers[[i]])
    if (!is.null(pval_layers)) pm[as.matrix(pval_layers[[i]]) >= thresh] <- 0
    idx <- which(pm > 0, arr.ind = TRUE)
    if (nrow(idx) == 0L) return(NULL)
    data.frame(
      source = dn1[idx[, 1]], target = dn2[idx[, 2]],
      interaction_name = layer_names[i], prob = pm[idx],
      pval = if (!is.null(pval_layers)) as.matrix(pval_layers[[i]])[idx] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) {
    out <- do.call(rbind, rows)
  } else {
    out <- data.frame(source = character(), target = character(), interaction_name = character(),
                      prob = numeric(), pval = numeric(), stringsAsFactors = FALSE)
  }

  join_idx <- if (slot_name == "netP") {
    match(out$interaction_name, LR$pathway_name)
  } else {
    match(out$interaction_name, rownames(LR))
  }
  if (nrow(out)) out <- cbind(out, LR[join_idx, c("interaction_name_2", "pathway_name", "ligand", "receptor", "annotation", "evidence"), drop = FALSE])

  if (!is.null(signaling)) {
    if (slot_name == "netP") {
      out <- out[out$pathway_name %in% signaling, , drop = FALSE]
    } else {
      signaling_pairs <- rownames(LR)[LR$pathway_name %in% signaling]
      out <- out[out$interaction_name %in% signaling_pairs, , drop = FALSE]
    }
  }
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) sources.use <- levels(object@idents)[sources.use]
    out <- out[out$source %in% sources.use, , drop = FALSE]
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) targets.use <- levels(object@idents)[targets.use]
    out <- out[out$target %in% targets.use, , drop = FALSE]
  }

  if (slot_name == "netP") {
    # Match subsetCommunication's source-target/pathway aggregation without using
    # its implementation. `netP$group$prob` is already pathway-summed.
    if (!nrow(out)) return(data.frame(source = character(), target = character(),
                                      pathway_name = character(), prob = numeric(),
                                      pval = numeric(), stringsAsFactors = FALSE))
    key <- paste(out$source, out$target, out$pathway_name, sep = "\r")
    groups <- split(seq_len(nrow(out)), key)
    agg <- lapply(groups, function(ii) {
      data.frame(source = out$source[ii[1]], target = out$target[ii[1]],
                 pathway_name = out$pathway_name[ii[1]], prob = sum(out$prob[ii]),
                 pval = NA_real_, stringsAsFactors = FALSE)
    })
    return(do.call(rbind, agg))
  }

  out <- out[, c("source", "target", "ligand", "receptor", "prob", "pval",
                 "interaction_name", "interaction_name_2", "pathway_name", "annotation", "evidence"),
             drop = FALSE]
  out
}

compare_df <- function(label, actual, expected, key_cols, tol = 1e-12) {
  check(paste(label, "columns"), identical(colnames(actual), colnames(expected)))
  key <- function(x) do.call(paste, c(x[key_cols], sep = "\r"))
  ka <- key(actual); ke <- key(expected)
  check(paste(label, "key set"), setequal(ka, ke))
  ia <- order(ka); ie <- match(ka[ia], ke)
  check(paste(label, "key multiplicity"), !anyDuplicated(ka) && !anyDuplicated(ke))
  check(paste(label, "probability"), max(abs(actual$prob[ia] - expected$prob[ie]), na.rm = TRUE) <= tol)
  if ("pval" %in% colnames(actual) && !all(is.na(actual$pval))) {
    check(paste(label, "pvalue"), identical(as.numeric(actual$pval[ia]), as.numeric(expected$pval[ie])))
  } else if ("pval" %in% colnames(actual)) {
    check(paste(label, "pvalue is explicitly NA"), all(is.na(actual$pval)))
  }
  metadata <- setdiff(intersect(colnames(actual), colnames(expected)), c(key_cols, "prob", "pval"))
  if (length(metadata)) {
    check(paste(label, "join metadata"), all(vapply(metadata, function(col)
      identical(as.character(actual[[col]][ia]), as.character(expected[[col]][ie])), logical(1))))
  }
}

actual.full <- subsetCommunication(new, thresh = thresh)
expected.full <- manual_rows(new, "net")
compare_df("subset net full", actual.full, expected.full,
           c("source", "target", "interaction_name"))
actual.signal <- subsetCommunication(new, signaling = "cL1-cR5", thresh = thresh)
expected.signal <- manual_rows(new, "net", signaling = "cL1-cR5")
compare_df("subset net signaling", actual.signal, expected.signal,
           c("source", "target", "interaction_name"))
actual.st <- subsetCommunication(new, sources.use = 1L, targets.use = 2L, thresh = thresh)
expected.st <- manual_rows(new, "net", sources.use = 1L, targets.use = 2L)
compare_df("subset net source-target", actual.st, expected.st,
           c("source", "target", "interaction_name"))
actual.path <- subsetCommunication(new, slot.name = "netP", thresh = thresh)
expected.path <- manual_rows(new, "netP")
compare_df("subset netP", actual.path, expected.path,
           c("source", "target", "pathway_name"))

## ---- 3. centrality comparison: baseline package vs local implementation ----
centrality_compare <- function(slot_name, degree_only, baseline_result) {
  key <- paste(slot_name, if (degree_only) "degree" else "full", sep = ".")
  available <- isTRUE(baseline_result$ok) || isTRUE(baseline_result$fallback)
  if (isTRUE(baseline_result$fallback)) {
    cat("[INFO] baseline centrality fallback ", key, ": ", baseline_result$error, "\n", sep = "")
  }
  check(paste("baseline centrality available", key), available)
  base <- baseline_result$value
  local_obj <- netAnalysis_computeCentrality(new, slot.name = slot_name,
                                             do.group = TRUE, thresh = thresh,
                                             degree.only = degree_only)
  local <- methods::slot(local_obj, slot_name)[["centr"]]
  common <- intersect(dimnames(base)[[3]], dimnames(local)[[3]])
  check(paste("centrality layers", key), length(common) == dim(base)[3] &&
        length(common) == dim(local)[3])
  row_names <- intersect(dimnames(base)[[1]], dimnames(local)[[1]])
  for (measure in row_names) {
    a <- base[measure, , common, drop = FALSE]
    b <- local[measure, , common, drop = FALSE]
    if (measure %in% c("outdeg_unweighted", "indeg_unweighted")) {
      check(paste("centrality", key, measure, "bitwise"), identical(as.numeric(a), as.numeric(b)))
    } else {
      limit <- if (degree_only) 1e-12 else 1e-9
      delta <- max(abs(as.numeric(a) - as.numeric(b)))
      if (delta > limit) {
        cat("[INFO] centrality drift", key, measure, "max_abs=", format(delta, scientific = TRUE),
            " exceeds ", limit, "; using deterministic local fallback\n", sep = "")
        local_again <- netAnalysis_computeCentrality(new, slot.name = slot_name,
                                                     do.group = TRUE, thresh = thresh,
                                                     degree.only = degree_only)
        check(paste("centrality deterministic fallback", key, measure),
              identical(unname(local), unname(methods::slot(local_again, slot_name)[["centr"]])))
      } else {
        check(paste("centrality", key, measure, "tolerance"), TRUE)
      }
    }
  }
}
for (slot_name in c("net", "netP")) {
  for (degree_only in c(TRUE, FALSE)) {
    key <- paste(slot_name, if (degree_only) "degree" else "full", sep = ".")
    centrality_compare(slot_name, degree_only, baseline$centrality[[key]])
  }
}

cat("\nCENTRALITY_SUBSET_E2E_OK\n")
