# Instrument a private copy of the installed function for the legacy fixture audit.
# Upstream worker results carry their trace as an attribute back to the caller.
trace_group_permutations <- function(average, engine, trace, reference.seed = NULL) {
  lines <- deparse(body(average), width.cutoff = 500L)
  insert_before <- function(anchor, code) {
    location <- which(trimws(lines) == anchor)
    if (length(location) != 1L) {
      stop("Cannot locate a unique trace point in ", engine, ": ", anchor)
    }
    lines <<- append(lines, deparse(code, width.cutoff = 500L), after = location - 1L)
  }

  insert_before("permutation <- replicate(nboot, sample.int(nC, size = nC))", quote({
    if (!is.null(.legacy_trace_seed)) {
      assign(".Random.seed", .legacy_trace_seed, envir = .GlobalEnv)
    }
    trace.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }))
  if (identical(engine, "upstream")) {
    insert_before("return(Pval.i)", quote({
      attr(Pval.i, "permutation.trace") <- list(Pboot = Pboot, Pnull = Pnull)
    }))
    insert_before("Pval[Prob == 0] <- 1", quote({
      trace.layers <- setNames(lapply(Pval_, attr, which = "permutation.trace"),
                              dimnames(Prob)[[3L]][LRsig.use.idx])
    }))
  } else {
    insert_before("run_perm <- function(i) {", quote({
      trace.layers <- setNames(vector("list", length(LRsig.use.idx.perm)),
                              dimnames(Prob)[[3L]][LRsig.use.idx.perm])
    }))
    insert_before("nReject <- rowSums(Pboot - Pnull > 0)", quote({
      trace.layers[[dimnames(Prob)[[3L]][i]]] <<- list(Pboot = Pboot, Pnull = Pnull)
    }))
  }
  insert_before("Pval[Prob == 0] <- 1", quote({
    .legacy_trace$value <- list(seed.before.drawing = trace.seed, permutation = permutation,
                                nboot = nboot, seed.use = seed.use, min.percent = min.percent,
                                min.cells.sr = min.cells.sr, axes = dimnames(Prob),
                                group = group, layers = trace.layers)
  }))
  body(average) <- parse(text = lines, keep.source = FALSE)[[1L]]
  environment(average) <- list2env(
    list(.legacy_trace = trace, .legacy_trace_seed = reference.seed),
    parent = environment(average)
  )
  average
}
