# ====== computeCommunProb v2: per-LR expression helpers (internal) ======
# Operate directly on the SPARSE scaled signaling layer (genes x cells).
# Semantics mirror baseline computeExpr_LR/_complex/_coreceptor/_agonist/_antagonist;
# per-LR rows are computed on demand with a per-symbol cache (genes repeat across
# LR pairs, distinct names are far fewer than LR pairs).

# Baseline scaling, kept per matrix class: a sparse layer is scaled through @x
# (implicit zeros untouched), a dense layer by dividing the whole matrix. Both
# branches agree whenever every stored value is non-negative, which holds for
# signaling and for the PPI projection.
.sc_scale_expression <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    if (!length(x@x))
      stop("expression layer has no non-zero values; cannot scale it",
           call. = FALSE)
    maximum <- max(x@x)
    if (!is.finite(maximum) || maximum <= 0)
      stop("expression layer has no positive values; cannot scale it", call. = FALSE)
    x@x <- x@x / maximum
    return(x)
  }
  maximum <- max(x)
  if (!is.finite(maximum) || maximum <= 0)
    stop("expression layer has no positive values; cannot scale it", call. = FALSE)
  x / maximum
}

.sc_expr_row <- function(gene, sig, complex_subunits, cache) {
  row <- cache[[gene]]
  if (!is.null(row)) return(row)
  if (gene %in% rownames(sig)) {
    row <- as.numeric(sig[gene, ])
  } else if (!is.null(sub <- complex_subunits[[gene]])) {
    # Upstream .sc_lr_feature_state guarantees all subunits of an available
    # complex are in gene.use; keep an explicit guard (stop vs baseline crash).
    sub <- sub[sub %in% rownames(sig)]
    if (!length(sub) || length(sub) != length(complex_subunits[[gene]]))
      stop("complex '", gene, "' has subunits missing from assay$signaling; ",
           "rerun `subsetData()` and `identifyOverExpressedInteractions()`",
           call. = FALSE)
    # baseline geometricMean semantics: exp(colMeans(log(x), na.rm=TRUE));
    # log(0) = -Inf -> mean -Inf -> exp = 0, i.e. any zero subunit zeroes the complex
    row <- geometricMean(as.matrix(sig[sub, , drop = FALSE]))
  } else {
    stop("gene/complex '", gene, "' not found in assay$signaling rows; ",
         "please rerun `subsetData()`", call. = FALSE)
  }
  cache[[gene]] <- row
  row
}

# Per-cell modulation factor for one cofactor entry (name row in DB$cofactor).
# type: "A" (co-activation, 1+e), "I" (co-inhibition, 1+e), "agonist" (1+Hill(e)),
#       "antagonist" (Hill^-1 = Kh^n/(Kh^n+e^n)). NULL name -> NULL (neutral skip).
# Multi-gene member branch uses apply(..., 2, prod); the single-gene branch of the
# baseline equals the length-1 product exactly (bitwise identical).
.sc_cofactor_factor <- function(cof.name, cofactor_input, sig,
                                type = c("A", "I", "agonist", "antagonist"),
                                Kh = NULL, n = NULL, cache) {
  type <- match.arg(type)
  if (is.null(cof.name) || length(cof.name) != 1L || is.na(cof.name) || !nzchar(cof.name))
    return(NULL)
  key <- paste0(type, ".", cof.name, ".", format(Kh, digits = 17), ".", format(n, digits = 17))
  f <- cache[[key]]
  if (!is.null(f)) return(f)
  cols <- grep("cofactor", colnames(cofactor_input))
  ind <- cofactor_input[cof.name, cols]
  genes <- as.character(unlist(ind, use.names = FALSE))
  # presence filter only: preserve member order and duplicates (baseline semantics)
  genes <- genes[genes != "" & genes %in% rownames(sig)]
  if (!length(genes)) {
    f <- rep(1, ncol(sig))
  } else {
    submat <- as.matrix(sig[genes, , drop = FALSE])
    if (type == "A" || type == "I") {
      f <- apply(1 + submat, 2, prod)
    } else if (type == "agonist") {
      f <- apply(1 + submat^n / (Kh^n + submat^n), 2, prod)
    } else {
      f <- apply(Kh^n / (Kh^n + submat^n), 2, prod)
    }
  }
  cache[[key]] <- f
  f
}
# v2 内部：手工 CSR（dgC -> 行指针 rp + 行内升序排列 ord + 条目列索引 jc）。
# Matrix 无 dgC->dgR coerce 方法；radix order 稳定，行内保持列主序升序。
.sc_csr_build <- function(m) {
  ord <- order(m@i, method = "radix")
  rp <- as.integer(c(0L, cumsum(tabulate(m@i + 1L, nbins = nrow(m)))))
  jc <- rep.int(seq_len(ncol(m)), diff(m@p))
  list(ord = ord, rp = rp, jc = jc)
}

# v2 内部：d 侧 pattern 行索引（等价 which(d[ci, ] > 0)：含显式零排除、含对角、升序）
.sc_pattern_row <- function(csr, m, ci, cache) {
  key <- as.character(ci)
  v <- cache[[key]]
  if (is.null(v)) {
    s <- csr$rp[ci]; e <- csr$rp[ci + 1L]
    seg <- csr$ord[(s + 1L):e]
    v <- csr$jc[seg][m@x[seg] > 0]
    cache[[key]] <- v
  }
  v
}

# Plan A 内部：正式的 numeric min.percent 门控使用一位有效数字舍入。
# contract = signif(mean(...), 1) >= min.percent；这里的 1 是 digits 参数，
# 不是比例阈值。组大小在标签置换下不变，因此每组可预计算一个整数阈值，
# kernel 内退化为 cnt >= thr；单调性在 tests_dev/test-computeAvgCommunProb.R
# 中对 fixture 组大小穷举验证。该契约有意不复现 pinned upstream 的
# format() 字符比较行为，后者仅作为历史对拍参考。
.sc_gate_threshold <- function(n, min.percent) {
  decision <- function(cnt)
    isTRUE(signif(mean(rep(c(1, 0), c(cnt, n - cnt))), 1) >= min.percent)
  if (!decision(n)) return(n + 1L)  # 即使整组全部表达也不通过门控
  lo <- 0L; hi <- n
  while (lo < hi) {
    mid <- lo + (hi - lo) %/% 2L
    if (decision(mid)) hi <- mid else lo <- mid + 1L
  }
  lo
}

#' computeCommunProb
#'
#' @description
#' Compute the communication probability/strength between any interacting individual
#' cells with spatial distance constraints, on the final 11-slot object schema.
#' Sparse optimized: per-LR work is O(nnz(P.spatial)) via the `cpp_prob_layer` kernel
#' (single fused pass, OpenMP-capable); cell-level layers are bitwise identical to the
#' baseline dense-crossprod implementation (see agent note 2026-09-10, Wave 1b tests).
#'
#' @param object A SpatialCellChat object with `assay$signaling` (from
#' \link{subsetData}) and `LR$LRsig` (from \link{identifyOverExpressedInteractions}).
#' @param LR.use Optional custom ligand-receptor table (same columns as `LR$LRsig`).
#' @param raw.use `TRUE` (default) reads `assay$signaling`; `FALSE` reads
#' `assay$smooth`, the PPI projection written by \link{projectData}. The two
#' must cover the same genes, and `assay$smooth` is cleared by `subsetData()`,
#' so `projectData()` has to be re-run after any re-subsetting.
#' @param Kh Hill-function EC50 parameter (input is globally max-scaled to [0, 1]).
#' @param n Hill coefficient (n = 1 uses the division fast path in the kernel).
#' @param distance.use Whether to weight communication by 1/(scaled distance).
#' @param tol Distance tolerance (defaults to `images$spatial.factors$tol`).
#' @param interaction.range Maximum diffusion length in microns.
#' @param scale.distance Distance unit normalization; the check requires
#' min(scaled distance) >= 1 so that SCL <= 1 keeps probabilities in [0, 1].
#' The thesis formula SCL = 1/d corresponds to scale.distance = 1.
#' @param use.AGAN Whether to include agonist/antagonist cofactor modulation.
#' @param contact.dependent Whether `Cell-Cell Contact` signaling is restricted
#' to the contact adjacency pattern.
#' @param contact.range Contact threshold in the working physical distance unit;
#' `tol` is added when building the contact adjacency pattern.
#' @param contact.dependent.forced Force all LRs to be contact-dependent.
#' @param spill.dir Optional directory to write per-LR layers to disk instead of
#' retaining them in the final object (for cell counts where the output tensor
#' exceeds RAM; `filterProbability` consumes and clears the spilled layers). NULL
#' keeps all layers in memory. In parallel mode, workers finish first and the
#' caller then writes the completed layers; this is correct but not streaming
#' during worker execution.
#' @param spill.backend Spill storage backend; currently only "rds" (per-layer
#' saveRDS, bitwise lossless). Further backends (e.g. BPCells) are reserved.
#' @param nthreads OpenMP threads for the kernel (default 1; results are bitwise
#' identical for any thread count).
#' @param parallel Whether to distribute LR iterations via `my_future_lapply`
#' (default FALSE: after sparse optimization the sequential kernel is fastest;
#' measured multisession scheduling overhead dominates at 4.6 ms/LR task size).
#' @param future.scheduling Optional future.scheduling value passed through when
#' parallel = TRUE (default: chunked to the number of workers).
#' @param verbose Logical. FALSE disables ordinary CLI messages and progress.
#' TRUE follows `options(SpatialCellChat.verbose = 0L)`, `1L`, or `2L` for
#' quiet, normal, or detailed output and never overrides global quiet.
#' Progress also respects progressr's enable policy; noninteractive scripts
#' must opt in, for example with `options(progressr.enable = TRUE)`.
#'
#' @return A SpatialCellChat object with `net$cell$prob` as a SparseChatArray
#' (nC x nC x nLR, one named dgCMatrix layer per LR pair) and parameters in
#' `misc$.param$communication`.
#' @export
#'
computeCommunProb <- function(
    object,
    LR.use = NULL,
    raw.use = TRUE,
    Kh = 0.5,
    n = 1,
    distance.use = TRUE,
    tol = NULL,
    interaction.range = 250,
    scale.distance = 0.01,
    use.AGAN = TRUE,
    contact.dependent = TRUE,
    contact.range = 10,
    contact.dependent.forced = FALSE,
    spill.dir = NULL,
    spill.backend = c("rds"),
    nthreads = 1L,
    parallel = FALSE,
    future.scheduling = NULL,
    verbose = TRUE) {
  t0 <- Sys.time()
  object <- .sc_assert_spatial_cell_chat(object)
  if (isTRUE(verbose)) .cli("computeCommunProb", .type = "subheader")
  if (length(raw.use) != 1L || !is.logical(raw.use) || is.na(raw.use))
    stop("`raw.use` must be TRUE or FALSE", call. = FALSE)
  spill.backend <- match.arg(spill.backend)
  if (!is.null(spill.dir)) {
    if (!dir.exists(spill.dir)) dir.create(spill.dir, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(verbose)) .cli("Spill mode: layers stream to {.val {spill.dir}} (backend {.val {spill.backend}})",
                              .type = "text", .verbose = 2L)
  }

  # ---- LR table (baseline semantics preserved) ----
  if (is.null(LR.use)) {
    pairLRsig <- object@LR$LRsig
    if (is.null(pairLRsig) || !nrow(pairLRsig))
      stop("`LR$LRsig` is empty; run `identifyOverExpressedInteractions()` first",
           call. = FALSE)
  } else {
    if (length(unique(LR.use$annotation)) > 1) {
      LR.use$annotation <- factor(LR.use$annotation,
        levels = c("Secreted Signaling", "ECM-Receptor",
                   "Non-protein Signaling", "Cell-Cell Contact"))
      LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
      LR.use$annotation <- as.character(LR.use$annotation)
    }
    pairLRsig <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group)))
    stop("Please check `unique(object@idents)` and drop unused levels with ",
         "`droplevels()` before running `computeCommunProb()`.", call. = FALSE)

  # raw.use = TRUE reads the subsetData gene subset; FALSE reads the PPI
  # projection written by projectData(), which must cover the same space.
  layer <- if (raw.use) "signaling" else "smooth"
  sig <- assay(object, layer)
  if (is.null(sig) || sum(dim(sig)) == 0) {
    if (raw.use)
      stop("assay$signaling is empty; run `subsetData()` first", call. = FALSE)
    stop("assay$smooth is empty; run `projectData()` first, or use raw.use = TRUE",
         call. = FALSE)
  }
  if (!raw.use) {
    signaling <- assay(object, "signaling")
    if (is.null(signaling) || !identical(rownames(sig), rownames(signaling)) ||
        !identical(colnames(sig), colnames(signaling)))
      stop("assay$smooth must carry the same genes and cells as assay$signaling; ",
           "rerun `projectData()` after the last `subsetData()`", call. = FALSE)
  }
  cell.names <- colnames(sig)
  nC <- ncol(sig)
  # 全局 max 缩放（基线语义：sparse 用 @x / max(@x)，dense 用 data / max(data)）；副本，不改对象
  sig.scaled <- .sc_scale_expression(sig)
  if (isTRUE(verbose)) .cli("Input: assay${layer} ({.val {nrow(sig)}} genes x {.val {nC}} cells); {.val {nLR}} LR pairs",
                            .type = "info")

  # ---- distance cache: validate or compute (baseline semantics) ----
  ratio <- object@images$spatial.factors[["ratio"]]
  if (is.null(tol)) tol <- object@images$spatial.factors[["tol"]]
  if (isTRUE(verbose)) .cli("Distance limits: interaction.range = {.val {interaction.range}}, contact.range = {.val {contact.range}}, tol = {.val {if (is.null(tol)) 0 else tol}}; ratio = {.val {if (is.null(ratio)) 'NULL (coordinates already in working units)' else ratio}}. Ranges and tolerance use the working distance unit after ratio conversion.",
                            .type = "text", .verbose = 2L)
  data.spatial <- as.matrix(object@images$coordinates)
  cache.d <- object@images$.distance
  expected <- list(interaction.range = interaction.range,
                   contact.range = contact.range, ratio = ratio, tol = tol)
  valid.cache <- is.list(cache.d) &&
    all(c("d.spatial", "adj.contact", ".parameters") %in% names(cache.d)) &&
    isTRUE(all.equal(cache.d$.parameters, expected)) &&
    ncol(cache.d$d.spatial) == nC
  if (!valid.cache) {
    if (isTRUE(verbose)) .cli("Distance cache absent or stale; recomputing with `computeCellDistance()`.",
                              .type = "text", .verbose = 2L)
    res <- computeCellDistance(coordinates = data.spatial, ratio = ratio,
                               interaction.range = interaction.range,
                               contact.range = contact.range, tol = tol)
    images(object, ".distance") <- res
    cache.d <- res
  }
  if (isTRUE(verbose)) .cli("Distance cache: {.val {length(cache.d$d.spatial@x)}} distance entries; {.val {length(cache.d$adj.contact@x)}} contact entries",
                            .type = "text", .verbose = 2L)
  res <- cache.d

  # ---- P.spatial + loop invariants ----
  if (distance.use) {
    d.work <- cache.d$d.spatial
    d.work@x <- d.work@x * scale.distance
    d.min <- min(d.work@x)
    if (d.min < 1)
      stop(sprintf(paste0("The minimum scaled distance is %.3g but must be >= 1 ",
                          "(SCL = 1/d must stay <= 1 for probability semantics). ",
                          "Increase `scale.distance` to at least %.4g (current: %g)."),
                   d.min, scale.distance / d.min, scale.distance), call. = FALSE)
    P.spatial <- createPspatialFrom_dspatial(d.work, distance.use = TRUE)
  } else {
    P.spatial <- createPspatialFrom_dspatial(cache.d$d.spatial, distance.use = FALSE)
  }
  x0 <- P.spatial@x
  row0 <- P.spatial@i
  col0 <- rep.int(0L:(nC - 1L), diff(P.spatial@p))  # 0-based，与 row0/kernel C++ 索引一致
  nnz.P <- length(x0)

  # contact adjacency mask aligned to the P.spatial pattern (computed once)
  adj.ts <- methods::as(cache.d$adj.contact, "TsparseMatrix")
  idx <- match(adj.ts@i + adj.ts@j * nC, row0 + col0 * nC)
  if (anyNA(idx)) stop("adj.contact pattern must be a subset of the distance pattern",
                       call. = FALSE)
  mask.one <- rep(1, nnz.P)
  # 接触 LR 的掩码：仅 adj 位置保留（1），其余位置清 0（基线 P1_Pspatial * adj.contact
  # 的元素乘语义 = 交集 pattern）
  mask.contact <- numeric(nnz.P)
  mask.contact[idx] <- adj.ts@x

  # contact branch flags (baseline semantics)
  all.contact.dependent <- FALSE
  all.diffusible <- FALSE
  if (contact.dependent.forced == TRUE) {
    mask.one <- mask.contact
    nLR1 <- nLR
    all.contact.dependent <- TRUE
  } else if (contact.dependent == TRUE && length(unique(pairLRsig$annotation)) > 0) {
    if (all(unique(pairLRsig$annotation) == "Cell-Cell Contact")) {
      mask.one <- mask.contact
      nLR1 <- nLR
      all.contact.dependent <- TRUE
    } else if (all(unique(pairLRsig$annotation) %in%
                   c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
      nLR1 <- nLR
      all.diffusible <- TRUE
    } else {
      nLR1 <- max(which(pairLRsig$annotation %in%
                          c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
    }
  } else {
    nLR1 <- nLR
  }
  if (isTRUE(verbose)) {
    nLR.diffusible <- if (all.contact.dependent) 0L else nLR1
    .cli("LR split: {.val {nLR.diffusible}} diffusible / {.val {nLR - nLR.diffusible}} contact-dependent (forced = {.val {contact.dependent.forced}})",
         .type = "text", .verbose = 2L)
    .cli("Execution: {.val {if (isTRUE(parallel)) 'future' else 'serial'}}; {.val {nthreads}} kernel threads; distance weighting = {.val {distance.use}}; agonists/antagonists = {.val {use.AGAN}}",
         .type = "text", .verbose = 2L)
  }

  # ---- per-LR loop ----
  complex_subunits <- .sc_complex_subunits(complex_input)
  expr.cache <- new.env(parent = emptyenv())
  cof.cache <- new.env(parent = emptyenv())
  ag.names <- as.character(pairLRsig$agonist)
  an.names <- as.character(pairLRsig$antagonist)
  has.ag <- !is.na(ag.names) & nzchar(ag.names)
  has.an <- !is.na(an.names) & nzchar(an.names)
  has.AGAN.lr <- use.AGAN & (has.ag | has.an)

  layers <- vector("list", nLR)
  names(layers) <- rownames(pairLRsig)
  nnz.cum <- 0
  guard.checked <- FALSE
  guard.nLR <- max(5L, ceiling(nLR * 0.05))
  GUARD.BYTES <- 16e9
  spill.i <- 0L

  run_one <- function(i) {
    Li <- .sc_expr_row(geneL[i], sig.scaled, complex_subunits, expr.cache)
    Ri <- .sc_expr_row(geneR[i], sig.scaled, complex_subunits, expr.cache)
    fA <- .sc_cofactor_factor(pairLRsig$co_A_receptor[i], cofactor_input, sig.scaled,
                              type = "A", cache = cof.cache)
    fI <- .sc_cofactor_factor(pairLRsig$co_I_receptor[i], cofactor_input, sig.scaled,
                              type = "I", cache = cof.cache)
    if (!is.null(fA) && !is.null(fI)) Ri <- Ri * fA / fI
    else if (!is.null(fA)) Ri <- Ri * fA
    else if (!is.null(fI)) Ri <- Ri / fI
    cmask <- if (i > nLR1) mask.contact else mask.one
    out <- cpp_prob_layer(x0, row0, col0, Li, Ri, Kh = Kh, n = n,
                          contact_mask = cmask,
                          fAG = if (has.AGAN.lr[i] && has.ag[i])
                            .sc_cofactor_factor(ag.names[i], cofactor_input, sig.scaled,
                                                type = "agonist", Kh = Kh, n = n,
                                                cache = cof.cache) else NULL,
                          fAN = if (has.AGAN.lr[i] && has.an[i])
                            .sc_cofactor_factor(an.names[i], cofactor_input, sig.scaled,
                                                type = "antagonist", Kh = Kh, n = n,
                                                cache = cof.cache) else NULL,
                          nC = nC, nthreads = nthreads)
    new("dgCMatrix", i = out$i, p = out$p, x = out$x, Dim = c(nC, nC))
  }

  if (isTRUE(parallel)) {
    sched <- if (is.null(future.scheduling)) ceiling(nLR / max(future::nbrOfWorkers(), 1L)) else future.scheduling
    Prob.cell_ <- my_future_lapply(seq_len(nLR), run_one,
                                   future.scheduling = sched, future.seed = TRUE,
                                   .progress = verbose)
    names(Prob.cell_) <- rownames(pairLRsig)
    nnz.cum <- sum(vapply(Prob.cell_, function(m) length(m@x), numeric(1)))
    if (nLR >= guard.nLR) {
      guard.checked <- TRUE
      proj <- nnz.cum / nLR * nLR * 12
      if (proj > GUARD.BYTES && is.null(spill.dir))
        stop(sprintf(paste0("Projected cell-level output tensor is ~%.1f GB ",
                            "(> %.0f GB in-memory threshold). Set `spill.dir` to ",
                            "stream layers to disk, or reduce broadly expressed ",
                            "LR pairs."), proj / 1e9, GUARD.BYTES / 1e9),
             call. = FALSE)
    }
    layers <- Prob.cell_
    if (!is.null(spill.dir)) {
      for (i in seq_len(nLR)) {
        saveRDS(layers[[i]], file.path(spill.dir, paste0(names(layers)[i], ".rds")),
                compress = FALSE)
        spill.i <- spill.i + 1L
      }
      layers <- NULL
      rm(Prob.cell_)
    }
  } else {
    progressr::with_progress({
      pr <- progressr::progressor(along = seq_len(nLR),
                                 message = "Computing LR probabilities",
                                 enable = .sc_progress_enabled(verbose) && nLR > 1L &&
                                   getOption("progressr.enable", TRUE))
      lapply(seq_len(nLR), function(i) {
        layer <- run_one(i)
        nnz.cum <<- nnz.cum + length(layer@x)
        if (!guard.checked && i >= guard.nLR) {
          guard.checked <<- TRUE
          proj <- nnz.cum / i * nLR * 12
          if (proj > GUARD.BYTES && is.null(spill.dir))
            stop(sprintf(paste0("Projected cell-level output tensor is ~%.1f GB ",
                                "(> %.0f GB in-memory threshold). Set `spill.dir` to ",
                                "stream layers to disk, or reduce broadly expressed ",
                                "LR pairs."), proj / 1e9, GUARD.BYTES / 1e9),
                 call. = FALSE)
        }
        if (!is.null(spill.dir)) {
          saveRDS(layer, file.path(spill.dir, paste0(names(layers)[i], ".rds")),
                  compress = FALSE)
          spill.i <<- spill.i + 1L
        } else {
          layers[[i]] <<- layer
        }
        pr()
        NULL
      })
    })
  }
  if (is.null(spill.dir)) {
    prob.array <- SparseChatArray(layers)
    dimnames(prob.array) <- list(cell.names, cell.names, rownames(pairLRsig))
    object@net$cell$prob <- prob.array
  } else {
    # spill 模式：层在磁盘上，net$cell$prob 由 filterProbability 流式回填（validator 容忍 NULL）
    object@net$cell$prob <- NULL
  }
  object@misc$.param$communication <- list(
    raw.use = raw.use, Kh = Kh, n = n, nLR = nLR, nLR1 = nLR1,
    layer.names = rownames(pairLRsig),
    scale.distance = scale.distance, use.AGAN = use.AGAN,
    distance.use = distance.use, interaction.range = interaction.range,
    contact.range = contact.range, contact.dependent = contact.dependent,
    contact.dependent.forced = contact.dependent.forced,
    all.contact.dependent = all.contact.dependent, all.diffusible = all.diffusible,
    spill.dir = if (is.null(spill.dir)) NULL else normalizePath(spill.dir),
    spill.backend = if (is.null(spill.dir)) NULL else spill.backend,
    spill.count = if (is.null(spill.dir)) NULL else spill.i,
    run.time = as.numeric(Sys.time() - t0, units = "secs")
  )
  object <- .log_operation(object, "computeCommunProb", params = list(
    layer = "signaling", nLR = nLR, nthreads = nthreads,
    spill.backend = spill.backend,
    spill.dir = if (is.null(spill.dir)) NULL else normalizePath(spill.dir),
    run.time = as.numeric(Sys.time() - t0, units = "secs")))
  if (isTRUE(verbose)) .cli("computeCommunProb done: {.val {nC}} x {.val {nC}} x {.val {nLR}} layers; nnz = {.val {nnz.cum}}; {.val {round(as.numeric(Sys.time() - t0, units = 'secs'), 2)}}s",
                            .type = "success")
  object
}


#' create P.spatial from d.spatial object
#'
#' @description
#' A function defined for \link{computeCommunProb}. To model ligands' diffusion attenuation, we take the reciprocal of cell-cell distance
#' as a factor (named as P.spatial) determining the cell-cell communication probability. Here we design
#' the function to help take reciprocal of a `d.spatial` matrix.
#'
#' @param d.spatial CsparseMatrix (dsC Matrix) object. A cell-cell distance matrix which has been scaled, see in \link{computeCommunProb}
#' @param distance.use Boolean. Whether to use distance constraints to compute communication probability, passed from \link{computeCommunProb}'s same parameter
#' Setting `distance.use = TRUE` indicates that the cell-cell communication probability is inversely proportional to the computed distance.
#' Otherwise, the cell-cell communication probability share the same computation factor 1 within the maximum interaction/diffusion length of ligands
#' @return CsparseMatrix object in `Matrix` Package
#' @export
#'
#' @examples
createPspatialFrom_dspatial <- function(d.spatial,distance.use=T){
  # d.spatial is a sparse matrix, note if dividing it directly, many 1/0 => inf will arise
  P.spatial <- d.spatial
  P.spatial@x <- 1/P.spatial@x
  # Fill up the diagonal: autocrine has the largest spatial CCC probability
  Matrix::diag(P.spatial) <- max(P.spatial@x)

  # make P.spatial binary
  if (distance.use) {
    return(P.spatial)
  } else {
    P.spatial@x <- rep.int(1,times = length(P.spatial@x))
    return(P.spatial)
  }
}






#' @title relabelSpatialCellChat
#' @description
#' Merge some cell groups into a new cell group, and recalculate cell-cell communication
#' probability/strength at the cell-group level.
#'
#' @param object SpatialCellChat object.
#' @param labelSet Named character vector. Make sure the names are cell group labels in
#' present `object@idents`, and values are new cell group labels respectively.
#' For example,`labelSet = c("Group1"="NewGroup1","Group2"="NewGroup1","Group3"="NewGroup1")`.
#' @param group.by Set a column name in `object@meta` if not use `object@idents` as cell group labels
#' to relabel.
#' @param nboot Integer. See the same parameter in \link{computeAvgCommunProb} or \link{computeAvgCommunProb_Visium}.
#' @param min.cells.sr Integer. See the same parameter in \link{computeAvgCommunProb}.
#' @param min.percent Numeric. See the same parameter in \link{computeAvgCommunProb}.
#' @param do.permutation Boolean. Integer. See the same parameter in \link{computeAvgCommunProb} or \link{computeAvgCommunProb_Visium}.
#' @param min.cells Integer. See the same parameter in \link{filterCommunication}.
#' @param thresh Numeric. See the same parameter in \link{computeCommunProbPathway}.
#'
#' @return SpatialCellChat object.
#' @export
relabelSpatialCellChat <- function(
    object,
    labelSet,
    group.by=NULL,
    nboot = 100,
    min.cells.sr = 10,
    min.percent = 0.1,
    do.permutation = T,
    min.cells = 10,
    thresh = 0.05
){
  # object <- Chat
  #
  # labelSet <- c(
  #   "Basal" = "Keratinocyte",
  #   "Spinous" = "Keratinocyte",
  #   "Supraspinous" = "Keratinocyte",
  #   "TC" = "Immune",
  #   "MYL" = "Immune"
  # )
  if(is.null(group.by)){
    group.lab <- as.character(object@idents)
    group.level <- levels(object@idents)
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `meta` slot, which will be used for cell grouping.")
    } else {
      group.lab <- as.character(object@meta[[group.by]])
      # group.lab <- factor(group.lab,levels = unique(group.lab))
      # object@meta[[group.by]] <- group.lab
      group.level <- unique(group.lab)
    }
  }

  if(is.null(names(labelSet))){
    stop("Please provide the names of cell groups which will be dropped
    by setting the `names` of labelSet.\n")
  }
  group.drop <- names(labelSet)
  if(!all(group.drop%in%group.level)) {
    stop("Please check the `labelSet`, make sure the cell groups
    to be dropped are in the `group.by` column.\n")
  }
  .cli("Relabeling {.val {length(group.drop)}} of {.val {length(group.level)}} cell groups; recomputing group and pathway summaries",
       .type = "text", .verbose = 2L)

  tmp.df <- data.frame(
    pre = group.lab,
    new = labelSet[group.lab]
  )
  tmp.df <- tmp.df %>% mutate(
    new = if_else(is.na(new),pre,new)
  )

  object@meta[["new.ident"]] <- tmp.df[["new"]]
  # spatialDimPlot(object,group.by = "new.ident")
  # spatialDimPlot(object,group.by = "ident")

  decomp <- object@misc$.param$averaging$cell.type.decomposition
  if(!is.null(decomp)){
    cell.type.decomposition <- as.matrix(decomp)
    cell.type.decomposition.left <- cell.type.decomposition[,base::setdiff(group.level,group.drop),drop=F]

    cell.type.decomposition.new <- lapply(
      X = unique(labelSet),
      FUN = function(i){
        group.drop.i <- names(labelSet[labelSet==i])
        cell.type.decomposition.i <-
          Matrix::rowSums(cell.type.decomposition[,group.drop.i,drop=F])
        return(cell.type.decomposition.i)
      }
    )
    names(cell.type.decomposition.new) <- unique(labelSet)
    cell.type.decomposition.new <- as.data.frame(cell.type.decomposition.new)
    cell.type.decomposition.new <- cbind(cell.type.decomposition.left,as.matrix(cell.type.decomposition.new))
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]],levels = colnames(cell.type.decomposition.new))
    object@idents <- object@meta[["new.ident"]]
    names(object@idents) <- colnames(object@assay$norm)


    object <- computeAvgCommunProb_Visium(
      object,
      cell.type.decomposition = cell.type.decomposition.new,
      nboot = nboot,
      # min.cells.sr = min.cells.sr,
      # min.percent = min.percent,
      do.permutation = do.permutation
    )
  } else {
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]],levels = unique(object@meta[["new.ident"]]))
    object@idents <- object@meta[["new.ident"]]
    names(object@idents) <- colnames(object@assay$norm)
    object <- computeAvgCommunProb(
      object,
      nboot = nboot,
      min.cells.sr = min.cells.sr,
      min.percent = min.percent,
      do.permutation = do.permutation
    )
  }


  object <- filterCommunication(
    object,
    min.cells = min.cells,
    min.links = NULL,
    min.cells.sr = NULL
  )

  object <- computeCommunProbPathway(
    object,
    thresh = thresh
  )

  object <- aggregateNet(object)
  .cli("Relabeling done: {.val {length(group.level)}} original groups to {.val {nlevels(object@idents)}} groups",
       .type = "text", .verbose = 2L)

  return(object)
}



#' Compute group-level cell-cell communication for Visium Low Definition (LD)
#'
#' @description
#' Consume cell-level communication probabilities from the final 11-slot schema
#' and aggregate them with a cell-type decomposition. For Visium HD ST datasets,
#' please use \link{computeAvgCommunProb} instead.
#'
#' @param object SpatialCellChat object with communication probabilities in
#'   `net$cell$prob`.
#' @param cell.type.decomposition NCells x NTypes matrix of cell-type proportions.
#' @param group.by Cell-group metadata column; NULL uses `object@idents`.
#' @param avg.type Methods for integrating communication probabilities per cell group.
#' @param do.permutation Whether to perform the group-column permutation test.
#' @param nboot Number of permutations when `do.permutation = TRUE`.
#' @param seed.use Random seed.
#' @param colocalization.use Whether to gate interactions by colocalization.
#' @param thresh.colo Remove interactions whose colocalization p-value exceeds this threshold.
#'
#' @return A SpatialCellChat object with updated final 11-slot `net`:
#'   `net$group$prob` is a `SparseChatArray` of group-level LR probabilities and
#'   `net$group$pval` is the corresponding `SparseChatArray`, or NULL without
#'   permutation. Runtime counts and parameters, including the decomposition, are
#'   stored in `net$tmp` and `misc$.param$averaging`, respectively.
#'
#' @export
computeAvgCommunProb_Visium <- function(
    object,
    cell.type.decomposition,
    avg.type = c("avg","sum"),
    group.by = NULL,
    do.permutation = TRUE,
    nboot = 100,
    seed.use = 1L,
    colocalization.use = FALSE,
    thresh.colo = 0.05
){
  t0 <- Sys.time()
  object <- .sc_assert_spatial_cell_chat(object)

  if (is.null(group.by)) {
    group <- object@idents
  } else {
    if (length(group.by) != 1L || is.na(group.by) ||
        !(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `object@meta`, which will be used for cell grouping.",
           call. = FALSE)
    }
    group <- object@meta[[group.by]]
  }
  if (is.null(group) || length(group) == 0L)
    stop("The selected group labels must not be empty.", call. = FALSE)
  if (anyNA(group) || any(!nzchar(trimws(as.character(group)))))
    stop("The selected group labels must be non-missing and non-empty.", call. = FALSE)
  if (!is.factor(group)) group <- factor(group)
  group <- droplevels(group)
  if (nlevels(group) == 0L)
    stop("The selected group labels must contain at least one group.", call. = FALSE)

  if (length(do.permutation) != 1L || !is.logical(do.permutation) ||
      is.na(do.permutation))
    stop("do.permutation must be one non-missing logical value.", call. = FALSE)
  do.permutation <- isTRUE(do.permutation)
  if (do.permutation &&
      (length(nboot) != 1L || !is.numeric(nboot) || !is.finite(nboot) ||
       nboot < 1 || nboot != floor(nboot) || nboot > .Machine$integer.max))
    stop("nboot must be a positive integer when do.permutation = TRUE.", call. = FALSE)
  if (do.permutation) nboot <- as.integer(nboot)

  prob.array <- object@net$cell$prob
  if (is.null(prob.array) || !inherits(prob.array, "SparseChatArray"))
    stop("`object@net$cell$prob` must be a SparseChatArray; run `computeCommunProb` first.",
         call. = FALSE)
  prob.dim <- dim(prob.array)
  prob.dimnames <- dimnames(prob.array)
  if (length(prob.dim) != 3L || length(prob.dimnames) != 3L ||
      is.null(prob.dimnames[[1L]]) || is.null(prob.dimnames[[2L]]) ||
      is.null(prob.dimnames[[3L]]) ||
      !identical(prob.dimnames[[1L]], prob.dimnames[[2L]]) ||
      prob.dim[[1L]] != prob.dim[[2L]] ||
      length(prob.dimnames[[1L]]) != prob.dim[[1L]] ||
      length(prob.dimnames[[3L]]) != prob.dim[[3L]] ||
      anyDuplicated(prob.dimnames[[1L]]) || anyDuplicated(prob.dimnames[[3L]]))
    stop("`object@net$cell$prob` must have named, identically ordered square cell axes and named LR layers.",
         call. = FALSE)
  cell.names <- prob.dimnames[[1L]]
  LRsig <- prob.dimnames[[3L]]
  prob.cell_ <- unclass(prob.array)
  if (!all(vapply(prob.cell_, inherits, logical(1L), "dgCMatrix")))
    stop("Every layer of `object@net$cell$prob` must be a dgCMatrix.", call. = FALSE)
  nC <- prob.dim[[1L]]

  if (length(group) != nC)
    stop("Group label count (", length(group), ") != cell count (", nC,
         "); `group.by` must cover exactly the cells of `net$cell$prob`.", call. = FALSE)
  if (!is.null(names(group))) {
    if (anyDuplicated(names(group)) || !setequal(names(group), cell.names))
      stop("Named group labels must contain exactly the cell names of `net$cell$prob`.",
           call. = FALSE)
    group <- group[match(cell.names, names(group))]
    names(group) <- cell.names
  }
  numCluster <- nlevels(group)
  if (numCluster < 1L || length(unique(group)) != numCluster)
    stop("The selected group labels must contain at least one observed group.", call. = FALSE)

  if (!(is.matrix(cell.type.decomposition) || inherits(cell.type.decomposition, "Matrix")) ||
      length(dim(cell.type.decomposition)) != 2L)
    stop("`cell.type.decomposition` must be a matrix or Matrix object.", call. = FALSE)
  if (!identical(as.integer(dim(cell.type.decomposition)), c(nC, numCluster)) ||
      !identical(rownames(cell.type.decomposition), cell.names) ||
      !identical(colnames(cell.type.decomposition), levels(group)))
    stop("`cell.type.decomposition` dimensions and dimnames must exactly match the cell axis and group levels.",
         call. = FALSE)
  decomposition.values <- .sc_matrix_values(cell.type.decomposition)
  if (!is.numeric(decomposition.values) || any(!is.finite(decomposition.values)) ||
      any(decomposition.values < 0))
    stop("`cell.type.decomposition` must contain only finite, non-negative values.",
         call. = FALSE)
  cell.type.onehot <- 1 * (cell.type.decomposition > 0)
  avg.type <- match.arg(avg.type)

  if (isTRUE(colocalization.use)) {
    pval.colo <- computeColocalization(
      coordinates = object@images$coordinates,
      group = group,
      nboot = nboot,
      seed.use = seed.use
    )
  } else {
    pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster,
                        dimnames = list(levels(group), levels(group)))
  }

  computeAvgCommunProb_LR_Visium <- if (identical(avg.type, "avg")) {
    computeAvgCommunProb_LR_Visium_Avg
  } else {
    computeAvgCommunProb_LR_Visium_Sum
  }
  cat(paste0(cli.symbol(), "Compute group-level cell-cell communication... <<< [",
             Sys.time(), "]\n"))

  ccc.counts <- vapply(prob.cell_, function(mat) length(mat@x), numeric(1L))
  names(ccc.counts) <- LRsig
  LRsig.use.idx <- which(ccc.counts > 0)
  if (length(LRsig.use.idx) < 1L)
    stop("Each LR pair does not have any spot-level links/interactions.", call. = FALSE)

  prob.results <- lapply(seq_len(length(LRsig)), function(i)
    matrix(0, nrow = numCluster, ncol = numCluster,
           dimnames = list(levels(group), levels(group))))
  names(prob.results) <- LRsig
  for (i in LRsig.use.idx) {
    prob.avg <- computeAvgCommunProb_LR_Visium(
      prob = prob.cell_[[i]],
      cell.type.proportion = cell.type.decomposition,
      cell.type.onehot = cell.type.onehot
    )
    dimnames(prob.avg) <- list(levels(group), levels(group))
    prob.avg[pval.colo > thresh.colo] <- 0
    prob.results[[i]] <- prob.avg
  }

  ggc.counts <- vapply(prob.results, function(mat) sum(mat > 0), numeric(1L))
  names(ggc.counts) <- LRsig
  LRsig.use.idx.perm <- which(ggc.counts > 0)
  pval.results <- NULL

  if (do.permutation) {
    cat(paste0(cli.symbol(), "Perform permutation test for group-level communication... <<< [",
               Sys.time(), "]\n"))
    set.seed(seed.use)
    permutation <- matrix(
      replicate(nboot, sample.int(numCluster, size = numCluster)),
      nrow = numCluster, ncol = nboot
    )
    pval.results <- lapply(seq_len(length(LRsig)), function(i)
      matrix(1, nrow = numCluster, ncol = numCluster,
             dimnames = list(levels(group), levels(group))))
    names(pval.results) <- LRsig
    for (i in LRsig.use.idx.perm) {
      pnull <- as.vector(prob.results[[i]])
      pboot <- vapply(seq_len(nboot), function(nE) {
        decomposition.boot <- cell.type.decomposition[, permutation[, nE], drop = FALSE]
        onehot.boot <- 1 * (decomposition.boot > 0)
        as.vector(computeAvgCommunProb_LR_Visium(
          prob = prob.cell_[[i]],
          cell.type.proportion = decomposition.boot,
          cell.type.onehot = onehot.boot
        ))
      }, numeric(numCluster * numCluster))
      nReject <- rowSums(pboot - pnull > 0)
      pval.i <- matrix(nReject / nboot, nrow = numCluster, ncol = numCluster,
                       byrow = FALSE,
                       dimnames = list(levels(group), levels(group)))
      pval.i[prob.results[[i]] == 0] <- 1
      pval.results[[i]] <- pval.i
    }
  }

  group.layers <- lapply(prob.results, function(mat) as(mat, "dgCMatrix"))
  names(group.layers) <- LRsig
  if (is.null(object@net$group)) object@net$group <- list()
  object@net$group$prob <- SparseChatArray(
    group.layers,
    dimnames = list(levels(group), levels(group), LRsig)
  )
  if (do.permutation) {
    pval.layers <- lapply(pval.results, function(mat) as(mat, "dgCMatrix"))
    names(pval.layers) <- LRsig
    object@net$group$pval <- SparseChatArray(
      pval.layers,
      dimnames = list(levels(group), levels(group), LRsig)
    )
  } else {
    object@net$group$pval <- NULL
  }

  if (is.null(object@net$tmp) || !is.list(object@net$tmp)) object@net$tmp <- list()
  object@net$tmp$prob.cell <- NULL
  object@net$tmp$cell.type.decomposition <- NULL
  object@net$tmp$LRsig.CCC.counts <- ccc.counts
  object@net$tmp$LRsig.GGC.counts <- ggc.counts
  object@net$tmp$LRsig.use.idx <- LRsig.use.idx
  object@net$tmp$LRsig.use.idx.perm <- LRsig.use.idx.perm

  run.time <- as.numeric(Sys.time() - t0, units = "secs")
  averaging.params <- list(
    avg.type = avg.type,
    do.permutation = do.permutation,
    nboot = nboot,
    seed.use = seed.use,
    colocalization.use = colocalization.use,
    thresh.colo = thresh.colo,
    group.by = group.by,
    nLR = length(LRsig),
    numCluster = numCluster,
    LRsig.CCC.counts = ccc.counts,
    LRsig.GGC.counts = ggc.counts,
    LRsig.use.idx = LRsig.use.idx,
    LRsig.use.idx.perm = LRsig.use.idx.perm,
    cell.type.decomposition = cell.type.decomposition,
    run.time = run.time
  )
  if (is.null(object@misc$.param) || !is.list(object@misc$.param))
    object@misc$.param <- list()
  object@misc$.param$averaging <- averaging.params
  if (isTRUE(colocalization.use)) object@images$colocalization <- pval.colo
  object <- .log_operation(object, "computeAvgCommunProb_Visium", params = averaging.params)
  .cli("Inference of group-level cell-cell communication is done. Parameters are stored in `object@misc$.param$averaging`",
       .type = "success")
  object
}


#' @title computeAvgCommunProb_LR_Visium_Sum
#'
#' @description
#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#' for 10X Visium Low Definition.
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param cell.type.proportion Ncells x NTypes Matrix, which comes from any kind of cell type decomposition method.
#' @param cell.type.onehot NULL. Placeholder.
#'
#' @return A NTypes x NTypes matrix of communication probabilities for pairwise individual cell types for one particular ligand-receptor pair/signaling pathway
#' @export
computeAvgCommunProb_LR_Visium_Sum  <- function (
    prob,
    cell.type.proportion,
    cell.type.onehot = NULL
){
  # calculate communication strength across cell type pairs
  Prob.avg <- Matrix::crossprod(x = cell.type.proportion,
                                y = prob %*% cell.type.proportion)
  Prob.avg <- as.matrix(Prob.avg)
  return(Prob.avg)
}


#' @title computeAvgCommunProb_LR_Visium_Avg
#'
#' @description
#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#' for 10X Visium Low Definition.
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param cell.type.proportion Ncells x NTypes Matrix, which comes from any kind of cell type decomposition method.
#' @param cell.type.onehot Ncells x NTypes Matrix, which differs from \code{cell.type.proportion}, can get through \code{cell.type.onehot <- 1*(cell.type.proportion>0)}
#'
#' @return A NTypes x NTypes matrix of communication probabilities for pairwise individual cell types for one particular ligand-receptor pair/signaling pathway
#' @export
computeAvgCommunProb_LR_Visium_Avg  <- function (
    prob,
    cell.type.proportion,
    cell.type.onehot
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5
  # cell.type.proportion <- cell.type.decomposition

  # calculate communication strength across cell type pairs
  Prob.avg <- Matrix::crossprod(x = cell.type.proportion,
                                y = prob %*% cell.type.proportion)

  # count communication links  across cell type pairs
  prob@x <- rep.int(1, times = length(prob@x))
  Prob.scale.factor <- Matrix::crossprod(x = cell.type.onehot,
                                         y = prob %*% cell.type.onehot)
  Prob.avg <- Prob.avg/Prob.scale.factor
  Prob.avg <- as.matrix(Prob.avg)
  Prob.avg[is.nan(Prob.avg)] <- 0
  return(Prob.avg)
}



# #' @title computeAvgCommunProb_LR_Visium
# #'
# #' @description
# #' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
# #' for 10X Visium Low Definition.
# #'
# #' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
# #' @param group cell group information used for computing averaged communication probabilities
# #' @param cell.type.decomposition Ncells x NTypes Matrix, which comes from any kind of cell type decomposition.
# #' @param dataLR a nCell*2 expression matrix of a given pair of ligand-receptor.
# #' @param do.filter Boolean. Whether to filter out some interactions of a given pair of ligand-receptor
# #' between several cell groups according to the parameters `min.percent` and `min.cells.sr`
# #' when computing the average communication probability of a given pair of ligand-receptor between cell groups.
# #' Default is TRUE.
# #' @param min.percent Numeric from 0 to 1. The expressed-ligand/receptor proportion is rounded to one significant digit numerically and must meet this threshold for group-level signaling
# #' Default is 0.1.
# #' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
# #' Default is 5.
# #'
# #' @return
# #' @export
# computeAvgCommunProb_LR_Visium  <- function (
#     prob,
#     group,
#     cell.type.decomposition,
#     dataLR = NULL,
#     do.filter = T,
#     min.percent = 0.1,
#     min.cells.sr = 5
# ){
#   # code for test
#   # prob <- prob.cell.i
#   # group <- group
#   # dataLR = dataLR_temp
#   # min.percent = 0.1
#   # min.cells.sr = 5
#
#
#   # group is vector/factor
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   cell.type.mat <- cell.type.decomposition
#
#   if(do.filter){
#     # min.percent
#     dataLR_temp <- 1 * (dataLR > 0)
#     dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#     dataLR_percent <- 1 * (signif(dataLR_temp[, -1], 1) >= min.percent)
#     Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#
#     if(sum(Prob_percent)==0){
#       Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
#     } else {
#
#       Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                     y = prob %*% cell.type.mat)
#
#       # min.cells.sr
#       prob@x <- rep.int(1, times = length(prob@x))
#       sender.counts <- Matrix::rowSums(prob)
#       receptor.counts <- Matrix::colSums(prob)
#       cells.sr <- cbind(sender.counts, receptor.counts)
#       cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
#       cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
#       cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))
#
#       Prob.avg <- Prob.avg * Prob_percent * cells.sr
#     }
#   } else {
#     Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                   y = prob %*% cell.type.mat)
#   }
#
#   # return
#   Prob.avg <- as.matrix(Prob.avg)
#   return(Prob.avg)
# }




#' Compute group-level cell-cell communication
#'
#' @description
#' Aggregate cell-level communication probabilities (`net$cell$prob`, from
#' \link{computeCommunProb}, optionally \link{filterProbability}) into
#' group-level communication with a permutation test, on the final 11-slot
#' object schema.
#'
#' Plan A (2026-09-28): accelerated group aggregation under the declared numeric
#' `min.percent` gate. The per-boot gated aggregation runs in the
#' `cpp_group_avg_obs` / `cpp_group_avg_perm` kernels. The pinned upstream
#' `format()` character-comparison behavior is a historical comparator, not the
#' current numeric contract. See `docs/COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md`
#' for the compatibility boundary and `docs/PERMUTATION_TEST_AUDIT.md` for the
#' statistical audit of the test itself, which is intentionally out of scope.
#'
#' @param object SpatialCellChat object with `net$cell$prob`.
#' @param group.by cell group information used for computing average communication probabilities.
#' @param avg.type methods for integrating communication probabilities per cell group.
#' @param min.percent Numeric from 0 to 1. The expressed-ligand/receptor proportion is rounded to one significant digit numerically and must meet this threshold for group-level signaling.
#' Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' Default is 5.
#' @param do.permutation whether performing permutation test
#' @param nboot the number of permutations
#' @param seed.use set a random seed. By default, set the seed to 1.
#' @param colocalization.use whether filtering out spatially distant cell groups based on colocalization analysis between any cell groups
#' @param thresh.colo removal of cell-cell communication with no significant colocalizations (fdr < 0.05)
#' @param nthreads OpenMP threads for the permutation kernel (default 1; results are bitwise identical for any thread count).
#' @param verbose Logical. FALSE disables ordinary CLI messages and progress.
#' TRUE follows `options(SpatialCellChat.verbose = 0L)`, `1L`, or `2L` for
#' quiet, normal, or detailed output and never overrides global quiet.
#' Progress also respects progressr's enable policy; noninteractive scripts
#' must opt in, for example with `options(progressr.enable = TRUE)`.
#'
#' @return A SpatialCellChat object with updated slots:
#'
#' `net$group$prob` is the inferred group-level communication probability
#' (SparseChatArray, K x K x nLR; first/second dimensions = source/target group),
#' and `net$group$pval` the corresponding permutation p-values (NULL when
#' `do.permutation = FALSE`). Parameters are stored in `misc$.param$averaging`.
#'
#' USER can access all the inferred cell-cell communications using the function 'subsetCommunication(object)', which returns a data frame.
#' @export
#'
computeAvgCommunProb <- function(
    object,
    group.by = NULL,
    avg.type = c("avg","sum"),
    min.percent = 0.1,
    min.cells.sr = 5,
    do.permutation = TRUE,
    nboot = 100,
    seed.use = 1L,
    colocalization.use = FALSE,
    thresh.colo = 0.05,
    nthreads = 1L,
    verbose = TRUE
){
  t0 <- Sys.time()
  object <- .sc_assert_spatial_cell_chat(object)
  if (isTRUE(verbose)) .cli("computeAvgCommunProb", .type = "subheader")
  avg.type <- match.arg(avg.type)
  do_avg <- avg.type == "avg"

  # ---- group (baseline semantics preserved) ----
  if (is.null(group.by)) {
    group <- object@idents
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `object@meta`, which will be used for cell grouping.",
           call. = FALSE)
    }
    group <- object@meta[[group.by]]
    if (!is.factor(group)) group <- factor(group)
  }
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`",
         call. = FALSE)
  }
  if (isTRUE(verbose)) .cli("Groups: {.val {levels(group)}}",
                            .type = "text", .verbose = 2L)

  # ---- cell-level layers (11-slot input) ----
  if (is.null(object@net$cell$prob)) {
    stop("Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! ",
         call. = FALSE)
  }
  prob.array <- object@net$cell$prob
  prob.cell_ <- unclass(prob.array)
  nC <- dim(prob.array)[1]
  nLR <- dim(prob.array)[3]
  LRsig <- dimnames(prob.array)[[3]]
  if (length(group) != nC)
    stop("Group label count (", length(group), ") != cell count (", nC, "); ",
         "`group.by` must cover exactly the cells of `net$cell$prob`.", call. = FALSE)

  # ---- LR table & per-LR support rows ----
  # Baseline gate inputs are net$tmp$Lavg/Ravg (cached by the pre-migration
  # computeCommunProb), consumed only via `1 * (dataLR > 0)`. The cached Ravg is
  # cofactor-adjusted (Rexpr * coA / coI), but coA/coI factors equal 1 + expr >= 1
  # on non-negative scaled data, so the > 0 pattern is cofactor-invariant and is
  # taken directly from assay$signaling rows (geometricMean > 0 <=> all complex
  # subunits > 0). Gate inputs are therefore bitwise-equal to the baseline.
  LRtab <- object@LR$LRsig
  LRidx <- match(LRsig, rownames(LRtab))
  if (anyNA(LRidx))
    stop("LR pairs missing from `LR$LRsig`: ",
         paste(LRsig[is.na(LRidx)], collapse = ", "), call. = FALSE)
  geneL <- as.character(LRtab$ligand[LRidx])
  geneR <- as.character(LRtab$receptor[LRidx])
  sig <- assay(object, "signaling")
  complex_subunits <- .sc_complex_subunits(object@DB$complex)
  expr.cache <- new.env(parent = emptyenv())
  group_int <- as.integer(group) - 1L

  # min.percent thresholds: one integer threshold per group under the formal
  # numeric gate `signif(mean, 1) >= min.percent`. Group sizes (the mean
  # denominators) are invariant under label permutation, so a single per-group
  # threshold, probed against R's own mean()+signif(), reproduces the declared
  # numeric gate for every boot and LR (see .sc_gate_threshold).
  grp.size <- tabulate(as.integer(group), nbins = numCluster)
  thr <- vapply(grp.size, .sc_gate_threshold, integer(1), min.percent = min.percent)

  # ---- colocalization (baseline semantics, unchanged) ----
  if (colocalization.use) {
    data.spatial <- object@images$coordinates
    pval.colo <- computeColocalization(coordinates = data.spatial, group = group,
                                       nboot = nboot, seed.use = seed.use)
  } else {
    pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster)
  }

  if (isTRUE(verbose)) {
    .cli("Input: {.val {nC}} cells; {.val {numCluster}} groups; {.val {nLR}} LR pairs; method = {.val {avg.type}}",
         .type = "info")
    .cli("Gates: min.percent = {.val {min.percent}}, min.cells.sr = {.val {min.cells.sr}}; colocalization = {.val {colocalization.use}}; permutation = {.val {do.permutation}}",
         .type = "text", .verbose = 2L)
  }

  Prob <- array(0, dim = c(numCluster, numCluster, nLR),
                dimnames = list(levels(group), levels(group), LRsig))
  Pval <- array(1, dim = c(numCluster, numCluster, nLR), dimnames = dimnames(Prob))

  set.seed(seed.use)  # baseline position: identical RNG stream for the permutation draws

  # ---- observed averages (baseline my_future_sapply section) ----
  ccc.counts <- vapply(prob.cell_, function(Mat) length(Mat@x), numeric(1))
  names(ccc.counts) <- LRsig
  LRsig.use.idx <- which(ccc.counts > 0)
  if (length(LRsig.use.idx) < 1) {
    stop("Each LR pair does not have any cell-level links/interactions.", call. = FALSE)
  }

  supp_srin <- function(i) {
    m <- prob.cell_[[i]]
    suppL <- as.numeric(.sc_expr_row(geneL[i], sig, complex_subunits, expr.cache) > 0)
    suppR <- as.numeric(.sc_expr_row(geneR[i], sig, complex_subunits, expr.cache) > 0)
    sr_out <- tabulate(m@i + 1L, nbins = nC)                         # per-sender stored counts
    sr_in  <- tabulate(rep.int(seq_len(nC), diff(m@p)), nbins = nC)  # per-receiver stored counts
    list(suppL = suppL, suppR = suppR, sr_out = sr_out, sr_in = sr_in)
  }

  run_avg <- function(i) {
    s <- supp_srin(i)
    m <- prob.cell_[[i]]
    v <- cpp_group_avg_obs(m@x, m@i, m@p, s$suppL, s$suppR, s$sr_out, s$sr_in,
                           group_int, numCluster, thr, min.cells.sr,
                           if (do_avg) 0L else 1L, nC)
    Prob.avg <- matrix(v, nrow = numCluster, ncol = numCluster)
    Prob.avg[pval.colo > thresh.colo] <- 0
    Prob.avg
  }

  progressr::with_progress({
    pr <- progressr::progressor(along = seq_along(LRsig.use.idx),
                               message = "Averaging LR probabilities",
                               enable = .sc_progress_enabled(verbose) && length(LRsig.use.idx) > 1L &&
                                 getOption("progressr.enable", TRUE))
    for (x in seq_along(LRsig.use.idx)) {
      i <- LRsig.use.idx[[x]]
      Prob[, , i] <- run_avg(i)
      pr()
    }
  })

  # ---- permutation (baseline my_future_lapply section) ----
  ggc.counts <- apply(Prob > 0, 3, sum)
  LRsig.use.idx.perm <- which(ggc.counts > 0)
  if (do.permutation) {
    if (nboot < 1) stop("nboot must be >= 1 when do.permutation = TRUE", call. = FALSE)
    if (isTRUE(verbose)) .cli("Permutation test: {.val {nboot}} boots x {.val {length(LRsig.use.idx.perm)}} active LR pairs; {.val {nthreads}} kernel threads",
                              .type = "text", .verbose = 2L)
    # Baseline RNG: this single replicate() call after set.seed(seed.use) above.
    permutation <- replicate(nboot, sample.int(nC, size = nC))

    run_perm <- function(i) {
      s <- supp_srin(i)
      m <- prob.cell_[[i]]
      # k*k x nboot; column = one boot, sender-group fastest (== baseline
      # matrix(unlist(Pboot), nrow = k*k, ncol = nboot) layout)
      Pboot <- cpp_group_avg_perm(m@x, m@i, m@p, s$suppL, s$suppR, s$sr_out, s$sr_in,
                                  group_int, numCluster, thr, min.cells.sr,
                                  if (do_avg) 0L else 1L, permutation, nthreads)
      Pnull <- as.vector(Prob[, , i])
      nReject <- rowSums(Pboot - Pnull > 0)
      p <- nReject / nboot
      matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    }

    progressr::with_progress({
      pr <- progressr::progressor(along = seq_along(LRsig.use.idx.perm),
                                 message = "Testing LR permutations",
                                 enable = .sc_progress_enabled(verbose) && length(LRsig.use.idx.perm) > 1L &&
                                   getOption("progressr.enable", TRUE))
      for (x in seq_along(LRsig.use.idx.perm)) {
        i <- LRsig.use.idx.perm[[x]]
        Pval[, , i] <- run_perm(i)
        pr()
      }
    })
    Pval[Prob == 0] <- 1
  } else {
    Pval <- NULL
  }

  # ---- write (11-slot: net$group$prob / net$group$pval, params in misc$.param) ----
  group.layers <- lapply(seq_len(nLR), function(i)
    as(matrix(Prob[, , i], nrow = numCluster, ncol = numCluster), "dgCMatrix"))
  names(group.layers) <- LRsig
  prob.group <- SparseChatArray(group.layers)
  dimnames(prob.group) <- list(levels(group), levels(group), LRsig)
  object@net$group$prob <- prob.group

  if (do.permutation) {
    pval.layers <- lapply(seq_len(nLR), function(i)
      as(matrix(Pval[, , i], nrow = numCluster, ncol = numCluster), "dgCMatrix"))
    names(pval.layers) <- LRsig
    pval.group <- SparseChatArray(pval.layers)
    dimnames(pval.group) <- list(levels(group), levels(group), LRsig)
    object@net$group$pval <- pval.group
  } else {
    object@net$group$pval <- NULL
  }

  object@misc$.param$averaging <- list(
    avg.type = avg.type, min.percent = min.percent, min.cells.sr = min.cells.sr,
    do.permutation = do.permutation, nboot = nboot, seed.use = seed.use,
    colocalization.use = colocalization.use, thresh.colo = thresh.colo,
    group.by = group.by, nLR = nLR, numCluster = numCluster, nthreads = nthreads,
    LRsig.CCC.counts = ccc.counts, LRsig.use.idx = LRsig.use.idx,
    LRsig.GGC.counts = ggc.counts, LRsig.use.idx.perm = LRsig.use.idx.perm,
    run.time = as.numeric(Sys.time() - t0, units = "secs")
  )
  object <- .log_operation(object, "computeAvgCommunProb", params = list(
    nLR = nLR, numCluster = numCluster, nboot = nboot, avg.type = avg.type,
    do.permutation = do.permutation, nthreads = nthreads,
    run.time = as.numeric(Sys.time() - t0, units = "secs")))
  if (isTRUE(verbose)) .cli("computeAvgCommunProb done: {.val {numCluster}} x {.val {numCluster}} x {.val {nLR}} group-level layers; {.val {length(LRsig.use.idx.perm)}} active LR pairs; {.val {round(as.numeric(Sys.time() - t0, units = 'secs'), 2)}}s",
                            .type = "success")
  object
}







# #' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
# #'
# #' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
# #' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
# #' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
# #' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
# #' Default is 0.1.
# #' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
# #' Default is 5.
# #'
# #' @return Returns a matrix containing the interaction weights between any two cell groups.
# #' @export
# computeAvgCommunProb_LR <- function (
#     prob,
#     group,
#     dataLR = NULL,
#     min.percent = 0.1,
#     min.cells.sr = 5
# ){
#   # code for test
#   # prob <- prob.cell.i
#   # group <- group
#   # dataLR = dataLR_temp
#   # min.percent = 0.1
#   # min.cells.sr = 5
#
#   # group is vector/factor
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   cell.type.mat <- model.matrix(~group-1)
#
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#
#   if(sum(Prob_percent)==0){
#     Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
#   } else {
#
#     Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                   y = prob %*% cell.type.mat)
#
#     # min.cells.sr
#     prob@x <- rep.int(1, times = length(prob@x))
#     sender.counts <- Matrix::rowSums(prob)
#     receptor.counts <- Matrix::colSums(prob)
#     cells.sr <- cbind(sender.counts, receptor.counts)
#     cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
#     cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
#     cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))
#
#     Prob.avg <- Prob.avg * Prob_percent * cells.sr
#   }
#
#   # return
#   Prob.avg <- as.matrix(Prob.avg)
#   return(Prob.avg)
# }
#


# computeAvgCommunProb_LR <- function (prob, group, dataLR = NULL, min.percent = 0.1, min.cells.sr = 5)
# {
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   level.use <- levels(group)
#   level.use <- level.use[level.use %in% unique(group)]
#   numCluster <- length(level.use)
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#   melt_Prob_percent <- reshape2::melt(Prob_percent,value.name = "Binary.Prob")
#
#   Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
#   # min.cells.sr
#   if(NROW(melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]) == 0){
#     return(Prob.avg)
#   } else {
#     df_nonZero_Prob <- melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]
#     for (x in seq_len(NROW(df_nonZero_Prob))) {
#
#       ii <- df_nonZero_Prob[x,1,drop=T]
#       jj <- df_nonZero_Prob[x,2,drop=T]
#       prob.ij <- prob[group %in% level.use[ii], group %in% level.use[jj], drop = FALSE]
#       if ((sum(Matrix::rowSums(prob.ij > 0) > 0) >= min.cells.sr) &
#           (sum(Matrix::colSums(prob.ij > 0) > 0) >= min.cells.sr)) {
#
#         # Prob.avg[ii, jj] <- sum(prob.ij@x)/(cumprod(prob.ij@Dim)[[2]])
#         Prob.avg[ii, jj] <- sum(prob.ij@x)
#       }
#     } # forloop
#     return(Prob.avg)
#   }
# }


# Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
# for (ii in 1:numCluster) {
#   for (jj in 1:numCluster) {
#     Prob.temp <- as.vector(P.spatial[group %in% level.use[ii], group %in% level.use[jj]]) > 0
#     # Prob.sum <- sum(Prob.temp)
#     # Prob.avg[ii,jj] <- Prob.sum/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#     Prob.avg[ii,jj] <- sum(Prob.temp)/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#   }
# }


# Sum SparseChatArray layers (list of K x K / nC x nC dgCMatrix) into a single
# dgCMatrix.  Mirrors the delayed cpp lookup of marginSums.SparseChatArray.
.sc_sum_layers <- function(layers) {
  if (length(layers) == 1L) return(layers[[1L]])
  cpp_sum_layers <- get("cpp_sum_layers", inherits = TRUE)
  if (!is.function(cpp_sum_layers)) {
    stop("cpp_sum_layers is required to sum SparseChatArray layers", call. = FALSE)
  }
  cpp_sum_layers(layers)
}

#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param net A list from object@net; If net = NULL, net = object@net
#' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group whether to compute the group-level signaling based on the cell group information in `object@idents`
#' @param do.cell whether to compute the individual-cell signaling at signaling pathway level. This works when "prob.cell" exists in `object@net`.
#' @return A CellChat object with updated slot 'netP':
#'
#' object@netP$group$prob is the group-level communication probability array (SparseChatArray) on signaling pathway level; layers are pathways with significant communications, ordered by the descending total communication probability. USER can access all significant interactions using the function \code{\link{subsetCommunication}}
#'
#' object@netP$cell$prob is the individual-cell-level communication probability array (SparseChatArray) on signaling pathway level.
#'
#' object@netP$pathways (group level) and object@netP$pathways.cell (cell level) list all the signaling pathways with significant communications.
#'
#' From version >= 1.1.0, pathways are ordered based on the total communication probabilities. NB: pathways with small total communication probabilities might be also very important since they might be specifically activated between only few cell types.
#'
#' @export

computeCommunProbPathway <- function(
    object = NULL,
    net = NULL,
    pairLR.use = NULL,
    thresh = 0.05,
    do.group = TRUE,
    do.cell = TRUE
) {
  .cli("computeCommunProbPathway", .type = "subheader")
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  pathways <- unique(pairLR.use$pathway_name)
  .cli("Input: {.val {nrow(pairLR.use)}} LR annotations; {.val {length(pathways)}} candidate pathways",
       .type = "info")
  .cli("Levels: group = {.val {do.group}}, cell = {.val {do.cell}}; group p-value threshold = {.val {thresh}}",
       .type = "text", .verbose = 2L)

  if (do.group) {
    if ( !inherits(net$group$prob, "SparseChatArray") ) {
      stop("Please run `computeAvgCommunProb` to compute the group-level signaling! ")
    }
    .cli("Summing {.val {dim(net$group$prob)[3]}} group-level LR layers; significance filter = {.val {!is.null(net$group$pval)}}",
         .type = "text", .verbose = 2L)
    prob.layers <- unclass(net$group$prob)
    # layer names = dimnames[[3]] = rownames(LRsig) (Plan A invariant); layers
    # without a matching pairLR.use row belong to no pathway and are excluded
    layer.pathway <- pairLR.use$pathway_name[match(dimnames(net$group$prob)[[3]], rownames(pairLR.use))]

    # significance gating per layer; net$group$pval NULL = no significance filter
    if (!is.null(net$group$pval)) {
      pval.layers <- unclass(net$group$pval)
      prob.layers <- mapply(function(pm, pv) {
        pm <- as.matrix(pm)
        pm[as.matrix(pv) >= thresh] <- 0
        as(pm, "dgCMatrix")
      }, prob.layers, pval.layers, SIMPLIFY = FALSE)
    }

    # one summed matrix per pathway (keeps pairLR.use order; includes all-zero candidates)
    prob.pathways <- lapply(pathways, function(one_pathway) {
      idx <- which(layer.pathway == one_pathway)
      if (length(idx) == 0L) {
        return(as(matrix(0, nrow(prob.layers[[1]]), ncol(prob.layers[[1]])), "dgCMatrix"))
      }
      .sc_sum_layers(prob.layers[idx])
    })
    names(prob.pathways) <- pathways

    # `prob.sum != 0` means group-group signaling exists in the specific pathway
    prob.sum <- vapply(prob.pathways, function(m) sum(as.matrix(m)), numeric(1))
    pathways.sig <- pathways[prob.sum != 0]

    # sort `prob.pathways` according to group-group communication probability sum
    idx <- sort(prob.sum[pathways.sig], decreasing = TRUE, index.return = TRUE)$ix
    pathways.sig <- pathways.sig[idx]
    prob.pathways.sig <- prob.pathways[pathways.sig]
    group.prob <- SparseChatArray(
      unname(prob.pathways.sig),
      dimnames = list(dimnames(net$group$prob)[[1]],
                      dimnames(net$group$prob)[[2]],
                      pathways.sig)
    )
  } else {
    pathways.sig <- NULL
    group.prob <- NULL
  }


  if (do.cell) {
    if (!is.null(net$cell$prob)) {
      if ( !inherits(net$cell$prob, "SparseChatArray") ) {
        stop("Please run `computeCommunProb` to compute the individual cell-level signaling! ")
      }
      cell.layers <- unclass(net$cell$prob)
      layer.pathway.cell <- pairLR.use$pathway_name[match(dimnames(net$cell$prob)[[3]], rownames(pairLR.use))]

      .cli("Summing {.val {length(cell.layers)}} individual-cell LR layers",
           .type = "text", .verbose = 2L)
      cell.pathways <- lapply(pathways, function(one_pathway) {
        idx <- which(layer.pathway.cell == one_pathway)
        if (length(idx) == 0L) {
          return(as(matrix(0, nrow(cell.layers[[1]]), ncol(cell.layers[[1]])), "dgCMatrix"))
        }
        .sc_sum_layers(cell.layers[idx])
      })
      names(cell.pathways) <- pathways

      # `prob.sum != 0` means cell-cell signaling exists in the specific pathway
      prob.sum <- vapply(cell.pathways, function(m) sum(m@x), numeric(1))
      PathwaySig.use.idx <- which(prob.sum > 0)


      # sort according to the corresponding prob.sum value
      idx <- sort(prob.sum[PathwaySig.use.idx], decreasing = TRUE, index.return = TRUE)$ix
      pathways.sig.cell <- pathways[PathwaySig.use.idx][idx]

      # sort `cell.pathways`
      .cli("Retaining {.val {length(pathways.sig.cell)}} non-zero cell-level pathways in decreasing total-probability order",
           .type = "text", .verbose = 2L)
      cell.prob <- SparseChatArray(
        unname(cell.pathways[pathways.sig.cell]),
        dimnames = list(dimnames(net$cell$prob)[[1]],
                        dimnames(net$cell$prob)[[2]],
                        pathways.sig.cell)
      )
    } else {
      .cli("Cell-level pathway summary skipped: no cell-level probabilities",
           .type = "text", .verbose = 2L)
      pathways.sig.cell <- NULL
      cell.prob <- NULL
    }
  } else {
    pathways.sig.cell <- NULL
    cell.prob <- NULL
  }
  # group-level: pathways; group$prob
  # individual cell-level: pathways.cell; cell$prob
  .cli("computeCommunProbPathway done: {.val {length(pathways.sig)}} group-level and {.val {length(pathways.sig.cell)}} cell-level pathways",
       .type = "success")
  if (is.null(object)) {
    return(list(pathways = pathways.sig,
                group = list(prob = group.prob),
                cell = list(prob = cell.prob),
                pathways.cell = pathways.sig.cell))
  }
  # merge-update: the validator requires netP$cell / netP$group to remain non-NULL lists
  object@netP$pathways <- pathways.sig
  if (is.null(object@netP$group)) object@netP$group <- list()
  object@netP$group$prob <- group.prob
  if (is.null(object@netP$cell)) object@netP$cell <- list()
  object@netP$cell$prob <- cell.prob
  object@netP$pathways.cell <- pathways.sig.cell
  object@netP$tmp <- NULL
  return(object)
}

# computeAvgCommunProb_LR <- function (prob, group, dataLR = NULL, min.percent = 0.1, min.cells.sr = 5)
# {
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   level.use <- levels(group)
#   level.use <- level.use[level.use %in% unique(group)]
#   numCluster <- length(level.use)
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#   melt_Prob_percent <- reshape2::melt(Prob_percent,value.name = "Binary.Prob")
#
#   Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
#   # min.cells.sr
#   if(NROW(melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]) == 0){
#     return(Prob.avg)
#   } else {
#     df_nonZero_Prob <- melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]
#     for (x in seq_len(NROW(df_nonZero_Prob))) {
#
#       ii <- df_nonZero_Prob[x,1,drop=T]
#       jj <- df_nonZero_Prob[x,2,drop=T]
#       prob.ij <- prob[group %in% level.use[ii], group %in% level.use[jj], drop = FALSE]
#       if ((sum(Matrix::rowSums(prob.ij > 0) > 0) >= min.cells.sr) &
#           (sum(Matrix::colSums(prob.ij > 0) > 0) >= min.cells.sr)) {
#
#         # Prob.avg[ii, jj] <- sum(prob.ij@x)/(cumprod(prob.ij@Dim)[[2]])
#         Prob.avg[ii, jj] <- sum(prob.ij@x)
#       }
#     } # forloop
#     return(Prob.avg)
#   }
# }


# Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
# for (ii in 1:numCluster) {
#   for (jj in 1:numCluster) {
#     Prob.temp <- as.vector(P.spatial[group %in% level.use[ii], group %in% level.use[jj]]) > 0
#     # Prob.sum <- sum(Prob.temp)
#     # Prob.avg[ii,jj] <- Prob.sum/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#     Prob.avg[ii,jj] <- sum(Prob.temp)/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#   }
# }




#' Filter cell-cell communication if there are only few number of cells in certain cell groups or only few interactions
#'
#' @param object CellChat object
#' @param min.cells the minimum number of cells required in each cell group for filtering cell group-level communication
#' @param min.links the minimum number of links/interactions required in the ligand-receptor pair for filtering individual cell-level communication
#' @param min.cells.sr the minimum number of cells required as senders or receivers for filtering individual cell-level communication
#' @return CellChat object with an updated slot net
#' @export
#'
filterCommunication <- function(object, min.cells = 10, min.links = 5, min.cells.sr = 5) {
  .cli("filterCommunication", .type = "subheader")
  .cli("Input: {.val {length(object@idents)}} cells; {.val {nlevels(object@idents)}} groups; group filter = {.val {!is.null(min.cells)}}, cell filter = {.val {!is.null(min.links) || !is.null(min.cells.sr)}}",
       .type = "info")
  n.groups.excluded <- 0L
  n.lr.links.excluded <- 0L
  n.lr.cells.excluded <- 0L

  if (!is.null(min.cells)) {
    .cli("Group-level filter: min.cells = {.val {min.cells}}",
         .type = "text", .verbose = 2L)
    group.prob <- object@net$group$prob
    if (!inherits(group.prob, "SparseChatArray")) {
      stop(
        "Please run `computeAvgCommunProb` to compute the group-level communication!",
        call. = FALSE
      )
    }
    group.pval <- object@net$group$pval
    if (!is.null(group.pval) && !inherits(group.pval, "SparseChatArray")) {
      stop("net$group$pval must be a SparseChatArray when present", call. = FALSE)
    }
    if (!is.null(group.pval) &&
        !identical(as.integer(dim(group.pval)), as.integer(dim(group.prob)))) {
      stop("net$group$pval dimensions must match net$group$prob", call. = FALSE)
    }

    cell.excludes <- which(as.numeric(table(object@idents)) < min.cells)
    n.groups.excluded <- length(cell.excludes)
    if (length(cell.excludes) > 0L) {
      group.levels <- levels(object@idents)
      group.names <- dimnames(group.prob)[[1L]]
      excluded.names <- group.levels[cell.excludes]
      group.excludes <- match(excluded.names, group.names)
      if (anyNA(group.excludes)) {
        stop(
          "net$group$prob group names do not match levels(object@idents)",
          call. = FALSE
        )
      }

      .cli("Groups excluded for insufficient cells: {.val {excluded.names}}",
           .type = "text", .verbose = 2L)

      group.prob <- .sc_map_sparse_layers(group.prob, function(layer) {
        layer[group.excludes, ] <- 0
        layer[, group.excludes] <- 0
        Matrix::drop0(layer)
      })
      object@net$group$prob <- group.prob

      if (!is.null(group.pval)) {
        pval.layers <- Map(
          f = function(pval.layer, prob.layer) {
            pval.layer[prob.layer == 0] <- 1
            Matrix::drop0(pval.layer)
          },
          unclass(group.pval),
          unclass(group.prob)
        )
        object@net$group$pval <- .new_SparseChatArray(
          pval.layers,
          dimnames(group.pval)
        )
      }
    }
    rm(group.prob, group.pval)
    gc()
  }

  if (!is.null(min.links) || !is.null(min.cells.sr)) {
    .cli("Filtering individual-cell LR layers", .type = "text", .verbose = 2L)
    prob.cell <- object@net$cell$prob
    if (!inherits(prob.cell, "SparseChatArray")) {
      stop(
        "Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! ",
        call. = FALSE
      )
    }
    prob.cell_ <- unclass(prob.cell)

    prob.sum <- vapply(
      prob.cell_,
      function(mat) length(mat@x),
      numeric(1)
    )
    dimArr <- dim(prob.cell)
    .cli("Cell-level input: {.val {dimArr[3]}} LR layers; {.val {dimArr[1]}} cells",
         .type = "text", .verbose = 2L)

    AllzeroMat <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      repr = "C",
      dims = dimArr[c(1L, 2L)],
      index1 = TRUE
    )

    if (!is.null(min.links)) {
      .cli("Link filter: min.links = {.val {min.links}}", .type = "text", .verbose = 2L)
      idx.signaling.excludes <- which((prob.sum < min.links) & (prob.sum > 0))
      n.lr.links.excluded <- length(idx.signaling.excludes)
      if (length(idx.signaling.excludes) > 0L) {
        .cli("Excluded {.val {n.lr.links.excluded}} LR pairs for insufficient links",
             .type = "text", .verbose = 2L)
        prob.cell_[idx.signaling.excludes] <- rep(
          list(AllzeroMat),
          length(idx.signaling.excludes)
        )
      }
    }

    if (!is.null(min.cells.sr)) {
      .cli("Sender/receiver filter: min.cells.sr = {.val {min.cells.sr}}",
           .type = "text", .verbose = 2L)
      if (is.null(min.links)) {
        pathway.idx <- which(prob.sum > 0)
      } else {
        pathway.idx <- which(prob.sum >= max(1, min.links))
      }

      if (length(pathway.idx) > 0L) {
        pathways.remove.idx <- progressr::with_progress({
          pr <- progressr::progressor(along = pathway.idx,
                                     message = "Filtering sender/receiver support",
                                     enable = .sc_progress_enabled() && length(pathway.idx) > 1L &&
                                       getOption("progressr.enable", TRUE))
          lapply(
            X = pathway.idx,
            FUN = function(x) {
              prob.cell.i <- prob.cell_[[x]] > 0
              excluded <-
                (sum(Matrix::rowSums(prob.cell.i) > 0) < min.cells.sr) |
                (sum(Matrix::colSums(prob.cell.i) > 0) < min.cells.sr)
              result <- if (excluded) x else NULL
              pr()
              result
            }
          )
        })
        pathways.remove.idx <- unlist(pathways.remove.idx, use.names = FALSE)
        n.lr.cells.excluded <- length(pathways.remove.idx)

        if (length(pathways.remove.idx) > 0L) {
          .cli("Excluded {.val {n.lr.cells.excluded}} LR pairs for insufficient sending/receiving cells",
               .type = "text", .verbose = 2L)
          prob.cell_[pathways.remove.idx] <- rep(
            list(AllzeroMat),
            length(pathways.remove.idx)
          )
        }
      }
    }

    object@net$cell$prob <- .new_SparseChatArray(
      prob.cell_,
      dimnames(prob.cell)
    )
  }

  .cli("filterCommunication done: excluded {.val {n.groups.excluded}} groups; {.val {n.lr.links.excluded}} LR pairs by link count and {.val {n.lr.cells.excluded}} by sender/receiver support",
       .type = "success")
  object
}



#' @title scMatrixTruncation
#' @description
#' Single cell data is usually saved as "dgCMatrix" or "dgTMatrix" in R,
#' use this function to cut off Expression Matrix's values to make
#' Mat'v >= cutoff.
#'
#' @param Mat dgCMatrix or dgTMatrix
#' @param cutoff Numeric. Cut off Expression Matrix's values to make
#' Mat' values >= cutoff.
#' @param remain.cutoff.v Boolean. Use `>=` or `>`.
#' @param repr Character string. One of "C", "T", specifying the representation of the sparse matrix result.
#'
#' @return dgCMatrix or dgTMatrix
#' @export
scMatrixTruncation <- function(Mat,cutoff=NULL,remain.cutoff.v=T,repr = c("C", "T")){
  if (!inherits(x = Mat, what = c("dgCMatrix","dgTMatrix"))) {
    stop("Please use `dgCMatrix` or `dgTMatrix` as a Matrix input!")
  }

  if((!is.null(cutoff)) & is.numeric(cutoff)){
    if(length(Mat@x) > 0){
      Mat <- as(Mat,Class = "TsparseMatrix")

      if(remain.cutoff.v){
        remain.idx <- which(Mat@x >= cutoff)
      } else{
        remain.idx <- which(Mat@x > cutoff)
      }

      # update
      Mat@i <- Mat@i[remain.idx]
      Mat@j <- Mat@j[remain.idx]
      Mat@x <- Mat@x[remain.idx]
    }
  }

  repr <- match.arg(repr)
  if(repr=="C"){
    Mat <- as(Mat,Class = "CsparseMatrix")
  } else if (repr=="T"){
    Mat <- as(Mat,Class = "TsparseMatrix")
  }
  return(Mat)
}



#' @title filterProbability
#' @description
#' Filter out statistically non-significant communication probability at the level of individual cells after running \link{computeCommunProb}
#'
#' @param object CellChat object
#' @param nboot Numeric. The number of bootstrap samples, 100 by default.
#' @param seed.use Integer. The random seed used when taking a sample from the communication probabilities
#' @param thresh Numeric. The threshold for defining significant individual cell-cell communication at (1-thresh) of a shuffled distribution of each L-R pair.
#'
#' @return CellChat object
#' @export
filterProbability <- function (
    object,
    nboot = 100,
    seed.use = 666L,
    thresh = 0.05
){
  .cli("filterProbability", .type = "subheader")
  quantile.prob <- 1-thresh
  if(quantile.prob==0){
    .cli("filterProbability skipped: quantile is zero; probabilities unchanged", .type = "success")
    return(object)
  } else { # quantile.prob<1
    comm.param <- object@misc$.param$communication
    spill.dir <- if (is.null(comm.param)) NULL else comm.param$spill.dir
    if (is.null(object@net$cell$prob) && is.null(spill.dir)) {
      stop(
        "Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! "
      )
    } else {
      # 层来源：内存模式读 net$cell$prob（SparseChatArray 即命名 list）；spill 模式逐层读盘并消费
      if (!is.null(spill.dir)) {
        pair.LR.use <- comm.param$layer.names
        if (is.null(pair.LR.use))
          stop("Spill mode requires `layer.names` in `misc$.param$communication`", call. = FALSE)
        missing.files <- pair.LR.use[!file.exists(file.path(spill.dir, paste0(pair.LR.use, ".rds")))]
        if (length(missing.files))
          stop("Spilled layer files missing for LR pairs: ",
               paste(missing.files, collapse = ", "), call. = FALSE)
        prob.cell_ <- NULL
      } else {
        prob.cell_ <- unclass(object@net$cell$prob)
        pair.LR.use <- dimnames(object@net$cell$prob)[[3]]
      }
      cell.names <- colnames(assay(object, "signaling"))

      d.spatial <- object@images$.distance$d.spatial
      Matrix::diag(d.spatial) <- 1
      adj.contact <- object@images$.distance$adj.contact

      nLR <- length(pair.LR.use)
      nLR1 <- comm.param$nLR1
      nC <- NROW(d.spatial)
      .cli("Input: {.val {nC}} cells; {.val {nLR}} LR pairs; {.val {nboot}} sampled cells per LR; quantile = {.val {quantile.prob}}",
           .type = "info")
      .cli("Layer source: {.val {if (is.null(spill.dir)) 'memory' else spill.dir}}; seed = {.val {seed.use}}",
           .type = "text", .verbose = 2L)

      set.seed(seed.use)
      if (comm.param$all.contact.dependent == TRUE) {
        d.spatial <- adj.contact
      }
      # v2 优化①②：模式行计数与 d 侧模式索引每支只构建一次（原逐层 rowSums/单行提取）。
      # tabulate(@i[x != 0]) 与 rowSums(x != 0) 计数逐位一致；radix 稳定 order 保证行内列升序
      rows.branch1 <- tabulate(d.spatial@i[d.spatial@x != 0] + 1L, nbins = nC)
      rows.branch2 <- tabulate(adj.contact@i[adj.contact@x != 0] + 1L, nbins = nC)
      csr.branch1 <- .sc_csr_build(d.spatial)
      csr.branch2 <- .sc_csr_build(adj.contact)
      pat.cache.1 <- new.env(hash = TRUE)
      pat.cache.2 <- new.env(hash = TRUE)

      # dim(permutation) = nboot x nLR
      permutation <- replicate(nLR, base::sample(x = 1:nC, size = nboot, replace = F))

      prob.cell_ <- my_future_lapply(X = 1:nLR, FUN = function(i) {

        if (i <= nLR1) {
          d_spatial <- d.spatial
        } else {
          d_spatial <- adj.contact
        }
        sample.cells <- permutation[ ,i,drop=T]
        Prob.cell.i <- if (!is.null(spill.dir)) {
          readRDS(file.path(spill.dir, paste0(pair.LR.use[i], ".rds")))
        } else {
          prob.cell_[[i]]
        }

        # 高效分位捷径：基线把各采样源在 pattern 列上的稠密行拼接后取 type-7 分位。
        # 该分位为 0 当且仅当全零前缀长度 n0 满足 n0 >= ceiling(1 + (|S|-1) * quantile.prob)。
        # 零值数下界 = |S| - 层非零值数（O(1) 判定）；不满足时再行级精确判定，
        # 两者均命中则直接返回原层（与基线 nboot.quantile == 0 分支逐位等价），
        # 未命中才构造完整采样向量走基线 quantile 路径。
        pattern.rows <- if (i <= nLR1) rows.branch1 else rows.branch2
        total.sampled <- sum(pattern.rows[sample.cells])
        cutoff.zero <- FALSE
        if (total.sampled > 0) {
          h.cut <- ceiling(1 + (total.sampled - 1) * quantile.prob)
          if (total.sampled - sum(Prob.cell.i@x != 0) >= h.cut) {
            cutoff.zero <- TRUE
          } else {
            sub.x <- Prob.cell.i[sample.cells, , drop = FALSE]@x
            if (total.sampled - sum(sub.x != 0) >= h.cut) cutoff.zero <- TRUE
          }
        }
        if (cutoff.zero) {
          # nboot.quantile == 0：基线直接返回原层
          return(Prob.cell.i)
        }

        # v2 优化③：P 侧手工 CSR + d 侧模式缓存。与基线 purrr::map 逐位同构：
        # 同置换顺序 × prob.index 升序（which 语义）×（列,值）对放置，缺失列补 0
        cP <- .sc_csr_build(Prob.cell.i)
        cD <- if (i <= nLR1) csr.branch1 else csr.branch2
        pats <- if (i <= nLR1) pat.cache.1 else pat.cache.2
        sample.prob.cell.i <- unlist(lapply(seq_along(sample.cells), function(jj) {
          cell.index <- sample.cells[jj]
          prob.index <- .sc_pattern_row(cD, d_spatial, cell.index, pats)
          vals <- numeric(length(prob.index))
          s <- cP$rp[cell.index]; e <- cP$rp[cell.index + 1L]
          if (e > s) {
            seg <- cP$ord[(s + 1L):e]
            hit <- match(cP$jc[seg], prob.index)
            keep <- !is.na(hit)
            vals[hit[keep]] <- Prob.cell.i@x[seg][keep]
          }
          vals
        }), use.names = FALSE)
        nboot.quantile <- quantile(sample.prob.cell.i, probs = quantile.prob)
        gc()

        if(nboot.quantile == 0){
          # sparse enough, return directly
          return(Prob.cell.i)
        } else if (nboot.quantile > 0) {
          # filter `Prob.cell.i` to make it sparse enough
          Prob.cell.i <- scMatrixTruncation(Prob.cell.i,cutoff = nboot.quantile,remain.cutoff.v = T,repr = "C")
        }

        return(Prob.cell.i)
      }, future.globals = future::nbrOfWorkers() != 1L, simplify = F)
    names(prob.cell_) <- pair.LR.use
    prob.cell <- SparseChatArray(prob.cell_)
    dimnames(prob.cell) <- list(cell.names, cell.names, pair.LR.use)
    object@net$cell$prob <- prob.cell

    if (!is.null(spill.dir)) {
      rds.files <- list.files(spill.dir, pattern = "\\.rds$", full.names = TRUE)
      unlink(rds.files)
      if (length(list.files(spill.dir, all.files = TRUE, no.. = TRUE)) == 0)
        unlink(spill.dir, recursive = TRUE)
      comm.param$spill.dir <- NULL
      comm.param$spill.backend <- NULL
      comm.param$spill.count <- NULL
      comm.param$layer.names <- NULL
      object@misc$.param$communication <- comm.param
      .cli("Consumed spilled LR layers and cleared spill metadata", .type = "text", .verbose = 2L)
    }
    object <- .log_operation(object, "filterProbability", params = list(
      nLR = nLR, nboot = nboot, thresh = thresh,
      spill.consumed = !is.null(spill.dir)))

    .cli("filterProbability done: {.val {nLR}} LR layers retained at quantile {.val {quantile.prob}}",
         .type = "success")
    return(object)
    } # whether probability computed
  } # whether to filter out
}


#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param sources.use,targets.use,signaling,pairLR.use Please check the description in function \code{\link{subsetCommunication}}
#' @param remove.isolate whether removing the isolate cell groups without any interactions when applying \code{\link{subsetCommunication}}
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.object whether return an updated CellChat object
#' @importFrom  dplyr group_by summarize groups
#' @importFrom stringr str_split
#'
#' @return Return an updated CellChat object:
#'
#' `object@net$group$count` is a matrix: rows and columns are sources and targets respectively, and elements are the number of significant ligand-receptor interactions between any two cell groups. USER can convert a matrix to a data frame using the function `reshape2::melt()`
#'
#' `object@net$group$weight` is also a matrix containing the interaction weights between any two cell groups
#'
#' `object@net$group$LR.sig` lists the ligand-receptor pairs with significant communications
#'
#' `object@net$cell$count` / `object@net$cell$weight` are the individual-cell-level count and weight matrices; `object@net$cell$LR.sig` lists the significant ligand-receptor pairs
#'
#' @export
#'
aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
  .cli("aggregateNet", .type = "subheader")
  net <- object@net
  .cli("Input: {.val {dim(net$group$prob)[1]}} groups; {.val {dim(net$group$prob)[3]}} group-level LR layers; p-value threshold = {.val {thresh}}",
       .type = "info")
  if (is.null(sources.use) && is.null(targets.use) && is.null(signaling) && is.null(pairLR.use)) {
    if ( !inherits(net$group$prob, "SparseChatArray") ) {
      stop("Please run `computeAvgCommunProb` to compute the group-level signaling! ")
    }
    prob.layers <- unclass(net$group$prob)
    pval.layers <- if (!is.null(net$group$pval)) unclass(net$group$pval) else NULL
    K <- nrow(prob.layers[[1]])
    weight.acc <- matrix(0, nrow = K, ncol = K)
    count.acc <- matrix(0, nrow = K, ncol = K)
    layer.sum <- numeric(length(prob.layers))
    for (k in seq_along(prob.layers)) {
      pm <- as.matrix(prob.layers[[k]])
      if (!is.null(pval.layers)) {
        vm <- as.matrix(pval.layers[[k]])
        vm[pm == 0] <- 1
        pm[vm >= thresh] <- 0
      }
      weight.acc <- weight.acc + pm
      count.acc <- count.acc + (pm != 0)
      layer.sum[k] <- sum(pm)
    }
    dn.group <- dimnames(net$group$prob)[1:2]
    dimnames(weight.acc) <- dn.group
    dimnames(count.acc) <- dn.group
    net$group$count <- as(count.acc, "dgCMatrix")
    net$group$weight <- as(weight.acc, "dgCMatrix")
    net$group$LR.sig <- dimnames(net$group$prob)[[3]][layer.sum > 0]
  } else {
    if ( !inherits(net$group$prob, "SparseChatArray") ) {
      stop("Please run `computeAvgCommunProb` to compute the group-level signaling! ")
    }
    .cli("Subsetting group-level sources, targets, or signaling before aggregation",
         .type = "text", .verbose = 2L)
    prob.layers <- unclass(net$group$prob)
    pval.layers <- if (!is.null(net$group$pval)) unclass(net$group$pval) else NULL
    LR <- object@LR$LRsig
    dn.group <- dimnames(net$group$prob)[1:2]
    layer.names <- dimnames(net$group$prob)[[3]]
    # natively replicate subsetCommunication(slot.name = "net") on the SparseChatArray
    rows <- lapply(seq_along(prob.layers), function(k) {
      lname <- layer.names[k]
      pm <- as.matrix(prob.layers[[k]])
      if (!is.null(pval.layers)) {
        pm[as.matrix(pval.layers[[k]]) >= thresh] <- 0
      }
      idx <- which(pm > 0, arr.ind = TRUE)
      if (nrow(idx) == 0L) return(NULL)
      data.frame(
        source = dn.group[[1]][idx[, 1]],
        target = dn.group[[2]][idx[, 2]],
        interaction_name = lname,
        pathway_name = LR$pathway_name[match(lname, rownames(LR))],
        prob = pm[idx],
        pval = if (!is.null(pval.layers)) as.matrix(pval.layers[[k]])[idx] else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    df.net <- do.call(rbind, rows)
    if (is.null(df.net)) {
      df.net <- data.frame(source = character(), target = character(),
                           interaction_name = character(), pathway_name = character(),
                           prob = numeric(), pval = numeric(),
                           stringsAsFactors = FALSE)
    }
    if (!is.null(signaling)) {
      df.net <- subset(df.net, pathway_name %in% signaling)
    }
    if (!is.null(pairLR.use)) {
      if ("interaction_name" %in% colnames(pairLR.use)) {
        df.net <- subset(df.net, interaction_name %in% pairLR.use$interaction_name)
      } else if ("pathway_name" %in% colnames(pairLR.use)) {
        df.net <- subset(df.net, pathway_name %in% pairLR.use$pathway_name)
      }
    }
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- levels(object@idents)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- levels(object@idents)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (nrow(df.net) == 0) {
      stop("No significant signaling interactions are inferred based on the input!")
    }
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      .cli("Removing isolated groups; set `remove.isolate = FALSE` to retain them",
           .type = "text", .verbose = 2L)
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }

    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$group$count <- count
    net$group$weight <- prob
    net$group$weight[is.na(net$group$weight)] <- 0
    net$group$count[is.na(net$group$count)] <- 0
  }

  if (!is.null(net$cell$prob)) {
    if ( !inherits(net$cell$prob, "SparseChatArray") ) {
      stop("Please run `computeCommunProb` to compute the individual cell-level signaling! ")
    }
    cell.layers <- unclass(net$cell$prob)
    # binarize the stored entries: explicit zeros still count as links (baseline quirk)
    net$cell$count <- .sc_sum_layers(lapply(cell.layers, function(m) { m@x <- rep.int(1, length(m@x)); m }))
    net$cell$weight <- .sc_sum_layers(cell.layers)
    dn.cell <- dimnames(net$cell$prob)[1:2]
    dimnames(net$cell$count) <- dn.cell
    dimnames(net$cell$weight) <- dn.cell
    net$cell$LR.sig <- dimnames(net$cell$prob)[[3]][vapply(cell.layers, function(m) sum(m@x), numeric(1)) > 0]
    if (!is.null(sources.use) | !is.null(targets.use) | !is.null(signaling) | !is.null(pairLR.use)) {
      .cli("Source, target, and signaling subsets apply only to group-level aggregation; individual-cell layers remain unsubsetted",
           .type = "text", .verbose = 2L)
    }
  }

  .cli("aggregateNet done: {.val {nrow(net$group$count)}} source x {.val {ncol(net$group$count)}} target groups; cell-level network = {.val {!is.null(net$cell$prob)}}",
       .type = "success")

  if (return.object) {
    object@net <- net
    # object@net$tmp <- NULL
    return(object)
  } else {
    return(net)
  }
}




#' Compute averaged expression values for each cell group
#'
#' @param object CellChat object
#' @param features a char vector giving the used features. default use all features
#' @param group.by cell group information; default is `object@idents` when input is a single object and `object@idents` when input is a merged object; otherwise it should be one of the column names of the meta slot
#' @param type methods for computing the average gene expression per cell group.
#'
#' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
#'
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function \code{\link[base]{mean}}.
#'
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param slot.name the data in the slot.name to use
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the 'slot.name' is used
#'
#' @return Returns a matrix with genes as rows, cell groups as columns.

#' @export
#'
computeAveExpr <- function(object, features = NULL, group.by = NULL, type = c("triMean", "truncatedMean", "median"), trim = NULL,
                           slot.name = c("data.signaling", "data"), data.use = NULL) {
  type <- match.arg(type)
  slot.name <- match.arg(slot.name)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  if (is.null(data.use)) {
    data.use <- slot(object, slot.name)
  }
  if (is.null(features)) {
    features.use <- row.names(data.use)
  } else {
    features.use <- intersect(features, row.names(data.use))
  }
  data.use <- data.use[features.use, , drop = FALSE]
  data.use <- as.matrix(data.use)

  if (is.null(group.by)) {
    labels <- object@idents
    if (!is.factor(labels)) {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents
    }
  } else {
    labels <- object@meta[[group.by]]
  }
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }
  # compute the average expression per group
  data.use.avg <- aggregate(t(data.use), list(labels), FUN = FunMean)
  data.use.avg <- t(data.use.avg[,-1])
  rownames(data.use.avg) <- features.use
  colnames(data.use.avg) <- levels(labels)
  return(data.use.avg)
}



#' Compute the expression of complex in individual cells using geometric mean
#' @param complex_input the complex_input from CellChatDB
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex the names of complex
#' @return
#' @importFrom dplyr select starts_with
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_complex <- function(complex_input, data.use, complex) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))

  data.complex = my_future_sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      return(geometricMean(data.use[RsubunitsV,,drop=F]))
    },
  )
  data.complex <- t(data.complex)
  return(data.complex)
}


#' Compute the expression of ligands or receptors using geometric mean
#' @param geneLR a char vector giving a set of ligands or receptors
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex_input the complex_input from CellChatDB
# #' @param group a factor defining the cell groups; If NULL, compute the expression of ligands or receptors in individual cells; otherwise, compute the average expression of ligands or receptors per cell group
# #' @param FunMean the function for computing average expression per cell group
#' @return
#' @export
computeExpr_LR <- function(geneLR, data.use, complex_input){
  nLR <- length(geneLR)
  numCluster <- ncol(data.use)
  index.singleL <- which(geneLR %in% rownames(data.use))
  dataL1avg <- data.use[geneLR[index.singleL],]
  dataLavg <- matrix(nrow = nLR, ncol = numCluster)
  dataLavg[index.singleL,] <- dataL1avg
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneLR[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataLavg[index.complexL,] <- data.complex
  }
  return(dataLavg)
}










#' Compute the geometric mean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return Numeric.
#' @export
geometricMean <- function(x,na.rm=TRUE){
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=na.rm))
  } else {
    exp(apply(log(x),2L,mean,na.rm=na.rm))
  }
}


# arithmeticMean <- function(x,na.rm=TRUE){
#   if (is.null(nrow(x))) {
#     BiocGenerics::mean(x,na.rm = na.rm)
#   } else {
#     apply(x,2, BiocGenerics::mean ,na.rm=na.rm)
#   }
# }


#' Compute the Tukey's trimean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @importFrom stats quantile
#' @export
triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}

#' Compute the average expression per cell group when the percent of expressing cells per cell group larger than a threshold
#' @param x a numeric vector
#' @param trim the percent of expressing cells per cell group to be considered as zero
#' @param na.rm whether remove na
#' @return
#' @importFrom Matrix nnzero
#' @export
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}

#' Identify all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' @param object CellChat object
#' @param from a vector giving the index or the name of source cell groups
#' @param to a corresponding vector giving the index or the name of target cell groups. Note: The length of 'from' and 'to' must be the same, giving the corresponding pair of cell groups for communication.
#' @param bidirection whether show the bidirectional communication, i.e., both 'from'->'to' and 'to'->'from'.
#' @param pair.only whether only return ligand-receptor pairs without pathway names and communication strength
#' @param pairLR.use0 ligand-receptor pairs to use; default is all the significant interactions
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
identifyEnrichedInteractions <- function(object, from, to, bidirection = FALSE, pair.only = TRUE, pairLR.use0 = NULL, thresh = 0.05){
  pairwiseLR <- object@net$pairwiseRank
  if (is.null(pairwiseLR)) {
    stop("The interactions between pairwise cell groups have not been extracted!
         Please first run `object <- rankNetPairwise(object)`")
  }
  group.names.all <- names(pairwiseLR)
  if (!is.numeric(from)) {
    from <- match(from, group.names.all)
    if (sum(is.na(from)) > 0) {
      message("Some input cell group names in 'from' do not exist!")
      from <- from[!is.na(from)]
    }
  }
  if (!is.numeric(to)) {
    to <- match(to, group.names.all)
    if (sum(is.na(to)) > 0) {
      message("Some input cell group names in 'to' do not exist!")
      to <- to[!is.na(to)]
    }
  }
  if (length(from) != length(to)) {
    stop("The length of 'from' and 'to' must be the same!")
  }
  if (bidirection) {
    from2 <- c(from, to)
    to <- c(to, from)
    from <- from2
  }
  if (is.null(pairLR.use0)) {
    k <- 0
    pairLR.use0 <- list()
    for (i in 1:length(from)){
      pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
      idx <- pairwiseLR_ij$pval < thresh
      if (length(idx) > 0) {
        k <- k +1
        pairLR.use0[[k]] <- pairwiseLR_ij[idx,]
      }
    }
    pairLR.use0 <- do.call(rbind, pairLR.use0)
  }

  k <- 0
  pval <- matrix(nrow = length(rownames(pairLR.use0)), ncol = length(from))
  prob <- pval
  group.names <- c()
  for (i in 1:length(from)) {
    k <- k+1
    pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
    pairwiseLR_ij <- pairwiseLR_ij[rownames(pairLR.use0),]
    pval_ij <- pairwiseLR_ij$pval
    prob_ij <- pairwiseLR_ij$prob
    pval_ij[pval_ij > 0.05] = 1
    pval_ij[pval_ij > 0.01 & pval_ij <= 0.05] = 2
    pval_ij[pval_ij <= 0.01] = 3
    prob_ij[pval_ij ==1] <- 0
    pval[,k] <- pval_ij
    prob[,k] <- prob_ij
    group.names <- c(group.names, paste(group.names.all[from[i]], group.names.all[to[i]], sep = " - "))
  }
  prob[which(prob == 0)] <- NA
  # remove rows that are entirely NA
  pval <- pval[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  pairLR.use0 <- pairLR.use0[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  prob <- prob[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  if (pair.only) {
    pairLR.use0 <- dplyr::select(pairLR.use0, ligand, receptor)
  }
  return(pairLR.use0)
}


#' Compute the region distance based on the spatial locations of each splot/cell of the spatial transcriptomics
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param group a factor vector defining the regions/labels of each cell/spot
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before computing the average distance per cell group.
#' @param interaction.length The maximum interaction/diffusion length of ligands. This hard threshold is used to filter out the connections between spatially distant regions
#' @param spot.size theoretical spot size; e.g., 10x Visium (spot.size = 65 microns)
#' @param spot.size.fullres The number of pixels that span the diameter of a theoretical spot size in the original,full-resolution image.
#' @param k.min the minimum number of interacting cell pairs required for defining adjacent cell groups
# #' @param k.spatial Number of neighbors in a knn graph, which is used to filter out the connections between spatially distant regions that do not share many neighbor spots/cells
#' @importFrom BiocNeighbors queryKNN KmknnParam
#' @return A square matrix giving the pairwise region distance
#'
#' @export
computeRegionDistance <- function(coordinates, group, trim = 0.1,
                                  interaction.length = NULL, spot.size = NULL, spot.size.fullres = NULL, k.min = 10
) {
  if (ncol(coordinates) != 2) {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  if (!is.factor(group)) {
    stop("Please input the `group` as a factor!")
  }
  # type <- match.arg(type)
  type <- "truncatedMean"
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use <- level.use[level.use %in% unique(group)]
  d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
  adj.spatial <- matrix(0, nrow = numCluster, ncol = numCluster)
  for (i in 1:numCluster) {
    for (j in 1:numCluster) {
      data.spatial.i <- coordinates[group %in% level.use[i], , drop = FALSE]
      data.spatial.j <- coordinates[group %in% level.use[j], , drop = FALSE]
      qout <- suppressWarnings(BiocNeighbors::queryKNN(data.spatial.j, data.spatial.i, k = 1, BNPARAM = BiocNeighbors::KmknnParam(), get.index = TRUE))
      if (!is.null(spot.size) & !is.null(spot.size.fullres)) {
        qout$distance <- qout$distance*spot.size/spot.size.fullres
        idx <- qout$distance - interaction.length < spot.size/2
        adj.spatial[i,j] <- (length(unique(qout$index[idx])) >= k.min) * 1
      }
      d.spatial[i,j] <- FunMean(qout$distance) # since distances are positive values, different ways for computing the mean have little effects.

    }
  }
  d.spatial <- (d.spatial + t(d.spatial))/2
  if (!is.null(spot.size) & !is.null(spot.size.fullres)) {
    adj.spatial <- adj.spatial * t(adj.spatial) # if one is zero, then both are zeros.
    adj.spatial[adj.spatial == 0] <- NaN
    d.spatial <- d.spatial * adj.spatial
  }

  rownames(d.spatial) <- levels(group); colnames(d.spatial) <- levels(group)
  return(d.spatial)

}

