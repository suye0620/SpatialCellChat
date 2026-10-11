# Agent Note: centrality/subsetCommunication 11-slot migration + reticulate removal

Status: implemented
Governance: v1
Date: 2026-10-10
Decision type: migration
Scope: R/analysis.R netAnalysis_computeCentrality + subsetCommunication(_internal); consumers in visualization.R, utilities.R updateClusterLabels, netMappingDEG
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

Wave 3a of the 11-slot migration. Two analysis-layer functions still speak the legacy
object schema and are broken on real project objects:

1. `netAnalysis_computeCentrality` (R/analysis.R L82-201) reads legacy slots
   `slot$prob`/`slot$pval`/`slot$tmp$prob.cell`; the final 11-slot schema stores
   `net$group$prob`/`net$group$pval`/`net$cell$prob` and `netP$group$prob`/`netP$cell$prob`
   (SparseChatArray, see implemented Note 2026-08-15-final-11-slot-schema). It also
   contains the only reticulate call site in analysis.R (L184,
   `reticulate::array_reshape(centr.all, ..., order = "F")`); reticulate init measured
   146-175 s on this workstation with one Python-discovery halt. `order = "F"` is exactly
   R column-major order, so the call is replaceable by `array(as.vector(x), dim)`.
2. `subsetCommunication`/`subsetCommunication_internal` (L2266-2365, L2370-2530):
- BUG-B (crash): wrapper accessed nonexistent `object@options$mode`; new-schema objects store mode in `misc$.mode`, so the old access errored before dispatch.
   - BUG-A (crash): internal join does `dplyr::select(LR, ..., "evidence")` but LRsig has
     no `evidence` column -> every net/netP path on real objects errors
     ("Column `evidence` doesn't exist").
   - BUG-C (crash): netP path melts a NULL pval (`reshape2::melt(NULL)`) -> names-attribute error.
   - BUG-D (wrong join): netP join matches layer names against `rownames(LR)` (interaction
     names); pathway-level layers are pathway names, so multi-LR pathways yield all-NA rows.
   - BUG-E (misleading error): out-of-range numeric `sources.use`/`targets.use` silently
     index to NA and surface as "No significant signaling interactions are inferred!".
   - netP gating: legacy netP objects carry no pval; `pval[prob==0] <- 1` on NULL is a
     silent no-op (empirically verified), so legacy effectively performs no gating at netP
     level. `computeCommunProbPathway` already gates `netP$group$prob` by thresh.

## Decision

Rewrite both functions against the final 11-slot schema, preserving observable behavior
on legal inputs (discrete outputs bitwise; subset probability comparisons retain the
Plan A <= 1e-12 contract; centrality weighted values retain the existing group-level
probability reduction drift documented in Evidence), and remove the Python dependency:

- `netAnalysis_computeCentrality`: signature unchanged. Entry assert
  `.sc_assert_spatial_cell_chat` when object is non-NULL; `slot.name <- match.arg(slot.name, c("net","netP"))`;
  do.group reads `unclass(slot$group$prob)` (+ `$group$pval`, NULL = no gating, explicit
  for netP) with the legacy two-step gate `vm[pm==0] <- 1; pm[vm>=thresh] <- 0` kept verbatim
  for net; do.cell reads `unclass(slot$cell$prob)` (replaces `$tmp$prob.cell` and
  spatstat.sparse dimnames helper); `net` argument (2D/3D dense) and `signaling.name`
  filtering preserved verbatim; the two duplicated `my_future_sapply` branches merged via a
  single `centr.name` vector; L184 replaced by
  `array(as.vector(centr.all), c(nrow(centr.all)/N, N, nrun))`; slot write-back via
  copy-out/whole-slot assignment (`netobj <- slot(object, slot.name); netobj[["centr"]] <- ...; slot(object, slot.name) <- netobj`);
  `.log_operation` before returning the object.
- `subsetCommunication`: wrapper reads `misc$.mode` (missing mode defaults to single), entry assert; no `.log_operation` (pure getter). pairLR.use pathway_name flip and signaling-conflict stop preserved (`&` → `&&`, scalar operands, behavior identical).
- `subsetCommunication_internal`: replace only the `if (!is.data.frame(net))` block —
  per-layer dense extraction from the SparseChatArray, single-step `pm[pv >= thresh] <- 0`
  for net (pval NULL = no gating), netP layers ungated with `pval` column = `NA_real_`
  (BUG-C); join keys: net uses `match(net$interaction_name, rownames(pairLR))` verbatim,
  netP uses `match(net$interaction_name, pairLR$pathway_name)` (BUG-D, recorded deviation);
  join columns missing from LR are back-filled with `NA_character_` before
  `dplyr::select` (BUG-A); numeric sources/targets out of range -> explicit stop
  (BUG-E). All downstream filters, guards, stop/warning texts, netP summarize step and
  final column selections preserved verbatim. `reshape2::melt` leaves this function.
- NAMESPACE/DESCRIPTION untouched (reticulate remains a dependency for utilities.R
  umap/geosketch).

## Constraints and invariants

- Discrete results (support sets, counts, measure names, column names, error texts)
  bitwise identical to legacy semantics on legal inputs; subset probability floats
  retain the Plan A group-level prob contract within 1e-12, while centrality weighted
  values use the E2E fallback/tolerance documented in Evidence for inherited reduction
  drift.
- Gate order/step count preserved exactly: net centrality = two-step
  (`pval[prob==0] <- 1` then `prob[pval>=thresh] <- 0`); net subsetCommunication = single
  step (`prob[pval>=thresh] <- 0`); netP paths never gate (already gated upstream by
  computeCommunProbPathway).
- `computeCentralityLocal`, `my_future_sapply`, `searchPair` unchanged.
- centr outputs keep shape `(n measures, node.names, signaling.name)`; downstream readers
  (netAnalysis_signalingRole_*, visualization.R L3119/L3166, netMappingDEG,
  updateClusterLabels) consume `slot[["centr"]]`/`slot[["centr.cell"]]` 3D arrays only —
  zero-change requirement.
- No nested slot replacement; whole-slot write-back only.
- Wave 2b aggregateNet subset branch (modeling.R) untouched; internal carries its own row
  construction.

## Alternatives considered

- Keep reticulate::array_reshape: rejected — 146-175 s init cost, discovery halt observed,
  and `order="F"` is provably identical to R column-major `array(as.vector(x), dim)`.
- Drop the `evidence` column when absent instead of back-filling NA: rejected by default —
  output schema must stay identical to upstream CellChat documentation (evidence is a
  constant metadata column there); back-fill keeps all 3+1 final-column variants valid.
- Refactor a shared helper with Wave 2b aggregateNet subset branch: rejected — Wave 2b code
  is tested and frozen; duplicating the small row-construction avoids cross-wave coupling.
- Keep `options$mode` strict check plus a producer-side fix writing `options$.mode`: rejected — the final class has no `options` slot; its validator owns `misc$.mode`, so the wrapper must follow the shipped schema.

## Consumer impact

- Producers: `computeCommunProb` (net$cell$prob), `computeAvgCommunProb` (net$group$prob/pval),
  `computeCommunProbPathway` (netP$group$prob, netP$cell$prob).
- Callers: visualization.R L3119/L3166 (subsetCommunication net slot), analysis.R L3231
  netMappingDEG (thresh only), utilities.R L1229 updateClusterLabels
  (netAnalysis_computeCentrality slot.name="netP"), centr readers L2562-3178
  (netAnalysis_signalingRole_network/scatter/heatmap etc.).
- Tests: new tests_dev/test-centrality-subset.R; regression suites test-netpathway-aggregate.R,
  test-computeAvgCommunProb.R, test-filterCommunication.R, test-SparseChatArray.R,
  test-computeCommunProb.R; e2e tests_dev/simulated_data/compare_centrality_subset_s3_s2.R.
- Docs: docs/COMPUTATION_OPTIMIZATION.md §9.4 appends the Wave 3a row. roxygen updates for
  the two functions. NAMESPACE/DESCRIPTION unchanged.

## Consequences

- New-schema objects gain working centrality and communication-subset paths; legacy
  CellChat objects now fail fast with a clear assert message instead of cryptic NULL errors.
- Python/reticulate removed from the analysis path: no interpreter init in this pipeline
  (saves ~150 s per session here); risk of Python environment drift eliminated.
- BUG-A~E are crash-path/misleading-error fixes, not semantic changes; the only content
  difference on legal inputs is the netP `pval` column being explicitly NA_real_ instead of
  a crash, and evidence = NA when LRsig lacks the column.
- Known deviation recorded: netP join by `LR$pathway_name` fixes multi-LR pathways (legacy
  matched interaction-name rownames and produced all-NA rows); on the current fixture
  (rownames(LRsig) == pathway_name) output is bitwise identical to legacy intent.
- merged-mode branch of the wrapper becomes unreachable (kept; cleanup deferred to Wave 3b).

## Evidence

- Legacy source: R/analysis.R L82-201 (centrality), L2266-2365 (wrapper), L2370-2530
  (internal); reticulate call L184.
- 11-slot writers: R/modeling.R L395 (net$cell$prob), L1292-1308 (net$group$prob/pval),
  L1750-1753 (netP$group$prob / netP$cell$prob as SparseChatArray).
- Schema/mode: R/SpatialCellChat_class.R L1240 (misc$.mode), L1334
  (.sc_assert_spatial_cell_chat); validator L683-688.
- Baseline package probe (`compare_centrality_subset_s3_s2.R`) on `SpatialChat_2.rds` records the expected evidence-column errors for both net and netP. Baseline centrality calls were attempted for net/netP × degree/full; the frozen reticulate path reproducibly failed with `should not be called with handlers on the stack`, so the oracle used the baseline package's pure-R `computeCentralityLocal` over the original arrays as a documented fallback.
- `Rscript tests_dev/test-centrality-subset.R` passed all checks, including 5/11 measure names, gate order and pval boundary, native reshape, cell/group write-back, netP no-gate, pathway filtering/join, evidence NA, pval NA, BUG-A~E, assert/log/validObject, and determinism.
- `Rscript tests_dev/simulated_data/compare_centrality_subset_s3_s2.R` ended with `CENTRALITY_SUBSET_E2E_OK`: four independent subset oracles passed; unweighted degrees matched bitwise; weighted degree max absolute drift was 5.684342e-12 / 4.831691e-12 from the already-recorded group-probability reduction difference; full measures stayed within 1e-9 and deterministic fallback checks passed.
- `Rscript` parse check passed after the final source edits. Documentation was updated in `docs/COMPUTATION_OPTIMIZATION.md` §9.4.

## Acceptance criteria

1. `Rscript tests_dev/test-centrality-subset.R` ended with "All checks passed."; zero FAIL. It covers reshape equivalence, netP no-gating, BUG-A~E, assert/log/write-back/validObject, filters and determinism.
2. `Rscript tests_dev/simulated_data/compare_centrality_subset_s3_s2.R` ended with `CENTRALITY_SUBSET_E2E_OK`; baseline evidence crashes were recorded; four subset calls matched independent manual rows (keys, pval, prob tolerance); unweighted centrality was bitwise and full-measure comparisons were within 1e-9. The frozen reticulate centrality call was unavailable in this environment and is explicitly recorded with a pure-R baseline fallback; weighted degree drift was 5.684342e-12 / 4.831691e-12, attributable to the existing group-probability reduction difference.
3. Regression suites `test-netpathway-aggregate.R`, `test-computeAvgCommunProb.R`, `test-filterCommunication.R`, `test-SparseChatArray.R`, and `test-computeCommunProb.R` were run after implementation and passed.
4. The note validator reports zero errors; no scratch/scaffold files remain in the final tree.
