# Agent Note: netpathway-aggregate-migration

- **Date**: 2026-10-08
- **Status**: implemented
- **Governance**: v1
- **Wave**: 2b (inference chain: computeCommunProbPathway / aggregateNet / relabelSpatialCellChat)

## Problem

The inference chain functions up to `filterCommunication` run on the 11-slot schema (`net$cell$prob` / `net$group$prob` SparseChatArray). The next three chain functions still read legacy slots and break on migrated objects:

1. `computeCommunProbPathway` (R/modeling.R L1616-1745) reads `net$prob` / `net$pval` (do.group) and `net$prob.cell` + `net$tmp$prob.cell` (do.cell, via deprecated `my_as_sparse3Darray` + `spatstat.sparse::marginSumsSparse`), writes `netP$prob` / `netP$prob.cell` / `netP$tmp`.
2. `aggregateNet` (R/modeling.R L2247-2305) reads `net$prob` / `net$pval` / `net$prob.cell` (sparse3Darray), writes `net$count` / `net$weight` / `net$LR.sig` / `net$count.cell` / `net$weight.cell` / `net$LR.sig.cell`; its subset branch delegates to legacy `subsetCommunication(slot.name="net")`.
3. `relabelSpatialCellChat` (R/modeling.R L616-617) dispatches on `object@net[["tmp"]]$cell.type.decomposition` (written by legacy `computeAvgCommunProb_Visium`).

## Schema I/O Contract

- `computeCommunProbPathway` reads `net$group$prob` / `net$group$pval` (SparseChatArray; pval may be NULL), `net$cell$prob` (SparseChatArray), `object@LR$LRsig` (layer names = dimnames[[3]] = rownames). Writes `netP$group$prob` (SparseChatArray K×K×nPathways.sig), `netP$cell$prob` (SparseChatArray nC×nC×nPathways.sig), `netP$pathways` / `netP$pathways.cell` (char, descending by strength). `netP$tmp` is never written (explicitly removed). Merge-update of `object@netP` branches (validator requires netP$cell / netP$group to be non-NULL lists; wholesale replacement forbidden). `object=NULL` returns an isomorphic list.
- `aggregateNet` writes `net$group$count` / `net$group$weight` (dgCMatrix, dimnames = group names), `net$group$LR.sig` (char), `net$cell$count` / `net$cell$weight` (dgCMatrix, dimnames = cell names), `net$cell$LR.sig` (char). Legacy `net$count` / `net$weight` / `net$LR.sig` / `net$count.cell` / `net$weight.cell` / `net$LR.sig.cell` are no longer written.
- `relabelSpatialCellChat` dispatch reads `object@misc$.param$averaging$cell.type.decomposition`.
- Layer names are read from `dimnames(...)[[3]]`, never `names(unclass(...))` — the SparseChatArray constructor strips internal list names.

## Equivalence Contract

- Discrete results bitwise identical to baseline: pathways.sig membership and ordering, integer counts, LR.sig sets, error paths and messages.
- Cross-layer float sums (pathway aggregation, weight): same real-number multiset summed in different order → ≤ 1 ulp (Plan A precedent); tests assert ≤ 1e-12. In practice the fixture comparisons are bitwise identical (`.sc_sum_layers` preserves layer order per pathway).
- Baseline quirks preserved verbatim: two-step gating order in aggregateNet default branch (`pval[prob==0] <- 1` then `prob[pval >= thresh] <- 0`); single-step `>=` gating in computeCommunProbPathway and subset paths; cell-level count binarizes `m@x <- 1` (explicit zeros count as links); subset-branch `count = n()` rows, `weight = sum(prob)`; stop message "No significant signaling interactions are inferred based on the input!".
- Documented divergence: legacy `subsetCommunication` silently flipped slot.name to "netP" when pairLR.use contained pathway_name — quirk not reproduced in the native subset branch (no internal caller passes pairLR.use to aggregateNet).
- pval NULL (do.permutation=FALSE objects): new explicit semantics = skip significance gating (count/weight over all layers). Old behavior on this path was an accidental error/no-op.
- Visium branch of relabelSpatialCellChat unchanged; unreachable on 11-slot objects until computeAvgCommunProb_Visium migrates (dispatch field absent → standard branch, same as old dispatch).

## Implementation Notes

- New helper `.sc_sum_layers` (R/modeling.R, @noRd): single layer passthrough, else delayed `get("cpp_sum_layers", inherits = TRUE)` (mirrors marginSums.SparseChatArray).
- R/List `$` partial matching pitfall: `object@netP$pathways` on an object whose netP lacks `pathways` silently matched `pathways.cell`. External readers of `netP$pathways` keep working (element exists after do.group=TRUE); tests use exact `[[`.
- relabelSpatialCellChat: `object@meta[["new.ident"]] <- factor(...)` strips names via data.frame `[[<-`; both branches now restore `names(object@idents) <- colnames(object@assay$norm)` (validator requires idents names to match cell names).
- Cross-layer sums that exceed the 1-ulp contract in principle (pathway aggregation across gated layers) are computed with cpp_sum_layers; the test fixture confirms bitwise equality at fixture scale, the 1e-12 assertion remains the general guarantee.

## Acceptance Checklist

- [x] New suite `tests_dev/test-netpathway-aggregate.R` all green (38 checks: hand-built fixtures + subset replication + relabel e2e + determinism).
- [x] Regression suites green: test-computeAvgCommunProb.R ("All computeAvgCommunProb Plan A checks passed", 0 FAIL), test-filterCommunication.R (all PASS), test-SparseChatArray.R (41 passed / 0 failed), test-computeCommunProb.R (Wave 1c + 1d all passed).
- [x] validObject passes on migrated objects; no `net$tmp` / `netP$tmp` written anywhere in the three functions.
- [x] Callers adapted: relabelSpatialCellChat dispatch read (`misc$.param$averaging`); analysis.R L2674 `$count` → `$group$count`.
- [x] Roxygen @return updated; docs/COMPUTATION_OPTIMIZATION.md §9.4 row added; obsolete pbsapply risk row updated.
- [ ] Commit (pending user approval, per session convention).

## Evidence

- New suite: 38/38 `[PASS]`, final line "All checks passed." (`Rscript tests_dev/test-netpathway-aggregate.R`, ~21 s).
- Regressions: test-computeAvgCommunProb.R → "All computeAvgCommunProb Plan A checks passed." (0 `[FAIL]`); test-filterCommunication.R → all PASS incl. rejection paths; test-SparseChatArray.R → "41 passed, 0 failed"; test-computeCommunProb.R → "All Wave 1c main-body checks passed." + "All Wave 1d filterProbability checks passed.".
- Gating anchor (case 1c): P1 (1,2) equals g1[1,2] exactly — LR2's contribution is zeroed by `pval == thresh` under baseline `>=` semantics (bitwise `identical`).
- Determinism (case 11): same-seed relabel chain reproduces `net$group$weight`, `netP$group$prob`, `net$group$count`, `netP$pathways`, `net$group$LR.sig` bitwise (`identical`).
