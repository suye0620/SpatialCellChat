# Agent Note: archive the superseded expression / aggregation toolchain

Status: implemented
Governance: v1
Date: 2026-10-26
Decision type: cleanup
Scope: R/modeling.R, R/utilities.R, R/legacy/, NAMESPACE, man/legacy/, tests_dev reference oracles
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

The 11-slot migration replaced the per-pair expression gating with the Rcpp kernels
(`cpp_prob_layer`, `cpp_group_avg_obs`, `cpp_group_avg_perm`) and replaced the
`sparse3Darray` container with `SparseChatArray`. Fifteen functions that implemented the
old path survived in `R/modeling.R` and `R/utilities.R` with **zero internal callers**,
and eight of them are not runnable at all on the new schema.

They keep the package's dependency surface and `man/` in a state that no longer matches
the code, and they make the function inventory unreadable: a reader cannot tell a live
helper from a superseded one.

## Evidence

Reference scan across `R/`, `tests_dev/`, `man/` (2026-10-26):

- All 15 have zero callers inside `R/` except `my_as_sparse3Darray`, called once by
  `netVisual_CommunField` (R/visualization.R L5516) on the retired `net$field` slot.
- No test calls them except the three frozen-comparison scripts listed under Consumer impact.
- `computeAvgCommunProb_LR_Avg/Sum` are superseded by `computeAvgCommunProb`; the old
  `format(digits = 1)` gate they implemented is no longer the acceptance contract
  (see `docs/COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md`).
- Inventory: `docs/FUNCTION_INVENTORY.md` §0 and §4, rows labelled 未动 / 疑似.

Original locations at archive time:

| File | Line | Function |
|---|---:|---|
| R/modeling.R | 483 | `HillFunctionFordataLR` |
| R/modeling.R | 505 | `elementwiseProductForSparseMatandDenseMat` |
| R/modeling.R | 536 | `myElementwiseProduct` |
| R/modeling.R | 1378 | `computeAvgCommunProb_LR_Sum` |
| R/modeling.R | 1449 | `computeAvgCommunProb_LR_Avg` |
| R/modeling.R | 2046 | `modifySparse3Darray` |
| R/modeling.R | 2565 | `.computeExprGroup_complex` |
| R/modeling.R | 2631 | `computeExpr_coreceptor` |
| R/modeling.R | 2678 | `.computeExprGroup_coreceptor` |
| R/modeling.R | 2737 | `computeExprGroup_agonist` |
| R/modeling.R | 2769 | `computeExprGroup_antagonist` |
| R/modeling.R | 2803 | `computeExpr_agonist` |
| R/modeling.R | 2838 | `computeExpr_antagonist` |
| R/utilities.R | 361 | `==.sparse3Darray` |
| R/utilities.R | 491 | `my_as_sparse3Darray` |

## Decision

Move all fifteen definitions verbatim into `R/legacy/expression-toolchain.R`, together
with their roxygen blocks and a header recording the mapping and the restore procedure.
Because `R CMD INSTALL` only sources top-level `.R` files under `R/`, the subdirectory is
not compiled: the code stays in the repository for reference but leaves the package.

Concretely:

1. 540 lines removed from `R/modeling.R`, 131 from `R/utilities.R`.
2. 12 `export()` lines removed from `NAMESPACE`; the other three
   (`.computeExprGroup_complex`, `.computeExprGroup_coreceptor`, `==.sparse3Darray`) were
   never exported. `S3method("==",sparse3Darray)` removed as well.
3. 12 matching `man/*.Rd` moved to `man/legacy/` so that `R CMD check` does not report
   usage entries without code.
4. Three frozen-comparison scripts gained an idempotent guard before use:

```r
if (!exists("<archived fn>", mode = "function")) source("R/legacy/expression-toolchain.R")
```

## Consumer impact

Breaking change for twelve exported functions. Downstream callers of
`HillFunctionFordataLR`, `myElementwiseProduct`, `computeExpr_*`, `computeExprGroup_*`,
`computeAvgCommunProb_LR_Avg/Sum`, `modifySparse3Darray`, `my_as_sparse3Darray` and
`elementwiseProductForSparseMatandDenseMat` will get "could not find function".

Replacement mapping recorded for users:

| Archived | Use instead |
|---|---|
| `HillFunctionFordataLR`, `elementwiseProduct*`, `myElementwiseProduct` | `cpp_prob_layer` |
| `computeExpr_*`, `computeExprGroup_*`, `.computeExprGroup_*` | `.sc_expr_row` + `.sc_complex_subunits` + `.sc_cofactor_factor` |
| `computeAvgCommunProb_LR_Avg/Sum` | `computeAvgCommunProb` |
| `modifySparse3Darray`, `my_as_sparse3Darray`, `==.sparse3Darray` | `SparseChatArray` and its `[` / `[[` / `dim` methods |

Two known residues, both deliberate:

- `netVisual_CommunField` (R/visualization.R L5516) still calls `my_as_sparse3Darray` on
  `net$field`. That path was already unreachable before this archive (nothing writes
  `net$field`); the function is queued for the next round together with
  `netVisual_CommunFieldZoomIn`.
- `tests_dev/test_data/{benchmark_mbm05_upstream_vs_refactored.R,debug_upstream_filterProbability_rscript.R}`
  call `SpatialCellChat::my_as_sparse3Darray`, i.e. through the **installed baseline
  package**, not this source tree. They remain valid as long as they run against the
  pinned `be88f30` build and must not be re-run against a future build of the current source.

## Constraints and invariants

- `R/legacy/` must stay out of `Collate`; adding it back would silently re-register the
  archived functions in the namespace.
- `man/legacy/` must not be moved back to `man/` without also restoring the functions.
- The archived bodies are frozen: they document what the pre-migration code did. Do not
  "fix" them.
- `spatstat.sparse` is now only referenced by the legacy file; its `importFrom()` lines in
  `NAMESPACE` and the `Imports:` entry in `DESCRIPTION` are candidates for a later cleanup,
  not part of this change.

## Consequences

- `R/modeling.R` 3071 → 2531 lines, `R/utilities.R` 2205 → 2074 lines.
- `R CMD INSTALL` still fails, but only with the four pre-existing undefined exports
  (`liftSpatialCellChat`, `mergeSpatialCellChat`, `subsetSpatialCellChat`,
  `updateSpatialCellChat`). This archive neither causes nor fixes that.
- The frozen-comparison oracles keep working, so the numeric-gate evidence remains
  reproducible.

## Alternatives considered

- **Leave them in place with `@deprecated`.** Rejected: the package would keep 670 lines of
  dead code, `man/` would keep 12 stale pages, and `R CMD check` would keep checking them.
- **`git rm` outright.** Rejected: the bodies are the only readable record of how the old
  gate and the old 3-D array handling worked, and several are still quoted in
  `docs/COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md`.
- **Archive the whole `疑似` list (40 functions).** Rejected for now: that list mixes
  genuinely superseded helpers with exported functions that merely lost their callers
  (`computeAveExpr`, `communicationDistPlot`). Group B (broken by removed slots) is
  handled separately.

## Acceptance criteria

Verified on 2026-10-26:

- All 13 `.R` files under `R/` (including `R/legacy/`) parse.
- `grep` over `R/` excluding `R/legacy/` finds no call to any archived function except
  the documented `netVisual_CommunField` site.
- `tests_dev/test-computeCommunProb-helpers.R` 16 checks pass,
  `test-cpp-prob-layer.R` 20 pass, `test-computeCommunProb.R` 54 pass,
  `test-computeAvgCommunProb.R` 162 pass.
- `R CMD INSTALL` reports only the four undefined exports and no new warnings.
