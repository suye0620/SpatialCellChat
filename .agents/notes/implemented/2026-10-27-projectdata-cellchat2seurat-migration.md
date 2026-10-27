# Agent Note: migrate projectData / CellChat2Seurat to the 11-slot schema

Status: implemented
Governance: v1
Date: 2026-10-27
Decision type: migration
Scope: R/utilities.R (`projectData`, projection chain, `CellChat2Seurat`, `subsetData`), R/modeling.R (`computeCommunProb` raw.use branch), R/SpatialCellChat_class.R (`.sc_validate_assay`, `.sc_assay_layer`), tests_dev
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

Two exported functions still speak the pre-migration schema and cannot run on a current object:

1. `projectData()` (R/utilities.R L1689) reads `object@data.signaling` and writes
   `object@data.project`. Neither slot exists in the 11-slot class, so any call fails
   with "no slot of name". It is the documented entry point of the optional PPI-smoothing
   branch (`docs/COMPUTE_COMMUN_PROB_PREREQUISITES.md` §8) and the only producer of the
   matrix that `computeCommunProb(raw.use = FALSE)` is supposed to consume.
2. `CellChat2Seurat()` (R/utilities.R L581) reads `object@meta`, `object@idents`,
   `object@data` and `object@images$coordinates`. All four moved. It is the prerequisite
   of the topic/community chain (PP1 batch B5) and of `mergeChatIntoSeurat()`
   (R/utilities.R L123), which calls it twice.

`computeCommunProb()` currently short-circuits this whole branch with a `stop()`
(R/modeling.R L182), so the regression is masked: nothing else fails, and the documented
`raw.use = FALSE` path is simply unavailable.

## Evidence

- `R CMD INSTALL` and `R CMD check` are already blocked by four undefined exports; these
  two functions are not part of that failure, they fail at call time instead.
- `docs/PP1_DOWNSTREAM_REFACTORING_PLAN.md` L252 records the actual failure:
  "CellChat2Seurat 显式传 counts 仍因旧 @data 失败".
- `docs/ANALYSIS_SEMANTICS.md` A-17 flags that the Seurat-side scale (raw counts placed
  into `data`, AverageExpression behaviour) is still unresolved and must be pinned down
  here, not assumed.
- Algorithm audit against netSmooth (Ronen & Akalin 2018, GPL-3), from which these six
  helpers were vendored: the four core functions (`projectOnNetwork`,
  `randomWalkBySolve`, `smoothAndRecombine`/`projectAndRecombine`,
  `projectFromNetworkRecombine`) are faithful line-by-line copies. Six defects listed
  below were found by comparing against the upstream sources.
- `assay<-` calls `.sc_validate_after_update()` → `validObject()`, so any written layer is
  validated immediately. A validator change and the writer change must land together.

## Decision

### D1. Slot mapping

| Old slot | New |
|---|---|
| `@data` | `assay(object, "norm")` |
| `@data.signaling` | `assay(object, "signaling")` |
| `@data.project` | `assay(object, "smooth")` |
| `@meta` | `meta(object)` |
| `@idents` | `idents(object)` |
| `@images$coordinates` | `images(object, "coordinates")` |
| `@options$parameter$raw.use` | `misc$.param$communication$raw.use` (already written) |

The `signaling → smooth` mapping is already the documented contract
(`docs/COMPUTE_COMMUN_PROB_PREREQUISITES.md` §8), so this note only closes the gap between
documentation and code.

### D2. `assay$smooth` lives in the signaling gene space — needs sign-off

`projectFromNetworkRecombine()` writes the smoothed values back into the *original*
matrix, which is `assay$signaling`. The result therefore has signaling's rows, not the
full gene set.

`.sc_validate_assay()` (R/SpatialCellChat_class.R L362) currently allows only `signaling`
to be a gene subset; every other layer must be full `genes × cells`. Writing the
projection as-is would fail `validObject()` immediately.

**Decision: put `smooth` into the same subset branch as `signaling`.** Rules become:

- columns must equal the object's cells, in order;
- rownames must be a subset of `assay$norm` rownames;
- row order must follow `assay$norm`.

Rationale: matches baseline semantics, keeps the stored matrix small (a 1.5k × 30k
signaling space instead of a 20k × 30k dense full-gene matrix), and no writer currently
produces a full-dimension `smooth`, so nothing regresses.

Rejected alternative: store a full-gene matrix with zeros outside the network. It would
avoid the validator change but would falsely imply "these genes were smoothed to zero".

### D3. Six implementation defects, fixed in the same change

**P1 (blocker at scale) — dense-only path; the projected matrix is sized by the PPI.**
`projectOnNetwork()` allocates `matrix(0, nrow = length(rownames(adj_matrix)), ncol = n_cells)`.
The row count comes from the **network**, not from signaling: a 1.5k-gene signaling set
against an 18k-gene PPI still allocates 18k × n_cells. With 30k cells that is ~4.3 GB,
plus `diag(18000)` (2.6 GB), plus `AA` (2.6 GB), plus an O(n³) `solve()`.
netSmooth ships `randomWalkByIterations()` as an escape hatch
(`f ← α·A_norm·f + (1−α)·f₀`, cost `O(nnz(A)·m·iters)`, no identity matrix, same fixed
point). We only vendored the `solve()` path.

Plan: add internal `randomWalkByIterations()`; add `solver = c("solve", "iterative")` to
`projectData()`, **default `"solve"`** so the numeric result stays identical to the
baseline and to any previously produced `smooth`. `"iterative"` is opt-in with
`tol = 1e-8`, `max.iter = 500`, and a message stating that it is not bit-identical.

**P2 (silent wrong answer) — no gene-name overlap check.**
If PPI rownames and signaling rownames do not intersect (Ensembl vs symbol, case,
version suffixes), `projectOnNetwork()` yields an all-zero network space, the walk of a
zero vector is zero, and `projectFromNetworkRecombine()` writes nothing back: the
function returns approximately its input, with no warning.

Plan: after resolving `new_features`, compute
`overlap <- intersect(rownames(sig), rownames(adj))`; `stop()` on zero overlap;
`warning()` when `length(overlap) < 0.1 * nrow(adj)`, reporting both counts.

**P3 (stale state) — `subsetData()` does not invalidate `assay$smooth`.**
`normalizeData()` clears scale/smooth/signaling (R/utilities.R L818); `subsetData()` only
writes signaling (L1199). After this migration a user could smooth, re-subset, and then
run `raw.use = FALSE` against a smooth matrix from a different gene set.

Plan: `assay(object, "smooth") <- NULL` in `subsetData()` right before writing signaling,
matching `normalizeData()`.

**P4 (inherited, doc/code conflict) — dead alpha check.**
`stopifnot(is.numeric(alpha) && (alpha > 0 && alpha < 1))` already excludes 0 and 1; the
following `if (alpha < 0 | alpha > 1) stop(...)` can never fire. Meanwhile the roxygen
block advertises `alpha = 0: no smoothing`, so a documented value produces
`alpha > 0 && alpha < 1 is not TRUE`.

Plan: keep rejecting 0 (no numeric surprises), delete the dead branch, and rewrite the
parameter description as "numeric in (0, 1); 0 is not supported".

**P5 (silent wrong answer) — adjacency matrix is not validated.**
Projection uses `rownames(A)`; the walk uses `A`'s actual row/column order. If
`colnames(A)` is a different permutation, the two spaces disagree and the result is
silently scrambled.

Plan: `stopifnot(nrow(A) == ncol(A), !is.null(rownames(A)), identical(rownames(A), colnames(A)))`.

**P6 — `alpha = "auto"` is deliberately not supported.**
netSmooth's default is `alpha = "auto"`: it runs the full smoothing once per candidate in
`0.1 * (1:9)`, scores every result by cluster robustness or 2D-embedding entropy, and
keeps the best. Three reasons to keep it out of this package:

1. **Dependency cost.** `scoreSmoothing()`'s robustness branch needs `clusterExperiment`;
   the entropy branch needs a 2D PCA embedding plus `pickDimReduction()`. Those are
   exactly the `Depends` (`scater`, `clusterExperiment`) that vendoring the six helpers
   was meant to avoid.
2. **Cost multiplier on the wrong axis.** `auto` performs nine full smoothings. With P1
   unresolved each one allocates a dense `n_net × n_cells` matrix and an `O(n³)` solve;
   adding a 9× multiplier before fixing the memory profile is the wrong order.
3. **Misaligned objective.** `auto` optimises "are the cells better separated", not "is
   the ligand/receptor expression estimate more trustworthy". Here the smoothed matrix is
   only an input to `computeCommunProb()`, so an alpha that separates clusters better can
   flatten a group-specific ligand.

Decision: `projectData()` accepts a single numeric `alpha` only; anything else stops with
a message pointing at this note. The roxygen block states that netSmooth's
`alpha = "auto"` usage is not reachable. Revisit only if (a) P1 is solved, (b) a scoring
rule aligned with the communication objective is defined, and (c) the candidate range
becomes an explicit argument.

### D4. `CellChat2Seurat()`

Straight substitution, no logic change:

```r
meta.df <- meta(object)
meta.df$group.cellchat <- idents(object)
norm <- assay(object, "norm")
if (is.null(norm)) stop("assay$norm is empty", call. = FALSE)
data.use <- if (is.null(counts)) norm else counts
obj.new <- Seurat::CreateSeuratObject(counts = data.use, meta.data = meta.df, assay = assay)
if (identical(assay, "Spatial")) {
  coords <- images(object, "coordinates")
  if (!is.null(coords)) {
    obj.new@images$image <- new("SlideSeq", assay = "Spatial", key = "image_",
                                coordinates = coords)
  }
}
Seurat::Idents(obj.new) <- if (is.null(group.by)) "group.cellchat" else group.by
```

Verified facts that shape this code:

- `Seurat::Idents(obj) <- "group.cellchat"` resolves the string to the metadata column
  (probed on the pinned stack: idents become `A,A,B,B,B`, `nlevels = 2`). It is not a bug
  and must not be "fixed".
- `SlideSeq` still exists in SeuratObject 5.0.2 with slots `coordinates / assay / key`.
- `obj.new@assays[[...]]@misc[["features.info"]]` stays as a direct slot write: this
  project targets the v3 assay layout that Seurat 4.4.0 creates; do not abstract it for
  Assay5.
- `images$coordinates` may be NULL; the baseline assumed non-NULL and crashed.

Optional (P2 of this migration, separate commit): add `layer = c("norm", "raw", "signaling")`
so the community chain can export signaling/raw without going through `counts`.

## Consumer impact

- `mergeChatIntoSeurat()` becomes functional again for free (it only calls
  `CellChat2Seurat()`).
- `computeCommunProb(raw.use = FALSE)` stops erroring; its input changes from
  `@data.project` to `assay$smooth`, same values, same gene space.
- Both functions stay exported; no NAMESPACE change.
- `makeGridSpatialCellChat(data.slot = "smooth")` would now receive a gene subset and
  build a grid object whose `norm` is a subset. Either document it as unsupported for
  `smooth`, or aggregate to the full gene space. Decide in that commit; not a blocker here.

## Constraints and invariants

- `assay$smooth` must never exist without a matching `assay$signaling`; any writer of
  signaling clears smooth (P3).
- The `solve` path must stay numerically identical to the baseline for a fixed
  `(signaling, adjMatrix, alpha)`. Changing the default solver is out of scope.
- The 11-slot class keeps five assay layers: `raw, norm, scale, smooth, signaling`.
  Do not add a sixth.
- Seurat interop targets the v3 assay layout produced by Seurat 4.4.0; no Assay5 branches.

## Consequences

- `raw.use = FALSE` becomes usable, and `projectData()` regains a documented path.
- Large networks are still slow by default; users must opt into `solver = "iterative"`.
  This is recorded rather than silently changed, because the default affects numbers.
- `subsetData()` gains one line and one new invalidation rule.
- The validator accepts a gene-subset `smooth` from now on; anything that assumed
  full-dimension smooth must be re-checked (currently nothing writes it).

## Alternatives considered

- **Archive both functions (group B treatment).** Rejected: both sit on documented user
  paths, and the community/topic chain (PP1 B5) needs `CellChat2Seurat`. Archiving would
  delete functionality the project still plans to ship.
- **Store the projection in `misc$` instead of `assay$smooth`.** Rejected: breaks the
  "all expression lives in assay" invariant and contradicts the already-published
  documentation.
- **Switch the default solver to iterative.** Rejected: the result would differ from the
  baseline at ~1e-8 relative, which is exactly the kind of silent numeric drift this
  project has been documenting against.
- **Support `alpha = "auto"`.** Rejected, see P6: it needs `clusterExperiment`/`scater`,
  multiplies the most expensive step by nine, and optimises a criterion that does not
  match this package's objective.

## Acceptance criteria

1. `R CMD INSTALL` shows no new failure beyond the four pre-existing undefined exports.
2. `tests_dev/test-projectData.R` (new):
   - empty signaling, alpha outside the range, PPI with a zero row/column, non-square
     adjacency, mismatched `rownames`/`colnames`, zero gene overlap → all rejected with
     the intended message;
   - output row names equal `rownames(assay$signaling)`, column names equal the cells;
   - `validateSpatialCellChat()` passes after the write (proves D2);
   - on a 3-gene chain PPI, `solve` and `iterative` agree to `<= 1e-8` relative and both
     match a hand-computed `(I − αA)⁻¹(1 − α)E`;
   - `normalizeData()` and `subsetData()` both clear `smooth`.
3. `tests_dev/test-seurat-bridge.R` (new):
   - cell count, `Idents` levels and the `group.cellchat` column match the source object;
   - `counts=` overrides `assay$norm`;
   - `assay = "Spatial"` creates one image, `assay = "RNA"` creates none;
   - matrix round trip `CellChat2Seurat()` → `LayerData()` equals `assay$norm`;
   - `check.object = TRUE` is never exercised in tests (it prints).
4. End-to-end: on one fixture, `projectData(alpha = 0.5)` then
   `computeCommunProb(raw.use = FALSE)` runs and `misc$.param$communication$raw.use` is
   `FALSE`; results differ from `raw.use = TRUE` only on network genes.

## Rollout

Four commits, each independently reviewable:

1. Validator (`smooth` subset rules) + `subsetData` invalidation + adjacency/gene-overlap
   validation helpers.
2. `projectData` slot migration + `solver` parameter + `randomWalkByIterations()`.
3. `computeCommunProb` `raw.use = FALSE` branch + end-to-end test.
4. `CellChat2Seurat` slot migration + Seurat bridge tests.

Open item before step 1: confirm D2 (smooth in the signaling gene space).

## Implementation outcome

D2 was confirmed; the four commits landed as planned:

| Commit | Content |
|---|---|
| `4e26ad6` | `.sc_validate_subset_layer()` + smooth in the subset branch; `subsetData()` clears smooth |
| `8f2b131` | `projectData()` slot migration, `.sc_validate_adjacency()`, overlap stop/warning, `randomWalkByIterations()`, `solver` argument; fixes the `normalizeAjdMatrix` typo |
| `eeb4312` | `computeCommunProb(raw.use = FALSE)`, `.sc_scale_expression()` |
| `cd8035b` | `CellChat2Seurat()` accessor migration, `.sc_assay_layer()` |

Deviations from the plan, all recorded deliberately:

1. The validation helpers live next to their only caller (`projectData`) instead of
   landing in step 1, so no commit carries unused code.
2. Two defects surfaced only while porting `CellChat2Seurat`, so they were fixed there
   rather than in a separate change:
   - `SlideSeq` validates `coordinates` as a **data.frame** while the object stores a
     matrix, so `new("SlideSeq", ...)` failed on every call;
   - the function's own `assay` argument shadows the `assay()` accessor, so the layer is
     read through `.sc_assay_layer()`. A test fixture that defines a global `meta` object
     hits the same shadowing and fails with "could not find function"; the test renames
     its fixture instead of weakening the accessor.
3. `computeCommunProb` scales by matrix class (`.sc_scale_expression()`): the sparse path
   keeps `@x / max(@x)`, the dense projection uses `data / max(data)`, matching the two
   baseline branches. An all-zero layer now errors instead of yielding NaN.

Verification on 2026-10-27: `test-projectData.R` 19 checks, `test-seurat-bridge.R` 14,
`test-computeCommunProb.R` 61 (up from 54), `test-computeAvgCommunProb.R` 162,
`test-filterCommunication.R` 22, `test-cpp-prob-layer.R` 20,
`test-computeCommunProb-helpers.R` 16, `test-parallel.R` 147 passed,
`test-normalizeData.R` 9, `test-pre-inference-new-schema.R` 15,
`test-SpatialCellChat-class.R` / `-accessors.R` 22 each.
`R CMD INSTALL` still reports only the four pre-existing undefined exports.
`test-SparseChatArray.R` cannot run here: it requires the uninstalled `microbenchmark`.

Not done, by decision: `alpha = "auto"` (P6). Still open elsewhere: `makeGridSpatialCellChat()`
with `data.slot = "smooth"` now receives a gene subset and builds a grid object whose norm
is a subset; document or guard it when that path is touched.
