# Agent Note: consolidate the Seurat bridge into `asSeurat()`

Status: implemented
Governance: v1
Date: 2026-10-28
Decision type: api
Scope: R/utilities.R (`asSeurat`, `.sc_as_seurat_*`, `.sc_seurat_assemble`, `CellChat2Seurat`, `mergeChatIntoSeurat`), R/analysis.R (topic/community call sites), NAMESPACE, man/, tests_dev/test-seurat-bridge.R
Owner: maintainer
Impact: medium (public API, one breaking behaviour change in a deprecated wrapper)
Supersedes: none
Superseded by: none

## Problem

Two exported functions put SpatialCellChat content into a Seurat object, with no shared
vocabulary:

| Function | Axes exposed |
|---|---|
| `CellChat2Seurat()` | expression only |
| `mergeChatIntoSeurat()` | `slot.name` (`net`/`netP`) × `pattern` (incoming/outgoing) × `merge.chat` |

Four defects came with that split:

1. **The list branch only merged two objects.** `merge(seu.list[[1]], seu.list[[2]])`
   silently dropped every element past the second, and `names(seu@images) <- names(object.list)`
   then failed on the length mismatch.
2. **A dead `slot.name` branch.** `slot.name %in% c("data", "data.signaling",
   "data.scale", "data.project", "data.raw")` cannot be reached on an 11-slot object:
   those slots no longer exist.
3. **`netP` was annotated with the wrong table.** `features.info$interaction_name` was set
   to the centrality layer names — LR pair names for `net`, but pathway names for `netP` —
   and then left-joined to `DB$interaction`, producing all-NA columns on the `netP` path.
4. **`pattern` never selected anything.** Both directions were always computed; the value
   only ordered the `rbind` and assigned `-O`/`-I`. The default assay was hard-wired to
   `incoming` regardless of what the caller passed.

## Decision

One entry point, two axes:

```r
asSeurat(object,
         slot.name = c("assay", "net", "netP"),
         what = c("incoming", "outgoing", "merge"),
         counts = NULL, normalize.chat = TRUE,
         normalization.method = "LogNormalize",
         assay.name = "Spatial", group.by = NULL, check.object = FALSE)
```

`slot.name` selects the source slot and therefore the feature space; `what` selects which
part of that slot is written. `what` is ignored for `slot.name = "assay"` and `counts` is
ignored for the centrality modes; passing either in the wrong mode warns instead of
failing silently.

**`slot.name = "assay"`** maps the layers the way Seurat expects them:

| SpatialCellChat | Seurat assay |
|---|---|
| `assay$raw` (else `assay$norm`) | `counts` |
| `assay$norm` | `data` |
| `assay$scale` | `scale.data` |
| `assay$signaling`, `assay$smooth` | `misc$subsets` (gene subsets have no Seurat layer) |

`data` and `scale.data` are only attached when their dimnames match the counts layer, so a
user-supplied `counts` with other genes is never mixed with `norm`.

**`slot.name = "net"` / `"netP"`** reads `centr.cell` (requires
`netAnalysis_computeCentrality(do.group = FALSE)`), turns `indeg`/`outdeg` into a
features × cells matrix, and for `what = "merge"` stacks both with `-I` / `-O` suffixes.

`assay` (the target assay name) was renamed to `assay.name`: the old name shadowed the
`assay()` accessor and forced every internal read through `.sc_assay_layer()`.

## Alternatives rejected

- **Keep them separate and only fix the bugs.** Leaves two vocabularies for one operation
  and does not give `expression` a list form.
- **`what = c("expression", "centrality", "pathway.centrality")` with no `slot.name`.**
  Considered first, then dropped: `slot.name` is what the caller already knows from
  `netAnalysis_computeCentrality()`, and keeping it lets `what` stay a pure direction
  selector instead of a three-way enum that mixes source and direction.
- **Name the parameter `features`.** More precise about the result, but it reads like the
  gene vector taken by `subsetData(features = )`.
- **Keep `pattern` + `merge.chat`.** `what = "merge"` expresses the same thing with one
  value instead of a cross product.

## Consequences

- **Breaking for `mergeChatIntoSeurat()` callers**: `merge.chat = FALSE` used to return an
  object carrying `Spatial` = incoming plus `incoming`/`outgoing` aliases; it now returns
  only the requested direction. `what = "merge"` gives both in one assay.
- **`CellChat2Seurat()` default `check.object` stays `TRUE`** so existing calls behave the
  same; `asSeurat()` defaults to `FALSE`.
- Both old names remain exported as `.Deprecated()` wrappers for one release. The internal
  callers in `identifyCellTopics()` / `identifyTopicCommunities()` were moved to
  `asSeurat()` so package code never emits its own deprecation warning.
- Seurat rewrites `_` to `-` in feature names, so a centrality feature `L1_R1` arrives as
  `L1-R1`. `features.info$orig.features` keeps the untouched name; `srt.features` is now
  read back from the object, so it matches the actual rownames (it did not before).
- `asSeurat()` accepts a list in every mode, including `expression`, which the old API did
  not support at all. Elements are merged with `Condition` set to the list names.

## Verification

`tests_dev/test-seurat-bridge.R`, 41 checks: layer mapping (`raw`→counts, `norm`→data,
`scale`→`scale.data`, subsets→misc), SlideSeq image hand-off, counts override with other
genes staying out of the data layer, argument gating warnings, `net` incoming/outgoing/
merge with values checked against `indeg`/`outdeg`, `normalize.chat` both ways, `netP`
feature and `pathway_name` mapping, three-object list merge (`Condition` per cell), and
both deprecated wrappers still working.

Regression: projectData 19, computeCommunProb 61, computeAvgCommunProb 162,
filterCommunication 22, normalizeData 9, pre-inference 15, class 22, accessors 22,
cpp-prob-layer 20, helpers 16 — all green. All 13 files under `R/` parse.
`R CMD INSTALL` still stops only on the four pre-existing undefined exports.
