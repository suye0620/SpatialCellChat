# Agent Note: merged SpatialCellChat producer and dataset snapshot contract

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: schema
Scope: `mergeSpatialCellChat()`, `SpatialCellChat@misc$.datasets`, merged comparison consumers in `R/analysis.R` and `R/visualization.R`, tests_dev fixtures
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

Before this decision was implemented, `NAMESPACE` exported `mergeSpatialCellChat` and comparison consumers expected merged objects, but no producer existed under `R/`. Hand-built mode changes could make dataset, group, LR, and cell axes disagree while still producing a plot. The 11-slot schema required a producer that preserved sparse communication layers and made name alignment explicit.

## Decision

Implement `mergeSpatialCellChat()` as the only supported producer for merged `SpatialCellChat` objects in this migration. Accept a named list of at least two valid single objects, require one common datatype and identical group levels, prefix cell identifiers by dataset by default, and store per-dataset snapshots under `misc$.datasets[[dataset]]`. Keep the top-level object structurally valid with empty merged communication branches rather than fabricating a joint probability tensor. Align LR/pathway layers by exact names and synchronize cell-level sparse matrix dimnames with any prefixed cell identifiers.

## Constraints and invariants

- Inputs must be valid `SpatialCellChat` objects in `single` mode; dataset names must be unique and non-empty.
- `datatype` must be identical across inputs. Spatial coordinate systems, coordinates, rasters, and spatial factors are preserved only inside per-dataset snapshots unless a common coordinate frame is explicitly available; the merged top-level object uses `images$coordinate.scope = "dataset-local"` and does not expose concatenated coordinates.
- Spatial plotting entry points reject a merged object with dataset-local coordinates; users plot the original single-dataset objects separately unless a shared transform is explicitly established.
- Group labels are aligned by name, not position. This first implementation requires identical group levels and fails clearly otherwise.
- Cell identifiers are unique in the merged object. If `cell.prefix = TRUE`, every snapshot cell-level matrix uses the same prefixed names as the merged assay and metadata; group-level matrices retain group names.
- Sparse communication layers remain `SparseChatArray` objects. No full cell-by-cell-by-LR dense tensor may be materialized.
- LR/pathway layer names are preserved exactly within snapshots. Missing optional branches are allowed only when the source object already validates them; no zero-valued fake communication is inserted.
- The returned object has `misc$.mode = "merged"`, named `misc$.datasets`, and passes `validateSpatialCellChat()`/`validObject()`.

## Alternatives considered

1. Reuse or recreate the historical `mergeCellChat()` layout. Rejected because it is not present in the current package and conflicts with the 11-slot schema and sparse layer contract.
2. Allow arbitrary group unions with zero padding. Deferred because silently changing group composition changes comparison denominators and requires an explicit scientific decision about absent groups.
3. Combine all dataset networks into a fabricated top-level joint network. Rejected because the top-level network has no single-condition meaning and would invite consumers to compare the wrong axis.

## Consumer impact

- Producers: `mergeSpatialCellChat()` and future subset/update helpers.
- Consumers: `rankNet()` snapshots, merged `netVisual_diffInteraction()`, heatmap/barplot/bubble comparison, pairwise embedding, and merged analysis functions.
- Tests: `tests_dev/test-analysis-network-migration.R` and a new merged producer regression script using named snapshots and prefixed sparse matrices.
- Documentation: `docs/DATA_STRUCTURE.md`, `man/mergeSpatialCellChat.Rd` through roxygen, and P2 migration issue/review records.
- Dynamic lookups: `object@misc$.mode`, `object@misc$.datasets`, `meta$datasets`, and all comparison index mapping must use the same dataset order.

## Consequences

The producer enables executable comparison consumers through dataset-local snapshots. Requiring identical gene order and group levels narrows the API but prevents silent axis corruption; group union or gene intersection policies require a separate scientific decision. The existing Visium PP1 and Slide-seq MBM05 fixtures remain unsuitable for independent biological comparison because their full gene axes differ and their group labels have no overlap. A Visium self-comparison validates a zero-difference technical path only.

## Evidence

- `R/merging.R::mergeSpatialCellChat()` creates a valid 11-slot merged container with named snapshots, aligned assay/meta/idents, cell-level sparse layer names, dataset-local spatial state, and safe snapshot embeddings.
- `R/analysis.R` and `R/visualization.R` consume those snapshots for rank comparisons, count/weight differences, pairwise embedding, communication extraction, and signaling-role comparisons.
- `docs/DATA_STRUCTURE.md` and `.agents/notes/implemented/2026-08-15-final-11-slot-schema.md` define the current slot and sparse-array contracts.
- `tests_dev/test-merge-spatialcellchat.R` covers producer validation errors, assay and embedding axes, sparse LR layers, count/weight group aggregation, pairwise consumers, and rejection of missing merged measures rather than falling back to ordinary matrices.
- `tests_dev/audit-merge-real-fixtures.R` documents why the available independent Visium and Slide-seq fixtures cannot be merged under the named-axis contract. `tests_dev/render-merged-visium-self-control.R` exercises ten comparison consumers through the real producer as a same-sample technical control.

## Acceptance criteria

1. `tests_dev/test-merge-spatialcellchat.R` rejects invalid lists, duplicate dataset names, datatype/group mismatches and duplicate unprefixed cell IDs; it verifies aligned metadata, idents, assay layers, safe embeddings and named snapshots.
2. The same test checks prefixed cell-level `SparseChatArray` axes, unchanged LR layer names, named group axes, spatial snapshot-local coordinates, `validObject()`, and explicit spatial-plot rejection.
3. Real producer consumers pass for rankNet comparison, pairwise similarity/embedding/zoom, merged communication extraction, barplot/diffInteraction/bubble/heatmap and signaling-role comparisons.
4. `mergeInteractions()` count and weight matrices preserve mapped axes and values; missing `count.merged` or `weight.merged` is an explicit error rather than fallback.
5. `tests_dev/render-merged-visium-self-control.R` passes all ten comparison consumers through the real producer on a real 905-spot object. This verifies a technical zero-difference control, not independent biological replication.
6. Independent-sample comparison remains unverified because the existing real fixtures fail the exact gene/group axis contract; do not pad or relabel them implicitly.
7. The current tarball builds, excludes development paths and compiled artifacts, installs into an isolated library, and passes source-package check with its remaining warnings/notes and unavailable Suggests reported separately. Agent Notes validator results are reported separately; historical validator errors are not presented as a clean gate.

Verified on 2026-11-11: `test-merge-spatialcellchat.R`, `test-remaining-migration.R`, `test-spatial-p1-migration.R`, and `test-analysis-network-migration.R` pass. The source tarball contains 223 entries and excludes development paths and local binaries. `R CMD check --no-manual --no-build-vignettes` exits 0 with 2 WARNINGs and 6 NOTEs; isolated installation, `test-streamline-installed.R`, and `test-installed-visualization-smoke.R` pass. The governance validator reports 20 pre-existing errors and 12 legacy warnings across the repository; this Note validates successfully.
