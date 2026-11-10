# Agent Note: merged SpatialCellChat producer and dataset snapshot contract

Status: proposed
Governance: v1
Date: 2026-11-10
Decision type: schema
Scope: `mergeSpatialCellChat()`, `SpatialCellChat@misc$.datasets`, merged comparison consumers in `R/analysis.R` and `R/visualization.R`, tests_dev fixtures
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

`NAMESPACE` exports `mergeSpatialCellChat`, and comparison consumers already expect a merged object, but no producer exists under `R/`. The current comparison code therefore has no safe way to construct `misc$.mode = "merged"` and `misc$.datasets`. Hand-built mode changes can make dataset, group, LR, and cell axes disagree while still producing a plot. The 11-slot schema requires a producer that preserves sparse communication layers and makes name alignment explicit.

## Proposal

Implement `mergeSpatialCellChat()` as the only supported producer for merged `SpatialCellChat` objects in this migration. Accept a named list of at least two valid single objects, require one common datatype and identical group levels, prefix cell identifiers by dataset by default, and store per-dataset snapshots under `misc$.datasets[[dataset]]`. Keep the top-level object structurally valid with empty merged communication branches rather than fabricating a joint probability tensor. Align LR/pathway layers by exact names and synchronize cell-level sparse matrix dimnames with any prefixed cell identifiers.

## Constraints and invariants

- Inputs must be valid `SpatialCellChat` objects in `single` mode; dataset names must be unique and non-empty.
- `datatype` must be identical across inputs. Spatial coordinate systems and rasters are preserved only inside per-dataset snapshots unless a common coordinate frame is explicitly available; no silent geometric concatenation is allowed.
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

The merged producer makes comparison paths constructible and gives downstream code one authoritative snapshot location. It intentionally does not solve all comparison rendering in the same change. Requiring identical group levels narrows the initial API but prevents silent axis corruption; a future union policy can be added as a separate decision. Prefixing cell names avoids collisions but requires synchronized renaming in every cell-level assay, embedding, raster coordinate, and sparse communication layer.

## Evidence

- `R/SpatialCellChat_class.R` defines the 11-slot schema and validates `misc$.mode`/`misc$.datasets`, but no producer is present.
- `R/analysis.R` and `R/visualization.R` contain comparison consumers that need dataset snapshots.
- `docs/DATA_STRUCTURE.md` and `.agents/notes/implemented/2026-08-15-final-11-slot-schema.md` define the current slot and sparse-array contracts.
- `tests_dev/test-analysis-network-migration.R` currently uses a temporary snapshot mirror because a real producer is absent; that test is not sufficient producer evidence.

## Acceptance criteria

1. `mergeSpatialCellChat()` rejects invalid object lists, duplicate dataset names, datatype mismatches, group-level mismatches, and duplicate resulting cell identifiers with actionable errors.
2. Two or more valid objects produce a valid merged object with aligned metadata, idents, assay layers, embeddings where safe, `meta$datasets`, and named `misc$.datasets` snapshots.
3. Cell-level `SparseChatArray` layers have synchronized prefixed dimnames without dense expansion; group-level layers retain group names; LR layer names are unchanged.
4. A producer-specific regression test passes on simulated objects and exercises the real `rankNet(..., mode = "comparison")` path.
5. `validObject()` and `validateSpatialCellChat()` pass on the returned object; the package source defines every exported function needed for installation.
6. The Note remains `proposed` until comparison consumers and the documented acceptance checks have been implemented and verified; it must not be described as a completed merged-comparison migration before then.