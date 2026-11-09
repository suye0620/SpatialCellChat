# Agent Note: Spatial feature plot schema migration

Status: implemented
Governance: v1
Date: 2026-11-09
Decision type: migration
Scope: spatialFeaturePlot expression and LR selection
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The function reads deleted expression slots and old network fields, and binds cells to coordinates by position.

## Decision

Use the current assay and network schema while retaining the existing `slot.name` values as aliases for `signaling` and `norm`; also accept the current layer names. Align input cells by name and preserve the distinction between pathway cache selection and explicit LR selection.

## Constraints and invariants

- Default to `assay$signaling`; `data` selects `assay$norm`.
- Explicit LR selection requires a positive probability; exclude group p-values above `thresh` when available. Cell selection does not use group p-values.
- Pathway selection continues using the existing `extractEnrichedLR()` cache semantics.
- Each coordinate receives the value of the same named cell; no probability array is densified.

## Alternatives considered

Removing old argument values breaks existing calls. Using cached `LR.sig` for explicit pairs changes the current-threshold behavior. Densifying the probability array scales poorly with cell count.

## Consumer impact

Producers: `createSpatialCellChat()`, `subsetData()`, `computeCommunProb()`, `computeAvgCommunProb()` and `aggregateNet()`. Consumers: users and `tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.Rmd`; repository search found no other direct R callers. Update roxygen, `man/spatialFeaturePlot.Rd`, semantics documentation and regression tests. No interface removal, configuration or dynamic lookup identified.

## Consequences

Old argument values remain accepted without accessing removed slots. Explicit LR selection iterates sparse matrices instead of materializing the cell-by-cell-by-LR array. Pathway cache semantics and expression meanings remain unchanged.

## Evidence

`R/SpatialCellChat_class.R` defines the schema; `R/utilities.R` creates signaling from norm; `R/visualization.R` contains the plot logic; `R/analysis.R` defines pathway cache selection. `tests_dev/test-spatialFeaturePlot-schema.R` passes the expression, selection, alignment, complex and plotting checks. `test-extractEnrichedLR.R`, `test-plotGeneExpression-schema.R`, `test-SpatialCellChat-class.R`, `test-netpathway-aggregate.R`, and `test-comm-distance.R` pass. `tools::checkRd()` and R parsing pass.

## Acceptance criteria

`tests_dev/test-spatialFeaturePlot-schema.R` verifies expression mappings, coordinate alignment, group p-value boundary, sparse explicit-zero and cell selection, pathway cache selection, binary and continuous plots, and missing prerequisites. Existing extraction, plotting, class, aggregation and distance regressions pass; `git diff --check` passes.
