# Agent Note: Merged LR consumers use dataset-local metadata

Status: implemented
Governance: v1
Date: 2026-11-10
Decision type: migration
Scope: R/analysis.R extractEnrichedLR and merged visualization consumers
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

Merged `extractEnrichedLR()` applied the top-level union `LR$LRsig` to every dataset snapshot. This loses the dataset-local candidate boundary even though the producer retains each snapshot's own LR table. Further, `extractEnrichedLR_internal()` delegated exact pathway selection to `searchPair()`, whose `future_sapply()` simplification is unstable for a one-interaction result. A pathway present only in some datasets could error or yield incorrect/empty candidates instead of producing the union of eligible snapshot results.

## Decision

Use each merged snapshot's own `LR$LRsig`, skip a snapshot when it does not contain the requested pathway, and error only when a requested pathway is absent from the merged LR union. Replace the helper's exact pathway lookup with direct data-frame filtering. Keep single-object behavior, `LR.sig` intersection, `thresh` semantics, gene expansion and merged union return shape unchanged.

## Constraints and invariants

- Snapshot-local LR metadata determines which LR candidates belong to that dataset.
- `enriched.only=TRUE` remains the existing candidate/prob-layer/`LR.sig` intersection; this change does not turn `thresh` into a fresh p-value gate.
- `enriched.only=FALSE` remains all LRsig candidates for the requested pathway within each snapshot.
- The merged result remains the union across snapshots, with genes derived only from returned LR candidates.
- A one-row pathway table must return the same stable data-frame schema as a multi-row table.

## Alternatives considered

Continue using the top-level union for every dataset: rejected because it erases dataset-local metadata and makes the snapshot LR tables irrelevant to pathway extraction. Keep `searchPair()` and special-case one-row output: rejected because it adds branching around an asynchronous helper for a task that is a simple exact column filter.

## Consumer impact

Producer: `mergeSpatialCellChat()` in `R/merging.R`, which stores dataset-local `LR$LRsig`. Consumer: `extractEnrichedLR()` in `R/analysis.R`; downstream consumers include `spatialFeaturePlot()` and `plotGeneExpression()` in `R/visualization.R`. Tests: `tests_dev/test-extractEnrichedLR.R` and `tests_dev/test-merge-spatialcellchat.R`. Documentation: `docs/ANALYSIS_SEMANTICS.md` and `docs/remaining-migration-issues-2026-11-10.md`. No generated public API, configuration, or dynamic lookup is changed.

## Consequences

Merged results now follow source dataset metadata and tolerate pathway absence in individual snapshots. Direct filtering removes `future.apply` dispatch and shape simplification from the helper's candidate selection. Public function signatures and single-object return types are unchanged. Existing cached-selection behavior and its threshold limitations remain; callers still receive a cross-dataset union rather than dataset-tagged results.

## Evidence

The producer contract and snapshot-local LR tables are implemented in `R/merging.R::mergeSpatialCellChat()`. The prior consumer used `object@LR` in each merged iteration in `R/analysis.R::extractEnrichedLR()`. Current regressions cover one-row selection, empty cache intersection, heterogeneous dataset LR tables and missing pathway snapshots.

## Acceptance criteria

`Rscript tests_dev/test-extractEnrichedLR.R` passes the single-object and heterogeneous snapshot checks. `Rscript tests_dev/test-merge-spatialcellchat.R` passes merged LR extraction through a real producer, including a dataset that lacks the requested pathway. `docs/ANALYSIS_SEMANTICS.md` describes the verified cache and threshold behavior without claiming a fresh significance test.
