# Agent Note: Dataset-aligned DEG mapping for merged communication results

Status: proposed
Governance: v1
Date: 2026-11-10
Decision type: migration
Scope: netMappingDEG(), identifyOverExpressedGenes(), merged misc$.datasets snapshots, tests_dev, and API documentation
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`mergeSpatialCellChat()` preserves each source object's `misc$.var.features` inside its named dataset snapshot but initializes the merged object's top-level `misc$.var.features` as empty (`R/merging.R`). `netMappingDEG()` previously read only the top-level table. Consequently, DEG mapping on a newly produced merged object failed even when each source dataset had its marker table, and using one dataset's table for every dataset would misattribute results.

## Proposal

For merged objects, map each dataset's communication table only with DEG rows aligned to that dataset: prefer an explicit top-level DEG table carrying a `datasets` column and filter by dataset; otherwise read that dataset's snapshot-local `misc$.var.features[[paste0(features.name, ".info")]]`. Do not silently reuse an unstratified top-level marker table across datasets. Single-object mapping retains the existing marker table and mapping algorithm.

## Constraints and invariants

- Communication edges remain sourced from each named `misc$.datasets[[dataset]]` snapshot.
- A DEG row may annotate only the dataset whose identity is recorded in the row or whose snapshot contains the marker table.
- Ligand/receptor matching, complex-subunit aggregation, and p-value/threshold calculations remain unchanged.
- Missing dataset-aligned marker evidence produces an actionable error rather than silently returning unannotated or cross-dataset annotations.

## Alternatives considered

1. Reuse one top-level unstratified DEG table for every dataset. Rejected because merged network rows would appear dataset-specific while their annotations could come from a pooled or unrelated comparison.
2. Require users to call `netMappingDEG()` separately on each source object and omit merged support. Rejected because the merged producer already preserves source marker tables and the function can safely map named snapshots independently.

## Consumer impact

Producer: `mergeSpatialCellChat()` snapshots source misc; `identifyOverExpressedGenes()` may also produce a merged top-level table with `datasets` when testing a positive dataset. Consumer: `netMappingDEG()`; downstream DEG filters in `subsetCommunication()`. Tests: `tests_dev/test-remaining-migration.R`. Docs: `man/netMappingDEG.Rd`, `docs/ANALYSIS_SEMANTICS.md`, and `docs/remaining-migration-issues-2026-11-10.md`. No schema or serialization change is required.

## Consequences

Single-object behavior is preserved. Merged mapping becomes supported when per-dataset marker evidence exists; merged calls with an unstratified pooled table fail explicitly and require dataset-specific evidence. Runtime remains linear in the combined communication table and uses existing data-frame operations.

## Evidence

`R/merging.R` stores source misc in `misc$.datasets`; `R/analysis.R::netMappingDEG()` reads the merged marker list; `R/utilities.R::identifyOverExpressedGenes()` can attach a `datasets` column when testing a positive dataset; `subsetCommunication()` returns one named table per merged dataset. Existing single-object mapping regression is in `tests_dev/test-remaining-migration.R`.

## Acceptance criteria

1. A merged object with distinct marker tables in its snapshots returns communication rows annotated with the matching dataset's distinct values.
2. A merged top-level marker table with a `datasets` column is filtered independently for each dataset.
3. A merged call with neither dataset-local tables nor a dataset-stratified top-level table fails with a message naming the missing dataset and required producer.
4. Existing single-object marker mapping remains unchanged.
5. Run the focused remaining migration test, the merged producer regression, package installation, and the Agent Notes validator. Keep this Note proposed until all criteria pass and maintainer review adopts the contract.
