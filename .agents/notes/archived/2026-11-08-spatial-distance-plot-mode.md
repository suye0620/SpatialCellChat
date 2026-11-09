# Agent Note: Name distance plot selection modes explicitly

Status: superseded
Governance: v1
Date: 2026-11-08
Decision type: api
Scope: spatialC2CDistPlot public API, callers, docs and tests
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: .agents/notes/implemented/2026-11-08-undirected-cell-pair-distance-density.md

## Problem

`R/spatial.R:spatialC2CDistPlot()` calls its record selection `enriched.only`, although this is not an enrichment test. The name suggests a statistical result that the implementation does not establish.

## Decision

Replace the public flag with `mode = c("geometry", "communication")`, in that order, before `split.by`. Do not change the underlying record selection or KDE.

## Constraints and invariants

- Geometry remains the default and needs no inference. Communication retains its provenance checks and cell-pair-by-layer records. Both splits retain their behavior.
- Explicit `split.by` in geometry mode still warns and is ignored.
- This decision does not introduce probability weighting or fix explicit zeros.

## Alternatives considered

Keeping the Boolean retains ambiguous terminology. A compatibility alias would extend the ambiguous public API; this refactor is deliberately breaking.

## Consumer impact

Interface: `R/spatial.R:spatialC2CDistPlot`; help: `man/spatialC2CDistPlot.Rd`; callers: `tests_dev/render_viz_figures.R` and `tests_dev/audit_visualization_findings.R`; tests: `tests_dev/test-comm-distance.R` and `tests_dev/test-SpatialCellChat-spatial-image.R`. Repository text search found no other live calls or dynamic lookups. Historical Notes preserve prior decisions; other functions' `enriched.only` are distinct.

## Consequences

Old named and positional calls must migrate; no compatibility alias is provided. No performance or density semantics change is intended.

## Evidence

`R/spatial.R:.sc_comm_distance_frame()` assembles records for stored layer entries and distances. `R/spatial.R:.sc_c2c_inference_cache()` checks provenance. On R 4.4.3, `tests_dev/test-comm-distance.R` and `tests_dev/test-SpatialCellChat-spatial-image.R` passed; the optional scatterpie check was skipped. `tools::checkRd("man/spatialC2CDistPlot.Rd")` and `git diff --check` passed. The governed-note validator accepts this Note; it reports 20 pre-existing errors and 13 warnings elsewhere.

## Acceptance criteria

- Both modes, default, invalid mode, old argument rejection, splits and stale inference are checked in `tests_dev/test-comm-distance.R`.
- Spatial image integration checks pass and active scripts use the new API.
- `tools::checkRd()` accepts the help file; governed-note validation is run.
