# Agent Note: Use undirected cell pairs as distance observations

Status: implemented
Governance: v1
Date: 2026-11-08
Decision type: algorithm
Scope: spatialC2CDistPlot communication and geometry distance observation construction
Owner: maintainer
Impact: high
Supersedes: .agents/notes/archived/2026-11-08-spatial-distance-plot-mode.md
Superseded by: none

## Problem

The distance plot currently counts both matrix orientations of a geometric pair and repeats a pair once per stored LR layer. That makes the curve depend on matrix orientation and the number of LR records rather than on the intended observation unit: one spatial cell pair.

## Decision

Use one undirected cell pair as one observation in every panel. Geometry mode takes the strict upper triangle of the cached distance matrix. Communication mode first aggregates LR layers by annotation, tests aggregated probability `> 0`, merges the two communication directions with logical OR, and then takes the strict upper triangle. The same directional merge is used for the all-communication distance split.

## Constraints and invariants

- Diagonal entries are excluded in every mode.
- A and B contribute once per panel even when both A→B and B→A have communication or when multiple LR layers are positive.
- Annotation panels aggregate only layers in their own class; the all communication panel aggregates all layers.
- Distances remain bounded by the cached interaction range and contact filtering uses `adj.contact`.
- Probability magnitudes do not become KDE weights.

## Alternatives considered

Keeping directed pairs duplicates symmetric geometry. Keeping one row per LR layer makes the result depend on LR database coverage. Summing probabilities as KDE weights changes the scientific estimand to a probability-mass distance distribution, so it is not used for this plot.

## Consumer impact

Producer and plot helpers: `R/spatial.R` (`.sc_comm_distance_frame`, `.sc_comm_dist_panels`, `.sc_c2c_panels`). Consumers are the public `spatialC2CDistPlot`, its geometry/communication callers, and the distance tests in `tests_dev/test-comm-distance.R` and `tests_dev/test-SpatialCellChat-spatial-image.R`. The public help and visualization semantics must describe undirected cell-pair observations. The prior mode Note is superseded by this decision.

## Consequences

Communication output changes from directed LR records to undirected cell-pair observations, so density values and panel sample counts change intentionally. Per-class `marginSums()` uses the existing C++ layer aggregation and avoids materializing all LR-distance rows. The implementation must validate explicit zeros through `> 0` and retain sparse operations.

## Evidence

`R/spatial.R:computeCellDistance()` constructs symmetric neighbor distances and excludes self/zero-distance edges. `R/SpatialCellChat_class.R:marginSums.SparseChatArray()` delegates multi-layer sums to `src/SpatialChat_Rcpp.cpp:cpp_sum_layers()`. The previous frame builder expanded every stored LR layer entry, which was the source of duplicate observations. The implementation now uses per-class `marginSums()`, a positive-probability mask, logical direction merging, and a strict upper-triangle distance frame. On R 4.4.3, `tests_dev/test-comm-distance.R`, `tests_dev/test-computeCellDistance.R`, and `tests_dev/test-SpatialCellChat-spatial-image.R` passed; the optional scatterpie check was skipped. `tools::checkRd("man/spatialC2CDistPlot.Rd")` and `git diff --check` passed.

## Acceptance criteria

- Geometry frame contains one strict-upper-triangle distance per spatial pair and a logical contact flag.
- Communication frame contains at most one row per unordered pair per annotation class; all-distance mode contains at most one row per unordered pair overall.
- A pair positive in either direction is retained; diagonal and explicit-zero-only entries are excluded.
- Existing cache provenance tests continue to pass, and new fixtures verify duplicate LR layers, reverse directions, contact filtering, and annotation-specific deduplication.
- R documentation, `git diff --check`, and the governance validator pass (apart from known legacy repository errors).
