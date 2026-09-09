# Agent Note: Grid 操作可扩展性与诊断契约

Status: implemented
Governance: v1
Date: 2026-09-02
Decision type: performance
Scope: `R/spatial.R`, `R/visualization.R`, `R/utilities.R` CLI diagnostics, grid regression tests, generated help, and `docs/COMPUTATION_OPTIMIZATION.md`
Owner: project maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

Before this decision, the grid paths duplicated default-size calculations based on dense `Rfast::Dist()` matrices. `computeGridSize()` assumed legacy coordinate columns in its plotting path, despite the canonical object schema storing analysis coordinates in `images$coordinates` with `x`/`y` axes. Dense point-grid membership (`st_intersects(..., sparse = FALSE)`) created an avoidable `n × G` logical allocation. The affected paths were `computeGridSize()`, `makeGridSpatialCellChat()`, `netVisual_CommunFieldGrid()`, and `netVisual_CommunFlow()`.

Baseline evidence: `R/spatial.R` contained `computeGridSize()` and the sparse aggregation path of `makeGridSpatialCellChat()`; `R/visualization.R` contained two duplicated `Rfast::Dist()` blocks; `R/utilities.R` routed `identifyOverExpressedGenes(do.grid = TRUE)` through `makeGridSpatialCellChat()`; existing regression coverage was in `tests_dev/test-makeGridSpatialCellChat.R` and the class/image scripts. The performance audit recorded the unresolved dense grid cost in `docs/COMPUTATION_OPTIMIZATION.md`.

## Decision

SpatialCellChat now uses a private grid-size resolver and sparse membership helper across all four grid consumers without changing their public contracts.

1. `.sc_resolve_grid_size()` in `R/spatial.R` accepts finite 2D/3D coordinate matrices, validated explicit scalar or length-two `cellsize`, and positive `grid.resolution`. When `cellsize = NULL`, it derives the minimum non-self Euclidean distance without materializing a dense distance matrix and records base/effective size, resolution, ratio, source, and method metadata.
2. `.sc_grid_membership(points, grid)` uses `sf::st_intersects(..., sparse = TRUE)` and returns hit lists, point/grid hit counts, and point/grid indices while preserving edge/corner multi-grid membership and grid ordering.
3. `computeGridSize()` operates on canonical `images$coordinates`, retains its plot/non-plot return behavior, uses sparse membership diagnostics, and emits `.cli()` diagnostics through `SpatialCellChat.verbose` levels 0–3.
4. `makeGridSpatialCellChat()` consumes the shared membership result while retaining occupied-grid filtering, `Grid<id>` naming, identity tie handling, counts, centroids, expression `rowMeans`, `images$.grid$within.nGrid`, and `recommended.contact.range`.
5. `netVisual_CommunFieldGrid()` and `netVisual_CommunFlow()` consume the shared size resolver and sparse membership mechanics while retaining coordinate preparation, interpolation, aggregation, empty-grid behavior, vector direction, and plot class.
6. Roxygen/generated help and the performance audit describe the implemented behavior. The implementation leaves `renv.lock` unchanged and adds no dependency.

Implementation followed the recorded order: baseline/contract tests; resolver; sparse membership and `computeGridSize()`; `makeGridSpatialCellChat()` migration; visualization migration; CLI contract; documentation; final exactness and 10k/100k scale checks. This sequence is retained as implementation provenance, not as a second authoritative plan.

## Constraints and invariants

- Do not modify, move, or replace any existing exact `sf::st_crs(df.sf) <- 3857` statement.
- `makeGridSpatialCellChat()` output and `identifyOverExpressedGenes(do.grid = TRUE)` behavior remain unchanged.
- `computeGridSize(do.plot = TRUE)` returns a `ggplot`; `do.plot = FALSE` remains an invisible `NULL` diagnostic path.
- Default sizing remains the minimum non-self Euclidean distance, multiplied component-wise by `grid.resolution`; the implementation changes, not the statistic.
- Duplicate coordinates preserve zero nearest-neighbor distance semantics; fewer than two points fail clearly rather than propagating `Inf`.
- Canonical grid geometry uses the first two coordinate axes (`x`, `y`) without a new silent axis swap or y-direction transformation.
- `sf::st_intersects()` edge/corner behavior is preserved, including multi-grid hits.
- No dense `n × n` distance matrix or dense `n × G` point-grid matrix is introduced.
- `cellsize` rejects zero, negative, non-finite, and unsupported lengths; `grid.resolution` is one positive finite number.
- `ratio = NULL` or invalid calibration is reported as uncalibrated; no fabricated physical unit is printed.
- Existing `.cli()` is the only verbosity gateway: level 0 warnings/errors, level 1 concise summary, level 2 grid diagnostics, level 3 algorithm/allocation details. Warnings/errors remain visible at level 0.
- No dependency or `renv.lock` changes.

## Alternatives considered

- Keep separate `Rfast::Dist()` and aggregation logic in each function: rejected because it repeats dense allocations and lets default-size semantics drift.
- Replace all grid geometry construction with one large geometry abstraction: rejected because it risks moving tested CRS assignments and changing visualization-specific coordinate behavior; only sizing and membership are shared.
- Return a structured diagnostic object from `computeGridSize()`: rejected for this change because existing callers and the public contract require `ggplot` for `do.plot = TRUE` and invisible `NULL` otherwise.
- Use a dense logical intersection matrix for simpler code: rejected because its memory cost scales as `n × G` and is the defect being removed.

## Consumer impact

- Producers: `computeGridSize()` and `makeGridSpatialCellChat()` in `R/spatial.R`; field/flow visualization functions in `R/visualization.R`.
- Consumers: `identifyOverExpressedGenes(do.grid = TRUE)` in `R/utilities.R`; visualization callers and users of `images$.grid` metadata.
- Tests: existing `tests_dev/test-makeGridSpatialCellChat.R`, `tests_dev/test-pre-inference-new-schema.R`, `tests_dev/test-SpatialCellChat-class.R`, and `tests_dev/test-SpatialCellChat-spatial-image.R`; direct contracts are covered by `tests_dev/test-computeGridSize.R`, `tests_dev/test-grid-visualization.R`, and `tests_dev/test-grid-scale-100k.R`.
- Documentation/generated artifacts: roxygen in `R/spatial.R` and `R/visualization.R`, generated `man/computeGridSize.Rd`, `man/makeGridSpatialCellChat.Rd`, relevant visualization `.Rd` files, and `docs/COMPUTATION_OPTIMIZATION.md`.
- Configuration: `SpatialCellChat.verbose` is reused; no new global option. `DESCRIPTION`/renv dependencies remain unchanged. No other runtime callers of `computeGridSize()` were found in the pre-change repository search.

## Consequences

The grid paths avoid dense distance and point-grid allocations and share one size/membership contract. Small inputs should remain numerically equivalent, including boundary and duplicate-coordinate behavior. Migration touches implementation, tests, generated help, and performance evidence; it must retain the existing geometry and coordinate conventions rather than silently normalize them. Large polygon-grid construction and downstream vector-field interpolation may remain material costs. Windows may not provide portable true peak RSS evidence, so performance reporting must distinguish elapsed time and R-level allocation deltas from unsupported process peak memory claims.

## Evidence

- `R/spatial.R` current `computeGridSize()` and `makeGridSpatialCellChat()` implementations, including exact CRS assignments and existing sparse aggregation.
- `R/visualization.R` current `netVisual_CommunFieldGrid()` and `netVisual_CommunFlow()` duplicated default-size paths.
- `R/utilities.R` `identifyOverExpressedGenes()` grid call and `.cli()` wrapper.
- `tests_dev/test-makeGridSpatialCellChat.R` existing sparse/dense equivalence, boundary, simulated-data, and stress contracts.
- `docs/COMPUTATION_OPTIMIZATION.md` outstanding grid allocation findings.
- `DESCRIPTION`, `renv.lock`, and the approved dependency baseline note confirm required packages are already declared; this decision adds no dependency.
- `tests_dev/test-computeGridSize.R` passed deterministic 2D/3D nearest-neighbor equality, duplicate-coordinate zero distance, explicit-size single-point behavior, validation errors, canonical plotting, and verbosity 0--3 diagnostics.
- `tests_dev/test-makeGridSpatialCellChat.R` passed square, hexagonal, boundary, simulated, pre-inference-compatible, and 10,000-point sparse equivalence checks; the 10k run reported 3.980 s and `Ncells_delta=2.6 MB`.
- `tests_dev/test-grid-scale-100k.R` passed with 100,172 points, 99,540 occupied grids, and 398,160 sparse hits in 168.070 s with `Ncells_delta=18.1 MB`; the result validated and retained a `dgCMatrix` assay.
- `tests_dev/test-SpatialCellChat-class.R`, `tests_dev/test-SpatialCellChat-spatial-image.R`, `tests_dev/test-normalizeData.R`, `tests_dev/test-SpatialCellChat-accessors.R`, `tests_dev/test-computeCellDistance.R`, and `tests_dev/test-pre-inference-new-schema.R` all passed on the final tree. Source and target `.Rd` parsing passed, and repository searches found no grid-path `Rfast::Dist()` or dense point-grid intersection.
- `tests_dev/test-grid-visualization.R` passed the explicit/implicit cellsize Field/Flow rendering smoke cases, adaptive arrow/streamline density, magnitude-to-width mapping, and input immutability checks. A current-source reconstruction of `tests_dev/simulated_data/SpatialChat_2.rds` rendered both optimized functions for `eL5-eR6` successfully (`1050 x 1050` PNG outputs); the render script also confirmed `input_unchanged=TRUE`.

## Acceptance criteria

- Repository search shows no `Rfast::Dist(` in the four grid-related functions and no dense point-grid intersection in `computeGridSize()`.
- The nearest-neighbor default equals `min(as.vector(dist(coords)))` within `1e-12` on deterministic 2D and 3D fixtures, including duplicate-coordinate zero distance behavior; invalid inputs fail clearly.
- `computeGridSize()` works with canonical `x`/`y` coordinates, returns `ggplot` when plotting, invisible `NULL` otherwise, and preserves exact CRS assignment statements.
- Square, hexagonal, boundary, simulated, and 10k-cell `makeGridSpatialCellChat()` results remain numerically/schema equivalent to the existing dense reference; `identifyOverExpressedGenes(do.grid = TRUE)` remains valid.
- Field/flow visualization smoke cases with explicit/implicit cellsize and integer/double resolution preserve returned plot class and existing coordinate/vector semantics.
- Captured CLI output at verbosity 0–3 reports size source, effective size, resolution, geometry, generated/occupied/empty grids, assigned/unassigned/multi-hit points, calibration state, and next-step usage guidance; warnings remain visible at level 0.
- Roxygen and generated help describe canonical coordinates, size validation, nearest-neighbor defaults, resolution, geometry, return behavior, verbosity, and uncalibrated ratios; stale old-schema terms are removed.
- Deterministic 10k and 100k sparse checks complete without an `n × n` or `n × G` dense allocation and record elapsed/allocation evidence. Unsupported peak-RSS claims are not made.

## Validation gate

The targeted exactness, class/image, visualization, documentation, 10k, and 100k validation checks recorded above completed successfully. The decision is therefore implemented.
