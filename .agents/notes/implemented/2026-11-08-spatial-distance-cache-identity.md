# Agent Note: 空间距离缓存必须绑定坐标身份

Status: implemented
Governance: v1
Date: 2026-11-08
Decision type: architecture
Scope: R/spatial.R, R/modeling.R, R/SpatialCellChat_class.R, R/visualization.R, distance-cache tests
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

`images$.distance` is generated from `images$coordinates`, but the cache matcher
currently compares only `interaction.range`, `contact.range`, `ratio`, and `tol`.
Updating coordinates through `images<-` leaves the old cache in place, so both
`spatialC2CDistPlot(enriched.only = FALSE)` and `computeCommunProb()` can reuse
distances and contact edges for a different coordinate frame. The same defect is
observable when a caller mutates the coordinate slot directly because the cache
matcher has no coordinate identity check.

Evidence: `R/spatial.R:.spatial_distance_cache_matches()` and
`.sc_c2c_geometry_cache()`, `R/modeling.R:265-285`, and the 15-cell reproduction
where moving one cell outside the interaction range leaves the matcher returning
`TRUE` and the old edge present.

## Decision

Coordinate identity is part of the distance-cache contract. `computeCellDistance()`
normalizes and stores the coordinates used to build the cache. The single
`.spatial_distance_cache_matches()` helper requires matching coordinates,
coordinate row order/names, matrix dimensions, and distance parameters. Both
`spatialC2CDistPlot()` geometry mode and `computeCommunProb()` use that helper.
The `images<-` replacement resets `images$.distance` to an empty cache when
`coordinates` or `spatial.factors` changes, while the matcher remains the defense
against direct slot mutation and legacy caches without coordinate metadata.
`computeCommunProb()` records the coordinate matrix, distance parameters, and
spatial factors used for inference. Enriched distance plots require these to match
the current cache and spatial data, preventing a freshly rebuilt distance cache from
being combined with probability layers inferred from older coordinates.

## Constraints and invariants

- A cache is reusable only when its `d.spatial` and `adj.contact` matrices were
  built from the current coordinate values in the current cell order.
- `interaction.range`, `contact.range`, `ratio`, and normalized `tol` must still
  match exactly under the existing numeric comparison semantics.
- A missing coordinate identity or an old cache without the new metadata is stale
  and must be recomputed, never silently accepted.
- Enriched plots reject missing or mismatched inference provenance, including when
  only the distance cache has been rebuilt after coordinates change.
- Coordinate or spatial-factor updates through `images<-` must not leave a valid
  looking distance cache behind.
- `computeCellDistance()` remains sparse and must not materialize a dense `n x n`
  distance matrix.
- Existing public distance outputs (`d.spatial`, `adj.contact`, `.parameters`) and
  distance values remain unchanged for unchanged inputs.

## Alternatives considered

- Clearing the cache only in `images<-`: rejected as insufficient because direct
  slot mutation remains possible in this S4 codebase and enriched consumers need a
  final validity check.
- Recomputing on every plot and inference call: rejected because it removes the
  correctness bug at a substantial cost for large spatial datasets and defeats the
  existing sparse cache design.
- Adding a new hashing dependency: rejected because the exact normalized coordinate
  matrix is small relative to the distance/probability data and avoids a new package
  dependency or collision semantics.

## Consumer impact

Producer: `R/spatial.R:computeCellDistance()`.

Consumers: `R/spatial.R:.sc_c2c_geometry_cache()` and
`.sc_c2c_inference_cache()`, `R/modeling.R:computeCommunProb()`, and
`R/visualization.R:spatialLeePlot()`.

Validation and regression coverage: `tests_dev/test-computeCellDistance.R`,
`tests_dev/test-comm-distance.R`, and spatial object/accessor tests that exercise
`images<-`. Existing distance-cache fixtures and RDS objects without coordinate
metadata are intentionally treated as stale and rebuilt.

## Consequences

Distance caches gain one hidden metadata field containing the normalized coordinate
matrix. This adds only O(n * dimensions) storage and an O(n * dimensions) comparison
when validating a cache, which is small compared with the sparse distance and
communication layers. Old in-memory or serialized caches without the metadata are
not backward-compatible for reuse and are rebuilt automatically by geometry and
inference. Enriched plots reject legacy inference without provenance and require
rerunning `computeCommunProb()`. Changing spatial coordinates likewise requires
rerunning inference before enriched plots can be trusted.

## Evidence

 - `R/spatial.R` cache producer, matcher, and geometry consumer.
 - `R/modeling.R` independent inference cache check.
 - `R/SpatialCellChat_class.R` image setter and cache validation.
 - `R/visualization.R` Lee-plot cache consumer.
 - Reproduction: old cache had 6 entries while recomputation after moving one cell
   had 2; the old matcher still returned `TRUE`.
 - R 4.4.3: `tests_dev/test-computeCellDistance.R`,
   `tests_dev/test-comm-distance.R`, `tests_dev/test-computeCommunProb.R`,
   `tests_dev/test-SpatialCellChat-class.R`, and
   `tests_dev/test-SpatialCellChat-accessors.R` passed.
 - `tests_dev/test-SpatialCellChat-spatial-image.R` passed its executed checks;
   its existing optional `scatterpie` check was skipped because it was unavailable.

## Acceptance criteria

- The matcher accepts an unchanged coordinate matrix and rejects changed values,
  changed row order/names, changed dimensions, changed parameters, and caches that
  lack coordinate metadata.
- Updating coordinates or spatial factors through `images<-` leaves an empty
  `.distance` cache and the next geometry request rebuilds it.
- `computeCommunProb()` rebuilds rather than reuses a cache after coordinates change.
- A fresh distance cache cannot make an enriched plot accept probability layers
  inferred using earlier coordinates.
- `spatialLeePlot()` passes the original analysis coordinate frame to the same matcher.
- `tests_dev/test-computeCellDistance.R` and `tests_dev/test-comm-distance.R` pass
  under R 4.4.3. The governed-note validator accepts this Note; repository-wide it
  still reports pre-existing errors in six other Notes.
