# Agent Note: Canvas htmlwidget 通讯流可视化

Status: implemented
Governance: v1
Date: 2026-09-02
Decision type: architecture
Scope: `netVisual_CommunFlowWidget()`、通讯向量场数据准备层、`inst/htmlwidgets/` 本地前端资产及其测试
Owner: SpatialCellChat maintainers
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`computeCommunField()` stores one signaling layer's cell-level `outgoing` and `incoming` vector fields, but the package currently exposes only static ggplot2 visualizations for the grid/flow paths. The stored field is a steady field derived from one communication inference result, not a time series: an animation can show particle advection over that fixed field, but must not be described as biological temporal evolution or a measured propagation speed. Existing `netVisual_CommunFlow()` also owns useful grid aggregation and Barnes interpolation semantics that a new renderer must preserve.

## Decision

Add `netVisual_CommunFlowWidget()` as an additive public API. Its R side validates the precomputed cell field, applies the same axis swap and y-direction convention as `netVisual_CommunFlow()`, aggregates cell vectors on an `sf` grid, interpolates `u` and `v` with `oce::interpBarnes()`, masks grid locations outside occupied tissue cells, and serializes one signaling layer plus cell points. Its browser side is a self-contained local `htmlwidget` using Canvas 2D, `requestAnimationFrame`, bilinear field sampling, deterministic particle seeding, trail fading, play/pause, visual-speed and particle-count controls, cell labels/colors, and hover metadata.

The widget is explicitly a **particle animation over a fixed communication vector field**. Particle brightness/width encodes local net-current magnitude. It does not infer time, velocity units, or biological transport. The initial implementation carries one `incoming` or `outgoing` layer per widget and does not serialize the original cell-to-cell edge matrix.

## Constraints and invariants

- Existing files under `R/` are not modified; the implementation is an additive R source file plus local widget assets and tests.
- `object` must be a `SpatialCellChat` with a precomputed cell-level `SparseChatArray` field in the requested `net` or `netP` slot; the function never recomputes or mutates the object.
- `signaling` must identify exactly one stored field layer and `pattern` must be exactly `incoming` or `outgoing`.
- Cell IDs, point coordinates, and group labels remain one-to-one and in canonical cell order.
- The visualization coordinate transform remains compatible with `netVisual_CommunFlow()`: the first two canonical axes are swapped for the plot frame, the plotted y coordinate is negated, and the vector y component is negated with it.
- The browser payload contains a bounded regular grid, validity mask, and point table, not an `n x n` communication matrix.
- Default grid cap is `128 x 128`; the public hard cap is `256 x 256`. Particle count is bounded independently of cell count.
- Invalid, empty, all-zero, duplicate-coordinate, and out-of-domain states have explicit R errors or browser reseeding behavior; no silent biological interpretation is added.
- `particle.speed` is a visual coordinate-units-per-animation-second multiplier only. No `m/s` or other physical transport unit is exposed.
- No runtime CDN or network request is required by the widget.

## Alternatives considered

1. **Plotly scatter frames**: useful for a low-cost animation prototype, but requires precomputing positions for every particle and frame, scales payload/repaint cost with particle count times frame count, and has no native two-dimensional continuous particle integrator. It loses to Canvas for the production Earth-style trail renderer.
2. **Leaflet Velocity / `leaflet.extras2::addVelocity()`**: technically accepts regular U/V grids, but current coordinates are tissue pixel/arbitrary coordinates rather than WGS84, and the package is not a current direct dependency. Its map/weather conventions risk implying geographic direction and physical speed; it is not the default tissue renderer.
3. **Modify `netVisual_CommunFlow()` to share a helper**: would reduce duplication, but violates the requested additive-only boundary and creates a larger compatibility surface. The new function mirrors the established extraction contract without changing existing return values or plots.

## Consumer impact

The new exported function is additive. Existing `netVisual_CommunFlow()`, `netVisual_CommunFieldGrid()`, `computeCommunField()`, their callers, existing tests, and serialized object schemas remain unchanged. New consumers are R users who have already run `computeCommunField()`, htmlwidgets hosts (RStudio Viewer, HTML documents, Shiny, or `htmltools::save_html()`), and the new R/browser contract tests. `htmlwidgets` becomes a direct package import; `oce` remains an optional runtime requirement already used by the static flow path. No existing configuration, generated artifact, or dynamic lookup is removed.

## Consequences

The package gains a browser-native interactive renderer and a small local JavaScript/CSS asset surface. R-side preparation repeats the current Flow extraction logic until a later additive helper extraction is approved. The payload is bounded by grid and particle caps, but Barnes interpolation and sf grid construction remain non-trivial for large cell sets. The widget is portable to arbitrary tissue coordinate frames and can be embedded offline, but it does not include a histology raster in this first API. A future multi-layer/Shiny selector must load layers on demand rather than serializing all signaling layers into one default widget.

## Evidence

- `R/analysis.R:3942-4004` computes each cell's probability-weighted net displacement; opposite directions can cancel, so vector magnitude is net direction current.
- `R/analysis.R:4041-4094` requires precomputed cell probabilities and stores `cell$field$outgoing` and `cell$field$incoming` as `SparseChatArray` values without mutating the input object.
- `R/visualization.R:6153-6448` defines the current Flow extraction, grid aggregation, `oce::interpBarnes()` interpolation, tissue mask, coordinate swap, and y/vector negation semantics.
- `R/spatial.R:187-283` provides the existing grid-size resolver and sf membership helper used by the static visualizations.
- `DESCRIPTION:12-19` declares `htmlwidgets` as a direct runtime import; `oce` remains an optional runtime dependency for Barnes interpolation and `tidyr` remains optional.
- `DESCRIPTION:12-19` declares `htmlwidgets` as a direct runtime import; `inst/htmlwidgets/spatialcellchat-commun-flow.yaml` registers the local JavaScript and CSS assets.
- `tests_dev/test-commun-flow-widget.R` passes all contract checks, including immutability, coordinate transform, bounded grid, mass/current separation, rejection paths, and all-zero behavior.
- `node --check inst/htmlwidgets/spatialcellchat-commun-flow.js` passes; browser smoke observed automatic DOM rendering with two canvases, moving particles, stable pause, resumed motion, bounded controls, re-render cancellation, destroy cleanup, and no remote resources. A screenshot confirmed the rendered control surface.
- Official htmlwidgets development guidance: <https://www.htmlwidgets.org/develop_intro.html>.
- Plotly animation reference: <https://plotly.com/javascript/animations/>.
- Leaflet Velocity interface: <https://rdrr.io/cran/leaflet.extras2/man/addVelocity.html>.

## Acceptance criteria

- `netVisual_CommunFlowWidget()` returns an `htmlwidget` for a valid precomputed field and leaves the input object identical.
- The R payload has one cell record per communication cell, canonical labels/colors, a regular `u`/`v` grid, and a same-length logical validity mask; it has no cell-to-cell edge table.
- The function rejects missing fields, unknown signaling layers, invalid patterns, invalid particle/grid caps, and missing optional runtime packages with actionable errors.
- Coordinate/vector direction and tissue masking are numerically equivalent to the existing Flow extraction semantics for a deterministic fixture.
- The local browser binding renders cells and moving particles, pauses without advancing positions, resumes, and cancels its animation loop on re-render/destroy.
- The browser smoke check observes the widget DOM/canvas and changed particle positions without network requests.
- Status transitioned to `implemented` after the R contract test, JavaScript syntax check, package-load target check, browser lifecycle smoke, and Agent Notes validation passed; existing legacy Note warnings remain unchanged.

## Revision 2026-09-03 — static streamline removal

After visual review, the widget's static interpolated streamline overlay was
removed. The static canvas now contains only identity-colored cell points with
white outlines; the trail canvas remains the only animated layer and moving
particles are the only flow-direction cue. This avoids presenting sampled
field curves as an additional biological or edge-level communication path.

The legend and header now say only `particles = signal-flow direction` and
`particle motion = signal flow`. Obsolete `drawStreamlines()`, `sampleLoose()`,
streamline legend CSS, and related documentation were removed from the source
assets and public widget description. Existing cell vectors and net-current
payload values remain available for hover/status behavior and are unchanged.

Evidence:

- `inst/htmlwidgets/spatialcellchat-commun-flow.js` contains no
  `drawStreamlines`, `sampleLoose`, or streamline legend implementation;
  `drawCells()` is called during resize/render and the animation frame only
  calls `drawParticles()`.
- `R/visualization_widget.R` and `man/netVisual_CommunFlowWidget.Rd` describe
  identity-colored cell points and moving particles without static streamlines
  or line-width semantics.
- `tests_dev/test-commun-flow-widget.R` passes all contract checks; the R and
  Rd sources parse successfully; `node --check
  inst/htmlwidgets/spatialcellchat-commun-flow.js` passes.
- The page generated from `tests_dev/simulated_data/SpatialChat_2.rds` selected
  `dL7-dR5`, rendered 2,000 cells and a 46 x 46 grid, and exposed 900
  particles. Headless Chromium observed a stable cell canvas (0 changed bytes
  over 700 ms) and a changing particle canvas (46,942 changed bytes), with no
  streamline code in the page.


## Revision 2026-09-03 — SpatialChat_2 LR validation

The widget was exercised with the real `tests_dev/simulated_data/SpatialChat_2.rds` fixture rather than the 9-cell contract fixture. The selected non-zero layer was `dL7-dR5` (ligand `dL7`, receptor `dR5`), which has 200,322 cell-level non-zero probabilities and total probability mass 824.7772188523128. The fixture contains 2,000 cells, 1,110 genes, and five identity groups.

Because the RDS is serialized in the legacy 14-slot object schema while the current source uses the final 11-slot schema, the test first extracted the fixture's normalized expression, metadata, coordinates, spatial factors, and selected sparse probability layer, then reconstructed a valid current-schema object. No probabilities, coordinates, or identities were synthesized or recomputed.

Evidence:

- `computeCommunField(slot.name = "net", signaling.name = "dL7-dR5", top = 0.8, sparse = TRUE)` completed with 1,500 non-zero outgoing cell vectors.
- `netVisual_CommunFlowWidget(pattern = "outgoing")` produced a bounded 46 x 46 payload grid with 1,323 valid cells and 600 particles.
- `commun-flow-spatialchat2-dL7-dR5.html` rendered two 1498 x 1000 Canvas layers, 2,000 colored cell points, five identity colors, `Sources` legend, and the selected LR title with no browser errors.
- Browser lifecycle smoke confirmed particle positions remained identical while paused and changed after resume.

Remaining risk: `SpatialChat_2.rds` contains simulated labels and simulated LR names, so this validates real fixture scale, sparse LR extraction, coordinate transform, field computation, and rendering—not biological interpretation of `dL7-dR5`. The legacy-to-current extraction bridge is test-only and is not a public migration API.

## Revision 2026-09-03 — final browser verification and docs alignment

The current widget build was re-opened from `commun-flow-current.html` after
embedding the latest JS/CSS and confirmed in headless Chromium:

- two canvases with z-order `sc-flow-cells` below `sc-flow-trails`
- control surface exposes Visual speed, Particles, Point opacity, and Zoom
- pause/resume toggles animation state; particle snapshots remain stable while
  paused and advance after resume
- wheel zoom and slider zoom both update the payload and redraw without console
  errors
- screenshot `C:\Users\Administrator\AppData\Local\Temp\omp-sshots-1575c28f0d0100ba.webp`
  shows the final surface

Residual risk remains the same: the real field may be visually sparse when the
interpolated net current is weak or canceling, but this is faithful to the
stored field rather than a fabricated edge display.

The Rd usage block now matches the implementation defaults for
`image.alpha = 0.32`, `particle.speed = 3`, `trail.length = 24L`, and
`zoom = 1`.

## Revision 2026-09-03 — direct SpatialChat_2 RDS LR rerun

After user review, the real-data validation was rerun directly from
`tests_dev/simulated_data/SpatialChat_2.rds` rather than relying on the small
9-cell contract fixture. The script selected the highest total-mass LR layer in
the RDS, `dL7-dR5`, reconstructed only the current-schema shell needed because
the serialized object is an older schema, and preserved the original expression,
metadata, coordinates, identities, and selected sparse probability layer.

Evidence from `C:/Users/Administrator/AppData/Local/Temp/render-spatialchat2-real-lr.R`:

- source RDS: `tests_dev/simulated_data/SpatialChat_2.rds`
- selected real LR: `dL7-dR5`, rank 1 by probability mass among 35 LR layers
- object scale: 2,000 cells, 1,110 genes, five identity groups
- selected probability layer: 200,322 non-zero entries, mass
  `824.7772188523128`
- `computeCommunField(slot.name = "net", signaling.name = "dL7-dR5", top = 0.8,
  sparse = TRUE)` produced 1,500 non-zero outgoing cell vectors
- widget payload: 2,000 cells, 46 x 46 interpolated grid, 1,323 valid non-zero
  grid points, 900 particles
- generated page: `commun-flow-spatialchat2-real-dL7-dR5.html`
- browser verification: the page rendered a `sc-flow-root`, title
  `Outgoing communication flow of dL7-dR5 (SpatialChat_2 real LR)`, two canvases
  (`sc-flow-cells` under `sc-flow-trails`), Visual speed / Particles / Point
  opacity / Zoom controls, and moving particles with no console error
- screenshot: `C:\Users\Administrator\AppData\Local\Temp\omp-sshots-1575d78fb936ebfa.webp`

This confirms the visualization implementation was tested on a real LR layer
from `SpatialChat_2.rds`, not on the 9-cell synthetic contract fixture.
