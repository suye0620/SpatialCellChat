# Agent Note: Explicit arrow display scale for field grids

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: api
Scope: netVisual_CommunFieldGrid arrow drawing
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

On the full Slide-seq LAMININ field, a majority of drawn nonzero arrows are too short to distinguish on the default canvas. The three-object comparison in `tests_dev/compare-arrow-display-scale.R` shows a median shaft length of about 0.075 mm at scale 1; a global increase improves visibility of some arrows but creates overlap among strong signals. Existing `arrow.length` changes arrowheads rather than the shaft length.

## Decision

Add a named-only positive finite `arrow.scale` to `netVisual_CommunFieldGrid()` after `...`, defaulting to 1. It multiplies only the rendered arrow shaft magnitude after the metR stat has applied its threshold and sampling. For non-default scales, label the figure with the chosen display multiplier so it is not mistaken for a change in the underlying field. Do not automatically choose a data-dependent multiplier or change the default.

## Constraints and invariants

- At scale 1 the existing mapping, layer data, threshold, sampling, and resulting plot remain unchanged.
- Do not modify stored cell fields, aggregation, interpolation, coordinates, vector directions, `min.mag`, arrow width, opacity, or arrowhead settings.
- Multiply all finite drawn arrow lengths by the same constant; preserve zero vectors and relative magnitudes. Missing/nonfinite values remain subject to the existing renderer.
- Make the changed scale explicit on the plot and in help. Different display multipliers must not be mistaken for cross-plot comparisons of absolute magnitude.
- Keep earlier positional arguments unchanged and accept only one finite positive numeric value.

## Alternatives considered

- Automatically normalize/clamp arrow lengths: changes relative visual magnitudes and could mask strong-vs-weak differences.
- Multiply input vectors: changes stat input and can change threshold decisions; rejected.
- Set a larger universal default: can create severe overlaps, so remains the caller's explicit choice.
- Expose only private ggplot layer edits: technically possible but fragile and easy to misinterpret.

## Consumer impact

The named-only `arrow.scale` parameter affects `netVisual_CommunFieldGrid()`, its ggplot layer mapping and optional caption, help, serialization, source-loaded and installed-package usage, and tests/rendered artifacts. `computeCommunField()`, `netVisual_CommunFlow()`, spatial data, grouping, and other visualization entry points remain unchanged.

## Consequences

The user can improve visibility without rewriting data, but choosing a large multiplier can cause overlapping arrows. A caption communicates the chosen display scaling; it is not a calibration legend. Multiplier 1 preserves existing output. Any altered multiplier should be chosen and reviewed in context.

## Evidence

- `docs/visualization-migration-P3-issues-2026-11-10.md` P3-018.
- `tests_dev/compare-arrow-display-scale.R` compares each actual producer-derived field at scales 1 and 4 while checking all other built-layer columns; PDFs and summaries are in `tests_dev/figures/P3/arrow-display-scale-2026-11-11/`.
- `tests_dev/test-arrow-display-scale.R` checks threshold preservation, custom geom arguments, actual arrow grob shaft/head lengths, invalid inputs, zero input, serialization and default behavior; it also runs the grid smoke suite.
- `tests_dev/test-streamline-compat.R`, `tests_dev/test-streamline-installed.R` with a temporary installed package, and `tools::checkRd("man/netVisual_CommunFieldGrid.Rd")` passed. The installed test reloads a scaled plot in a fresh R process.

## Acceptance criteria

- Default output has unchanged mapping, data, and built-layer values; scale 1 matches the existing drawing.
- At non-default scales only shaft magnitude changes, by the specified factor; raw layer inputs, locations, vector direction, normalized magnitude, colors, alpha, linewidth, threshold selection and background remain unchanged, with an explicit figure caption.
- Reject invalid values; preserve zero fields and test `min.mag` selection, explicit custom arguments, serialization and actual rendered output.
- Compare the real producer-derived fields for simulated, Visium and full Slide-seq at scale 1 and a non-default scale; report and retain the images without claiming visual or biological acceptance.
- Run streamline/grid regression and `tools::checkRd`; validate this Note and report unrelated existing validator failures.
