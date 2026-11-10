# Agent Note: All-zero communication field rendering

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: algorithm
Scope: netVisual_CommunFlow default streamline layer
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

An all-zero field through `netVisual_CommunFlow()` results in an all-zero interpolated grid. The installed metR 0.19.0 computes an undefined integration step count and the public plot fails at build time; with an explicit step count it may instead issue an empty-data warning. See `tests_dev/diagnose-zero-streamline-public.R` and P3-020.

## Decision

For the default metR stat only, if every interpolated vector component is finite and exactly zero, render an empty streamline data layer, retaining the background, mapping, and original stored field. Do not invoke the integrator for this case. Explicit custom stats and any nonzero or nonfinite input keep the existing path.

## Constraints and invariants

- Do not change field computation, interpolation, input object, probability, metadata, or the values of any nonzero vector.
- Only an exactly zero, entirely finite and nonempty grid is eligible; do not hide missing or infinite data.
- Keep custom stat data unchanged and let it decide its own zero-field behavior.
- Preserve background positions and labels, ggplot return and serialized plot rendering.
- Do not modify the metR namespace, calculate invented streamlines, or silence unrelated warnings.

## Alternatives considered

- Patch metR's setup and internal empty integration cases: unnecessary dependence on multiple private dependency methods.
- Suppress the error or replace NA values: masks invalid data and could alter nonzero fields.
- Omit the layer entirely: changes the plot structure and complicates existing consumers; use an empty layer instead.

## Consumer impact

`netVisual_CommunFlow()` default returned ggplot, its build and save paths, source and installed-package callers, `tests_dev/test-streamline-compat.R`, `tests_dev/test-grid-visualization.R`, help and P3 issue log. The data of the default streamline layer is empty only for this zero-grid case; stored object fields and background layers remain available.

## Consequences

A default all-zero plot displays the original background with no streamline and no metR zero-integrator error. Callers inspecting the default layer's raw data for a zero field see zero rows; to inspect the original field use the object. A zero grid avoids the metR compatibility probe, while any nonzero grid continues using it.

## Evidence

- `tests_dev/diagnose-zero-streamline-public.R`: synthetic public-entry failure and a prototype empty layer that renders the background.
- `docs/visualization-migration-P3-issues-2026-11-10.md`: P3-020 and P3-021.
- `tests_dev/test-streamline-compat.R`: public zero-field and custom-stat cases, controls, and real behavioral passthrough.
- `tests_dev/test-streamline-installed.R`: installed package and fresh-process reload for both zero and nonzero plots.
- `tests_dev/test-grid-visualization.R`, `tests_dev/test-viz-p3-migration.R`, and `tools::checkRd("man/netVisual_CommunFlow.Rd")` passed. The full Note validator has 20 pre-existing errors and 12 legacy warnings; this Note validates.

## Acceptance criteria

- Public zero-field plot builds and saves with no streamline integration error or warning, retaining all background points and leaving the source object unchanged.
- Mixed/nonzero, NA, infinite, and empty input must not be classified as all-zero; custom stat receives its original data.
- Three-object producer-field comparisons remain unchanged because their grids are nonzero; retain earlier comparison evidence and run focused controls.
- Run streamline, grid and P3 regressions, check updated help, installation/reload as appropriate, and the Note validator; document existing unrelated warnings.
