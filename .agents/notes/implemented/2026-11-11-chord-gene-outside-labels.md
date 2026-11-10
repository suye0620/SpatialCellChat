# Agent Note: Optional outside labels for gene chord diagrams

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: api
Scope: netVisual_chord_gene plotting API and graphics devices
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

Dense gene chord diagrams overlap sector labels, notably 111 sectors in the Visium COLLAGEN example. The candidate outside layout avoids strict leader crossings in 12 tested cases but can overflow labels or collide with the existing legend on a small device (P2-024–027).

## Decision

Named-only `label.layout` defaults to the existing radial layout and provides an opt-in outside layout. The outside layout spreads displayed sector labels into ordered left and right columns, connects each to its sector, and warns when measured label boxes overlap, overflow the plotting canvas, or overlap the measured legend. It uses the graphics device and existing `lab.cex`, legend position and colors. Scientific computations and links remain unchanged.

## Constraints and invariants

- Existing calls, including positional and `...` arguments, preserve the radial path and output.
- No filtering, collapsing, weighting, renaming or reordering of network edges or sectors is introduced.
- Each displayed sector has one label and leader; preserve group colors, label text, and the recordedplot return contract.
- Never silently hide text to make a dense graph appear readable; warn about measurable collisions and document unmeasured ones.
- Clear circlize state even if the new branch fails. Do not change user objects.

## Alternatives considered

Always change the layout: rejects established defaults and cannot guarantee adequate space for all figures. Hide or subset labels: loses identifiers. Existing `circos.labels` on newly created tracks placed labels inside the chords in the recorded experiment.

## Consumer impact

The exported `netVisual_chord_gene` API, its roxygen source and Rd documentation change; existing `tests_dev/render_viz_figures.R`, `tests_dev/render_chord_review.R`, experiment scripts and downstream callers continue to use the radial default. New tests and example captures cover explicit outside use. No stored schema, producer, dependency, or migration changes.

## Consequences

Outside mode is an additional rendering option that requires sufficient device area; detection is based on measured bounding boxes and does not guarantee readability or cover intersections with arbitrary custom graphical additions. Dense graphs remain visually complex. Displayed layout changes but calculations do not.

## Evidence

`docs/visualization-migration-P2-issues-2026-11-10.md` P2-024–027; `tests_dev/diagnose-chord-external-unfold.R`; `tests_dev/diagnose-chord-external-boundaries.R`; `tests_dev/audit-chord-boundaries.py`; three-object experiment outputs under `tests_dev/figures/P2/chord-external-boundaries-2026-11-11/`.

## Acceptance criteria

- `tests_dev/test-chord-outside-layout.R` compares the three fixture chord inputs with `f643c28` and checks simulated default pixel identity with positional parameters and `link.sort` forwarded through `...`.
- The same test records actual outside paths, verifies labels, sector anchors, edge inputs and colors for three fixtures on large and small devices, and checks measured capacity and legend warnings. `tests_dev/audit-chord-production-paths.py` checks strict crossings for the resulting nine outside plots.
- The test covers empty and invalid inputs, long labels, and error cleanup; a temporary-library install and actual installed-package plot plus `tools::checkRd` passed. The previous implementation comparison needs Git history and the fixture objects.
- Reviewable outputs and limits are documented in P2-028; visual and biological review remain separate.
