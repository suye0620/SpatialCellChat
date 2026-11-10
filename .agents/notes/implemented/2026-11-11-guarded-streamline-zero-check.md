# Agent Note: Guarded correction of streamline zero-vector filtering

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: algorithm
Scope: `netVisual_CommunFlow()` default static streamline layer and its metR compatibility adapter; no widget or stored-field changes
Owner: project maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

metR 0.19.0 filters integration points with `dx + dy != 0`, removing nonzero vectors whose components cancel. `tests_dev/diagnose-streamline-cancellation.R` reproduces the resulting missing line through the public SpatialCellChat entry point. A local, non-global experiment in `tests_dev/diagnose-streamline-local-patch.R` corrects only the two forward/backward predicates; extended controls and the full Slide-seq comparison show the scope and remaining limitations.

## Decision

Use a private, locally scoped, compatibility-guarded adapter for the default `metR::StatStreamline` layer. Probe the installed metR behavior once per session. If it handles a constant nonzero cancelling field, use it unchanged. If it drops that field, support the pinned metR 0.19.0 implementation only after checking the expected internal structure and exactly two predicate occurrences; clone its stat and integrator, replace only these predicates with a true component-wise zero test, and verify the resulting local stat against the probe and an unaffected control. Do not assign into any dependency namespace. Fail with an actionable error for an affected implementation that cannot be checked safely. Explicit caller-supplied custom stats remain untouched.

## Constraints and invariants

- Do not modify stored communication probabilities, computed fields, coordinates, grid aggregation, interpolation, masking, weighting, thresholds, seeds, or the original metR namespace.
- Preserve metR's numerical integration, plotting parameters, return class, and numerical output for unaffected fields. Only retention of nonzero cancelling vectors may differ.
- Cache verified adapter selection per R session; preserve the caller's random seed during probes, including when no seed existed.
- Preserve serialization of returned plots, explicit custom stats, and optional dependency failure behavior.
- Future metR versions that pass the behavioral probe run unmodified; affected unknown versions must fail clearly rather than silently drop signals.

## Alternatives considered

- Wait for an upstream release: least maintenance, but leaves a demonstrated silent false negative in current supported environments.
- Replace or vendor the whole integrator: removes private-interface dependence but changes a much larger algorithmic surface and creates licensing, equivalence, and maintenance costs.
- Perturb vector values or rotate fields: changes input and can move the cancellation to another direction, violating the field invariant.
- Modify metR's namespace or installed package: affects unrelated consumers globally and is not acceptable.

## Consumer impact

The adapter is called only for the default streamline stat constructed by `R/visualization.R:netVisual_CommunFlow()`. Consumers are ggplot rendering, `ggsave`, saved/reloaded ggplot objects, source-loaded development scripts, package namespace loading, user-supplied stat overrides, and static field tests. Existing `computeCommunField`, `netVisual_CommunFieldGrid`, HTML widget, class schemas, network producers, thresholds, and stored results are unchanged. `memoise` is made an explicit dependency because the local cloned integrator uses it; lockfile metadata must be checked. Help and P3 issue/review records must explain the guarded compatibility behavior and any outstanding visual limitations.

## Consequences

The adapter intentionally depends on vetted metR 0.19.0 internals and must reject changed affected internals rather than guessing. First use incurs a small deterministic probe; a future corrected metR version requires no adapter. Maintaining three separate validation paths (fixed upstream, supported affected, unsupported affected) is necessary. Default display density and integration resolution are not changed.

## Evidence

- `docs/visualization-migration-P3-issues-2026-11-10.md` P3-017 through P3-019.
- `tests_dev/test-streamline-compat.R`: production adapter branch, zero/axis/cancellation/near/curved/heterogeneous/boundary/periodic/jitter, random state, custom stat, and fresh-process serialization checks.
- `tests_dev/test-streamline-installed.R`: temporary package installation and fresh-process plot reload and rendering.
- `tests_dev/compare-streamline-production-fields.R` and `tests_dev/figures/P3/streamline-production-comparison-2026-11-11/`: producer-derived fields for three objects at res 0.5 and 8; one outgoing pathway per object, layer data unchanged in these cases.
- `tests_dev/test-grid-visualization.R`, `tests_dev/test-viz-p3-migration.R`, and `tools::checkRd("man/netVisual_CommunFlow.Rd")` passed. Full-zero warnings and Slide-seq visual limits remain recorded in P3-020.

## Acceptance criteria

- Exact cancelling nonzero synthetic vectors render; near, axis, curved, heterogeneous, boundary, periodic, and jitter controls match original metR exactly where unaffected; zero-field warnings/errors remain honestly reported.
- Full simulated, Visium, and Slide-seq flow inputs use actual producer-derived fields; compare old and adapted layer data at default and non-default resolution, report any differences, and render PDFs.
- Test corrected upstream passthrough, unsupported/mismatched affected versions, custom stats, preserved random seed for the probe, no metR namespace mutation, dependency declaration, and plot serialization/reload.
- Run existing grid visualization and P3 migration regression tests, `tools::checkRd` for changed help, and note validator; record pre-existing validator failures separately.
- Document the unresolved default visual quality, upstream zero-field issue, and limited coverage; do not mark visual or biological acceptance as passed.
