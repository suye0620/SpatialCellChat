# Agent Note: Render valid zero-difference comparison results

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: api
Scope: Merged comparison heatmap, signaling changes scatter and network
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

Merging two copies of the real Visium fixture yields identical group count/weight and group centralities. The comparison heatmap errors with a misleading claim that communication was not inferred; the signaling-changes scatter errors assigning a scalar to a zero-row data frame. The differential network renders but warns while processing an empty edge set. The prior heatmap regression intentionally expected errors for all-zero matrices, so the behavior change must be explicit.

## Decision

Merged comparisons with a nonempty all-zero displayed matrix render the actual zero matrix with an explicit no-difference-in-selection title and a neutral color, without a misleading numerical color scale or zero-bar axis. When no signaling rows pass the existing directional-change filter, the function returns an empty ggplot explicitly stating that no pathways passed that filter. The differential network skips edge attributes for an empty edge set. Single-object zero heatmaps and zero groups after remove-isolate retain their prior error behavior.

## Constraints and invariants

- Preserve subtraction order, all existing selection, p-value, threshold, centrality and `rowSums` change-filter semantics; a no-rows scatter does not assert that the underlying directional changes are biologically zero.
- Never fabricate a nonzero differential value, pathway, edge or comparison sample.
- Preserve the return classes, group axes and options for nonempty plots; keep existing no-data error when no groups remain after removal.
- Self-comparison of a single fixture is a technical control, not evidence from two independent biological samples.

## Alternatives considered

Retain errors: prevents displaying an ordinary valid null comparison and misstates why it is empty. Use synthetic nonzero placeholders: misrepresents the result. Change the row-sum filter to test directions independently: alters legacy scientific semantics and requires a separate decision.

## Consumer impact

The exported `netVisual_heatmap`, `netVisual_diffInteraction` and `netAnalysis_signalingChanges_scatter` comparison branches, their documentation, `tests_dev/test-merged-heatmap-edgecases.R`, the merged producer and comparison consumers, and their tests and figure captures are affected. Single-object heatmaps and other drawing functions are unchanged. No schema or dependency changes.

## Consequences

Valid zero comparisons produce explicit neutral results; users who depended on the previous error for nonempty zero merged matrices will receive a plot. All-zero output is a selected-result statement, not statistical evidence of biological equivalence. A truly empty selection remains an error.

## Evidence

`tests_dev/render-merged-visium-self-control.R` and its `tests_dev/figures/merged/visium-self-control-2026-11-11/` output reproduce the zero-difference failures. `tests_dev/test-merged-heatmap-edgecases.R` documents the previous zero error and nonzero signed-cancellation cases. `docs/remaining-migration-issues-2026-11-10.md` P5-010 records the separately deferred row-sum semantics.

## Acceptance criteria

- `tests_dev/render-merged-visium-self-control.R` verified equal matrices and centralities, merged the real Visium fixture with itself, and rendered 10 comparison consumers including pairwise embedding. This is a technical control, not independent replication.
- `tests_dev/test-merged-zero-comparison.R` checks zero merged count/weight matrices and group names, single-object zero and empty-after-isolate errors, empty and cancelling-direction scatter cases including both label settings, and edgeless network warnings.
- `tests_dev/test-merged-heatmap-edgecases.R` covers nonzero signed cancellation, negative values, isolate removal and zero heatmaps. `tests_dev/test-remaining-entrypoints.R` covers the original row-sum filter for nonzero values and nonzero comparison consumers.
- The package installed in a temporary library, and installed-package zero plots and both Rd checks passed. Real-data figures and the visual and independent-sample limits are in P5-015–017.
