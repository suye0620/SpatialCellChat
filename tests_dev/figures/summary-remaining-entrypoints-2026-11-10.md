# Remaining visualization migration: verification

## Completed in this pass

- `test-remaining-entrypoints.R`: comparison entry points, merged pairwise
  producer/consumer alignment, missing ranking, embedding row permutation and
  missing dataset layer, invalid endpoint rejection, and merged spatial-coordinate guard.
- `test-viz-p4-migration.R`: 3D feature and Gi trace point counts, highlight
  partition, fixed z coordinate, stack input/trace/color-domain checks, plus
  existing interactive payload and Shiny server checks.
- `test-commun-flow-widget.R`, `test-remaining-migration.R`, and
  `test-merge-spatialcellchat.R`: focused widget, P5 and merged-chain regressions.
- `test-viz-p3-migration.R`: group pattern, topic and embedding producer/consumer
  tests, with the test-specific warnings reported in the command output.

## Three-object figures and limits

- P0 old-cache failures are preserved in `summary-*-P0.csv`; the separate
  `summary-*-P0-spatialC2CDistPlot-fresh.csv` results pass after real inference.
- P3 Slide-seq Grid and Flow remain `VISUAL-REVIEW`, not visual acceptance.
- `P4/{simulated,visium,slideseq}/summary.csv` records generated HTML; the
  15-case Chrome CDP result is in `P4/browser-summary-2026-11-10.csv` and the
  simulated Shiny mouse-selection result is in
  `P4/shiny-browser-summary-2026-11-10.csv`.
- `plotly_spatialStackPlot()` is a completed single-dataset public entry point;
  multi-slice tissue registration remains outside its explicit contract.

See `docs/remaining-migration-issues-2026-11-10.md` and the phase issue logs for
causes, retained scientific semantics, and reproducibility limitations.


## Merged comparison figure follow-up

Ten comparison consumers now have actual PDFs under `merged/producer-comparison/`, generated from the real merge producer using a simulated fixture and an intentionally modified copy. The separate `merged/heatmap-edgecases/` cases check signed cancellation and isolate handling. These are rendering and contract checks, not biological comparisons or manual visual acceptance.

## Streamline compatibility follow-up, 2026-11-11

The default Flow layer now has a guarded correction for the metR 0.19.0 nonzero-diagonal filtering defect. The synthetic and producer-field comparisons, their limits, the zero-field fix, and remaining Slide-seq display issues are recorded in P3-020 and P3-021. The six three-object/resolution comparisons and their original/adapted PDFs are in `P3/streamline-production-comparison-2026-11-11/`. These are not visual or biological acceptance.
