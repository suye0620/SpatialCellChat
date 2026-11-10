# Remaining migration audit and fixes: P5 analysis/expression chain

Date: 2026-11-10

## P5-001: `subsetCommunication(slot.name = "netP", signaling = ...)` returned no rows

Location: `R/analysis.R:subsetCommunication_internal()`.

Cause: In the new schema, the third axis of `netP$group$prob` contains pathway names. The migrated code still converted `signaling` into LR `interaction_name` values and then filtered the pathway-level table by those names. A valid pathway selection was therefore removed completely. The bug only appeared when `signaling` was supplied explicitly.

Impact: `subsetCommunication(..., slot.name = "netP", signaling = ...)`, and consumers using that pathway-level query, could report no communication. The merged snapshot path had the same defect.

Fix: Use `pathway_name` selectors for `netP`; retain LR `interaction_name` selectors for `net`. Do not change pathway probabilities or add a new p-value gate.

Acceptance: Check explicit pathway selection for single objects and merged snapshots, including exclusion of other pathways.

## P5-002: `selectK()` still read removed top-level `slot$prob`

Location: `R/analysis.R:selectK()`.

Cause: The 11-slot schema stores group probabilities under `slot$group$prob`, with one sparse matrix per layer. The old entry point read `slot$prob` and also retained an undefined `data_sr` variable. This is a broken producer/consumer storage contract, not an algorithm change.

Impact: Users could not select the number of communication patterns from a new object. The old code also ignored caller-supplied `k.range` and `nrun`.

Fix: Use `.sc_group_prob_array()` to materialize only the small group-level NMF input, correct the `data` variable, and honor `k.range`/`nrun`. Cell-level sparse communication is not densified.

Acceptance: Run with a small `k.range` and `nrun`; check the returned ggplot and x-axis range.

## P5-003: `netMappingDEG()` still read removed `@var.features`

Location: `R/analysis.R:netMappingDEG()`.

Cause: The DEG producer writes results to `misc$.var.features`, while the consumer still accessed a deleted S4 slot. A second bug treated every data frame as a merged dataset list because `is.list(data.frame)` is true in base R.

Impact: After `identifyOverExpressedGenes()`, `netMappingDEG()` still failed on a new object; downstream ligand/receptor DEG filters were unavailable.

Fix: Read `object@misc$.var.features[[paste0(features.name, ".info")]]`, distinguish a data frame from merged list output, and keep explicit type/missing-data errors. The ligand/receptor and complex-subunit mapping algorithm is unchanged.

Acceptance: Check new-schema marker data for ordinary ligand/receptor rows, complex subunits, and missing markers.

## P5-004: `barPlot()` reused the previous `df.use` for the last panel

Location: `R/visualization.R:barPlot()`.

Cause: The loop updated `df.use` only for non-final iterations. The final feature therefore reused the previous feature; a single-feature call referenced an undefined object. Separately, one selected feature dropped the matrix dimension before `aggregate()`, and the transposed result was not a data frame.

Impact: `plotGeneExpression(type = "bar")` could fail for one feature and show the wrong gene in the last panel of a multi-feature plot.

Fix: Preserve the feature dimension through aggregation, convert the transposed result to a data frame, and build `df.use` for the current feature at the start of every iteration. Existing summary methods, grouping, and patchwork layout are unchanged.

Acceptance: Render one-feature and two-feature bar plots and check the final panel label/data.

## Audited conclusions intentionally unchanged

- Centrality remains in `slot$centr` / `slot$centr.cell`. This is the current 11-slot compatibility contract; `group$centr` is not the target of this pass. `spatialVisual_scoring()` and signaling-role consumers follow this contract.
- `netVisual_aggregate()` and `netVisual_individual()` retain their existing pathway semantics: select LR group layers, apply the call-time `pval > thresh` gate, then aggregate or draw per LR. They must not mechanically read `netP`.
- `netVisual_chord_cell()` retains the legacy real-time LR aggregation semantics in its `netP` branch; changing it to the production-time `netP` layer requires a separate decision.
- `prob.cell` in `computeAvgCommunProb_Visium()` is a local layer-list variable or a new `net$cell$prob` field access, not an old top-level object field.
- `@meta.data` in `identifyTopicCommunities()` and `barPlot()` applies to temporary Seurat objects; SpatialCellChat itself uses `@meta`.

## P5-005: signaling-role heatmaps ignored `cluster.cols`

Location: `R/analysis.R:netAnalysis_signalingRole_network()` and `netAnalysis_signalingRole_heatmap()`.

Cause: Both constructors passed `cluster.rows` to the column clustering argument, so callers could not independently control row and column ordering.

Impact: Requested column clustering settings were silently ignored, changing presentation but not the underlying scores.

Fix: Pass `cluster.cols` to `cluster_columns`; no data or scoring semantics change.

Acceptance: The P5 render test exercises both signaling-role consumers after the correction.

## P5-006: merged `netMappingDEG()` needed dataset-aligned marker evidence

Location: `R/analysis.R::netMappingDEG()` and `R/merging.R::mergeSpatialCellChat()`.

Cause: merged construction intentionally leaves top-level `misc$.var.features` empty, while preserving each source object's marker tables inside `misc$.datasets[[dataset]]`. A single-table lookup therefore either failed on a valid merged object or risked applying one dataset's DEG evidence to every dataset.

Impact: DEG-based communication filtering could not be used safely after merged construction, even though the source snapshots retained the required evidence.

Fix: map each named communication table with its matching snapshot marker table; when a merged marker table has a `datasets` column, filter that table by dataset first. Unstratified merged marker tables are not silently reused. Missing evidence reports the dataset name and the required producer.

Acceptance: The P5 regression covers snapshot-local mapping with distinct values, dataset-tagged merged tables, missing-evidence failure, and the existing single-object path.

## P5-007: remaining entry-point and widget audit checkpoint

Checkpoint after the P5 migration pass: the broad old-slot scan no longer finds
SpatialCellChat consumers of removed `@data`, `@data.signaling`, `@data.project`,
`@net$prob`, or `@net$pval` fields. The remaining work is focused on omitted
comparison/analysis entrances and optional interactive/3D rendering paths.

The first focused rerun found one test-environment defect rather than a payload
defect: `tests_dev/test-commun-flow-widget.R` can construct the widget, but its
source-mode asset assertion is not valid when `renv` has an installed
`SpatialCellChat` package, because `system.file()` resolves to the installed
package and `htmlwidgets::saveWidget()` copies assets from that package. The
contract needs to accept both source-tree and installed-package execution while
still requiring the JS/CSS files to be present in the rendered dependency
directory.

The next audit must cover `rankSimilarity()`, `rankNetPairwise()`,
`compareInteractions()`, `suggestK()`/`suggestK1()`, the comparison scatter
functions, `spatialVisual_scoring()` merged behavior, topic/statistics plots,
and all available 3D/widget branches. Any behavioral change must preserve the
11-slot schema and be recorded here with a focused regression check.
