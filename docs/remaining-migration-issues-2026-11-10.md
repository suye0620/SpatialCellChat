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

The focused audit is now complete for the previously listed entrances:
`rankSimilarity()`, `rankNetPairwise()`, `compareInteractions()`, the
comparison scatter functions, `spatialVisual_scoring()` merged behavior,
topic/statistics plots, and the available 3D/widget branches. `suggestK()` is
an optional Liger-only path and does not consume SpatialCellChat slots;
`suggestK1()` is a matrix-only path and now has a real rank-one/rank-two
regression. The remaining boundaries are optional dependency availability,
manual visual review, and the explicitly preserved statistical semantics below.

## P5-018: rank-one NMF topic/rank selection failed under current Matrix

Location: `R/analysis.R::identifyCellTopics()` and `suggestK1()`.

Cause: `Matrix::diag(sqrt(d))` interprets a length-one numeric argument as the
order of a square matrix. For a rank-one NMF result, the resulting matrix can
have an unrelated dimension, so the factor rescaling fails with a
non-conformable matrix multiplication error. This was a dependency/API
compatibility defect in an omitted analysis entrance, not a change to the NMF
algorithm.

Fix: Use the shared `.sc_nmf_sqrt_diagonal()` helper backed by
`Matrix::Diagonal(x = sqrt(d))`, with finite/non-negative validation. The
existing factor orientation, reconstruction, topic names, and cross-validation
formula are unchanged.

Validation: `tests_dev/test-remaining-migration.R` runs `suggestK1()` on an 8 x
8 sparse matrix for ranks 1 and 2 with two replicates and verifies the ggplot
and four statistic rows. The previous dimension error is no longer
reproducible. The durable decision and consumer map are recorded in
`.agents/notes/implemented/2026-11-11-nmf-rank-one-diagonal.md`.

## P5-008: merged pairwise rankings were not consumable

Location: `R/modeling.R::identifyEnrichedInteractions()` and
`R/analysis.R::rankNetPairwise()`.

Cause and impact: the producer writes merged rankings into each dataset snapshot,
but the consumer read only `object@net$pairwiseRank`. Consequently the ordinary
producer-then-consumer sequence failed for merged objects.

Fix: consume each matching snapshot and return a dataset-named list, consistent
with merged `subsetCommunication()`. A single object still returns a data frame;
missing rankings identify the dataset and request `rankNetPairwise()`. The legacy
filtering algorithm is unchanged. The proposed decision record is
`.agents/notes/proposed/2026-11-10-merged-pairwise-enrichment.md`.

Validation: `tests_dev/test-remaining-entrypoints.R` checks the single-object
path, different dataset-local results and missing-ranking errors.

## P5-009: empty embedding layer selection did not trigger the intended error

Location: `R/analysis.R::rankSimilarity()` row index construction.

Cause: `paste0(character(0), "--", dataset)` can yield a nonempty string in R.
Matching that string produces `NA`, so a missing dataset appeared to have one
index. This did not supply a valid comparison, but the subsequent error misreported
the problem as no shared layers.

Fix and validation: keep empty index vectors empty; focused tests check absent
dataset layers and row-name permutation. No distance calculation is changed.

## P5-010: historical threshold and directional-change semantics remain open

`identifyEnrichedInteractions()` uses `pval < thresh` for initial candidates,
but subsequently applies a fixed 0.05 cutoff, including when `pairLR.use0` is
provided. Thus `thresh > 0.05` does not simply expand the final result. This
pre-existing behavior is documented, not silently changed in the migration;
resolving its intended statistical meaning requires a separate decision.

`netAnalysis_signalingChanges_scatter()` filters using `rowSums(mat.diff) != 0`.
An outgoing change of +2 and incoming change of -2 cancels and disappears even
though both changed. This is A-13 in `docs/ANALYSIS_SEMANTICS.md`. The current
pass preserves legacy behavior as requested; the focused regression verifies the
current behavior, not its scientific correctness. A later decision should test
each direction independently and address zero and missing-layer interpretation.

`netAnalysis_contribution()` is explicitly single-dataset. It now rejects a
merged container with an actionable error rather than reporting a misleading
missing-probability failure; comparison consumers remain `rankNet()` and the
dedicated merged comparison plots.

## P5-011: validation and dependency boundaries

`suggestK()` consumes a Liger object and `suggestK1()` a matrix; neither directly
reads SpatialCellChat slots, so they are not slot-migration consumers. Their
NMF/statistical interpretation and optional dependencies remain separate issues.

P4 HTML generation now has a 15-case Chrome CDP smoke result (Plotly
initialization, point/cell coverage, 3D scene, stack trace/colorbar/layer ticks,
and relayout). Shiny Chrome sessions on all three datasets cover real mouse
selection, replacement and reset. Slide-seq updates took about 11.6 and 9.2
seconds in one local run, so performance still needs investigation. P3 widget
browser checks cover animation and selected controls on all three datasets,
but do not establish biological correctness. P3 Slide-seq Grid and Flow remain
`VISUAL-REVIEW`, including the integration warning on the default Flow plot.
The 3D stack production entry point and its single-dataset input contract are
implemented and verified; multi-slice tissue registration remains outside that
contract.

The Agent Notes validator reports pre-existing errors in other records; the new
merged-pairwise and 3D-stack Notes validate. Do not claim that the full Note
set is clean until those independent records are repaired and revalidated.

## P5-012: invalid group names could silently change cell-pair identity

Location: `identifyEnrichedInteractions()` previously removed invalid `from`
and `to` names independently before checking the remaining lengths. For example,
`from = c("A", "missing")`, `to = c("missing", "B")` was silently converted
to A -> B even though A -> B was not one of the requested pairs. This affects
single and merged results and can yield plausible but scientifically wrong data.

Fix: reject invalid names, missing values, empty selections and non-integer or
out-of-range indices before pairing; retain valid pair order and bidirectional
logic. This deliberately replaces the unsafe historical partial-name behavior
with an error. Focused tests cover mismatched invalid names and invalid indices.

When no pair passes the requested candidate threshold, the same consumer now
returns a typed empty data frame instead of failing while binding an empty list.

## P5-013: merged LR extraction used union metadata for dataset-local caches

Location: `extractEnrichedLR()` and `extractEnrichedLR_internal()` in
`R/analysis.R`.

Cause: The merged wrapper passed the merged top-level LR union to every dataset
snapshot. A pathway absent from one dataset could therefore be searched against
the wrong candidate table; additionally, `searchPair()` can simplify a
single-interaction result through `future_sapply()` into an unstable shape.

Fix: Use each snapshot's own `LR$LRsig`, skip snapshots that do not contain a
requested pathway, reject a pathway absent from the merged union, and perform
the helper's exact pathway selection by direct data-frame filtering. The
existing `LR.sig` cache intersection and call-time `thresh` behavior remain
unchanged.

Validation: `tests_dev/test-extractEnrichedLR.R` checks a one-interaction
dataset-local candidate table, and `tests_dev/test-merge-spatialcellchat.R`
checks merged extraction through the real producer when one dataset lacks a
requested pathway.


## P5-014: real fixture pair does not satisfy the merged producer contract

Location: `.sc_merge_assert_inputs()` and `mergeSpatialCellChat()` in `R/merging.R`. `tests_dev/audit-merge-real-fixtures.R` reads the existing independent Visium PP1 and Slide-seq MBM05 fixtures without modifying them. The 905-spot and 29,536-cell objects have 23,718 shared gene names but their full gene axes and order differ; their group-level names have no intersection. The stored DB objects are identical. The real producer rejects the pair at the gene-axis check, before it reaches the group-axis check; summary: `tests_dev/figures/merged/real-fixture-contract-2026-11-11.csv`.

Cause and significance: the current contract intentionally requires identical assay gene order and identical group levels, and the downstream comparison paths interpret groups as corresponding named biological populations. A group missing from one dataset is not necessarily a measured zero. Automatically padding groups or relabeling these unrelated annotations could produce plausible but misleading comparisons. The current test evidence therefore does not establish an independent-sample biological comparison.

Next direction: use independently acquired compatible samples with a defensible common annotation and aligned gene/DB inputs; separately decide how comparisons should represent groups or genes missing from one dataset before relaxing the existing contract. Merely having an identical stored DB does not establish sample comparability.


## P5-015：真实 Visium 对象的零差异 merged 技术对照

位置：`mergeSpatialCellChat()` 和比较绘图入口。此前 10 张比较图主要来自模拟对象和人为修改的副本，无法证明真实 fixture 经生产者进入同一链条后可以绘图。本次 `tests_dev/render-merged-visium-self-control.R` 将真实 Visium PP1（905 个 spot、6 个群、160 条 pathway）与其自身分别作为 `control_a`、`control_b` 传入真实生产者；逐项核对群 count、weight 和群中心性一致后，生成 10 个比较入口的 PDF，包括成对相似性、嵌入及其两种绘图。结果及警告记录位于 `tests_dev/figures/merged/visium-self-control-2026-11-11/`。

性质与限制：这两个输入是同一份样本，结果是可验证的技术零差异，不能当作两个独立样本，也不能用于声称组间生物学等价。现有 Visium 和 Slide-seq 的原始基因轴及群注释不符合合并契约（P5-014），独立样本比较仍需合适数据。所有图像仍待人工视觉和生物学验收；依赖弃用等警告保存在 `warnings.csv`，PNG 转换提示缺少 ArialUnicode 字体。

## P5-016：零差异热图错误地报告未推断到通讯

位置：`R/visualization.R:netVisual_heatmap()`。原因：合并对象的矩阵是第二个数据集减第一个数据集；两者完全相同时，差异矩阵为零，但旧的全零检查把它解释成“未推断到通讯”并报错。本质：混淆了“有效比较的差异为零”和“单对象没有推断到通讯”。影响所有合并比较中显示范围恰好全零的 count/weight 热图。

修复：对非空的全零 merged 比较保留原矩阵和群名称，以中性色绘制并明确标注“显示范围内无差异”，不生成虚构色标或零值柱图的自动微小刻度；单对象全零以及移除所有孤立群后的空选择仍沿用错误行为。`tests_dev/test-merged-heatmap-edgecases.R` 检查正负抵消、纯负变化、零矩阵和孤立群；`tests_dev/test-merged-zero-comparison.R` 检查 count/weight、轴名及单对象边界。

## P5-017：零行信号变化和空连边网络

位置：`R/analysis.R:netAnalysis_signalingChanges_scatter()` 和 `R/visualization.R:netVisual_diffInteraction()`。原因：既有 `rowSums(mat.diff) != 0` 筛选可能留下零行，随后给零行数据框赋单个值而报错；全零差异网络虽然能绘图，却对空边集求最大值并设置边属性，产生无意义警告。影响相同数据或筛选后没有保留通路/连边的 merged 比较。

修复：零行散点图返回无点 ggplot，并说明“没有通路通过当前变化筛选”；不更改 `rowSums` 筛选或将此结果解释为所有方向均无变化。空边网络跳过边属性处理，保留原节点。聚焦测试包含真实零差异、非零相反方向抵消、`do.label` 开关，以及原有非零比较回归。该图提示不能替代后续对 P5-010 方向抵消语义的单独决策。

## P5-019: 全零信号层在模式归一化中产生 NaN

位置：`selectK()` 和 `identifyCommunicationPatterns()` 的 outgoing/incoming 群网络预处理。原逻辑将每个信号层按该列最大值归一化；全零层的最大值为 0，除法生成 NaN。`rowSums(data) != 0` 随后可能保留含 NA 的行，NMF 最终报 `Input matrix x contains at least one null or NA-filled row`。这在 Visium fixture 的 1,032 个 LR 层中可复现，影响所选方向存在全零层的 rank 选择和模式识别。

修复：新增共享 `.sc_normalize_pattern_scores()`，先检查聚合分数为有限值，只对最大值大于 0 的信号列归一化；全零列从 NMF 输入移除，空网络仍由现有调用方给出明确错误。非零列的按列最大值归一化定义不变。`tests_dev/test-remaining-migration.R` 额外构造全零层，分别验证 `selectK()` 与 `identifyCommunicationPatterns()`；`tests_dev/render_p5_remaining_figures.R` 在模拟、Visium 和完整 Slide-seq fixture 上均生成 selectK 图。

边界：该修复只排除没有群级信号的信息列，不改变任一非零信号的权重或阈值。NMF 与 ggplot 依赖警告仍会输出；绘图写出及 CDP smoke 不等同于人工视觉或生物学验收。

## P5-020：静态 plotFeatures 模式此前只有合成数据验收

位置：`tests_dev/render_viz_p4.R` 及 P4 fixture 验收。此前三个真实对象的 P4 覆盖了 `plotFeatures()` 的 2D/3D Plotly 页面，但静态 `method="normal"` 和 `show.group=TRUE` 双面板只在合成数据路径上验证，不能证明 Visium 与完整 Slide-seq 的坐标/细胞标识输入也能正常绘制。

处理：新增 `SC_P4_NORMAL_ONLY=true` 测试模式，避免为静态验收重生成已存在的 Plotly 页面；对普通图及分组图的两个面板逐一断言每个 fixture 细胞恰好出现一次，并验证分组面板保留全部注释类别。脚本退出码现在同时覆盖 HTML 与静态 PDF 失败。审计器要求三个 fixture 各自存在普通/分组两行 `PDF-WRITTEN` 汇总并检查对应 PDF。

验证：R 4.4.3 下 simulated（20 细胞）、Visium（905 spot）和完整 Slide-seq（29,536 细胞）均成功写出两种静态图，共 6 份 PDF；汇总位于 `tests_dev/figures/continued-2026-11-12/P4-normal/`。这不替代人工视觉或生物学验收。
