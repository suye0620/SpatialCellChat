# SpatialCellChat 可视化层迁移总表（P0 产出）

> 2026-11-04 · 由三个只读侦察 agent 的结论合并（visualization.R 全量 / widget+spatial.R / net 参数流与绘图包用法）。
> 输入：docs/VISUALIZATION_AUDIT.md、docs/VISUALIZATION_SEMANTICS.md（V-001…V-018）。
> 判定口径：**OK** = 11-slot 全新槽可跑；**partial** = 主路径新、分支旧；**broken** = 读已删字段；**helper** = 不直接碰对象。

> **2026-11-10 状态校正：** §1-§3 的 broken/partial 数量和部分行号是 2026-11-04 的迁移前快照，不是当前验收状态。当前已通过 P0-P5 的核心回归：group LR 网络、空间 aggregate/individual、analysis network、embedding、CommunField、topic、plotly/Shiny、merged producer，以及 P5 的 netP 筛选、selectK、netMappingDEG、barPlot 和 signaling-role 图。详见 `tests_dev/test-*.R`、`tests_dev/figures/summary-*.txt` 和 `docs/remaining-migration-issues-2026-11-10.md`。

> **2026-11-11 补充验收：** 修复全零信号层归一化产生 NaN 的问题（P5-019）；P5 收尾的 8 类图已在模拟、Visium 和完整 Slide-seq fixture 上逐项生成。P4 另补 `plotFeatures()` 的三对象 2D/3D 页面，21 个页面通过 Chrome/CDP smoke。所有结果仍待人工视觉和生物学验收，不能把出图或浏览器 smoke 当作科学正确性证明。

## 1. 规模

| 文件 | 顶层函数 | broken | partial | OK/keep | helper |
|---|---|---|---|---|---|
| `visualization.R` | 61 | **18** | **3** | 11 | 29 |
| `visualization_widget.R` | 9 | 0 | 0 | 9 | — |
| `spatial.R` | 18 | 0 | 0 | 6 入口 | 12（共享基础设施 5 个，见 §5） |

`analysis.R` 另有 **10 个** broken（上游，见 §3），其中若干被可视化直接依赖。

## 2. broken / partial 清单（visualization.R，21 个）

| 函数 | 行 | 状态 | 旧依赖 | 关联 | 批次 |
|---|---|---|---|---|---|
| `spatialFeaturePlot` | 350 | migrated | 读取 `assay` 与 `net$group/cell`，显式 LR 稀疏逐层筛选，坐标按 cell name 对齐；保留旧 `slot.name` 值映射 | V-001/003/004 | P1 |
| `netVisual` | 964 | migrated | 逐 LR 使用 `net$group$prob/pval` 的稀疏层；空间分支显式传入 coordinates/idents | V-001/006 | P2 |
| `netVisual_aggregate` | 1353 | migrated | 逐 LR 使用 `net$group$prob/pval`，保留本次 `pval > thresh` 聚合语义 | V-001/002 | P2 |
| `netVisual_individual` | 1517 | migrated | 同上；空间分支使用对象 `@idents` 与坐标 | V-001/006 | P2 |
| `netVisual_diffInteraction` | 2485 | migrated | `@net[[comparison]][[measure]]` 旧 merged 布局 | V-001/007 | P2 |
| `netVisual_heatmap` | 2647 | migrated | `@options$mode`、`net[[i]][[measure]]`、remove.isolate 重复删两次、cluster 参数错位 | V-001/007/008 | P2 |
| `netVisual_barplot` | 2844 | migrated | 旧 merged 检测与旧容器 | V-001 | P2 |
| `netVisual_bubble` | 2992 | migrated | 主路径走 `subsetCommunication`（已适配）；`prob.original` 已改为显式列引用 | V-001/009 | P2 |
| `netVisual_chord_cell` | 3360 | migrated | 复用稀疏 group LR 层；保留 netP 分支实时 LR 聚合语义 | V-001/006 | P2 |
| `netVisual_chord_gene` | 3651 | migrated | 逐层读取稀疏 group LR 网络，按 `pval > thresh` 门控并关联 LR 注释 | V-001 | P2 |
| `netVisual_embedding` | 4250 | migrated | similarity 新、点大小读旧 `netP$prob` | V-001/017 | P3 |
| `netVisual_embeddingZoomIn` | 4344 | migrated | 同上 | V-001/017 | P3 |
| `netVisual_embeddingPairwise` | 4431 | migrated | 旧 merged 嵌套 | V-001/017 | P3 |
| `netVisual_embeddingPairwiseZoomIn` | 4546 | migrated | 同上 | V-001/017 | P3 |
| `netVisual_CommunField` | 5390 | migrated | 旧 `slot$field` 路径（V-013） | V-013 | P3 |
| `netVisual_CommunFieldZoomIn` | 6373 | migrated | 随基础图失效，修 CommunField 后自动恢复 | V-013 | P3 |
| `plotGeneExpression` | 6545 | migrated | 读取 `assay$signaling`；bar 分支支持单/多 feature | V-003/005 | P3 |
| `plotly_spatialLRpairPlot` | 7202 | migrated | `@data` + `apply(prob>0,3,sum)` 稠密膨胀 | V-001/004 | P4 |
| `plotly_spatialLRpairPlot_shiny` | 7831 | migrated | 启动即崩 + observer 累积 | V-001 | P4 |
| `netVisual_TopicComposition` | 4679 | migrated | SCC meta 读 `@meta.data`（错）、remove.isolate 错轴 | V-014 | P3 |
| `netVisual_TopicSignaling` | 4985 | migrated | SCC 主路径可跑、Seurat 分支必崩 | V-015 | P3 |

## 3. 上游报告（约束 5：本轮不擅改，需你裁定）

`analysis.R` 有 **10 个 broken**，其中 4 个被可视化直接依赖：

| 上游函数 | 破损点 | 谁依赖它 |
|---|---|---|
| `identifyCommunicationPatterns` | `slot$prob` + `apply(prob,c(1,3),sum)` | `netAnalysis_river`/`dot` 的数据源（P3） |
| `extractEnrichedLR_internal` | 参数 `net` 直接接收旧槽位 | `spatialFeaturePlot` signaling 模式、`plotGeneExpression`（P1/P3） |
| `netAnalysis_contribution` | `net$prob/prob.cell` | P2/P3 若使用 contribution 图 |
| `rankNet` | single/comparison 两分支都读旧槽 | bubble 的 rankNet 模式 |

**历史基线（已由后续迁移修正）**：当时 `mergeCellChat` 尚未定义，comparison consumers 因此没有可验证的 producer。当前实现由 `mergeSpatialCellChat()` 作为正式 producer，并保留仅供旧内部调用点使用的 `mergeCellChat()` adapter；实际状态和剩余边界见 [P2 migration issue log](visualization-migration-P2-issues-2026-11-10.md#p2-015-merged-comparison-producer-与空间坐标语义)。

历史上另有两个确定性变量名 bug：`selectK` 的 `data_sr` 和 `netVisual_bubble` 的 `prob.original`。两处现已分别改为 group-level 稀疏层展开后的 `data` 与 `df.net$prob.original`；后续回归仍需保留这两个断言。

## 4. OK / keep 清单（回归基线，不重构主体）

`spatialDimPlot`、`spatialDimPlotPoints`、`spatialC2CDistPlot`、`spatialVisual_scoring`、`spatialTopicPlot`、`netAnalysis_river`、`netAnalysis_dot`、`spatialLeePlot`、`netVisual_hierarchy1/2`（renderer）、`netVisual_circle`、`netVisual_spatial`、`netVisual_chord_cell_internal`、`netVisual_CommunFieldGrid`、`netVisual_CommunFlow`、`plotly_spatialDimPlot`、`netVisual_CommunFlowWidget`、`vizDimLoadings`、`dotPlot`、`StackedVlnPlot`、`modify_vlnplot`、`barPlot`、`visualizePie`、`plotStatistics_Lee/Gi`（helper）、`doHeatmap`（语法修复）、其余纯 helper 29 个。

keep ≠ 不动：本批仍要做**包用法核查**（`element_line(size=)`、`aes_string`、`guides(colour=FALSE)`、`coord_fixed+aspect.ratio` 并存等，见 §6）与 `options(warn=-1)` 清理（8 处）。

## 5. 共享基础设施（改动需谨慎）

1. `computeCellDistance` — 推理核心 + Lee 图 + C2C 图三方消费；返回契约三 key；
2. 缓存契约校验有 **3 处**：`.spatial_distance_cache_matches`（绘图）+ `modeling.R:272` 内联 all.equal（建模）——改缓存结构必须同步；
3. `.sc_resolve_grid_size` / `.sc_grid_membership` — 跨 3 文件 4 个消费方（含 widget）；
4. `makeGridSpatialCellChat` — 被 `identifyOverExpressedGenes(do.grid=TRUE)` 用作分析前处理；
5. `computeColocalization` — 与可视化无关但被 `computeAvgCommunProb*` 调用；硬依赖未安装的 `raster::pointDistance`，建议换 base/terra。

## 6. 绘图包用法偏离（合并 agent 3 的 14 类，映射批次）

| 类别 | 位置 | 批次 |
|---|---|---|
| `aes_string` ×3 | barplot_internal、doHeatmap | P2 |
| `ggplot_build` panel internals ×2 | dotPlot、doHeatmap | P2 |
| `size=` → `linewidth=` ×12 行 | embedding×5、barPlot、signalingRole_scatter×3、selectK、river、modify_vlnplot、netVisual_spatial、visualizePie | P2/P3 |
| `coord_fixed`+`aspect.ratio` 并存 ×8 行 | diffInteraction、plotFeatures、Lee/Gi helper | P2/P4 |
| 重复内联 theme ×7 块 | spatialFeaturePlot×4、CommunFieldGrid/Flow/plotFeatures 的 scChat_theme_opts 局部定义 | P1/P2/P3/P4 统一收敛 |
| ComplexHeatmap 手写 viewport 拼版 | identifyCommunicationPatterns（上游） | 上游报告 |
| circlize 无 `circos.par`、自环未去、两版 chord 行为不一致 | chord_cell_internal、chord_gene | P2 |
| plotly trace 硬编码 ×3、colorbar×3、无 WebGL×2、hover 未转义×2、observer×1 | plotFeatures、Gi 3d、shiny | P4 |

## 7. 共享 helper 抽取清单（迁移时一次到位）

1. `.sc_select_significant_layers(net, pairLR, thresh)` — 5 份"门控+筛层"逻辑（迁移到 `net$group$prob` + marginSums）；
2. `.sc_prep_plot_coords()` — 10 处坐标交换 + `scale_y_reverse()`（含 raster 分支）；
3. `.sc_make_colormap(palette, n, direction, zero.color)` — 7 份色图构造；
4. `themeFeaturePlot()` 修复后统一接入 — 7 块内联 theme；
5. `.sc_draw_chord()` — 两份 chordDiagram 全套；
6. `.sc_plotly_colorbar()` / `.sc_plotly_hover()` — trace 后处理与转义；
7. `.sc_gridded_field()` — FieldGrid/Flow 的 120 行重复插值段；
8. netVisual 的 svg/png/pdf 三份绘制循环、rankNet 两分支归零块、enriched.only 过滤块。

## 8. 批次划分（更新）

- **P1-P4**：已完成对应核心入口迁移和三类 fixture 的图形/交互回归；
- **P5 收尾**：已完成 netP 显式筛选、selectK、netMappingDEG（含 merged dataset 对齐）、barPlot 和 signaling-role 消费者的回归；全零信号列由共享预处理排除，并对三个 fixture 输出收尾图；
- **下一步**：完成 package install、源码 parse、diff check，处理 issue log 中仍未覆盖的可选入口/科学语义边界。
