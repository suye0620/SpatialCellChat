# Agent Note: 可视化函数全量迁移 / 重构 / 优化 总规划

Status: proposed（规划，未动代码）
Governance: v1
Date: 2026-11-04
Decision type: plan
Scope: R/visualization.R（55 函数）、R/visualization_widget.R（9）、R/spatial.R 辅助；测试底座 `tests_dev/`
Owner: maintainer
Impact: high
Supersedes: none（继承 docs/VISUALIZATION_AUDIT.md 与 docs/VISUALIZATION_SEMANTICS.md 的全部结论）
Superseded by: none

## 0. 输入与约束

输入：`docs/VISUALIZATION_AUDIT.md`（2026-08-09 审计：崩溃级修复状态、ggplot2 4.0.3 语法表、plotly 问题清单、P0–P4 优先级）与 `docs/VISUALIZATION_SEMANTICS.md`（23 张八字段语义卡 + V-001…V-018 稳定问题表）。

约束（用户 2026-11-03 指定）：
1. 迁移/重构/优化后，每个函数用 **模拟数据 + Slide-seq + Visium** 三套数据测试，图以 **PDF** 存到 `tests_dev/` 下的合适位置；
2. 先静态、后交互，两大步；
3. 决策点用卡片；
4. 先规划后改动；已迁移的上游不动，本轮只碰可视化；
5. 确需改上游（bug/性能）先报告。

## 1. 现状基线（来自两份文档 + 本次扫描）

- `visualization.R` 55 函数、`visualization_widget.R` 9、`spatial.R` 18；`visualization.R` 9800 行是最大单文件。
- 审计口径下**已适配 11-slot**：`spatialDimPlot`、`spatialDimPlotPoints`、`spatialC2CDistPlot`（原 spatialCCCDistPlot）、`plotly_spatialDimPlot`。
- **确认阻塞**：`spatialFeaturePlot`、`plotly_spatialLRpairPlot`（读旧 `@data*` / `@net$prob.cell`）；`netVisual_heatmap`（`@options$mode` + `net$pval`）；`plotGeneExpression`（`@data.signaling`）；Gi/Lee 外层依赖旧通信结构。
- 审计已给出：P1 收敛重复代码、P2 schema 适配、P3 Gi/Lee/trace 决策、P4 plotly 现代化（3D stack 原型已验证）。
- V-001…V-018 是未决问题表，实施时逐条走卡片。

## 2. 测试底座（P0，先于一切函数改动）

### 2.1 三套数据

| 数据 | 来源 | 现状 |
|---|---|---|
| 模拟 | 现测试脚本内的合成 fixture（小网格 + 推理完成对象） | **缺固化 builder**，需抽成 `tests_dev/test_data/make_simulated_chat.R` |
| Slide-seq | `MBM05Chat_pre_inference.rds`（**全量 29,536 细胞**，用户 2026-11-04 指示不取子集） | 首选 benchmark 产物 `new_postfilter_object.rds`（computeCommunProb→filterProbability 已就绪）→ 补 avg/pathway/aggregate/centrality；缓存 `slideseq_full_chat.rds`（1.5 GB，8 组、1563 层）。需 `future.globals.maxSize=20GB`。P0 基线：C2C enriched 1035s、geometry 28s、DimPlot 21s，3/3 PASS |
| Visium | `visium_human_psoriasis_PP1_...RData` | 审计已验证：905 spot、1 raster、600×595×3 图像，`validateSpatialCellChat()` TRUE |

### 2.2 figure 输出底座

- 目录：`tests_dev/figures/<批次>/<数据集>/<函数>.pdf`；
- runner：`tests_dev/render_viz_figures.R`——逐函数渲染 PDF，单函数失败不中断，输出 PASS/FAIL 清单；
- 已知依赖缺口：`scatterpie` 未安装（比例饼图跳过）；`ComplexHeatmap` 等在 devtools::check 缺失但本地已装（审计 §7）。

## 3. 批次划分（先静态后交互）

| 批次 | 内容 | 语义卡 | 主要 V-* |
|---|---|---|---|
| **P1 静态-空间点图** | `spatialFeaturePlot`（阻塞级）→ `spatialVisual_scoring` → `spatialTopicPlot` → `spatialGiPlot` → `spatialLeePlot`；`spatialDimPlot(Points)` 作已适配基线回归 | §1.1、§4、§6.1 | V-003/004、V-010、V-011、V-012 |
| **P2 静态-网络图** | `netVisual_aggregate`/`individual`（+circle/hierarchy/chord/spatial 四 renderer）→ `netVisual_bubble` → `netVisual_heatmap` → `netVisual_diffInteraction` → `netVisual_chord_gene`；`doHeatmap`/`extract_max` 语法修复 | §2、§3 | V-001/002/006/007/008/009 |
| **P3 静态-场与分解** | `netVisual_CommunField`/`ZoomIn`/`FieldGrid`/`Flow` → `netVisual_TopicComposition`/`TopicSignaling` → `netAnalysis_river`/`dot`/`embedding`×4 | §5.1/5.2、§6.2/6.3、§7 | V-013/014/015/016/017 |
| **P4 交互** | `plotly_spatialDimPlot`（基线）→ `plotFeatures` → `plotly_spatialLRpairPlot` + `_shiny`（observer 修复）→ Gi 3d → `netVisual_CommunFlowWidget` → 3D stack 生产化（审计 §5.2 原型转正） | §1.2、§5.3、审计§5 | V-004/018 |
| **P5 收尾** | 全量 PDF 回归、`VISUALIZATION_AUDIT.md` 状态刷新、V-* 表清账、man/ 与 roxygen 补齐 | — | — |

每个批次固定动作：**方案 → 卡片决策 → 实现 → 三数据集出图 → 回归 → 提交**。批内按"改一行 / 要重写 / 该归档"分档推进。

## 4. 上游越界预报（约束 5）

以下问题源在上游，本轮只报告不擅改：
- **V-002**（阈值等号 `>` vs `>=`）：生产者在 `modeling.R`（pathway producer）与 `analysis.R`（bubble 表）；统一属科学语义更正，需卡片；
- **V-003**（`extractEnrichedLR` 缓存与实时筛选不一致、空交集不清空）：在 `analysis.R`，`spatialFeaturePlot` 的 signaling 模式被它卡住，届时单独报告；
- **V-013**（field 读取路径 `slot$field` vs `slot$cell$field`）：生产者 `computeCommunField()` 在 `analysis.R`，切路径属于上游一行改动，届时报告；
- **`netP$centr` 误读**：`netVisual_spatial` 的 incoming/outgoing 权重选项读 `netP$centr`，而生产者写 `centr/centr.cell`——如需改，报告后卡片。

## 5. 并行策略

- **P0.1 清单分级**：3 个只读 agent 并行（visualization.R / visualization_widget.R+spatial.R / 全仓 `net` 参数流），产出"函数 → 状态（能跑/半迁移/坏）→ 依赖 V-* → 建议档位"表；我合并去重。
- **批内实施**：串行为主（同文件写入冲突），只有跨文件只读复查才并行。
- 测试 runner 可与实现并行准备。

## 6. 第一步行动

**P0：盘点 + 底座**（只读 + 新增测试文件，不碰任何可视化实现）：

1. 派 3 个只读 agent 出全量函数分级清单（我合并成 `docs/VISUALIZATION_MIGRATION_PLAN.md` 的总表：函数 / 类别 / 状态 / V-* / 批次归属）；
2. 固化三套数据 fixture 脚本（simulated builder、slide-seq 窗口、visium 导入）；
3. 建 `tests_dev/figures/` 与 `render_viz_figures.R` runner；
4. 用已适配的 `spatialC2CDistPlot`、`spatialDimPlot` 跑通 runner，验证底座（三数据集 × 2 函数 × PDF）。

P0 完成的标志：总表 + runner 跑通，P1 的卡片即可发。

## 8. 包用法核查（用户 2026-11-04 补充，纳入每个批次的标准动作）

四个官方文档源，实施时以当前安装版本核对（ggplot2 4.0.3、plotly 4.12.0，ComplexHeatmap/circlize 以 DESCRIPTION 为准）：

- ggplot2: https://ggplot2.tidyverse.org/
- ComplexHeatmap: https://jokergoo.github.io/ComplexHeatmap/
- circlize: https://jokergoo.github.io/circlize/
- plotly: https://plotly.com/r/

每批次"实现"步骤新增一个固定小节——**包用法核查**，逐函数确认调用是否符合官方建议与惯用写法，偏离项列出并修正：

| 包 | 核查要点（审计已发现的问题加粗） |
|---|---|
| ggplot2 | `aes()` + `after_stat()` 取代 **`aes_string`（deprecated）与 `..density..`**；不读 **`ggplot_build()` 的 panel internals**（改 `layer_scales()`）；`linewidth` 而非 `size`；`coord_fixed()` 与 `theme(aspect.ratio=)` **二选一**；主题收敛到统一 helper；`na.value` 显式 |
| ComplexHeatmap | `Heatmap()`/`HeatmapAnnotation()` 官方参数而非手写 grid；cluster 距离与方法显式传参；legend 合并用官方 `heatmap_legend` 机制；`draw()` 显式调用（PP1 已发现需显式 draw） |
| circlize | `chordDiagram()` 的 `grid.col`/`order`/透明度/去自环；`circos.par(gap.after=)` 显式；流向箭头用官方参数 |
| plotly | 单 trace 单 colorbar（长表 rbind）、**显式 `marker$colorbar` 取代 `pl$x$data[[N]]` 硬编码**；大数据 `toWebGL()`；hover HTML 统一转义；layout 全部走共享 helper；空 trace 条件添加 |

产出：每批次的实施结果小节附"包用法核查"记录——改了哪些偏离项、依据哪一节官方文档。审计 §4 与 §5 已列出的偏离项直接映射进对应批次，不重复开工单。

## 7. 待确认

- 批次顺序与范围是否认可（先 P0，然后 P1 静态点图）；
- figure 目录与 runner 位置（建议 `tests_dev/figures/` + `tests_dev/render_viz_figures.R`）；
- P0 是否现在启动。
