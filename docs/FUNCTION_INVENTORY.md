# SpatialCellChat 函数信息清单与重构进度
> 生成日期：2026-10-26　范围：`R/` 下 10 个脚本、300 个函数　基线：pinned upstream `be88f30`（`renv/library/windows/R-4.4/.../SpatialCellChat`）

本清单是**逐函数事实档案**，用于评估 11-slot 重构进度；不改变任何代码，也不代表迁移验收通过。

## 0. 关键结论
1. **包当前无法安装。** `R CMD INSTALL` 失败：`undefined exports: liftSpatialCellChat, mergeSpatialCellChat, subsetSpatialCellChat, updateSpatialCellChat`。NAMESPACE 仍导出这四个函数，但 `R/` 下已无定义。测试全部走 `source()` 直连源码，所以 `tests_dev/` 全绿也掩盖了这个问题。
2. **进度分布不均。** 300 个函数中，`迁移11slot` 24 个、`新增(11slot)` 51 个、`部分迁移` 9 个、`重构` 23 个、`新增` 54 个、`未动` 139 个。迁移集中在类定义、建模主链路与分析聚合层；**可视化层几乎没动**。
3. **计算链路已闭环，展示与下游未闭环。** `computeCommunProb` → `filterProbability` → `computeAvgCommunProb` → `computeCommunProbPathway` → `aggregateNet` 全部只读写新槽；但 `netVisual_*` 家族仍读 `net$prob` / `net$pval`，在新对象上直接失效。
4. **死代码与兼容层仍在。** 明确不可达 4 个，疑似 40 个，主要来自 `@data` / `@var.features` / `net$prob.cell` 等已删槽位的旧工具链。
5. **依赖可分三类。** 函数引用的外部包共 55 个，均已安装；未装的 Suggests 包 15 个（`RcppML`, `SNFtool`, `SingleCellExperiment`, `SummarizedExperiment`, `callr`, `gg.gap`, `metR`, `oce`, `pryr`, `raster`, `rliger`, `rstudioapi`, `scatterpie`, `spdep`, `wordcloud`），只影响少数分析/可视化的可选分支。

## 1. 口径
- **迁移11slot**：函数体只读写 11-slot 结构（`net$cell$prob`、`net$group$prob/pval`、`netP$group$prob`、`SparseChatArray`、`@images` 等），无旧字段。
- **新增(11slot)**：基线不存在、且属于 11-slot 基础设施的新函数（`SparseChatArray` 族、`.sc_*` 工具、访问器）。
- **部分迁移**：新旧字段混用，或仍依赖 `my_as_sparse3Darray` / `net$tmp$*` 兼容层。
- **重构**：相对基线有实质改动，但改动与槽位无关（性能、接口、日志）。
- **新增**：基线不存在，且不属于 11-slot 基础设施。
- **未动**：函数体与基线逐字一致（忽略空白与注释）。
- 标签由“基线对比”（R 解析函数体后逐字比对）与“人工读码”双重判定；两者冲突时以代码事实为准并在此修正。

## 2. 进度总览
| 文件 | 函数数 | 迁移11slot | 新增(11slot) | 部分迁移 | 重构 | 新增 | 未动 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `R/modeling.R` | 39 | 7 | 0 | 2 | 4 | 6 | 20 |
| `R/spatial.R` | 10 | 0 | 0 | 2 | 4 | 3 | 1 |
| `R/analysis.R` | 47 | 7 | 0 | 1 | 0 | 2 | 37 |
| `R/visualization.R` | 62 | 2 | 0 | 1 | 8 | 4 | 47 |
| `R/SpatialCellChat_class.R` | 70 | 3 | 51 | 0 | 1 | 14 | 1 |
| `R/utilities.R` | 45 | 5 | 0 | 3 | 4 | 8 | 25 |
| `R/database.R` | 9 | 0 | 0 | 0 | 2 | 0 | 7 |
| `R/visualization_widget.R` | 9 | 0 | 0 | 0 | 0 | 9 | 0 |
| `R/RcppExports.R` | 5 | 0 | 0 | 0 | 0 | 4 | 1 |
| `R/test_helpers.R` | 4 | 0 | 0 | 0 | 0 | 4 | 0 |
| **合计** | **300** | **24** | **51** | **9** | **23** | **54** | **139** |

## 3. 依赖核查（renv）
环境基准：R 4.4.3，`renv.lock` 237 个包，`renv::status()` 一致；本机库另有对拍用的上游 `SpatialCellChat@be88f30`，未写入锁。

**已安装、被函数引用的外部包**（按引用函数数排序，取前 20）：

| 包 | 函数数 | 代表函数 |
|---|---:|---|
| `Matrix` | 48 | `createPspatialFrom_dspatial`、`relabelSpatialCellChat`、`computeAvgCommunProb_LR_Visium_Sum` |
| `ggplot2` | 48 | `communicationDistPlot2`、`communicationDistPlot`、`computeGridSize` |
| `dplyr` | 33 | `relabelSpatialCellChat`、`aggregateNet`、`computeExpr_complex` |
| `Seurat` | 21 | `identifyCellTopics`、`identifyTopicCommunities`、`extractTopicSignaling` |
| `RColorBrewer` | 14 | `identifyCommunicationPatterns`、`netAnalysis_signalingRole_network`、`netAnalysis_signalingRole_heatmap` |
| `patchwork` | 13 | `communicationDistPlot2`、`communicationDistPlot`、`spatialCCCDistPlot` |
| `scales` | 12 | `spatialFeaturePlot`、`spatialVisual_scoring`、`netVisual_circle` |
| `ggrepel` | 11 | `netAnalysis_signalingRole_scatter`、`netAnalysis_diff_signalingRole_scatter`、`netAnalysis_signalingChanges_scatter` |
| `igraph` | 10 | `computeCentralityLocal`、`node_distance`、`alpha_centrality` |
| `methods` | 10 | `netAnalysis_computeCentrality`、`extractEnrichedLR`、`computeNetSimilarity` |
| `reshape2` | 10 | `netVisual_hierarchy1`、`netVisual_hierarchy2`、`netVisual_circle` |
| `stringr` | 9 | `aggregateNet`、`subsetCommunication_internal`、`identifyCellTopics` |
| `purrr` | 8 | `HillFunctionFordataLR`、`netVisual_CommunFieldGrid`、`netVisual_CommunFlow` |
| `sf` | 8 | `.sc_grid_membership`、`computeGridSize`、`makeGridSpatialCellChat` |
| `cowplot` | 8 | `netAnalysis_contribution`、`netAnalysis_river`、`netVisual_embeddingZoomIn` |
| `ComplexHeatmap` | 8 | `identifyCommunicationPatterns`、`netAnalysis_signalingRole_network`、`netAnalysis_signalingRole_heatmap` |
| `cli` | 7 | `suggestK`、`netVisual_CommunField`、`netVisual_CommunFieldGrid` |
| `progressr` | 6 | `computeCommunProb`、`computeAvgCommunProb`、`filterCommunication` |
| `future` | 5 | `computeCommunProb`、`filterProbability`、`.computeExprGroup_complex` |
| `plotly` | 5 | `plotFeatures`、`plotly_spatialLRpairPlot_shiny`、`plotly_spatialDimPlot` |

**未安装的 Suggests 包**（相关分支在本机不可执行）：

| 包 | 受影响函数 | 状态 |
|---|---|---|
| `RcppML` | `identifyCellTopics`、`suggestK1` | 未安装 |
| `SNFtool` | `identifyTopicCommunities` | 未安装 |
| `SingleCellExperiment` | `createSpatialCellChat` | 未安装 |
| `SummarizedExperiment` | `createSpatialCellChat` | 未安装 |
| `callr` | `start_session` | 未安装 |
| `gg.gap` | `rankNet` | 未安装 |
| `metR` | `netVisual_CommunField`、`netVisual_CommunFieldGrid`、`netVisual_CommunFlow`、`plotStatistics_Gi` | 未安装 |
| `oce` | `.sc_widget_interpolate`、`netVisual_CommunFieldGrid`、`netVisual_CommunFlow` | 未安装 |
| `pryr` | `my_obj_size` | 未安装 |
| `raster` | `computeColocalization` | 未安装 |
| `rliger` | `suggestK` | 未安装 |
| `rstudioapi` | `setpwd` | 未安装 |
| `scatterpie` | `spatialDimPlotPoints`、`visualizePie` | 未安装 |
| `spdep` | `plotStatistics_Gi`、`plotStatistics_Lee`、`spatialLeePlot` | 未安装 |
| `wordcloud` | `computeEnrichmentScore` | 未安装 |

## 4. 死代码与兼容层
判定依据为静态调用扫描（`R/` 内部调用关系）+ 槽位指纹 + 人工读码；「疑似」不等于可删，删除前需确认外部用户可能直接调用导出的部分。
### 4.1 明确不可达 / 重复定义

| 函数 | 文件:行 | 判定 | 依据 |
|---|---|---|---|
| `.H_to_conn_mat` | `R/analysis.R:1191` | 是 | 仅在注释中被引用；内部调用的 `max_factor_rcpp` 全包未定义 |
| `setIdent` | `R/SpatialCellChat_class.R:1314` | 是 | （重复定义，被 utilities.R:1122 覆盖） |
| `CellChat2Seurat` | `R/utilities.R:712` | 是 | 仍读已删除的 `@data` 槽，调用即报错 |
| `projectData` | `R/utilities.R:1820` | 是 | 仍读写已删除的 `@data.signaling` / `@data.project` 槽 |

### 4.2 疑似（无调用者、被取代、或依赖已删字段）

| 函数 | 文件:行 | 依据 |
|---|---|---|
| `HillFunctionFordataLR` | `R/modeling.R:497` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `elementwiseProductForSparseMatandDenseMat` | `R/modeling.R:518` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `myElementwiseProduct` | `R/modeling.R:547` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeAvgCommunProb_LR_Sum` | `R/modeling.R:1397` | 无内部调用者；相对基线有改动但已无调用方 |
| `computeAvgCommunProb_LR_Avg` | `R/modeling.R:1468` | 无内部调用者；相对基线有改动但已无调用方 |
| `modifySparse3Darray` | `R/modeling.R:2056` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeAveExpr` | `R/modeling.R:2498` | 无内部调用者；相对基线有改动但已无调用方 |
| `.computeExprGroup_complex` | `R/modeling.R:2577` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeExpr_coreceptor` | `R/modeling.R:2642` | 无内部调用者；相对基线有改动但已无调用方 |
| `.computeExprGroup_coreceptor` | `R/modeling.R:2691` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeExprGroup_agonist` | `R/modeling.R:2749` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeExprGroup_antagonist` | `R/modeling.R:2782` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeExpr_agonist` | `R/modeling.R:2815` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `computeExpr_antagonist` | `R/modeling.R:2851` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `communicationDistPlot2` | `R/spatial.R:11` | 无内部调用者；函数体仍命中旧槽位（net$prob / net$pval / net$prob.cell / net$tmp$）；相对基线有改动但已无调用方 |
| `communicationDistPlot` | `R/spatial.R:95` | 无内部调用者；函数体仍命中旧槽位（net$prob / net$pval / net$prob.cell / net$tmp$）；相对基线有改动但已无调用方 |
| `selectK` | `R/analysis.R:971` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `netVisual_CommunField` | `R/visualization.R:5477` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `netVisual_CommunFieldZoomIn` | `R/visualization.R:6460` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `.SparseChatArray_resolve_single_index` | `R/SpatialCellChat_class.R:282` | 非导出、无内部调用者 |
| `mergeChatIntoSeurat` | `R/utilities.R:123` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `==.sparse3Darray` | `R/utilities.R:368` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `my_obj_size` | `R/utilities.R:394` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `setpwd` | `R/utilities.R:406` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `my_as_sparse3Darray` | `R/utilities.R:501` | 无内部调用者；相对基线有改动但已无调用方 |
| `scaleData` | `R/utilities.R:982` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `sketchData` | `R/utilities.R:1044` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `setIdent` | `R/utilities.R:1122` | 与 `R/utilities.R:1122` 重复定义，按加载顺序被后者覆盖 |
| `updateClusterLabels` | `R/utilities.R:1172` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `identifyOverExpressedLigandReceptor` | `R/utilities.R:1656` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `projectAndRecombine` | `R/utilities.R:1864` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `projectOnNetwork` | `R/utilities.R:1884` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `randomWalkBySolve` | `R/utilities.R:1903` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `l1NormalizeColumns` | `R/utilities.R:1923` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `l1NormalizeRows` | `R/utilities.R:1934` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `projectFromNetworkRecombine` | `R/utilities.R:1946` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `runPCA` | `R/utilities.R:1967` | 无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `.error_if_no_Seurat` | `R/utilities.R:2072` | 非导出、无内部调用者；相对基线未改动，属被新实现取代后遗留 |
| `showDatabaseCategory` | `R/database.R:10` | 无内部调用者；相对基线有改动但已无调用方 |
| `pieChart` | `R/database.R:49` | 无内部调用者；相对基线有改动但已无调用方 |

## 5. 逐函数明细

### 5.1 `R/modeling.R`（39 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `.sc_expr_row` | 7 | 新增 | 取基因或复合体的单行表达 | 计算 | 子模块 | ↑ computeCommunProb, computeAvgCommunProb; ↓ geometricMean | — | — |
| `.sc_cofactor_factor` | 36 | 新增 | 计算辅因子逐细胞调制因子 | 计算 | 子模块 | ↑ computeCommunProb; ↓ — | — | — |
| `.sc_csr_build` | 67 | 新增 | 手工构建稀疏矩阵CSR索引 | 计算 | 子模块 | ↑ filterProbability; ↓ — | — | — |
| `.sc_pattern_row` | 75 | 新增 | 取距离矩阵的正值模式行 | 计算 | 子模块 | ↑ filterProbability; ↓ — | — | — |
| `.sc_gate_threshold` | 93 | 新增 | 预计算表达比例门控阈值 | 计算 | 子模块 | ↑ computeAvgCommunProb; ↓ — | — | — |
| `computeCommunProb` | 159 | 迁移11slot | 单细胞级通讯概率推断 | 计算 | 入口 | ↑ 用户; ↓ computeCellDistance, cpp_prob_layer | future, progressr | — |
| `createPspatialFrom_dspatial` | 465 | 未动 | 距离矩阵转空间权重矩阵 | 计算 | 子模块 | ↑ computeCommunProb; ↓ — | Matrix | — |
| `HillFunctionFordataLR` | 497 | 未动 | 对非零值施加Hill函数 | 计算 | 子模块 | ↑ —; ↓ — | purrr | 疑似 |
| `elementwiseProductForSparseMatandDenseMat` | 518 | 未动 | 稀疏与稠密矩阵逐元素乘 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `myElementwiseProduct` | 547 | 未动 | 稀疏矩阵与外积向量逐元素乘 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `relabelSpatialCellChat` | 586 | 迁移11slot | 合并细胞群并重算组级通讯 | 分析 | 入口 | ↑ 用户; ↓ computeAvgCommunProb, aggregateNet | Matrix, dplyr | — |
| `computeAvgCommunProb_Visium` | 734 | 部分迁移 | 反卷积数据组级通讯聚合 | 分析 | 主流程 | ↑ relabelSpatialCellChat; ↓ computeAvgCommunProb_LR_Visium_Avg, computeColocalization | — | — |
| `computeAvgCommunProb_LR_Visium_Sum` | 975 | 未动 | 单LR组级通讯求和聚合 | 分析 | 子模块 | ↑ computeAvgCommunProb_Visium; ↓ — | Matrix | — |
| `computeAvgCommunProb_LR_Visium_Avg` | 1000 | 未动 | 单LR组级通讯均值聚合 | 分析 | 子模块 | ↑ computeAvgCommunProb_Visium; ↓ — | Matrix | — |
| `computeAvgCommunProb` | 1156 | 迁移11slot | 组级通讯聚合与置换检验 | 分析 | 主流程 | ↑ relabelSpatialCellChat; ↓ cpp_group_avg_obs, cpp_group_avg_perm | progressr | — |
| `computeAvgCommunProb_LR_Sum` | 1397 | 重构 | 单LR组级通讯求和旧实现 | 分析 | 子模块 | ↑ —; ↓ — | Matrix | 疑似 |
| `computeAvgCommunProb_LR_Avg` | 1468 | 重构 | 单LR组级通讯均值旧实现 | 分析 | 子模块 | ↑ —; ↓ — | Matrix | 疑似 |
| `.sc_sum_layers` | 1649 | 新增 | 逐层求和通讯数组 | 分析 | 子模块 | ↑ computeCommunProbPathway, aggregateNet; ↓ cpp_sum_layers | — | — |
| `computeCommunProbPathway` | 1678 | 迁移11slot | 汇总通路级通讯概率 | 分析 | 主流程 | ↑ relabelSpatialCellChat, updateClusterLabels; ↓ .sc_sum_layers | — | — |
| `filterCommunication` | 1880 | 迁移11slot | 按细胞数与连接数过滤通讯 | 分析 | 主流程 | ↑ relabelSpatialCellChat; ↓ .sc_map_sparse_layers | Matrix, progressr | — |
| `modifySparse3Darray` | 2056 | 未动 | 稀疏三维数组裁剪与二值化 | 基础设施 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `scMatrixTruncation` | 2105 | 未动 | 稀疏矩阵按阈值截断 | 基础设施 | 子模块 | ↑ filterProbability; ↓ — | Matrix | — |
| `filterProbability` | 2149 | 迁移11slot | 细胞级通讯显著性过滤 | 计算 | 入口 | ↑ 用户; ↓ scMatrixTruncation, my_future_lapply | Matrix, future | — |
| `aggregateNet` | 2330 | 迁移11slot | 汇总网络权重与连接数 | 分析 | 主流程 | ↑ relabelSpatialCellChat, updateClusterLabels; ↓ .sc_sum_layers | dplyr, stringr | — |
| `computeAveExpr` | 2498 | 部分迁移 | 计算各细胞群平均表达 | 分析 | 入口 | ↑ —; ↓ triMean | — | 疑似 |
| `computeExpr_complex` | 2550 | 重构 | 几何均值算复合体表达 | 计算 | 子模块 | ↑ computeExpr_LR; ↓ geometricMean | dplyr | — |
| `.computeExprGroup_complex` | 2577 | 未动 | 组水平复合体平均表达 | 计算 | 子模块 | ↑ —; ↓ geometricMean | dplyr, future, future.apply | 疑似 |
| `computeExpr_LR` | 2614 | 未动 | 组装配体受体表达矩阵 | 计算 | 子模块 | ↑ spatialFeaturePlot; ↓ computeExpr_complex | — | — |
| `computeExpr_coreceptor` | 2642 | 重构 | 计算共受体调制因子 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `.computeExprGroup_coreceptor` | 2691 | 未动 | 组水平共受体调制因子 | 计算 | 子模块 | ↑ —; ↓ — | future, future.apply, pbapply | 疑似 |
| `computeExprGroup_agonist` | 2749 | 未动 | 组水平激动剂调制因子 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `computeExprGroup_antagonist` | 2782 | 未动 | 组水平拮抗剂调制因子 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `computeExpr_agonist` | 2815 | 未动 | 细胞水平激动剂调制因子 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `computeExpr_antagonist` | 2851 | 未动 | 细胞水平拮抗剂调制因子 | 计算 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `geometricMean` | 2880 | 未动 | 计算几何均值 | 基础设施 | 子模块 | ↑ .sc_expr_row, computeExpr_complex; ↓ — | — | — |
| `triMean` | 2904 | 未动 | 计算Tukey三均值 | 基础设施 | 子模块 | ↑ computeAveExpr; ↓ — | — | — |
| `thresholdedMean` | 2915 | 未动 | 表达比例达标后取均值 | 基础设施 | 子模块 | ↑ netVisual_TopicComposition; ↓ — | Matrix | — |
| `identifyEnrichedInteractions` | 2937 | 未动 | 提取指定群对显著LR | 分析 | 入口 | ↑ 用户; ↓ — | dplyr | — |
| `computeRegionDistance` | 3024 | 未动 | 计算区域间平均距离 | 计算 | 入口 | ↑ 用户; ↓ — | BiocNeighbors | — |

### 5.2 `R/spatial.R`（10 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `communicationDistPlot2` | 11 | 部分迁移 | 按信号类型绘距离密度 | 可视化 | 入口 | ↑ 用户; ↓ my_future_lapply | ggplot2, patchwork | 疑似 |
| `communicationDistPlot` | 95 | 部分迁移 | 绘通讯距离分布密度图 | 可视化 | 入口 | ↑ 用户; ↓ my_future_lapply | ggplot2, patchwork | 疑似 |
| `.sc_resolve_grid_size` | 187 | 新增 | 解析并校验网格尺寸 | 基础设施 | 子模块 | ↑ computeGridSize, makeGridSpatialCellChat; ↓ — | BiocNeighbors | — |
| `.sc_grid_membership` | 271 | 新增 | 计算点与网格稀疏归属 | 基础设施 | 子模块 | ↑ computeGridSize, makeGridSpatialCellChat; ↓ — | sf | — |
| `computeGridSize` | 313 | 重构 | 估计网格边长并预览 | 可视化 | 入口 | ↑ 用户; ↓ .sc_resolve_grid_size | ggplot2, sf | — |
| `makeGridSpatialCellChat` | 450 | 重构 | 按网格聚合细胞为对象 | 基础设施 | 入口 | ↑ identifyOverExpressedGenes; ↓ createSpatialCellChat | Matrix, sf | — |
| `.spatial_distance_cache_matches` | 574 | 新增 | 校验距离缓存参数 | 基础设施 | 子模块 | ↑ spatialCCCDistPlot; ↓ — | — | — |
| `computeCellDistance` | 628 | 重构 | 计算细胞距离与接触邻接 | 计算 | 入口 | ↑ computeCommunProb, spatialCCCDistPlot; ↓ createCellCellContactMatrixFrom_dspatial | BiocNeighbors, Matrix | — |
| `createCellCellContactMatrixFrom_dspatial` | 743 | 重构 | 由距离矩阵生成接触邻接 | 计算 | 子模块 | ↑ computeCellDistance; ↓ — | Matrix | — |
| `computeColocalization` | 782 | 未动 | 置换检验细胞群共定位 | 分析 | 子模块 | ↑ computeAvgCommunProb, computeAvgCommunProb_Visium; ↓ — | raster (!) | — |

### 5.3 `R/analysis.R`（47 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `computeCentralityLocal` | 9 | 未动 | 计算单个网络层的节点中心性指标 | 计算 | 子模块 | ↑ netAnalysis_computeCentrality; ↓ — | Matrix, igraph, sna | — |
| `netAnalysis_computeCentrality` | 83 | 迁移11slot | 计算各信号层细胞群或细胞中心性 | 分析 | 入口 | ↑ 用户, updateClusterLabels; ↓ computeCentralityLocal, my_future_sapply | methods | — |
| `netAnalysis_contribution` | 259 | 未动 | 绘制各LR对信号贡献条形图 | 可视化 | 入口 | ↑ 用户; ↓ searchPair | cowplot, ggplot2, dplyr | — |
| `rankNet` | 445 | 未动 | 按信息流或互作数排序信号网络 | 可视化 | 入口 | ↑ 用户; ↓ ggPalette, CellChat_theme_opts | ggplot2, gg.gap(!) | — |
| `extractEnrichedLR` | 856 | 未动 | 提取通路显著富集的LR对及基因 | 分析 | 入口 | ↑ plotGeneExpression, spatialFeaturePlot; ↓ extractEnrichedLR_internal | BiocGenerics, methods | — |
| `extractEnrichedLR_internal` | 913 | 未动 | 按通路筛显著LR对并取其基因 | 分析 | 主流程 | ↑ extractEnrichedLR; ↓ searchPair, extractGeneSubset | dplyr | — |
| `selectK` | 971 | 未动 | 用NMF稳定性曲线选择模式数 | 分析 | 入口 | ↑ 用户; ↓ — | NMF, ggplot2 | 疑似 |
| `suggestK` | 1044 | 未动 | 扫描k值评估iNMF目标误差 | 分析 | 入口 | ↑ 用户; ↓ .plotSuggestK | cli, rliger(!) | — |
| `.plotSuggestK` | 1109 | 未动 | 绘制k值误差与方差诊断图 | 可视化 | 子模块 | ↑ suggestK, suggestK1; ↓ CellChat_theme_opts | dplyr, ggplot2 | — |
| `.H_to_conn_mat` | 1191 | 未动 | 由H矩阵构建细胞连通矩阵 | 基础设施 | 子模块 | ↑ —（仅注释）; ↓ max_factor_rcpp（未定义） | Matrix | 是 |
| `identifyCommunicationPatterns` | 1249 | 未动 | NMF分解识别细胞通讯模式 | 分析 | 入口 | ↑ 用户; ↓ scPalette, scaleMat | NMF, RColorBrewer, ComplexHeatmap | — |
| `computeNetSimilarity` | 1345 | 未动 | 计算信号网络间功能或结构相似性 | 分析 | 入口 | ↑ 用户; ↓ buildSNN, computeNetD_structure | methods | — |
| `computeNetSimilarityPairwise` | 1421 | 未动 | 计算多样本网络两两相似性 | 分析 | 入口 | ↑ 用户; ↓ buildSNN, computeNetD_structure | methods | — |
| `netEmbedding` | 1533 | 未动 | 对网络相似性矩阵做UMAP降维 | 分析 | 入口 | ↑ 用户; ↓ runUMAP | uwot, methods | — |
| `netClustering` | 1594 | 未动 | 对网络嵌入结果聚类分型 | 分析 | 入口 | ↑ 用户; ↓ computeEigengap, my_future_sapply | Matrix, ggplot2 | — |
| `buildSNN` | 1671 | 未动 | 由kNN构建共享最近邻矩阵 | 计算 | 子模块 | ↑ computeNetSimilarity, computeNetSimilarityPairwise; ↓ ComputeSNN | FNN, Matrix | — |
| `computeEigengap` | 1699 | 未动 | 由拉普拉斯特征值间隙推断簇数 | 计算 | 子模块 | ↑ netClustering, identifyTopicCommunities; ↓ computeLaplacian | ggplot2 | — |
| `computeLaplacian` | 1750 | 未动 | 计算归一化图拉普拉斯特征值 | 计算 | 子模块 | ↑ computeEigengap; ↓ — | Matrix, RSpectra | — |
| `rankSimilarity` | 1785 | 未动 | 按嵌入距离排序共有通路 | 可视化 | 入口 | ↑ netVisual_embedding; ↓ — | ggplot2, methods | — |
| `compareInteractions` | 1876 | 未动 | 比较各数据集互作数量或强度 | 可视化 | 入口 | ↑ 用户; ↓ ggPalette | ggplot2 | — |
| `rankNetPairwise` | 1962 | 未动 | 排序任意两群间的LR互作 | 分析 | 入口 | ↑ 用户; ↓ — | — | — |
| `entropia` | 1996 | 未动 | 计算正向量香农熵 | 计算 | 子模块 | ↑ nnd, computeNetD_structure; ↓ — | — | — |
| `node_distance` | 2006 | 未动 | 计算图节点最短路径分布 | 计算 | 子模块 | ↑ nnd; ↓ — | Matrix, igraph | — |
| `nnd` | 2041 | 未动 | 计算节点距离分布及其熵 | 计算 | 子模块 | ↑ computeNetD_structure; ↓ entropia, node_distance | Matrix | — |
| `alpha_centrality` | 2059 | 未动 | 计算图alpha中心性向量 | 计算 | 子模块 | ↑ computeNetD_structure; ↓ — | igraph | — |
| `computeNetD_structure` | 2081 | 未动 | 计算两个网络的结构距离 | 计算 | 子模块 | ↑ computeNetSimilarity, computeNetSimilarityPairwise; ↓ nnd, alpha_centrality | igraph | — |
| `getMaxWeight` | 2189 | 未动 | 取各数据集互作指标最大值 | 分析 | 入口 | ↑ 用户; ↓ — | methods | — |
| `mergeInteractions` | 2222 | 未动 | 按合并分组汇总互作数与强度 | 分析 | 入口 | ↑ 用户; ↓ — | — | — |
| `subsetCommunication` | 2296 | 迁移11slot | 按细胞群或通路抽取通讯记录 | 分析 | 入口 | ↑ netMappingDEG, netVisual_bubble; ↓ subsetCommunication_internal | methods | — |
| `subsetCommunication_internal` | 2405 | 迁移11slot | 把稀疏通讯层展开为长表 | 分析 | 主流程 | ↑ subsetCommunication; ↓ searchPair | dplyr, stringr, BiocGenerics, magrittr | — |
| `netAnalysis_signalingRole_network` | 2633 | 未动 | 绘制单网络中心性热图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette | ComplexHeatmap, RColorBrewer, grid | — |
| `netAnalysis_signalingRole_scatter` | 2716 | 迁移11slot | 绘制细胞群收发角色散点图 | 可视化 | 入口 | ↑ 用户; ↓ aggregateNet, scPalette | ggplot2, ggrepel, grid | — |
| `netAnalysis_diff_signalingRole_scatter` | 2847 | 迁移11slot | 绘制两组收发角色差异散点图 | 可视化 | 入口 | ↑ 用户; ↓ mergeCellChat（未定义） | ggplot2, ggrepel, grid | — |
| `netAnalysis_signalingChanges_scatter` | 3011 | 迁移11slot | 绘制单群通路变化散点图 | 可视化 | 入口 | ↑ 用户; ↓ mergeCellChat（未定义） | ggplot2, ggrepel, plyr | — |
| `netAnalysis_signalingRole_heatmap` | 3192 | 未动 | 绘制通路信号强度热图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette | ComplexHeatmap, RColorBrewer, grid | — |
| `netMappingDEG` | 3296 | 未动 | 把DEG结果映射到通讯记录 | 分析 | 入口 | ↑ 用户; ↓ subsetCommunication | dplyr | — |
| `computeEnrichmentScore` | 3411 | 未动 | 计算LR富集分并绘制词云 | 分析 | 入口 | ↑ 用户; ↓ — | RColorBrewer, wordcloud(!), dplyr | — |
| `identifyCellTopics` | 3467 | 部分迁移 | 用NMF从细胞中心性提取主题 | 分析 | 入口 | ↑ 用户; ↓ suggestK1, CellChat2Seurat | Matrix, RcppML(!), Seurat, stringr | — |
| `identifyTopicCommunities` | 3595 | 未动 | 对细胞主题做空间聚类分群 | 分析 | 入口 | ↑ 用户; ↓ CellChat2Seurat, computeEigengap | Seurat, SNFtool(!), irlba, Matrix | — |
| `suggestK1` | 3812 | 未动 | 交叉验证选择NMF主题数 | 分析 | 入口 | ↑ identifyCellTopics; ↓ my_future_lapply, .plotSuggestK | RcppML(!), Matrix | — |
| `extractTopicSignaling` | 3943 | 未动 | 提取各主题代表性信号特征 | 分析 | 入口 | ↑ 用户; ↓ — | Seurat, dplyr | — |
| `.sc_communication_field_layer` | 4016 | 新增 | 计算单层通讯向量场 | 计算 | 子模块 | ↑ .sc_compute_communication_field; ↓ — | Matrix, methods | — |
| `.sc_compute_communication_field` | 4083 | 新增 | 汇总各层进出向量场 | 计算 | 子模块 | ↑ computeCommunField; ↓ .sc_communication_field_layer, SparseChatArray | — | — |
| `computeCommunField` | 4115 | 迁移11slot | 由细胞通讯概率算向量场 | 计算 | 入口 | ↑ 用户; ↓ .sc_compute_communication_field, .sc_validate_after_update | methods | — |
| `wkde2d` | 4211 | 未动 | 加权二维核密度估计 | 计算 | 子模块 | ↑ calculate_density; ↓ — | Matrix, ks | — |
| `get_dens` | 4271 | 未动 | 按网格取观测点密度值 | 计算 | 子模块 | ↑ calculate_density; ↓ — | — | — |
| `calculate_density` | 4312 | 未动 | 按方法估计并映射密度 | 计算 | 子模块 | ↑ plotStatistics_Gi; ↓ get_dens, wkde2d | ks | — |

### 5.4 `R/visualization.R`（62 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `generate_grid_nrows` | 9 | 新增 | 校验并生成 plotly 子图行数配置 | 可视化 | 子模块 | ↑ plotFeatures, plotly_spatialDimPlot; ↓ — | — | — |
| `generate_custom_scenes3d` | 16 | 新增 | 构建 plotly 三维场景布局参数 | 可视化 | 子模块 | ↑ plotFeatures, spatialGiPlot; ↓ — | — | — |
| `spatialCCCDistPlot` | 56 | 重构 | 绘制细胞间距离密度分布图 | 可视化 | 入口 | ↑ 用户; ↓ computeCellDistance | ggplot2, patchwork | — |
| `spatialDimPlotPoints` | 151 | 新增 | 按坐标绘制分组散点或饼图 | 可视化 | 子模块 | ↑ spatialDimPlot; ↓ scPalette | ggplot2, scatterpie (!) | — |
| `spatialDimPlot` | 326 | 重构 | 绘制空间分组图并叠加组织图像 | 可视化 | 入口 | ↑ 用户; ↓ spatialDimPlotPoints | — | — |
| `spatialFeaturePlot` | 437 | 未动 | 空间展示基因与配体受体表达 | 可视化 | 入口 | ↑ 用户; ↓ extractEnrichedLR, computeExpr_LR | ggplot2, RColorBrewer, patchwork, scales | — |
| `spatialVisual_scoring` | 724 | 未动 | 空间展示网络中心性得分 | 可视化 | 入口 | ↑ 用户; ↓ — | ggplot2, RColorBrewer, ggnewscale, patchwork, scales | — |
| `CellChat_theme_opts` | 944 | 未动 | CellChat 通用 ggplot 主题 | 可视化 | 子模块 | ↑ netVisual_embedding, netVisual_embeddingZoomIn; ↓ — | ggplot2 | — |
| `ggPalette` | 963 | 未动 | 生成 n 种 hcl 调色板颜色 | 可视化 | 子模块 | ↑ netVisual_bubble, netAnalysis_river; ↓ — | — | — |
| `scPalette` | 977 | 未动 | 生成单细胞可视化调色板 | 可视化 | 子模块 | ↑ netVisual_circle, netVisual_spatial; ↓ — | — | — |
| `netVisual` | 1051 | 未动 | 按布局绘制通路网络并导出图片 | 可视化 | 入口 | ↑ 用户; ↓ netVisual_circle, netVisual_hierarchy1 | svglite | — |
| `netVisual_aggregate` | 1440 | 未动 | 聚合配体受体绘制通路网络 | 可视化 | 入口 | ↑ 用户; ↓ netVisual_circle, netVisual_hierarchy1 | — | — |
| `netVisual_individual` | 1604 | 未动 | 逐个绘制配体受体对通讯网络 | 可视化 | 入口 | ↑ 用户; ↓ netVisual_circle, netVisual_hierarchy1 | — | — |
| `netVisual_hierarchy1` | 1795 | 未动 | 绘制层级图上半部通讯网络 | 可视化 | 入口 | ↑ netVisual, netVisual_aggregate; ↓ mycircle | Matrix, igraph, reshape2, shape | — |
| `netVisual_hierarchy2` | 1985 | 未动 | 绘制层级图下半部通讯网络 | 可视化 | 入口 | ↑ netVisual, netVisual_individual; ↓ mycircle | Matrix, igraph, reshape2, shape | — |
| `netVisual_circle` | 2176 | 未动 | 绘制环形细胞通讯网络图 | 可视化 | 入口 | ↑ netVisual, netVisual_individual; ↓ scPalette | Matrix, igraph, reshape2, scales | — |
| `mycircle` | 2310 | 未动 | igraph 自定义圆形顶点绘制回调 | 可视化 | 子模块 | ↑ netVisual_hierarchy1, netVisual_hierarchy2; ↓ — | — | — |
| `netVisual_spatial` | 2380 | 未动 | 在空间坐标上绘制通讯网络 | 可视化 | 入口 | ↑ netVisual, netVisual_individual; ↓ scPalette | Matrix, igraph, reshape2, ggnetwork, ggplot2 | — |
| `netVisual_diffInteraction` | 2572 | 未动 | 绘制两组间差异通讯环形图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette | Matrix, igraph, reshape2, scales | — |
| `netVisual_heatmap` | 2734 | 未动 | 绘制通讯网络热图 | 可视化 | 入口 | ↑ 用户; ↓ colorRamp3, scPalette | Matrix, RColorBrewer, ComplexHeatmap, reshape2 | — |
| `netVisual_barplot` | 2931 | 未动 | 绘制通讯数量条形图 | 可视化 | 入口 | ↑ 用户; ↓ barplot_internal, scPalette | reshape2 | — |
| `netVisual_bubble` | 3079 | 重构 | 绘制配体受体通讯气泡图 | 可视化 | 入口 | ↑ 用户; ↓ subsetCommunication, ggPalette | RColorBrewer, scales, stringr | — |
| `netVisual_chord_cell` | 3447 | 未动 | 绘制信号通路通讯弦图 | 可视化 | 入口 | ↑ 用户; ↓ netVisual_chord_cell_internal, searchPair | — | — |
| `netVisual_chord_cell_internal` | 3567 | 未动 | 由矩阵或数据框绘制弦图 | 可视化 | 入口 | ↑ netVisual, netVisual_chord_cell; ↓ scPalette | circlize, ComplexHeatmap, reshape2 | — |
| `netVisual_chord_gene` | 3738 | 未动 | 绘制配体受体基因弦图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette, searchPair | circlize, ComplexHeatmap, dplyr, reshape2, stringr | — |
| `netAnalysis_river` | 3978 | 未动 | 绘制潜在模式河流图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette, ggPalette | ggplot2, cowplot, ggalluvial | — |
| `netAnalysis_dot` | 4229 | 未动 | 绘制模式与细胞群关联点图 | 可视化 | 入口 | ↑ 用户; ↓ scPalette | ggplot2, dplyr | — |
| `netVisual_embedding` | 4337 | 未动 | 绘制信号网络二维流形图 | 可视化 | 入口 | ↑ 用户; ↓ CellChat_theme_opts, ggPalette | ggplot2, ggrepel | — |
| `netVisual_embeddingZoomIn` | 4431 | 未动 | 分簇放大绘制信号网络流形 | 可视化 | 入口 | ↑ 用户; ↓ CellChat_theme_opts, ggPalette | ggplot2, ggrepel, cowplot | — |
| `netVisual_embeddingPairwise` | 4518 | 未动 | 绘制两组联合流形图 | 可视化 | 入口 | ↑ 用户; ↓ CellChat_theme_opts, ggPalette | ggplot2, ggrepel | — |
| `netVisual_embeddingPairwiseZoomIn` | 4633 | 未动 | 分簇放大绘制联合流形图 | 可视化 | 入口 | ↑ 用户; ↓ CellChat_theme_opts, ggPalette | ggplot2, ggrepel, cowplot | — |
| `netVisual_TopicComposition` | 4766 | 未动 | 绘制 topic 与细胞群组成比例热图/点图 | 可视化 | 入口 | ↑ 用户; ↓ thresholdedMean, scPalette | ComplexHeatmap, RColorBrewer, Seurat, ggplot2, Matrix | — |
| `netVisual_TopicSignaling` | 5072 | 未动 | 绘制各 topic 的 signaling 载荷排名图 | 可视化 | 入口 | ↑ 用户; ↓ — | ggplot2, dplyr, ggrepel, patchwork, stringr | — |
| `vizDimLoadings` | 5323 | 未动 | 绘制 Seurat 主成分载荷排名图 | 可视化 | 入口 | ↑ 用户; ↓ — | Seurat, ggplot2, dplyr, ggrepel, patchwork | — |
| `netVisual_CommunField` | 5477 | 部分迁移 | 按细胞绘制通讯向量场箭头图 | 可视化 | 入口 | ↑ netVisual_CommunFieldZoomIn; ↓ my_as_sparse3Darray, scPalette | ggplot2, metR(!), dplyr, cli | 疑似 |
| `.sc_resolve_arrow_width` | 5743 | 新增 | 校验并规范化箭头线宽参数 | 基础设施 | 子模块 | ↑ netVisual_CommunFieldGrid, netVisual_CommunFlow; ↓ — | — | — |
| `netVisual_CommunFieldGrid` | 5810 | 迁移11slot | 网格聚合后绘制通讯向量场箭头图 | 可视化 | 入口 | ↑ 用户; ↓ .sc_resolve_grid_size, .sc_grid_membership | ggplot2, sf, metR(!), oce(!), dplyr, tidyr, purrr, cli | — |
| `netVisual_CommunFlow` | 6153 | 迁移11slot | 网格聚合后绘制通讯流线图 | 可视化 | 入口 | ↑ 用户; ↓ .sc_resolve_grid_size, .sc_grid_membership | ggplot2, sf, metR(!), oce(!), dplyr, tidyr, purrr, cli | — |
| `netVisual_CommunFieldZoomIn` | 6460 | 未动 | 并联绘制全场与局部放大通讯场 | 可视化 | 入口 | ↑ 用户; ↓ netVisual_CommunField | ggplot2, cowplot, stringr | 疑似 |
| `showDatabaseCategory` | 6549 | 重构 | 三张饼图汇总数据库构成（重复定义） | 可视化 | 入口 | ↑ 用户; ↓ pieChart | dplyr, cowplot | — |
| `pieChart` | 6587 | 重构 | 按分组绘制占比饼图（重复定义） | 可视化 | 子模块 | ↑ showDatabaseCategory; ↓ — | ggplot2, dplyr, scales, forcats, ggrepel | — |
| `plotGeneExpression` | 6632 | 重构 | 封装 Seurat 绘制信号基因表达 | 可视化 | 入口 | ↑ 用户; ↓ StackedVlnPlot, dotPlot | Seurat | — |
| `dotPlot` | 6703 | 未动 | 封装 Seurat 点图并调整主题 | 可视化 | 子模块 | ↑ plotGeneExpression; ↓ — | Seurat, ggplot2 | — |
| `StackedVlnPlot` | 6760 | 未动 | 纵向堆叠多基因小提琴图 | 可视化 | 主流程 | ↑ plotGeneExpression; ↓ modify_vlnplot, extract_max | Seurat, patchwork, purrr, ggplot2 | — |
| `modify_vlnplot` | 6818 | 未动 | 绘制单个基因小提琴图面板 | 可视化 | 子模块 | ↑ StackedVlnPlot; ↓ — | Seurat, ggplot2 | — |
| `extract_max` | 6858 | 未动 | 提取 ggplot 纵轴最大值并取整 | 可视化 | 子模块 | ↑ StackedVlnPlot; ↓ — | ggplot2 | — |
| `barPlot` | 6885 | 未动 | 按组绘制平均表达量柱状图 | 可视化 | 主流程 | ↑ plotGeneExpression; ↓ barplot_internal | Seurat, reshape2, patchwork | — |
| `barplot_internal` | 6969 | 未动 | 绘制单基因分组柱状图 | 可视化 | 子模块 | ↑ barPlot, netVisual_barplot; ↓ — | ggplot2, ggpubr, scales | — |
| `visualizePie` | 7029 | 未动 | scatterpie 绘制细胞比例饼图 | 可视化 | 入口 | ↑ 用户; ↓ — | ggplot2, scatterpie(!), gtools | — |
| `spatialTopicPlot` | 7123 | 未动 | 按空间坐标绘制 topic 得分图 | 可视化 | 入口 | ↑ 用户; ↓ — | ggplot2, Matrix, RColorBrewer, scales, patchwork | — |
| `plotly_spatialLRpairPlot` | 7289 | 未动 | 绘制 LR 特征 plotly 空间图（读旧槽） | 可视化 | 主流程 | ↑ plotly_spatialLRpairPlot_shiny; ↓ plotFeatures | dplyr | — |
| `plotFeatures` | 7457 | 未动 | 绘制单特征 2D/3D 空间得分图 | 可视化 | 子模块 | ↑ plotly_spatialLRpairPlot; ↓ generate_RedsBlues_colormap, generate_custom_scenes3d | ggplot2, plotly, patchwork, dplyr, RColorBrewer, scales | — |
| `plotly_spatialLRpairPlot_shiny` | 7918 | 未动 | 启动交互式选点 Shiny 应用 | 可视化 | 入口 | ↑ 用户; ↓ plotly_spatialLRpairPlot, plotly_spatialDimPlot | shiny, plotly, stringr | — |
| `plotly_spatialDimPlot` | 8139 | 重构 | 绘制可高亮的交互式空间分组图 | 可视化 | 入口 | ↑ plotly_spatialLRpairPlot_shiny; ↓ generate_grid_nrows, generate_custom_scenes3d | plotly | — |
| `themeFeaturePlot` | 8333 | 未动 | 提供特征图通用 ggplot 主题 | 可视化 | 入口 | ↑ 用户; ↓ — | ggplot2 | — |
| `spatialLeePlot` | 8368 | 重构 | 组织权重矩阵并绘制 Lee 统计图 | 可视化 | 主流程 | ↑ 用户; ↓ plotStatistics_Lee, computeCellDistance | spdep(!), Matrix | — |
| `plotStatistics_Lee` | 8580 | 未动 | 计算并绘制 Lee 空间共现统计 | 可视化 | 子模块 | ↑ spatialLeePlot; ↓ scPalette | ComplexHeatmap, RColorBrewer, spdep(!), purrr, Matrix, ggplot2 | — |
| `spatialGiPlot` | 8841 | 未动 | 绘制 Gi/Gi* 热点图并写回 meta | 可视化 | 入口 | ↑ 用户; ↓ plotStatistics_Gi, generate_grid_nrows | ggplot2, plotly, patchwork, Matrix, dplyr, purrr | — |
| `plotStatistics_Gi` | 9087 | 未动 | 计算并渲染 Gi/Gi* 多模式图 | 可视化 | 子模块 | ↑ spatialGiPlot; ↓ generate_RedsBlues_colormap, generate_custom_colorscale | spdep(!), metR(!), plotly, ggplot2, dplyr | — |
| `doHeatmap` | 9471 | 未动 | 绘制带分组色条的表达热图 | 可视化 | 入口 | ↑ 用户; ↓ generate_RedsBlues_colormap, scPalette | Seurat, ggplot2, patchwork, RColorBrewer, scales, purrr, rlang | — |
| `generate_custom_colorscale` | 9852 | 未动 | 生成调色板或 plotly 色阶 | 可视化 | 子模块 | ↑ plotStatistics_Gi; ↓ — | RColorBrewer, scales, BiocGenerics | — |
| `generate_RedsBlues_colormap` | 9901 | 未动 | 生成红蓝发散色阶函数 | 可视化 | 子模块 | ↑ doHeatmap, plotFeatures; ↓ colorRamp3 | — | — |

### 5.5 `R/SpatialCellChat_class.R`（70 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `SparseChatArray` | 13 | 新增(11slot) | 构造三维稀疏通讯数组容器 | 基础设施 | 主流程 | ↑ modeling/analysis 建模; ↓ .SparseChatArray_dimnames, .new_SparseChatArray | — | — |
| `[.SparseChatArray` | 26 | 新增(11slot) | 按行列层下标切分数组 | 基础设施 | 主流程 | ↑ 建模与可视化取值; ↓ .SparseChatArray_resolve_index, .SparseChatArray_subset_layer | — | — |
| `[[.SparseChatArray` | 58 | 新增(11slot) | 取单层或层子集矩阵 | 基础设施 | 子模块 | ↑ .sc_compute_communication_field; ↓ .SparseChatArray_resolve_index, .new_SparseChatArray | — | — |
| `dim.SparseChatArray` | 85 | 新增(11slot) | 返回行列与层数 | 基础设施 | 子模块 | ↑ 切片与校验; ↓ — | — | — |
| `dimnames.SparseChatArray` | 90 | 新增(11slot) | 读取三维名称 | 基础设施 | 子模块 | ↑ 切片与打印; ↓ .SparseChatArray_dimnames | — | — |
| `dimnames<-.SparseChatArray` | 96 | 新增(11slot) | 重设三维名称 | 基础设施 | 子模块 | ↑ 用户; ↓ .SparseChatArray_dimnames, .new_SparseChatArray | — | — |
| `length.SparseChatArray` | 102 | 新增(11slot) | 返回层数 | 基础设施 | 子模块 | ↑ dim 与切片; ↓ — | — | — |
| `names.SparseChatArray` | 105 | 新增(11slot) | 返回层名 | 基础设施 | 子模块 | ↑ as.data.frame.SparseChatArray; ↓ dimnames.SparseChatArray | — | — |
| `names<-.SparseChatArray` | 108 | 新增(11slot) | 重设层名 | 基础设施 | 子模块 | ↑ 用户; ↓ .SparseChatArray_validate_dimname | — | — |
| `print.SparseChatArray` | 115 | 新增(11slot) | 打印维度与稀疏度 | 基础设施 | 子模块 | ↑ tests_dev; ↓ dim, dimnames | — | — |
| `marginSums` | 135 | 新增 | 边缘求和的 S3 泛型 | 基础设施 | 子模块 | ↑ tests_dev; ↓ marginSums.SparseChatArray | — | — |
| `marginSums.SparseChatArray` | 138 | 新增(11slot) | 按行、列或层汇总稀疏矩阵 | 基础设施 | 子模块 | ↑ marginSums; ↓ cpp_sum_layers | Matrix | — |
| `as.data.frame.SparseChatArray` | 167 | 新增(11slot) | 三维数组转长表 | 基础设施 | 子模块 | ↑ tests_dev; ↓ names, dim | — | — |
| `t.SparseChatArray` | 198 | 新增(11slot) | 逐层转置三维数组 | 基础设施 | 子模块 | ↑ tests_dev; ↓ .new_SparseChatArray | pbapply | — |
| `.new_SparseChatArray` | 204 | 新增(11slot) | 低层装配容器并清层内名称 | 基础设施 | 主流程 | ↑ SparseChatArray 族; ↓ .SparseChatArray_drop_matrix_dimnames | — | — |
| `.sc_map_sparse_layers` | 215 | 新增(11slot) | 对每层独立应用函数 | 基础设施 | 子模块 | ↑ computeCommunProb; ↓ .new_SparseChatArray | — | — |
| `.SparseChatArray_dimnames` | 221 | 新增 | 归一化三维名称列表 | 基础设施 | 子模块 | ↑ SparseChatArray, dimnames<-; ↓ .SparseChatArray_common_matrix_dimname | — | — |
| `.SparseChatArray_common_matrix_dimname` | 242 | 新增 | 取各层一致的行或列名 | 基础设施 | 子模块 | ↑ .SparseChatArray_dimnames; ↓ — | — | — |
| `.SparseChatArray_validate_dimname` | 255 | 新增 | 校验名称长度 | 基础设施 | 子模块 | ↑ .SparseChatArray_dimnames, names<-; ↓ — | — | — |
| `.SparseChatArray_drop_matrix_dimnames` | 263 | 新增 | 清空单层矩阵自身名称 | 基础设施 | 子模块 | ↑ .new_SparseChatArray; ↓ — | — | — |
| `.SparseChatArray_resolve_index` | 268 | 新增 | 字符下标转位置 | 基础设施 | 子模块 | ↑ [.SparseChatArray, [[.SparseChatArray; ↓ — | — | — |
| `.SparseChatArray_resolve_single_index` | 282 | 新增 | 解析并校验单元素下标 | 基础设施 | 子模块 | ↑ —; ↓ .SparseChatArray_resolve_index | — | 疑似 |
| `.SparseChatArray_subset_dimname` | 298 | 新增 | 同步切片后的名称 | 基础设施 | 子模块 | ↑ [.SparseChatArray; ↓ — | — | — |
| `.SparseChatArray_subset_layer` | 304 | 新增 | 对单层做二维切片 | 基础设施 | 子模块 | ↑ [.SparseChatArray; ↓ — | — | — |
| `.sc_matrix_like` | 339 | 新增 | 判断是否为矩阵类对象 | 基础设施 | 子模块 | ↑ 各校验与归一化函数; ↓ — | — | — |
| `.sc_validate_matrix` | 343 | 新增(11slot) | 校验矩阵维度与名称 | 基础设施 | 子模块 | ↑ .sc_validate_assay, .sc_validate_net; ↓ .sc_matrix_like | — | — |
| `.sc_validate_assay` | 362 | 新增(11slot) | 校验五个表达层结构 | 基础设施 | 子模块 | ↑ .sc_validate_object; ↓ .sc_validate_matrix | — | — |
| `.sc_validate_idents` | 402 | 新增(11slot) | 校验分组因子对齐细胞 | 基础设施 | 子模块 | ↑ .sc_validate_object; ↓ — | — | — |
| `.sc_validate_coordinate_system` | 419 | 新增(11slot) | 校验坐标系元数据 | 基础设施 | 子模块 | ↑ .sc_validate_images, .sc_validate_rasters; ↓ — | — | — |
| `.sc_validate_spatial_factors` | 444 | 新增(11slot) | 校验比例与容差一致性 | 基础设施 | 子模块 | ↑ .sc_validate_images, createSpatialCellChat; ↓ — | — | — |
| `.sc_validate_transform` | 474 | 新增(11slot) | 校验三乘三栅格变换矩阵 | 基础设施 | 子模块 | ↑ .sc_validate_rasters; ↓ — | — | — |
| `.sc_validate_rasters` | 490 | 新增(11slot) | 校验栅格记录集合 | 基础设施 | 子模块 | ↑ .sc_validate_images; ↓ .sc_validate_transform | — | — |
| `.sc_validate_fov` | 536 | 新增(11slot) | 校验 FOV 几何记录 | 基础设施 | 子模块 | ↑ .sc_validate_images; ↓ — | — | — |
| `.sc_validate_distance_cache` | 567 | 新增(11slot) | 校验距离缓存矩阵 | 基础设施 | 子模块 | ↑ .sc_validate_images; ↓ .sc_validate_matrix | — | — |
| `.sc_validate_images` | 581 | 新增(11slot) | 校验 images 槽全部内容 | 基础设施 | 子模块 | ↑ .sc_validate_object; ↓ .sc_validate_rasters, .sc_validate_fov | — | — |
| `.sc_validate_net` | 621 | 新增(11slot) | 校验 net 与 netP 结构维度 | 基础设施 | 子模块 | ↑ .sc_validate_object; ↓ .sc_validate_matrix | — | — |
| `.sc_validate_object` | 660 | 新增(11slot) | 汇总校验十一个槽 | 基础设施 | 主流程 | ↑ validateSpatialCellChat, setValidity; ↓ .sc_validate_assay, .sc_validate_images | — | — |
| `.sc_default_images` | 709 | 新增(11slot) | 生成空 images 槽骨架 | 基础设施 | 子模块 | ↑ createSpatialCellChat; ↓ — | — | — |
| `.sc_normalize_coordinate_system` | 721 | 新增(11slot) | 补全坐标系默认字段 | 基础设施 | 子模块 | ↑ createSpatialCellChat; ↓ — | — | — |
| `.sc_normalize_spatial_factors` | 743 | 新增(11slot) | 补全空间因子默认值 | 基础设施 | 子模块 | ↑ createSpatialCellChat; ↓ — | — | — |
| `.sc_normalize_transform` | 754 | 新增(11slot) | 规范化栅格变换列表 | 基础设施 | 子模块 | ↑ .sc_normalize_raster; ↓ .sc_matrix_like | — | — |
| `.sc_normalize_raster` | 774 | 新增(11slot) | 把图像转为标准栅格记录 | 基础设施 | 子模块 | ↑ extractSeuratRasterImage, readSpatialImage; ↓ .sc_normalize_transform | SeuratObject | — |
| `readSpatialImage` | 865 | 新增(11slot) | 读取 Visium 组织图像 | 基础设施 | 子模块 | ↑ tests_dev; ↓ .sc_normalize_raster | Seurat | — |
| `extractSeuratRasterImage` | 878 | 新增(11slot) | 从 Seurat 对象提取栅格 | 基础设施 | 子模块 | ↑ createSpatialCellChat; ↓ .sc_normalize_raster | SeuratObject | — |
| `.sc_extract_seurat_coordinates` | 904 | 新增(11slot) | 从 Seurat 提取组织坐标 | 基础设施 | 子模块 | ↑ createSpatialCellChat; ↓ — | Seurat, SeuratObject | — |
| `.sc_matrix_values` | 929 | 新增 | 取出矩阵的数值向量 | 基础设施 | 子模块 | ↑ createSpatialCellChat, computeCommunProb; ↓ — | — | — |
| `validateSpatialCellChat` | 938 | 新增(11slot) | 对外校验对象并报告错误 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_validate_object | — | — |
| `createSpatialCellChat` | 1022 | 迁移11slot | 由矩阵、Seurat 或 SCE 建对象 | 基础设施 | 入口 | ↑ 用户; ↓ .sc_default_images, extractSeuratRasterImage | Seurat, SingleCellExperiment(!), SummarizedExperiment(!) | — |
| `idents` | 1295 | 新增(11slot) | 读取细胞分组标签 | 基础设施 | 子模块 | ↑ level.SpatialCellChat; ↓ — | — | — |
| `setIdent` | 1314 | 迁移11slot | 按元数据列重设分组 | 基础设施 | 入口 | ↑ —; ↓ validObject | — | 是（重复定义，被 utilities.R:1122 覆盖） |
| `level.SpatialCellChat` | 1325 | 重构 | 彩色打印细胞分组水平 | 基础设施 | 入口 | ↑ 用户; ↓ idents | cli | — |
| `.sc_assert_spatial_cell_chat` | 1334 | 新增 | 断言对象类型 | 基础设施 | 主流程 | ↑ 各访问器; ↓ — | — | — |
| `.sc_validate_after_update` | 1340 | 新增 | 写操作后触发整体校验 | 基础设施 | 主流程 | ↑ 各 setter; ↓ validObject | — | — |
| `assay` | 1349 | 新增(11slot) | 读取表达层 | 基础设施 | 主流程 | ↑ modeling/spatial/utilities; ↓ .sc_assert_spatial_cell_chat | — | — |
| `assay<-` | 1360 | 新增(11slot) | 替换表达层 | 基础设施 | 主流程 | ↑ utilities 预处理; ↓ .sc_validate_after_update | — | — |
| `meta` | 1377 | 新增(11slot) | 读取细胞元数据 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `meta<-` | 1387 | 新增(11slot) | 替换细胞元数据 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_validate_after_update | — | — |
| `images` | 1397 | 新增(11slot) | 读取空间数据与缓存 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `images<-` | 1407 | 新增(11slot) | 替换空间数据与缓存 | 基础设施 | 主流程 | ↑ computeCommunProb; ↓ .sc_validate_after_update | — | — |
| `LR` | 1423 | 新增(11slot) | 读取配体受体元数据 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `LR<-` | 1433 | 新增(11slot) | 替换配体受体元数据 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_validate_after_update | — | — |
| `DB` | 1449 | 新增(11slot) | 读取通讯数据库 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `DB<-` | 1459 | 新增(11slot) | 替换通讯数据库 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_validate_after_update | — | — |
| `communication` | 1477 | 新增(11slot) | 按槽、分辨率、指标取结果 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `communication<-` | 1497 | 新增(11slot) | 写入通讯结果 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_validate_after_update | — | — |
| `params` | 1526 | 新增(11slot) | 读取分析参数 | 基础设施 | 主流程 | ↑ tests_dev; ↓ .sc_assert_spatial_cell_chat | — | — |
| `params<-` | 1537 | 新增(11slot) | 写入分析参数 | 基础设施 | 主流程 | ↑ normalizeData, preProcessing; ↓ .sc_validate_after_update | — | — |
| `misc` | 1554 | 未动 | misc 的 S3 泛型 | 基础设施 | 入口 | ↑ 用户; ↓ misc.SpatialCellChat | — | — |
| `misc.SpatialCellChat` | 1567 | 迁移11slot | 按名读取 misc 条目 | 基础设施 | 主流程 | ↑ misc; ↓ .sc_assert_spatial_cell_chat | — | — |
| `%\|\|%` | 1587 | 新增 | 空值取默认的运算符 | 基础设施 | 子模块 | ↑ 全文多处; ↓ — | — | — |

### 5.6 `R/utilities.R`（45 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `.cli` | 12 | 新增 | 按详细级别统一输出 cli 消息 | 基础设施 | 子模块 | ↑ 各模块消息调用; ↓ cli | cli | — |
| `.sc_progress_enabled` | 41 | 新增 | 判断当前是否启用进度条 | 基础设施 | 子模块 | ↑ my_future_lapply, computeCommunProb; ↓ — | — | — |
| `.log_operation` | 47 | 新增 | 向 misc$.log 追加一条操作记录 | 基础设施 | 子模块 | ↑ preProcessing, netAnalysis_computeCentrality; ↓ — | — | — |
| `savePlotWithCow` | 68 | 未动 | 用 cowplot 保存当前绘图 | 可视化 | 入口 | ↑ 用户; ↓ cowplot::save_plot | cowplot | — |
| `mergeChatIntoSeurat` | 123 | 未动 | 把通讯中心性结果并入 Seurat 对象 | 基础设施 | 入口 | ↑ 用户; ↓ CellChat2Seurat | Seurat, dplyr, stringr, Matrix | 疑似 |
| `RunSeurat` | 339 | 未动 | 串联 Seurat 标准化、降维与聚类 | 计算 | 入口 | ↑ 用户; ↓ Seurat | Seurat | — |
| `==.sparse3Darray` | 368 | 未动 | 逐字段比较两个 sparse3Darray | 基础设施 | 子模块 | ↑ —; ↓ — | — | 疑似 |
| `my_obj_size` | 394 | 未动 | 估算对象内存占用并返回可读串 | 基础设施 | 入口 | ↑ 用户; ↓ pryr | pryr (!) | 疑似 |
| `setpwd` | 406 | 未动 | 把工作目录设为当前脚本目录 | 基础设施 | 入口 | ↑ 用户; ↓ rstudioapi | rstudioapi (!) | 疑似 |
| `cli.symbol` | 422 | 未动 | 返回统一的状态符号 | 基础设施 | 入口 | ↑ modeling, visualization; ↓ cli | cli | — |
| `setEnvironment` | 450 | 重构 | 配置并行后端、进度器与详细级别 | 基础设施 | 入口 | ↑ 用户; ↓ future::plan, progressr::handlers | future, progressr, reticulate | — |
| `my_as_sparse3Darray` | 501 | 重构 | 快速构造旧式 sparse3Darray | 基础设施 | 子模块 | ↑ netVisual_CommunField; ↓ spatstat.sparse | purrr, progressr, spatstat.sparse, spatstat.utils | 疑似 |
| `my_future_sapply` | 630 | 重构 | 带进度条的并行 sapply | 基础设施 | 子模块 | ↑ netAnalysis_computeCentrality, computeExpr 系列; ↓ my_future_lapply | — | — |
| `my_future_lapply` | 674 | 重构 | 带进度与稳定随机数的并行 lapply | 基础设施 | 子模块 | ↑ my_future_sapply, computeCommunProb; ↓ future.apply | future.apply, progressr | — |
| `CellChat2Seurat` | 712 | 部分迁移 | 旧版对象转 Seurat 并可视化校验 | 基础设施 | 入口 | ↑ mergeChatIntoSeurat, identifyCellTopics; ↓ Seurat | Seurat, SeuratObject | 是 |
| `preProcessing` | 785 | 迁移11slot | 用 ALRA 填补表达并写回指定层 | 计算 | 入口 | ↑ 用户; ↓ ALRA::alra, assay<- | ALRA, Matrix | — |
| `.sc_normalize_matrix` | 883 | 新增 | 稀疏矩阵按列归一化并取对数 | 计算 | 子模块 | ↑ normalizeData; ↓ Matrix::colSums | Matrix | — |
| `normalizeData` | 932 | 迁移11slot | 把 raw 计数归一化写入 norm 层 | 计算 | 主流程 | ↑ createSpatialCellChat; ↓ .sc_normalize_matrix | Matrix | — |
| `scaleData` | 982 | 未动 | 按基因中心化缩放表达矩阵 | 计算 | 子模块 | ↑ —; ↓ Matrix | Matrix | 疑似 |
| `scaleMat` | 998 | 未动 | 按行或列标准化、归一化矩阵 | 计算 | 子模块 | ↑ identifyCellTopics, identifyCommunicationPatterns; ↓ — | — | — |
| `sketchData` | 1044 | 未动 | 用 geosketch 下采样细胞 | 计算 | 入口 | ↑ —; ↓ runPCA, reticulate | reticulate, Seurat | 疑似 |
| `addMeta` | 1095 | 未动 | 整体替换对象的 meta 表 | 基础设施 | 入口 | ↑ —; ↓ — | — | — |
| `setIdent` | 1122 | 部分迁移 | 设置细胞分组并按旧网络重排 | 基础设施 | 子模块 | ↑ updateClusterLabels; ↓ — | — | 疑似 |
| `updateClusterLabels` | 1172 | 未动 | 重命名细胞分组并重算下游网络 | 分析 | 入口 | ↑ —; ↓ setIdent, computeCommunProbPathway | plyr | 疑似 |
| `.sc_clean_gene_names` | 1204 | 新增 | 清洗并去空基因名向量 | 基础设施 | 子模块 | ↑ subsetData, identifyOverExpressedInteractions; ↓ — | — | — |
| `.sc_complex_subunits` | 1209 | 新增 | 解析复合物亚基成员表 | 基础设施 | 子模块 | ↑ identifyOverExpressedInteractions, computeCommunProb; ↓ .sc_clean_gene_names | — | — |
| `.sc_db_signaling_genes` | 1225 | 新增 | 汇总数据库涉及的信号基因 | 基础设施 | 子模块 | ↑ subsetData; ↓ .sc_complex_subunits | — | — |
| `.sc_lr_feature_state` | 1260 | 新增 | 判断配受体可用与入选状态 | 基础设施 | 子模块 | ↑ identifyOverExpressedInteractions; ↓ — | — | — |
| `subsetData` | 1286 | 迁移11slot | 抽取信号基因写入 signaling 层 | 基础设施 | 入口 | ↑ 用户, identifyOverExpressedGenes; ↓ .sc_db_signaling_genes, assay<- | — | — |
| `identifyOverExpressedGenes` | 1372 | 迁移11slot | 鉴定组间高表达或空间可变基因 | 分析 | 入口 | ↑ 用户; ↓ subsetData, makeGridSpatialCellChat | Seurat, MERINGUE | — |
| `identifyOverExpressedLigandReceptor` | 1656 | 部分迁移 | 筛选各组过表达配体与受体 | 分析 | 入口 | ↑ —; ↓ my_future_lapply | dplyr | 疑似 |
| `identifyOverExpressedInteractions` | 1735 | 迁移11slot | 筛出过表达的配体受体对 | 分析 | 入口 | ↑ 用户; ↓ .sc_lr_feature_state, .sc_complex_subunits | — | — |
| `projectData` | 1820 | 未动 | 把表达投影到蛋白互作网络 | 计算 | 入口 | ↑ —; ↓ projectAndRecombine | Matrix | 是 |
| `projectAndRecombine` | 1864 | 未动 | 在网络上投影后回填原维度 | 计算 | 子模块 | ↑ projectData; ↓ projectOnNetwork, randomWalkBySolve | — | 疑似 |
| `projectOnNetwork` | 1884 | 未动 | 把表达投影到网络基因空间 | 计算 | 子模块 | ↑ projectAndRecombine; ↓ — | — | 疑似 |
| `randomWalkBySolve` | 1903 | 未动 | 解线性方程完成随机游走投影 | 计算 | 子模块 | ↑ projectAndRecombine; ↓ l1NormalizeRows, l1NormalizeColumns | — | 疑似 |
| `l1NormalizeColumns` | 1923 | 未动 | 按列做 L1 归一化 | 计算 | 子模块 | ↑ randomWalkBySolve; ↓ Matrix | Matrix | 疑似 |
| `l1NormalizeRows` | 1934 | 未动 | 按行做 L1 归一化 | 计算 | 子模块 | ↑ randomWalkBySolve; ↓ Matrix | Matrix | 疑似 |
| `projectFromNetworkRecombine` | 1946 | 未动 | 把投影值回填到原矩阵 | 计算 | 子模块 | ↑ projectAndRecombine; ↓ — | — | 疑似 |
| `runPCA` | 1967 | 未动 | 主成分分析降维 | 计算 | 子模块 | ↑ sketchData; ↓ irlba | irlba | 疑似 |
| `runUMAP` | 2022 | 未动 | 调用 umap-learn 做降维 | 计算 | 入口 | ↑ identifyCommunicationPatterns; ↓ reticulate | reticulate | — |
| `.error_if_no_Seurat` | 2072 | 未动 | 检查 Seurat 是否安装 | 基础设施 | 子模块 | ↑ —; ↓ Seurat | Seurat | 疑似 |
| `colorRamp3` | 2098 | 未动 | 按断点插值生成配色函数 | 可视化 | 入口 | ↑ netVisual_heatmap, generate_RedsBlues_colormap; ↓ .get_color | colorspace | — |
| `.restrict_in` | 2177 | 未动 | 把数值截断到给定区间 | 基础设施 | 子模块 | ↑ .get_color; ↓ — | — | — |
| `.get_color` | 2188 | 未动 | 在色彩空间内插值单点颜色 | 可视化 | 子模块 | ↑ colorRamp3; ↓ .restrict_in, colorspace | colorspace | — |

### 5.7 `R/database.R`（9 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `showDatabaseCategory` | 10 | 重构 | 数据库互作分类饼图组合；与 visualization.R 重复定义，被覆盖 | 可视化 | 入口 | ↑ 用户; ↓ pieChart | cowplot, dplyr | 疑似 |
| `pieChart` | 49 | 重构 | 分组占比饼图；与 visualization.R 重复定义，被覆盖 | 可视化 | 入口 | ↑ showDatabaseCategory; ↓ — | ggplot2, dplyr, scales, ggrepel | 疑似 |
| `searchPair` | 83 | 未动 | 按通路或配体检索配体受体对 | 分析 | 入口 | ↑ netVisual_aggregate, netVisual_chord_gene; ↓ — | dplyr, future.apply | — |
| `subsetDB` | 137 | 未动 | 按注释或关键词筛选数据库 | 基础设施 | 入口 | ↑ 用户; ↓ — | — | — |
| `extractGene` | 182 | 未动 | 提取数据库涉及的全部基因 | 基础设施 | 入口 | ↑ 用户; ↓ checkGeneSymbol, extractGeneSubset | dplyr | — |
| `extractGeneSubset` | 221 | 未动 | 将复合体展开为其亚基基因 | 基础设施 | 子模块 | ↑ extractGene, extractGeneSubsetFromPair; ↓ — | dplyr | — |
| `extractGeneSubsetFromPair` | 243 | 未动 | 从配体受体对提取信号基因 | 基础设施 | 入口 | ↑ 用户; ↓ extractGeneSubset | — | — |
| `checkGeneSymbol` | 278 | 未动 | 检查基因是否为官方符号 | 基础设施 | 子模块 | ↑ extractGene; ↓ — | — | — |
| `updateCellChatDB` | 311 | 未动 | 合并自定义配体受体库生成新数据库 | 基础设施 | 入口 | ↑ 用户; ↓ — | dplyr | — |

### 5.8 `R/visualization_widget.R`（9 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `.sc_widget_check_scalar` | 8 | 新增 | 校验单个有限数值参数 | 基础设施 | 子模块 | ↑ .sc_prepare_commun_flow_widget, netVisual_CommunFlowWidget; ↓ — | — | — |
| `.sc_widget_color_map` | 28 | 新增 | 生成细胞身份到颜色的映射 | 可视化 | 子模块 | ↑ .sc_prepare_commun_flow_widget; ↓ scPalette | — | — |
| `.sc_widget_assets` | 57 | 新增 | 定位 widget 前端 JS/CSS 资源 | 基础设施 | 子模块 | ↑ netVisual_CommunFlowWidget; ↓ — | — | — |
| `.sc_widget_dependency` | 84 | 新增 | 构造 htmlwidget 资源依赖对象 | 基础设施 | 子模块 | ↑ netVisual_CommunFlowWidget; ↓ — | htmltools | — |
| `.sc_widget_regular_grid` | 94 | 新增 | 抽稀矢量场为规则网格 | 可视化 | 子模块 | ↑ .sc_widget_interpolate; ↓ — | — | — |
| `.sc_widget_zero_grid` | 114 | 新增 | 零流场时生成退化空网格 | 可视化 | 子模块 | ↑ .sc_widget_interpolate; ↓ — | sf | — |
| `.sc_widget_interpolate` | 129 | 新增 | Barnes 插值生成平滑流场网格 | 可视化 | 子模块 | ↑ .sc_prepare_commun_flow_widget; ↓ .sc_widget_regular_grid, .sc_widget_zero_grid | oce(!), sf | — |
| `.sc_prepare_commun_flow_widget` | 185 | 新增 | 组装通讯流场 widget 数据载荷 | 可视化 | 主流程 | ↑ netVisual_CommunFlowWidget; ↓ .sc_widget_interpolate, .sc_resolve_grid_size | sf, Matrix | — |
| `netVisual_CommunFlowWidget` | 415 | 新增 | 渲染通讯流场粒子动画 widget | 可视化 | 入口 | ↑ 用户; ↓ .sc_prepare_commun_flow_widget, htmlwidgets::createWidget | htmlwidgets, htmltools | — |

### 5.9 `R/RcppExports.R`（5 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `ComputeSNN` | 4 | 未动 | Rcpp 导出包装：ComputeSNN | 基础设施 | 子模块 | ↑ buildSNN; ↓ — | — | — |
| `cpp_sum_layers` | 8 | 新增 | Rcpp 导出包装：cpp_sum_layers | 基础设施 | 子模块 | ↑ .sc_sum_layers, marginSums.SparseChatArray; ↓ — | — | — |
| `cpp_prob_layer` | 12 | 新增 | Rcpp 导出包装：cpp_prob_layer | 基础设施 | 子模块 | ↑ computeCommunProb; ↓ — | — | — |
| `cpp_group_avg_obs` | 16 | 新增 | Rcpp 导出包装：cpp_group_avg_obs | 基础设施 | 子模块 | ↑ computeAvgCommunProb; ↓ — | — | — |
| `cpp_group_avg_perm` | 20 | 新增 | Rcpp 导出包装：cpp_group_avg_perm | 基础设施 | 子模块 | ↑ computeAvgCommunProb; ↓ — | — | — |

### 5.10 `R/test_helpers.R`（4 个）

| 函数 | 行 | 标签 | 功能 | 流程 | 层级 | 上下游 | 依赖 | 死代码 |
|---|---|---|---|---|---|---|---|---|
| `start_session` | 17 | 新增 | 启动持久后台测试 R 会话 | 基础设施 | 子模块 | ↑ 开发者; ↓ — | callr(!) | — |
| `source_chat` | 29 | 新增 | 向测试会话加载全部 R 源码 | 基础设施 | 子模块 | ↑ 开发者; ↓ — | Seurat, Matrix, igraph | — |
| `init_test_data` | 57 | 新增 | 生成微型测试对象并缓存 | 基础设施 | 子模块 | ↑ 开发者; ↓ createSpatialCellChat | Matrix | — |
| `list_checkpoints` | 105 | 新增 | 列出测试缓存检查点文件 | 基础设施 | 子模块 | ↑ 开发者; ↓ — | — | — |

## 6. 读这份清单的方式
- **P0 只一条**：先让包能装。补回四个被导出的函数（或从 NAMESPACE 撤掉 export），否则 `R CMD INSTALL`、`load_all()`、任何 `library(SpatialCellChat)` 都不可用；现有测试全走 `source()`，覆盖不到这一层。
- **迁移顺序建议**：`R/visualization.R` 是最大缺口（62 个函数里 47 个未动，且 `netVisual_*` 主家族仍读 `net$prob`/`net$pval`）。这些函数在新对象上不是「结果不对」，而是取不到字段，适合按 `netVisual_aggregate` → `netVisual_circle`/`chord` → `heatmap`/`barplot` 的顺序成组迁移。
- **兼容层收尾**：`my_as_sparse3Darray` 仅剩 1 处调用（`netVisual_CommunField`），把它和 `computeAvgCommunProb_LR_Avg/Sum`、`modifySparse3Darray` 一起清掉。
- **死代码不要顺手删**：`是` 的 4 个可以删；`疑似` 的 40 个里有导出函数（如 `computeAveExpr`、`communicationDistPlot`），删除属 API 变更，建议先标 `@deprecated` 再删。
- **依赖缺口**：15 个未装包全部落在 Suggests 的可选分支，不影响主链路；若要跑 `identifyCellTopics` / `suggestK` / `computeEnrichmentScore` 这类高级分析，需单独 `renv::install()` 并重新快照。
