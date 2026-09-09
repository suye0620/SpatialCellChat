# SpatialCellChat 可视化代码审计报告

> 2026-08-09 · 针对 `spatialXXXPlot` 系列函数与 plotly 动图实现的系统审计
> 关联文档：[REFACTORING_PLAN.md](REFACTORING_PLAN.md)、[DATA_STRUCTURE.md](DATA_STRUCTURE.md)
> 本报告为重构 Phase 5（可视化）提供输入；本轮已落地距离计算、图像缓存迁移、plotly 共享布局辅助和 plotly 空间散点图矩阵适配。仍有依赖旧 slot 的 feature/LR 绘图函数，列为迁移阻塞而非隐藏兼容层。

---

## 1. 审计范围与方法

### 1.1 审计对象

| 文件 | 内容 | 规模 |
|---|---|---|
| `R/visualization.R` | spatialXXXPlot 系列 + 全部 plotly 路径 + 主题/色图辅助 | 9786 行 |
| `R/spatial.R` | `computeCellDistance`、`createCellCellContactMatrixFrom_dspatial` | — |
| `R/utilities.R` | `colorRamp3`（RedsBlues 色图底层） | — |
| `R/analysis.R` | `calculate_density`（Gi density 方法依赖） | — |

覆盖函数：

- `spatialCCCDistPlot`（L21）、`spatialDimPlotPoints`（L114）、`spatialDimPlot`（L252）、`spatialFeaturePlot`（L360）、`spatialVisual_scoring`（L647）、`spatialTopicPlot`（L6879）
- `spatialLeePlot`（L8226）、`plotStatistics_Lee`（L8429）、`spatialGiPlot`（L8690）、`plotStatistics_Gi`（L8936）
- `plotly_spatialLRpairPlot`（L7045）、`plotFeatures`（L7213）、`plotly_spatialLRpairPlot_shiny`（L7674）、`plotly_spatialDimPlot`（L7895）
- `generate_custom_colorscale`（L9701）、`generate_RedsBlues_colormap`（L9750）

### 1.2 验证环境与方法

- R 4.5.3 / ggplot2 4.0.3 / plotly 4.12.0 / spdep 1.x / ks 1.15.3 / circlize（均为实际安装版本）
- 方法：静态阅读 + 独立命名空间加载包内函数做运行时验证；审计脚本不修改包文件，也不在内存中注入缺失 helper。
- 所有"实测"结论均有运行输出佐证；可复现脚本：`tests_dev/audit_visualization_findings.R`
- 3D stack 原型演示：`tests_dev/audit_3dstack_prototype.html`

---

## 2. 崩溃级 Bug（当前状态）

### 2.1 `colorRamp3` 的颜色插值 helper 已存在

此前审计曾记录 `.get_color` 缺失。当前源码已在 `R/utilities.R:2032` 提供包内 `.get_color`，`colorRamp3()` 可直接运行，不再依赖 `circlize:::.get_color` 这样的内部符号。

运行证据：`tests_dev/audit_visualization_findings.R` 的 `colorRamp3 [package source]` 返回 `SUCCESS character`；`plotStatistics_Gi` 的 `normal`、`3d`、`density` 分支也均已成功构建。

维护约束：`.get_color` 是内部函数，不应导出；若后续拆分 `utilities.R`，必须把它与 `colorRamp3` 放在同一内部模块，避免重新引入裸符号丢失。

### 2.2 plotly 共享布局辅助已补回

当前定义位于 `R/visualization.R:1-34`：

- `custom_legend`：统一 legend 配置；
- `generate_grid_nrows()`：校验正整数并生成 subplot grid；
- `generate_custom_scenes3d()`：统一 scene、camera 和 aspect ratio。

运行证据：`plotFeatures` 的 2d/3d、`plotly_spatialDimPlot` 的 2d/3d 均在审计脚本中以包内定义直接返回 `plotly` 对象，不再依赖内存补丁。

### 2.3 `calculate_density(method = "wkde")` 已有本地实现

当前 `R/analysis.R:4199` 定义 `wkde2d()`，`calculate_density()` 的 wkde 分支调用该本地函数。此前关于“ks 中不存在 `wkde2d`”的结论只适用于外部 ks API，不适用于当前仓库源码。

运行证据：`calculate_density wkde [package source]` 返回 `SUCCESS numeric`。后续仍应补充带重复坐标、全零权重和单点输入的边界测试。


---

## 3. 指标计算原理审计

### 3.1 `spatialCCCDistPlot`：接触距离分布的零值与缓存路径

`createCellCellContactMatrixFrom_dspatial()` 会把接触矩阵对角设为 1，但新的 `computeCellDistance()` 明确排除自点与零距离边，`spatialCCCDistPlot()` 直接从稀疏距离矩阵的非零槽位生成密度数据。因此当前绘图路径不会把对角自距混入密度。

每个非零空间边仍以有向矩阵形式存储两次（`i → j` 与 `j → i`）。对密度图这只是等权重复；若未来从缓存计算总边数、均值或分位数，必须使用上三角去重或明确标注统计量为 directed-edge scale。

缓存位于 `object@images$.distance`，且现在携带 `.parameters`（interaction/contact range、ratio、tol）。绘图只有在缓存参数与本次请求一致时才复用；旧缓存或阈值变化会按 canonical coordinates 重算。旧的 `images$result.computeCellDistance` 不再作为绘图兼容路径。

### 3.2 `computeCellDistance`：从稠密距离矩阵改为精确范围查询

旧实现先用 `Rfast::Dist()` 构造 n×n 稠密矩阵，再逐列阈值化；大规模 Slide-seq 数据会产生不必要的 O(n²) 中间内存。

当前实现使用 `BiocNeighbors::queryNeighbors()` + `VptreeParam()` 做精确半径查询，只展开命中邻域并直接构造 `dgCMatrix`。同时修复了以下语义错误：支持 2D/3D 坐标；`ratio = NULL` 表示坐标已在工作单位；`tol = NULL` 按 0 处理；`ratio`、范围和坐标均做有限值校验。

证据：`tests_dev/test-computeCellDistance.R` 的 8 项回归检查通过；`tests_dev/benchmark-computeCellDistance.R` 以 dense reference 对照，1000 点基准本机测得树查询约 19.5× 加速（计时受机器影响），结果矩阵逐元素一致。该实现是树结构范围查询，不应在文档中误称为 FNN 的固定 k-NN。

### 3.3 `spatialGiPlot` / `plotStatistics_Gi`：统计正确，语义有坑

**统计正确性**：

- `spdep::knearneigh(k = n)`（不含自身）→ `knn2nb` → `localG(score, nb2listw(nb, style = "B"))` = Getis-Ord **Gi**，标准定义；
- `include.self()` 后 = Getis-Ord **Gi\***（局部和含自身）；
- 二进制权重 `style = "B"` 符合 Gi 统计定义；`longlat = FALSE`、`use_kd_tree = TRUE` 合理。

**问题清单**：

1. **默认 `do.binary = T`**：`centr[centr > 0] <- 1` 先把中心性得分二值化再算 Gi——这是"热点二值化"语义，文档未说明，且改变统计量含义（对 0/1 指标的 Gi 与对连续值的 Gi 不可直接比较）；
2. **`spot_highlight` 语义反转**（`visualization.R:8770` 附近）：高亮组得到 `F`，其余 `T`，变量名与含义相反；`sapply(spot_labels, function(x) if (x %in% group.highlight) F else T)` 应向量化为 `!(spot_labels %in% group.highlight)`；
3. **中心性聚合**：`apply(centr, c(1, 2), sum)` 跨全部通路求和得"incoming/outgoing 总分"——语义合理但需文档化（与 CellChat 的 per-pathway 展示不同）；
4. `class(gg[[1]])[[1]] %in% c('ggplot','gg',...)` 脆弱类型判断，改 `inherits(gg[[1]], "ggplot")`；
5. `scene.name = paste0("scene", ifelse(i == 1, "", i))` 与 `layout(scene=..., scene2=...)` 命名经核对恰好匹配 plotly.js 约定，无问题（勿"顺手修"）。

### 3.4 `spatialLeePlot` / `plotStatistics_Lee`：正确但默认行为需斟酌

- `spdep::lee` 对 k₁×k₂ 组别两两调用（每次 O(nnz(lw)) 空间滞后），正确但可向量化；
- `lw.type = "Itself"`：权重为恒等矩阵（`sparseMatrix(i=1:n, j=1:n, x=1)`）→ Lee's L 退化为 **Pearson 相关**。合理，但需文档化；
- `lw.type = "contact.range"` / `"interaction.range"`：分别用接触邻接（style "B"）与交互概率矩阵（style "W"），正确；
- **默认 `cutoff = TRUE, cutoff.Lee = 0`**：把所有 ≤ 0 的 Lee's L 置 NA——负空间关联是 Lee's L 的正常输出（负共现），默认隐藏会误导解读，建议默认 `cutoff = FALSE` 或文档说明；
- 数值型 `x.labels`/`y.labels` 输入会隐式强制 `type = "dot"`——隐式行为，需文档化；
- `%||%` 有定义（`SpatialCellChat_class.R:1130`），无问题。

### 3.5 `apply(prob > 0, 3, sum)`：稠密膨胀（必须迁移）

`spatialFeaturePlot`（L430 附近）与 `plotly_spatialLRpairPlot`（L7110 附近）的 cell-level 富集过滤：

```r
prob.cell <- object@net$prob.cell[, , pairLR.use.name, drop = FALSE]
prob.sum <- apply(prob.cell > 0, 3, sum)   # nC × nC × N 稠密逻辑数组
```

nC = 5 万时单层即 2.5 GB。这正是重构中 `SparseChatArray` / `marginSums`（Rcpp `cpp_sum_layers`）要消灭的模式——**所有此类调用统一迁移到 `marginSums`**。

---

## 4. ggplot2 语法规范审计（ggplot2 4.0.3 实测）

| 位置 | 问题 | 实测结果 | 建议 |
|---|---|---|---|
| `spatialCCCDistPlot`（当前 L58） | `geom_density(aes(x = x, y = ..density..))` | 已改为 `after_stat(density)` | 保持新语法，空/单值数据走 `geom_blank()` |
| `plotly_spatialDimPlot`（当前 L7967） | 直接给 matrix 使用 `$` 并触发 dplyr `mutate` | 已改为 data.frame，最终 schema 2d/3d 均可构建 | 保持坐标与 cell id 在同一 data.frame |
| `doHeatmap` | `aes_string` | deprecated | 改 `aes()` / `.data` 代词 |
| `extract_max` | `ggplot_build(p)$layout$panel_scales_y...` | 4.0.3 可用 | 这是 ggplot 内部 API，改用 `layer_scales()` 或捕获后缓存 |
| `doHeatmap` | `pbuild$layout$panel_params[[1]]$x.major` / `$x$get_breaks()` | 版本耦合 | 避免直接读 panel internals |
| `spatialDimPlotPoints` | `labels` 外部向量参与 aes | 已并入 `plot_data` | 保持 data.frame 显式映射 |
| `spatialDimPlot` | 直接修改 plot 对象的 layers/coordinates | 当前 raster wrapper 已集中处理 | 后续继续减少 ggplot 对象内部修改 |
| `spatialGiPlot` | `class(gg[[1]])[[1]]` 类型判断 | 可用但脆弱 | 改 `inherits(gg[[1]], "ggplot")` |
| `plotFeatures` normal | `coord_fixed()` 与 `theme(aspect.ratio = 1)` 并存 | 冗余 | 二选一 |
| 多处 | feature/topic/scoring 各自内联 theme | 重复且参数漂移 | 收敛到 `themeFeaturePlot()` |
| 多处 | x↔y 交换 + `scale_y_reverse()` 重复 | 语义隐含 | 统一坐标预处理 helper |

**良好实践确认**：`element_line(linewidth = ...)`（ggplot2 ≥ 3.4 规范，未用废弃的 `size`）；`scale_y_reverse() + coord_fixed()` 坐标约定在全系列一致；`na.value` 处理一致。

**坐标约定问题**：剩余绘图函数仍重复“x↔y 交换 + `scale_y_reverse()`”；raster 分支还引入 `spot.coordinates` 的图像轴方向。应把 analysis 坐标、raster pixel 坐标和绘图坐标的变换显式写入同一 accessor，避免未来图像来源切换时漏改。

---

## 5. plotly 动图审计与 3D Stack 方案

### 5.1 现状问题

1. **“3d”不是真 stack**：`plotFeatures` 3d 与 `plotly_spatialDimPlot` 3d 都把点放在固定 `z = three.dim.z` 平面；这适合空间位置的 3D 外壳，不等于多 feature/pathway 沿 z 轴堆叠。`plotly_spatialLRpairPlot` 仍按 feature 单独出图。
2. **trace 索引硬编码**：`plotFeatures`、`plotStatistics_Gi` 仍通过 `pl$x$data[[N]]` 修改 colorbar；highlight 是否为空会改变 trace 数量，必须改成按 trace 名称或构建时显式传 `marker$colorbar`。
3. **空 trace**：`plotly_spatialDimPlot` 已在无 highlight 时跳过空 trace；`plotFeatures` 与 `plotStatistics_Gi` 3d 仍需要同样的条件添加。
4. **标题与 hover 耦合**：`plotFeatures` 仍从 ggplot `labs()` 对象内部取标题；多个 plotly 函数直接拼接 HTML。新 `plotly_spatialDimPlot` 已对 cell id/type 做基本 HTML 转义，其余路径仍需统一。
5. **Shiny observer 泄漏**：`plotly_spatialLRpairPlot_shiny` 在 `plotly_selecting` observer 内注册 `input$update` observer，每次框选都会累积监听器；应改为 `reactiveVal` 保存选中 cell 集合。
6. **legend group/title 混用**：`legendgroup` 与 `legendgrouptitle` 当前仍共享同一个字符串参数，建议拆成机器 key 与显示标题两个参数。

### 5.2 3D Stack 正解（已验证）

**目标**：把 N 个 feature/通路画成 N 层水平"切片"沿 z 轴堆叠，共享一个颜色标尺，z 轴刻度直接标注 feature 名。

**关键设计点**（此前做不出来的常见原因）：

1. **单 trace 单 colorscale**：把 N 层 `rbind` 成一个长表，`color = ~value` 只映射一次 → 一根 colorbar。若每层一个 trace 各自 `color=`，plotly 会生成 N 根 colorbar 且色域各自归一；
2. **全局归一**：`vmin/vmax` 取全数据范围，0 锚定（`if (vmin >= 0) vmin <- 0`），保证层间可比；
3. **`aspectratio` 压低 z 轴**（如 `x:y:z = 1:1:0.35`），防止堆叠后呈长条形；
4. **z 轴 ticktext = feature 名**：`tickmode="array"` + `tickvals = seq(0, (nf-1)*z.space, by=z.space)` + `ticktext = feat`；
5. 层间间距 `z.axis.space` 可调；hover 显示层名 + cell_id + 值。

**原型实现**（`plotly_build` 实测 SUCCESS：3 层 × 300 点，z 轴标注 L1/L2/L3）：

```r
spatial3Dstack <- function(xy, layers, color.heatmap = "RdBu", point.size = 3, z.axis.space = 1) {
  feat <- colnames(layers)
  nf <- length(feat)
  vmin <- min(layers, na.rm = TRUE); vmax <- max(layers, na.rm = TRUE)
  if (vmin >= 0) vmin <- 0            # 0 锚定，保证层间可比
  if (vmax <= 0) vmax <- 0
  cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, color.heatmap)))(99)
  long <- do.call(rbind, lapply(seq_len(nf), function(k) data.frame(
    x = xy$x, y = xy$y, z = (k - 1) * z.axis.space,
    value = layers[[k]], layer = feat[k], cell_id = rownames(xy)
  )))
  plotly::plot_ly(long, x = ~x, y = ~y, z = ~z, color = ~value, colors = cols,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = point.size),
                  hoverinfo = "text",
                  text = ~paste("<b>", layer, "</b><br>", cell_id, "<br>score:", round(value, 3))) %>%
    plotly::layout(scene = list(
      xaxis = list(title = "x", showgrid = TRUE, zeroline = FALSE),
      yaxis = list(title = "y", showgrid = TRUE, zeroline = FALSE),
      zaxis = list(title = "", tickmode = "array",
                   tickvals = seq(0, (nf - 1) * z.axis.space, by = z.axis.space),
                   ticktext = feat, showgrid = TRUE, zeroline = FALSE),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.35),
      camera = list(eye = list(x = 1.6, y = 1.6, z = 0.9))),
      showlegend = FALSE)
}
```

**扩展方向**：`plotly::highlight` 框选联动（替代 shiny 轮询式 observer）；层间用值域归一化 z 间距；加层间半透明连接线展示"信号流"；大 spot 数时切换 `scattergl`/`toWebGL`。

---

## 6. 重构优化方向与优先级

| 优先级 | 项目 | 工作量 | 说明 |
|---|---|---|---|
| **P0（已闭环）** | 颜色插值、plotly 共享 helper、wkde 本地实现 | 小 | 审计脚本已验证 `colorRamp3`、Gi normal/3d/density、plotFeatures 2d/3d 和 plotly spatialDimPlot 2d/3d |
| **P1** | 收敛重复代码 | 中 | colormap、theme、坐标轴变换和 plotly trace 构建仍有重复 |
| **P2（进行中）** | 适配新 schema | 中/大 | `spatialDimPlot`、`spatialCCCDistPlot` 和 plotly spatialDimPlot 已适配；`spatialFeaturePlot`、`plotly_spatialLRpairPlot` 仍访问旧 `@data*`/旧通信结构 |
| **P3（部分完成）** | 指标与绘图细节 | 小/中 | 距离缓存、零值和树查询已修复；Gi `spot_highlight`、`do.binary`、Lee `cutoff`、trace 索引仍需决策或迁移 |
| **P4** | plotly 现代化 | 中 | 3D stack 单 trace 原型已验证；生产 API、动态框选、colorbar 显式配置和 Shiny observer 修复仍待落地 |

### 与现有重构计划的衔接

- Phase 5 应先完成 `assay()` / `communication()` / `LR()` 访问层迁移，再统一 theme 和坐标变换；不要在旧 slot 上继续加兼容分支。
- P1 收敛重复代码应和 `visualization.R` 拆分同步完成，避免迁移后再次复制。
- Gi/Lee 默认值属于接口决策：`do.binary`、`cutoff.Lee`、`spot_highlight` 必须在 API 文档和回归测试中固定。

---

## 7. 复现与验证指引

**证据脚本**：`tests_dev/audit_visualization_findings.R`。

- 直接 source 包内文件，无颜色/KDE/helper 内存补丁；
- `colorRamp3`、wkde、plotFeatures 2d/3d/normal、Gi normal/3d/density、plotly spatialDimPlot 2d/3d、spatialCCCDistPlot 和 3D stack 原型均成功；
- `spatialFeaturePlot` 与 `plotly_spatialLRpairPlot` 对最终 11-slot 对象仍分别报旧 `data.signaling` / `data` slot 错误，这是已确认的 schema 迁移阻塞；
- `tests_dev/test-SpatialCellChat-spatial-image.R`：全部图像/坐标/缓存检查通过，比例饼图因未安装 `scatterpie` 明确跳过；新增阈值变化时 CCC 距离缓存失配重算检查；
- 真实 psoriasis Seurat fixture 导入通过：`33538 × 905` norm 矩阵、1 个 `main` raster、`600 × 595 × 3` RGB 图像、905 个 spot 坐标，`validateSpatialCellChat()` 为 TRUE；
- `tests_dev/test-computeCellDistance.R`：8/8 通过；`tests_dev/benchmark-computeCellDistance.R`：树查询与 dense reference 逐元素一致，n=1000 的本机计时约 19.5×；n=250/500 的小样本计时受系统分辨率与启动成本影响（本次约 0.5×/NA），内存结论以不物化 n×n 矩阵为准。

**3D stack 演示**：`tests_dev/audit_3dstack_prototype.html`（3 层 × 300 点）。

**完整包检查**：已运行 `devtools::check('.', document = FALSE, args = c('--no-manual'))`；构建成功，但在依赖检查阶段因环境缺少 required packages（`ALRA`、`ggalluvial`、`svglite`、`MERINGUE`、`ComplexHeatmap`、`sna`、`forcats`、`ggpubr`、`ggnetwork`、`ggnewscale`）而停止。该结果不能作为源代码回归失败；本报告以独立 source 的 targeted checks 作为当前验证证据。feature/LR 绘图仍需最终 schema 的真实通信 fixture，不应把旧 slot 读取伪装成兼容成功。

---

## 8. 遗留与未验证项

- 真实 psoriasis Seurat 导入与 raster 自动提取已用 `tests_dev/test_data/visium_human_psoriasis_PP1_correctedLocations202401_annotated(1).RData` 验证；仍缺少更大规模坐标距离压力测试；
- `spatialFeaturePlot`、`plotly_spatialLRpairPlot` 和依赖旧通信结构的 Gi/Lee 外层函数仍未完成最终 schema 迁移；
- `plotStatistics_Gi` 的 trace 后处理、Shiny observer 累积和 HTML hover 转义仍需单独回归；
- 3D stack 当前是验证原型，不是已导出的生产函数；生产 API 需要确定 layer 输入、共享色域、z 轴单位和大数据渲染后再加入 NAMESPACE。
