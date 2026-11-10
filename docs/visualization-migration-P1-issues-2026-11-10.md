# 可视化迁移问题记录：P1 首轮验收

日期：2026-11-10

## 验收范围

对当前缓存的模拟对象执行 P1 函数：`spatialVisual_scoring()`、`spatialTopicPlot()`、`spatialGiPlot()`、`spatialLeePlot()`。

## 已观察问题

### P1-001：空间 topic 图的生产者与当前验收前置条件未闭合

`spatialTopicPlot()` 消费的是 `slot(object, slot.name)$topic[[pattern]]$cell`，其对应生产者是 `identifyCellTopics()`。此前把问题归因于 `identifyCommunicationPatterns()` 不准确：后者产生的是群级 `slot$pattern`，供 `netAnalysis_river()` / `netAnalysis_dot()` 使用，不能直接作为 cell × topic 的空间 loading。

当前缓存的三个验收对象都只有 `net` / `netP` 的 group centrality，没有 `centr.cell` 和 `topic`。因此直接调用 `spatialTopicPlot()` 会在 `NULL` topic 上失败。并且当前环境没有可选包 `RcppML`，`identifyCellTopics()` 的唯一 NMF 路径也不能作为验收前置。

影响范围：`spatialTopicPlot()`、`netVisual_TopicComposition()`、`netVisual_TopicSignaling()` 以及 `identifyTopicCommunities()` 的 topic 依赖链；`identifyCommunicationPatterns()` 属于另一条群级模式链，不在本问题的直接影响范围内。

处理方向：保持 `identifyCellTopics()` 读取已计算的 `centr.cell` 及原 `RcppML::nmf()` 的方法，测试时先准备中心性和 topic。已在本地 renv 安装 RcppML；尚需把验收依赖及版本固化，不能用另一种 NMF 算法替代原方法。

### P1-002：Gi 图单指标选择发生维度降维

`spatialGiPlot()` 构造 `df_sub <- df_score[c("x_cent", "y_cent", measure[i])]`。当 `df_score` 是 matrix 时，这是一维线性索引，而不是选择三列，结果成为向量；随后 `plotStatistics_Gi()` 把它当二维数据框设置三列名称，报 `attempt to set 'colnames' on an object with less than two dimensions`。

影响范围：`spatialGiPlot()` 的单指标调用，包括单 signaling、单 measure，以及 normal/scoring/density/contour/3d 的公共计算入口。

处理方向：选择 score 列时显式 `drop = FALSE`，并保持 `plotStatistics_Gi()` 的三列输入契约。

### P1-003：Lee/Gi 统计路径直接依赖当前环境缺失的 `spdep`

`spatialLeePlot()` 的 `Itself` 分支和 `plotStatistics_Lee()` 直接调用 `spdep::mat2listw()`、`spdep::lee()`；`plotStatistics_Gi()` 直接调用 `spdep::knearneigh()`、`spdep::knn2nb()`、`spdep::include.self()`、`spdep::nb2listw()`、`spdep::localG()`。当前项目 `DESCRIPTION` 将 `spdep` 放在 `Suggests`，当前 renv R 4.4.3 环境中没有安装它，因此默认 Lee/Gi 图无法完成三对象验收。

影响范围：`spatialLeePlot()`、`plotStatistics_Lee()`、`spatialGiPlot()`、`plotStatistics_Gi()` 的统计计算分支；绘图层本身并非根因。

处理方向：在本地 renv 安装 `spdep` 并继续使用原调用；尚需固化验收依赖。不得用未经严格等价验证的本地公式替代原统计实现。

### P1-004：Gi 全常数输入导致统计量未定义后色阶报错

模拟对象的部分单 pathway 在 `do.binary = TRUE` 时全部细胞得分均为 1。`spdep::localG()` 返回未定义统计量，现有配色函数对全非有限值取范围后报错。应保留未定义值并用 NA 颜色展示，不将其替换为 0 或伪造显著性。

### P1-005：Lee 热图仅有单一有限值时色阶报错

模拟对象用自身类别比较、`lw.type = "Itself"` 和默认 cutoff 时，只保留对角线上的 1；`ComplexHeatmap` 使用颜色向量自动构造色阶时没有两个不同断点而报错。应只扩展绘图色阶断点，保留原 Lee 数值、默认 cutoff 和矩阵方向。

### P1-006：源码直载的验收脚本缺少 colorspace 命名空间导入

Gi 的连续配色走 `colorRamp3()`；正式包的 `NAMESPACE` 导入 `colorspace::RGB` 等符号，但当前 runner 直接 `source()` 源文件，未加载 `colorspace`，因此报 `object 'RGB' not found`。这是验收脚本环境不完整，不应改动颜色算法。修复 runner 的加载环境后需重新保存并检查 PDF。

### P1-007：全量 Slide-seq 的中心性前置超过 future 默认上限

29,536 个细胞对象读取后，`netAnalysis_computeCentrality()` 的 future 闭包大小约 6.58 GiB，超过默认 500 MiB 上限，导致 scoring、topic、Gi 三图共同失败。`make_slideseq_chat()` 重建对象时已设置 `future.globals.maxSize = 20 * 1024^3`，但命中现有 RDS 缓存时提前返回，验收进程没有该设置。runner 应对全量 Slide-seq 显式设置相同上限并记录实际计算和绘图结果；提高限制不代表解决内存使用效率。

### P1-008：Gi 坐标顺序与统计结果对应

`spatialGiPlot()` 原来按位置拼接坐标和 centrality，也按位置写回 metadata；坐标行重排时会将细胞分数和统计结果分配给错误细胞，影响普通图、3D 图及返回对象。应核验坐标、centrality、identities、metadata 的细胞名并按名对齐，保留 Gi 计算和参数；增加同一对象重排坐标前后的结果一致性测试。已在不含等距邻居的模拟数据上验证重排一致性，并在最终代码上重新绘制三个对象的 Gi PDF。

### P1-009：Gi 和 Lee 测试中的邻接警告

最终代码上 Visium Gi 的 kNN 产生 2 个连通子图；Lee 的当前测试参数采用 `group.by1 = group.by2 = idents` 和 `lw.type = "Itself"`，Visium 提示 905 个子图，全量 Slide-seq 提示 29,536 个子图。模拟对象也有类似提示。它们没有阻止 PDF 绘制，但不能仅凭文件存在判断空间权重或结果的生物学合理性；需结合具体权重构造和真实分析目标核查，不能为了消除警告改变统计定义。

## 当前状态

- 三套对象的 P1 五个入口各 5/5 完成实际 PDF 绘制：`tests_dev/figures/P1/{simulated,visium,slideseq}/`；明细见 `tests_dev/figures/summary-*-P1.csv`。
- 两个真实 pathway 用于计算 cell centrality 和 outgoing topic，Slide-seq 使用全部 29,536 个细胞；这不构成全部 pathway 或参数组合的穷举验收。
- `tests_dev/test-spatial-p1-migration.R` 覆盖单指标 Gi、未定义 Gi、Lee 单值热图、topic 名称对齐与前缀。实际统计仍由 `spdep` 和 `RcppML` 完成。
- 四个本轮测试依赖已加入 `renv.lock`；项目的依赖扫描排除 Suggests，`renv::status()` 因此仍报告这些包已记录但未被扫描为使用，尚未完成独立环境恢复验证。

