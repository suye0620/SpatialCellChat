# 可视化迁移验收索引

日期：2026-11-10。此文件列出已生成的测试产物和已知限制，不代表全部迁移已完成。

## 测试对象

- `simulated`：模拟对象，20 个细胞。
- `visium`：Visium，905 个位置。
- `slideseq`：完整 Slide-seq，29,536 个细胞，未抽样。

## 产物位置

- `P0/<dataset>/`：距离图与空间维度图基线。原旧缓存对象的通信距离图会被来源校验拒绝，详见 `summary-<dataset>-P0.csv`；从原始输入重新推断的结果另存 `summary-<dataset>-P0-spatialC2CDistPlot-fresh.csv`，不会把旧缓存报错误称为通过。
- `P1/<dataset>/`：静态点图和统计图，Gi 在坐标对齐修复后重新绘制。早期 P1 全集汇总在 `summary-<dataset>-P1.csv`，Gi 最新单项汇总为 `summary-<dataset>-P1-spatialGiPlot.csv`。
- `P2/<dataset>/`：单对象网络图，最终代码重跑汇总为 `summary-<dataset>-P2.csv`。总入口导出采用只有一个 LR 的 pathway，以控制页面大小；individual 图采用第一条 pathway 的全部 LR，Visium 为 88 个，Slide-seq 为 16 个。
- `P3/<dataset>/`：场图、分解图以及通讯流 HTML。`netVisual_CommunFlowWidget-summary.csv` 标记 HTML 生成状态；浏览器测试见 `P3/commun-flow-widget-browser-summary-2026-11-10.csv`。全量 Slide-seq Grid/Flow 已输出，但显示效果待修复，另有预览 PNG。
- `P4/<dataset>/`：维度图、3 层 spatial stack、单基因表达图与单指标 Gi 的交互 HTML；每个目录的 `summary.csv` 记录生成状态，`P4/browser-summary-2026-11-10.csv` 记录 15 个 CDP smoke case 均通过。`P4/shiny-browser-summary-2026-11-10.csv` 记录 simulated Shiny 的真实鼠标框选、替换和 Reset 已通过；Visium 和 Slide-seq 浏览器交互结果见 `P4/shiny-browser-summary-visium-2026-11-10.csv`、`P4/shiny-browser-summary-slideseq-2026-11-10.csv`。

## 未通过与覆盖边界

- Visium `netVisual_chord_gene`：88 LR 的 sector 标签重叠，汇总记为 `VISUAL-REVIEW`。
- Slide-seq `netVisual_CommunFieldGrid`、`netVisual_CommunFlow`：默认背景遮挡向量，Flow 出现仅一步积分警告，记为 `VISUAL-REVIEW`。全量数据的参数对照改善了 Flow 显示并消除了该警告，但 Grid 箭头仍稀疏；两项对照也需视觉检查，不算通过。
- 比较/merged 图已由 `mergeSpatialCellChat()` 提供真实 producer；`tests_dev/test-merge-spatialcellchat.R` 覆盖 rankNet、pairwise embedding、差异图、bubble/heatmap、signaling-role 和 merged LR 提取。空间 merged 图仍因 dataset-local 坐标而明确拒绝，不能当作共享组织坐标比较。
- `netAnalysis_contribution()` 和 `rankNet()` 的单对象路径已迁移并完成三对象出图；`rankNet()` comparison 及 merged LR 提取已用真实 `mergeSpatialCellChat()` producer 验收。`netAnalysis_contribution()` 当前仍是单对象入口，merged 对象不应被误当作共享坐标或统一网络输入。
- 原始旧缓存对象缺少距离及通信概率来源快照，无法安全使用通信距离模式，需从原始输入重新推断。
- HTML 已完成 15 个页面的 CDP 初始化/覆盖/重布局 smoke，包含新 3D stack；simulated Shiny 已完成真实鼠标框选、连续替换和 Reset 浏览器验收。三对象的动态通讯流 widget 和 Shiny 已完成部分真实浏览器交互测试，但未穷举；Slide-seq 两次更新约 11.6 秒和 9.2 秒，性能仍需研究。
- topic 示例只覆盖两条 pathway；Seurat TopicSignaling 已用单个合成 Seurat 对象及显式 interaction 表实绘，但未覆盖三套真实 Seurat 数据。其他参数、三维堆叠及多样本比较并未穷举。
- `renv.lock` 已记录本轮额外所需的 Suggests 包及其依赖，但尚未在独立环境完成恢复；`renv::status()` 将新增 Suggests 包标为未被当前扫描设置使用。

## 本次复测

在 R 4.4.3 下再次运行 `tests_dev/test-netvisual-group-migration.R`、`tests_dev/test-viz-p3-migration.R`、`tests_dev/test-viz-p4-migration.R` 和 `tests_dev/test-analysis-network-migration.R`，四项退出码均为 0。P2 新增的 `netAnalysis_contribution()`、`rankNet-netP` 和 `rankNet-net` 均在三套对象上写出 PDF。P3 测试仍报告部分点被移除、UMAP 初始嵌入尺度提示以及依赖包弃用提示；群网络测试也有绘图库弃用提示。另已在浏览器中检查三对象的 Shiny 选择与重置，具体范围见上述汇总。复测环境启动时报告 locale 设置失败和 renv out-of-sync，不能据此认定环境已独立复现。

详细位置、原因、修复与待确认事项见 `docs/visualization-migration-P0-issues-2026-11-10.md` 至 `docs/visualization-migration-P4-issues-2026-11-10.md`。


## 网格聚合性能复测

Grid/Flow 的网格向量聚合已去除逐网格 `sf` 子集开销，并通过空网格、缺失值及两种聚合方法的等价性检查；完整 Slide-seq 两幅参数对照的优化前后渲染结果哈希一致。性能证据和比较范围见 `docs/visualization-migration-P3-issues-2026-11-10.md` 的 P3-014。该优化不意味着场图可读性问题已经解决。


## 零向量显示修正

`netVisual_CommunFieldGrid()` 的零长度箭头现在透明，避免绘图库把它们画成组织外的规则细点；非零箭头仍使用原位置、大小及用户指定透明度。三对象的结果分别在 `P3/<dataset>/zero-vector-display-check/`；完整 Slide-seq 是关闭箭头跳格、缩小背景点后的对照，并非默认参数结果。Grid 与 Flow 的整体可读性仍保留人工复核状态。


## Merged 比较实绘与热图边界

`merged/producer-comparison/summary.csv` 和同目录 10 张 PDF 使用真实 merge producer，但输入是模拟对象及其修改副本，不是独立生物样本。差异热图修复了正负变化求和抵消导致的误判和 `remove.isolate` 误删；边界图与三对象单图汇总分别见 `merged/heatmap-edgecases/` 和 `P2/heatmap-cancellation-regression-2026-11-10.csv`。图像仍需人工验收。


`merged/diff-isolate-edgecases/` 另含差异网络正负抵消、无孤立节点和单节点自环的实绘回归；全零网络按测试预期返回明确提示。


## 网络图孤立群回归

`P2/network-isolate-regression/` 含 17 个网络绘图边界情形；`P2/network-isolate-three-fixtures-2026-11-10.csv` 及各数据集的 `network-isolate-regression/` 含三对象 aggregate/individual 各三种布局的 18 张 PDF。每个数据集只选择一个 pathway 和一条 LR，全部使用 `remove.isolate = TRUE`，不能代表所有参数与信号。模拟数据产生绘图库弃用警告，其余本次图形调用没有捕获到警告。Slide-seq 空间网络点密、标签较小，hierarchy 部分长标签靠近中心，均需人工复核。问题定位及测试范围见 P2-021。


## 空间网络中心性节点大小

`P2/centrality-spatial-three-fixtures-2026-11-10.csv` 索引三对象各一条通路的 incoming/outgoing 结果；图像写出前将中心性群维顺序反转，并检查节点大小仍按群名对应。具体原因和验证范围见 P2-022，图像仍需人工验收。


## 差异网络删群属性回归

`merged/diff-isolate-edgecases/` 新增删除中间群的自定义与默认配色图，测试同时核对群名、权重、标签、形状和带自环的实际邻接矩阵；比较输入是模拟对象及其修改副本。详见 P2-023。

本索引列出的早期产物中，按主汇总表、补充入口和交互页面整理的文件已列入 `artifact-inventory-2026-11-10.csv` 并纳入 Git；工作区仍有其他未跟踪产物，不属于该清单。清单记录生成时的状态与 SHA-256，不证明历史图是用当前最新版代码重新绘制，也不代表视觉或生物学验收通过。


## 测试产物完整性与版本控制

`artifact-inventory-2026-11-10.csv` 收录 137 个 PDF、18 个交互 HTML、30 张预览或浏览器截图、147 个 HTML 直接依赖和 30 份相应依赖许可证；其中可视化结论沿用原汇总状态。`tests_dev/audit-visualization-artifacts.py --check-tracked` 校验 PDF 头尾、本地 HTML 依赖是否存在，以及清单文件的 Git 索引字节是否与工作区一致。该检查不验证 PDF 页面内容、所有外部请求、浏览器交互或生物学结论。完整对象和图像清单、各自生成与复核边界仍以各阶段汇总和问题记录为准。

2026-11-11：从提交 `bba1038` 导出到独立临时目录后重新运行产物审计，346 项文件通过结构与本地依赖检查，重生成的清单与提交版一致；该步骤没有重新运行 R 绘图或浏览器。


## 2026-11-11 通讯场显示诊断

`P3/slideseq/display-diagnostics-2026-11-11/` 包含完整 Slide-seq 默认及参数对照、计数和警告；默认跳格后 22 个非零箭头，不跳格 540 个，但后者非零箭头长度中位数仍约 0.075 mm。`P3/streamline-cancellation-2026-11-11/` 使用合成常量场复现 metR 非零对角向量被当作零向量而丢失流线，见 P3-017；这段是诊断当时的历史状态；后续局部兼容修复及其限制见 P3-020、P3-021，不能从 PDF 写出状态推定视觉或生物学验收通过。

P3 局部修复实验见 `P3/streamline-patch-experiment-2026-11-11/`：合成精确对角场恢复流线，两个控制场的绘图层逐项一致；该实验当时没有修改生产入口或已安装的 metR；后续已引入受版本保护的生产兼容层，见 P3-020，仍有依赖内部实现的维护边界。

## 2026-11-11 当前生产代码的独立复测

以生产提交 `449c181` 为基准，P0–P4 重跑结果及浏览器检查集中在 `recheck-449c181/`，没有覆盖本索引中的历史结果。`environment.txt` 记录生产代码版本和环境，`fixture-sha256.csv` 记录四个本地缓存或检查点的哈希；P0 的模拟与 Visium 使用 fresh fixture，Slide-seq 距离图使用经一致性检查的预过滤检查点。P1–P4 大多使用相应对象缓存，不能据此声称全部对象在本次重新推断。

测试脚本的隔离输出、待视觉复核的明确状态和浏览器测试边界见 `docs/visualization-migration-recheck-2026-11-11.md`。本批 `PASS` 仅表示指定调用成功出图或浏览器检查通过，仍需人工视觉和生物学验收；不代表所有函数参数、下游链或独立样本比较均已通过。

## 安装包复测

2026-11-11：当前生产提交 `449c181` 的安装包在临时库中完成安装；安装后流线兼容测试和真实模拟对象的 geometry、heatmap、表达、3D 堆叠 smoke 均通过。对应产物和日志位于 `recheck-449c181/installed-smoke/` 及其相邻日志。安装测试曾发现测试脚本错误使用未导出的 `assay()`，已改正并重跑；这不是生产 API 的改动。

## Source package 检查

补充：已修正 `.Rbuildignore`，排除 renv、tests_dev、调试 `.tmp_*` 和本地编译对象；干净 source tarball 的 `R CMD check` 在关闭强制 Suggests 后退出码为 0。检查仍有 2 个 WARNING 和 6 个 NOTE，主要是既有依赖声明、全局变量声明和 Rd 签名差异；本机缺少部分 Suggests 包，因此没有声称完整 CRAN 环境检查通过。

## 2026-11-11 迁移遗漏复核

当前源码复跑通过 `test-remaining-entrypoints.R`、`test-remaining-migration.R`、`test-analysis-network-migration.R`、`test-merge-spatialcellchat.R`、`test-merged-zero-comparison.R`、`test-viz-p3-migration.R`、`test-viz-p4-migration.R` 和 `test-plotly-spatial-stack.R`。其中新增 rank-one/rank-two `suggestK1()` 检查发现并修复了 Matrix 长度为 1 的 `diag()` 尺寸歧义；问题、原因和修复边界见 `docs/remaining-migration-issues-2026-11-10.md` P5-018。

P4 的 15 个页面浏览器检查记录在 `recheck-449c181/P4/browser-summary-2026-11-10.csv`，15/15 通过，覆盖全细胞覆盖、3D scene、stack 三层/单 colorbar 和 `Plotly.relayout`。同一生产提交的隔离重跑产物与已有 HTML 重复，未另行纳入版本控制。这仍是浏览器 smoke，不是人工视觉或所有浏览器交互验收。

merged 结构由真实 Visium 自身对照再次通过 `mergeSpatialCellChat()` producer 生成，并成功绘制 10 个 comparison consumers；摘要在 `merged/visium-self-control-2026-11-11/summary.csv`。该对照是同一样本的技术零差异，不是独立样本比较。完整独立样本比较仍受 P5-014 的基因轴和共同群注释契约阻塞。
