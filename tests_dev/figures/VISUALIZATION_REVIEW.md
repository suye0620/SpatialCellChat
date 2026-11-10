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
