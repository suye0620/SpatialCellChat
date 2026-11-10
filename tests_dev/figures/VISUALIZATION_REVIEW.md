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
- `P3/<dataset>/`：场图、分解图以及通讯流 HTML。`netVisual_CommunFlowWidget-summary.csv` 标记 HTML 和浏览器测试状态。全量 Slide-seq Grid/Flow 已输出，但显示效果待修复，另有预览 PNG。
- `P4/<dataset>/`：维度图、单基因表达图与单指标 Gi 的交互 HTML；每个目录的 `summary.csv` 记录细胞数与状态。`BROWSER-NOT-TESTED` 表示未完成对应图的交互验收，不代表文件不能打开。

## 未通过与覆盖边界

- Visium `netVisual_chord_gene`：88 LR 的 sector 标签重叠，汇总记为 `VISUAL-REVIEW`。
- Slide-seq `netVisual_CommunFieldGrid`、`netVisual_CommunFlow`：默认背景遮挡向量，Flow 出现仅一步积分警告，记为 `VISUAL-REVIEW`。全量数据的参数对照改善了 Flow 显示并消除了该警告，但 Grid 箭头仍稀疏；两项对照也需视觉检查，不算通过。
- 比较/merged 图及 pairwise embedding 缺少当前结构的真实对象构建入口和完整契约，尚未迁移验收。
- `netAnalysis_contribution()` 和 `rankNet()` 未注册进三对象绘图清单，仍读取旧概率字段，单对象与比较分支均未完成迁移验收，详见 P2-017。
- 原始旧缓存对象缺少距离及通信概率来源快照，无法安全使用通信距离模式，需从原始输入重新推断。
- HTML 浏览器截图仅为抽查，未完成全部图、动态 widget 和 Shiny 的浏览器交互测试。Shiny 模拟服务端会话测试不替代浏览器测试。
- topic 示例只覆盖两条 pathway；Seurat TopicSignaling 已用单个合成 Seurat 对象及显式 interaction 表实绘，但未覆盖三套真实 Seurat 数据。其他参数、三维堆叠及多样本比较并未穷举。
- `renv.lock` 已记录本轮额外所需的 Suggests 包及其依赖，但尚未在独立环境完成恢复；`renv::status()` 将新增 Suggests 包标为未被当前扫描设置使用。

## 本次复测

在 R 4.4.3 下再次运行 `tests_dev/test-netvisual-group-migration.R`、`tests_dev/test-viz-p3-migration.R` 和 `tests_dev/test-viz-p4-migration.R`，三项退出码均为 0。P3 测试仍报告部分点被移除、UMAP 初始嵌入尺度提示以及依赖包弃用提示；群网络测试也有绘图库弃用提示。Shiny 检查仅覆盖服务端模拟，不替代浏览器交互。复测环境启动时报告 locale 设置失败和 renv out-of-sync，不能据此认定环境已独立复现。

详细位置、原因、修复与待确认事项见 `docs/visualization-migration-P0-issues-2026-11-10.md` 至 `docs/visualization-migration-P4-issues-2026-11-10.md`。
