# Agent Note: main 分支分析与可视化 11-slot 迁移调研

Status: proposed
Governance: v1
Date: 2026-10-18
Decision type: migration
Scope: main branch R/analysis.R, R/visualization.R, R/visualization_widget.R, related tests and documentation
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

截至 2026-10-18，`main` 分支位于提交 `63d68fe`，最终 11-slot
`SpatialCellChat` schema 和核心通信计算已经投入使用，但分析层和可视化层仍处于
混合状态：一部分函数读取新结构，另一部分函数仍直接读取已废弃的旧 slot。

已完成或已有新 schema 路径的范围包括：

- `netAnalysis_computeCentrality()` 和 `subsetCommunication()` 已读取
  `net`/`netP` 的 `$cell`、`$group` 结构，并由 `tests_dev/test-centrality-subset.R`
  覆盖；
- `spatialDimPlot()`、`spatialCCCDistPlot()`、`computeCommunField()`、
  `netVisual_CommunField*()`、`netVisual_CommunFlow()` 和
  `netVisual_CommunFlowWidget()` 已使用 `images`、`idents` 或
  `SparseChatArray` 路径；
  这不表示场图全链已迁移：2026-10-20 语义核对发现基础 `netVisual_CommunField()`
  仍读 `slot$field` 并调用 `my_as_sparse3Darray()`，而生产者已写 `slot$cell$field`；
  `ZoomIn` 委托同一入口。见 [可视化语义卡 V-013](../../../docs/VISUALIZATION_SEMANTICS.md#9-稳定问题表证据消费者与未决边界)。
- `R/visualization.R` 中的 `spatialFeaturePlot()`、
  `plotly_spatialLRpairPlot()`、`netVisual_aggregate()` 和
  `netVisual_individual()` 仍访问 `object@data`、`object@data.signaling`、
  `object@net$prob` 或 `object@net$prob.cell`（见 `R/visualization.R:437-607`、
  `R/visualization.R:1440-1740`、`R/visualization.R:7289-7550`）；
- `R/analysis.R` 中的 `netAnalysis_contribution()`、`rankNet()`、相似性/模式
  分析以及 signaling-role 系列仍存在旧式 `object@options$mode`、裸
  `object@net`/`object@netP` 结构假设（见 `R/analysis.R:255-1990`、
  `R/analysis.R:2629-3280`）。

因此当前风险不是单个绘图 bug，而是最终 schema 的迁移边界没有在分析与可视化
消费者中闭合。现有审计也明确将 `spatialFeaturePlot()` 和
`plotly_spatialLRpairPlot()` 列为最终 schema 的迁移阻塞，见
`docs/VISUALIZATION_AUDIT.md` 第 3.5、6、7、8 节。

## Proposal

以最终 11-slot schema 作为 `main` 分支分析与可视化的唯一目标，按“共享访问契约
先行、函数分波迁移、每波独立验证”的顺序完成迁移；不在旧函数中继续叠加
`if (new schema)` 兼容分支。

后续语义调研将“共享契约先行”具体化为
[下游通讯分析与可视化语义契约](2026-10-20-downstream-communication-semantics-contract.md)：
在本批迁移前先核定科学问题、观察单位、约简/分解、视觉编码和统计边界，
区分存储迁移、文档/实现不一致与科学/API 更正。该子提案不替代本 Note 的
消费者范围和全量迁移验收，二者当前均保持 `proposed`。

[PP1 调用方增补](../../../docs/PP1_DOWNSTREAM_REFACTORING_PLAN.md) 覆盖真实脚本的
输入与 Visium 混合比例、Gi→metadata→Lee、topic→community→Seurat 及图/表导出。
它将语义门禁细化为 B0–B7 与链级验收，不替代下列四个优先入口或本 Note 的完整范围；
其中源码候选映射、已观察失败、局部运行成功与未执行场景分别标明。

具体顺序为：

1. 固定并复用 `assay()`、`communication()`、`LR()`、`DB()`、`images()`、
   `meta()` 和 `idents()` 的读取契约，统一 cell/group/pathway/LR 层的取数方式；
2. 优先迁移四个直接阻塞最终 schema 的可视化函数：
   `spatialFeaturePlot()`、`plotly_spatialLRpairPlot()`、
   `netVisual_aggregate()`、`netVisual_individual()`；
3. 再迁移 `netAnalysis_contribution()`、`rankNet()`、网络相似性/模式分析和
   signaling-role 消费者；
4. 对每一波补充最终 schema fixture、边界测试和文档，并在所有消费者迁移完成后
   删除对应旧 slot 读取路径。

拒绝“为每个函数单独保留旧/新两套分支”的方案，因为它会复制数据访问和坐标
   语义，延长废弃 schema 生命周期，并可能把旧 slot 读取误报为兼容成功。

## Constraints and invariants

- `SpatialCellChat` 保持最终 11 个顶层 slot；`.mode` 和 `.datatype` 位于
  `misc`，不得重新引入顶层 `options`、`data` 或 `data.signaling` 兼容 slot。
- `net` 保存 LR-level 结果，`netP` 保存 pathway-level 结果；两者的通信结果
  通过 `$cell`/`$group` 组织，概率和 p-value 层遵循 `SparseChatArray` 契约。
- `net` 的 group-level 显著性 gating 继续由其 p-value 语义控制；
  `netP$group$prob` 已在 pathway 生产阶段门控（有 group pval 时），
  `netP$cell$prob` 沿用 cell 结果按 pathway 求和，不新增 p-value 检验。
  消费者不得假设 netP 有 p-value 层；LR 实时聚合图也不能无条件替换为 netP 缓存。
- 空间分析坐标来自 `images$coordinates`；raster、距离缓存和坐标变换必须保持
  `docs/DATA_STRUCTURE.md` 与现有 `spatialDimPlot()` 的坐标语义。
- 迁移不得把 cell × cell × layer 的稠密数组重新物化；层级聚合应优先使用
  `SparseChatArray` 和 `marginSums()`。
- 已有函数的显著性过滤、返回列名、网络层名称、highlight 空值行为和错误边界
  必须由回归测试固定；若语义有意改变，必须先建立新的 Agent Note 或更新本 Note。
- 在代码、调用方、测试和文档完成并验证之前，本 Note 保持 `proposed`，不能作为
  “迁移已完成”的证据。

## Alternatives considered

1. **每个函数内增加新旧 schema 分支**：实现初期改动较小，但会复制访问逻辑、
   坐标变换和阈值语义；旧路径会继续被调用，难以证明最终 schema 的完整性，拒绝。
2. **先统一视觉主题和 plotly 外观，再迁移数据结构**：可以快速改善外观，但无法
   在旧数据访问仍存在时验证真实 API 契约，且主题重构会在迁移后重复修改，拒绝。
3. **一次性重写整个 `analysis.R`/`visualization.R`**：能减少过渡状态，但文件规模
   和函数耦合使回归范围过大，难以定位行为变化；采用按共享契约和消费者分波迁移。

## Consumer impact

**Schema 与生产者**

- `R/SpatialCellChat_class.R`：11-slot 类、校验和访问器；
- `R/modeling.R`：`net$cell`、`net$group`、`netP$cell`、`netP$group` 的通信结果
  生产者；
- `R/spatial.R`：空间坐标、距离和通信场相关结果；
- `R/utilities.R`：对象更新、并行和日志辅助。

**分析消费者**

- 已迁移基线：`netAnalysis_computeCentrality()`、
  `subsetCommunication()`/`subsetCommunication_internal()`；
- 待迁移重点：`netAnalysis_contribution()`、`rankNet()`、
  `identifyCommunicationPatterns()`、`computeNetSimilarity*()`、
  `netAnalysis_signalingRole_*()`、`netAnalysis_signalingChanges_scatter()`
  和 `netMappingDEG()` 的相关通信读取路径。

**可视化消费者**

- 已有新路径：`spatialDimPlot()`、`spatialCCCDistPlot()`、
  `plotly_spatialDimPlot()`、`netVisual_CommunField()`、
  `netVisual_CommunFieldGrid()`、`netVisual_CommunFlow()` 和
  `netVisual_CommunFlowWidget()`；
- 待迁移重点：`spatialFeaturePlot()`、`plotly_spatialLRpairPlot()`、
  `netVisual_aggregate()`、`netVisual_individual()`；
- 后续需复核：`plotGeneExpression()`、Gi/Lee 外层绘图、plotly trace/colorbar
  后处理和 Shiny selection observer 生命周期。

**测试、文档和记录**

- `tests_dev/test-centrality-subset.R`；
- `tests_dev/audit_visualization_findings.R`；
- `tests_dev/test-SpatialCellChat-spatial-image.R`；
- `tests_dev/test-grid-visualization.R`；
- `tests_dev/test-commun-flow-widget.R`；
- `tests_dev/test-netpathway-aggregate.R`、`test-filterCommunication.R`、
  `test-computeAvgCommunProb.R`、`test-computeCommunProb.R`；
- `docs/DATA_STRUCTURE.md`、`docs/REFACTORING_PLAN.md`、
  `docs/VISUALIZATION_AUDIT.md`；
- 相关既有记录：
  `.agents/notes/proposed/2026-08-15-spatial-plot-contract-migration.md`、
  `.agents/notes/proposed/2026-08-15-visualization-regression-baseline.md`、
  `.agents/notes/proposed/2026-08-15-spatial-statistics-and-plotly-hardening.md`、
  `.agents/notes/implemented/2026-10-10-centrality-subset-migration.md`。

本次调研未发现需要修改的配置文件或生成文件；NAMESPACE/API 导出只有在具体
函数签名发生变化时才需要同步更新。

## Consequences

- 迁移完成后，`main` 的分析与可视化路径将针对一个明确 schema 运行，旧
  `CellChat` 对象不会被隐式伪装成兼容对象，而应通过明确错误或显式升级流程处理。
- 共享访问契约和 `SparseChatArray` 聚合可以避免 `nC × nC × layer` 的不必要稠密
  膨胀，继续保留当前大规模空间数据的内存策略。
- 分波迁移会暂时保留部分函数不可用状态，但每一波可以用小型 fixture 隔离验证，
  比一次性重写更容易定位数值、列顺序、坐标和交互行为回归。
- `visualization.R` 的重复主题、坐标变换、plotly trace 和 hover 逻辑仍需在
  schema 迁移后收敛；这属于后续维护义务，不应通过兼容分支掩盖。
- 本次调研环境未提供 `Rscript`，因此本 Note 记录的是源代码审计和仓库已有的
  测试证据；不能把本次会话描述为重新执行了 R 回归测试。

## Evidence

- Git 基线：`main` 与 `origin/main` 同步，当前提交为 `63d68fe`
  （`Migrate centrality and communication subset to 11-slot schema`）；迁移前置
  提交包括 `dd0ed6a`、`f35cf34` 和 `bad5d33`。
- `R/analysis.R:83-224`：中心性分析已使用 `slot$group$prob`、
  `slot$group$pval` 和 `slot$cell$prob`；`R/analysis.R:2292-2530`：通信子集已
  使用新 schema 和 pathway join 规则。
- `R/visualization.R:56-435`：空间距离、坐标、raster 和 spatial dimension
  绘图的新路径；`R/visualization.R:5477-6548`：通信场与 flow 绘图路径。
- `R/visualization_widget.R:203-380`、`:415-490`：通信场 widget 使用
  cell-level `SparseChatArray`、`images$coordinates` 和明确参数校验。
- `R/visualization.R:437-607`、`:1440-1740`、`:7289-7550`：仍存在旧
  `@data*`/`net$prob`/`net$prob.cell` 访问，是当前迁移阻塞的直接证据。
- `R/analysis.R:870`、`:1532-1597`、`:1877-1983`、`:2629-3280`：仍存在
  `@options$mode`、旧 `net`/`netP` 结构或依赖旧字段的分析消费者。
- `docs/VISUALIZATION_AUDIT.md` 第 3.5、5、6、7、8 节：记录稠密膨胀风险、
  plotly trace/3D stack 问题、已适配函数、明确阻塞和现有验证边界。
- `docs/REFACTORING_PLAN.md` Phase 4-6：将下游分析、空间和可视化迁移列为
  后续阶段，并要求访问层先于业务函数迁移。
- 历史验证记录：`.agents/notes/implemented/2026-10-10-centrality-subset-migration.md`
  记录 `test-centrality-subset.R` 和 centrality/subset e2e 通过；
  `tests_dev/audit_visualization_findings.R` 记录 targeted visualization checks。
- 当前环境验证缺口：`Rscript --version` 无法执行，终端提示 `Rscript` 未被识别；
  因此本次没有重新运行 R 测试或 `devtools::check()`。

## Acceptance criteria

前置门禁：每批实现先完成关联语义卡与受影响争议的裁定；已知缺陷不自动作为
回归测试应保持的行为，具体更正另行记录。语义文档完成不能代替以下运行验收。

1. 对四个优先可视化函数建立最终 schema fixture；
   `spatialFeaturePlot()`、`plotly_spatialLRpairPlot()`、
   `netVisual_aggregate()` 和 `netVisual_individual()` 在合法 11-slot 对象上
   成功返回预期 ggplot/plotly/recorded plot 类型，并覆盖单层、多层、全零、NA、
   空 highlight 和 threshold 边界。
2. 上述函数及后续迁移的分析函数不再直接访问旧的 `@data`、
   `@data.signaling`、`net$prob`、`net$prob.cell` 或 `@options$mode`；静态扫描和
   运行时 fixture 都必须通过。允许外部 Seurat helper 的 `meta.data`/assay
   访问，但不得混入 `SpatialCellChat` 旧 slot 路径。
3. `net` 与 `netP` 的 cell/group 结果、p-value gating、pathway/LR join、返回列名
   和 `SparseChatArray` 稀疏性在迁移前后由针对性回归测试固定；不得重新引入稠密
   3D 通信数组。
4. 更新 `docs/VISUALIZATION_AUDIT.md` 和 `docs/REFACTORING_PLAN.md`，明确每个
   函数的迁移状态；对应代码、调用方和测试完成后，才允许将本 Note 从
   `proposed` 转为 `implemented`。
5. 在具备 R 环境后运行至少：
   `Rscript tests_dev/test-centrality-subset.R`、
   `Rscript tests_dev/audit_visualization_findings.R`、
   `Rscript tests_dev/test-SpatialCellChat-spatial-image.R`、
   `Rscript tests_dev/test-commun-flow-widget.R`，并运行
   `python .agents/skills/agent-notes-governance/scripts/validate_agent_notes.py --root .agents/notes`；
   所有目标检查通过后才可关闭本迁移边界。
