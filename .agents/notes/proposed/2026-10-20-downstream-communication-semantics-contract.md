# Agent Note: 下游通讯分析与可视化语义契约

Status: proposed
Governance: v1
Date: 2026-10-20
Decision type: migration
Scope: R/analysis.R, R/visualization.R, R/visualization_widget.R, R/modeling.R, R/utilities.R and downstream consumers, tests_dev, docs
Owner: maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

最终 11-slot 的存储迁移不能独自定义下游科学含义。当前源码存在以下可观察差异：

- `netVisual_aggregate()` 在绘图时按 pathway 选择 LR，清除 `pval > thresh` 后聚合；
  `computeCommunProbPathway()` 的 group 分支清除 `pval >= thresh`，结果还固定了生产时阈值。
  因图名含 pathway 而直接改读 `netP` 不是等价存储替换。
- `extractEnrichedLR_internal()` 最终与缓存 `LR.sig` 取交集，不从本次门控后的 `prob`
  重算 LR 集合；表达图的 `enriched.only` 不能一概解释为运行时重新检验。
- 中心性当前写 `slot$centr` / `slot$centr.cell`；部分 schema 示意中的
  `group$centr` / `cell$centr` 不是该生产者已实施的落点。
- `rankNet(do.stat=TRUE)` 以网络边配对检验，不能检验时写入零 p 值；
  `netAnalysis_dot()` 取 pattern 乘积的最大值，不是完整矩阵乘积。
  把所有现有输出当作正确基线会将科学解释错误固化为兼容性要求。

本提案细化 [迁移审计](2026-10-18-main-analysis-visualization-migration-audit.md) 的
“共享契约先行”门禁，不替代其全量迁移范围和验收，不改写既有 implemented 决策。

## Proposal

选择“逐函数语义卡 + 问题分流 + 有判别力的反例”作为下游迁移的前置条件，
因为同名参数和相同数组形状并不保证相同分析问题；拒绝仅按 slot 名机械改写，
也拒绝把统计更正、默认值修改与存储迁移混成一次不可归因的变更。

语义卡分别置于 [分析语义卡](../../../docs/ANALYSIS_SEMANTICS.md)、
[可视化语义卡](../../../docs/VISUALIZATION_SEMANTICS.md)，记录科学用途、观察单位、
计算、视觉编码、信息损失、统计边界、返回契约和反例。当前源码行为、预期用途和
迁移建议必须区分，静态推断标记为 `[INFERENCE]`；源码位置是事实依据，
不把语义卡当作函数已经完成迁移或端到端验证的证据。

每个迁移单元分为三类并单独审阅：

1. **存储迁移**：已核定语义不变，迁移所有生产者、消费者及名称对齐逻辑。
2. **实现/文档不一致**：明确错误位于哪一方；不能仅因存在旧输出就把错误写成新契约。
3. **科学/API 语义更正**：阈值、缺失值、归一化、检验单位、默认值、返回内容有变化时，
   先更新本提案或建立关联决策，再实施及验证；本提案不预先批准任何具体更正。

维护者负责争议裁定；实施者负责引用具体证据与受影响调用方。仅阻断受争议影响的
路径，不要求先修复全部下游统计问题才能迁移一个独立的表达展示路径。

[PP1 调用方方案](../../../docs/PP1_DOWNSTREAM_REFACTORING_PLAN.md) 将该门禁落实到
真实脚本的 IO、默认参数与二次消费链。采用生产者到消费者的完整契约检查，而不是
以绘图返回 class 或历史函数改名作为迁移完成依据：Visium 混合比例不得被普通硬标签
avg 静默替代，topic 社区不得被网络聚类替代，Seurat 桥接不得隐含改换通讯量纲。
该增补不批准新的群定义、归一化或统计方法；争议仍只阻断其依赖链。

## Constraints and invariants

- `net/netP`（LR/pathway）、`cell/group`（节点分辨率）、dataset（条件/样本）为独立轴。
  cell 到 group 的 `avg/sum`、支持度门控及分母沿用生产者定义，不在绘图层重新定义。
- 表达值、通讯权重、正边数、中心性、NMF loading、空间统计量不是可互换的量。
  每个图必须能说明颜色、大小、位置及缺失/零值各自代表什么。
- 逐函数记录筛选时点、`>`/`>=`、p 值缺失、缓存来源和归一化。
  当前差异是待核定事实，不意味着永久保留不同门限，也不意味着获准统一门限。
- 概率层复用 `SparseChatArray`；不得物化 cell × cell × layer 稠密张量。
  `marginSums(x,1/2)` 保留节点 × 层，`marginSums(x,3)` 跨层求和；
  后者不是 `apply(x,3,sum)` 的逐层总量。计数须区分正值边和存储条目。
- 中心性落点按当前已验证生产者使用 `centr/centr.cell`；若另行嵌套归位，必须作为
  独立的生产者/消费者切换，不可借本提案宣称已经归位。
- 比较分析必须明确同名节点与信号层的对齐、缺失层补零是否合理、比较方向和统计单位。
  `misc$.datasets` 是样本存档，顶层 joint 网络不可当作按 dataset 嵌套的列表。
- 不把表达共现、局部热点、网络差值或 DEG 注释描述为新证明的通讯因果效应。
  无法计算的 p 值不构成显著性证据；已定位错误单独裁定，不能通过兼容性测试固化。
- `.mode` 位于 `misc$.mode`，`params()` 访问的是 `misc$.param`；不新增旧 slot 或兼容 shim。
- 本次仅建立文档及迁移门禁，不更改函数签名、默认值、R 实现、生成 API 文档或统计方法。

## Alternatives considered

1. **只替换旧 slot 读取**：改动少，但实时阈值、缓存、节点单位、归一化和缺失值可能被
   静默改变；既不能证明正确，也难以定位科学解释回归，拒绝。
2. **一次性统一全部统计口径及 API**：可消除不一致，但扩大用户可见变化，无法区分迁移
   错误与新统计定义；裁定和回退成本高，拒绝。
3. **以所有旧输出为不可变黄金基线**：易对拍，但会将失败时 `p=0`、错误的轴对应等
   变为“兼容成功”；采用按语义分类的基线，已知缺陷不作为应保持的科学不变量。
4. **每个函数重复实现一套读取/聚合层**：局部灵活但容易漂移；复用已有访问器和稀疏
   运算，只有多消费者确有相同且已核定的语义时才提炼共享辅助，不预建通用抽象。

## Consumer impact

本次直接变更为两份语义卡、PP1 调用方方案、开发文档入口、重构计划、可视化审计
与父迁移提案的链接；下列代码、测试、教程和生成文档只作为消费者地图与后续迁移
检查范围，本次不修改。

- **生产者及支持层**：`R/SpatialCellChat_class.R` 的访问器/稀疏方法；
  `R/modeling.R` 的 cell/group 通讯、Visium 混合比例聚合、pathway、aggregate 和 complex expression；
  `R/database.R` 的 `searchPair()`/`extractGeneSubset()`；`R/utilities.R` 的 DEG、Seurat
  转换、NMF 缩放及 UMAP wrapper；`R/spatial.R` 的坐标和距离缓存。
- **分析消费者**：`extractEnrichedLR*`、`subsetCommunication*`、`netAnalysis_contribution`、
  `rankNet`、`netAnalysis_computeCentrality`、三个 `netAnalysis_signalingRole_*` 图、
  两类差异角色图、`selectK`、`identifyCommunicationPatterns`、`identifyCellTopics`、
  `identifyTopicCommunities`、`computeNetSimilarity*`、`netEmbedding`、`netClustering`、`rankSimilarity`、`netMappingDEG`。
- **绘图消费者**：四个优先入口及 circle/hierarchy/chord/spatial renderer；bubble、heatmap、
  diffInteraction、river/dot、embedding、topic、`spatialVisual_scoring`、Gi/Lee、
  `plotGeneExpression`；通信场/流及 `R/visualization_widget.R` 的 widget。
  `plotly_spatialLRpairPlot_shiny` 也消费表达图/选择数据；observer 生命周期不在本次改造。
- **动态选择**：`slot.name`、`do.group`、dataset comparison、层名/群名、缓存 `LR.sig`、
  `centr/centr.cell`、`pattern`、`topic`、`similarity` 与 `misc$.var.features`；PP1 的 Gi 统计列、
  ICR/OCR、Seurat reductions/assay/markers，以及重拟合后的关联状态。
  删除或更名任何接口之前仍须执行该批次完整引用检查，此清单不是删除授权。
- **现有测试**：`tests_dev/test-centrality-subset.R`、`test-netpathway-aggregate.R`、
  `test-SparseChatArray.R`、`test-pre-inference-new-schema.R`、
  `test-SpatialCellChat-spatial-image.R`、`test-grid-visualization.R`、
  `test-commun-flow-widget.R`、`audit_visualization_findings.R`；这些是已有测试入口，
  不是本次已通过结果，也不能替代尚缺的下游语义覆盖。
- **使用者与文档**：`tests_dev/test_data/visium_PP1.R` 的图、CSV、metadata 与 Seurat 消费链；
  `tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.Rmd`
  及生成 HTML；目标导出的 `man/*.Rd`、roxygen、`NAMESPACE`；`docs/DATA_STRUCTURE.md`、
  `REFACTORING_PLAN.md`、`VISUALIZATION_AUDIT.md`、`README.md`。
  未来签名/默认值/返回语义变化必须同步调用方、示例和帮助页；本次不重新生成 HTML/Rd。
- **配置与依赖**：本次没有配置或依赖变更；不把可选绘图/分解包缺失转换为算法回退，
  不修改 `renv.lock`、`DESCRIPTION` 或运行环境设置。

## Consequences

- 迁移批次增加一个可审阅的科学语义门禁，避免把“成功返回图对象”当作信息正确的证明。
- 后续文档可从同一语义卡生成完整的用途、计算、解释与限制描述，但当前帮助页不会因
  新卡出现而被视为已修正；语义卡也不替代源代码/API 文档的长期权威。
- 显著性、归一化与缺失值的争议保持显式开放；它们只阻断相关路径，不扩大为整包重写。
- 本记录与父审计保持 `proposed`。语义卡完成、validator 通过或独立审阅同意都不足以
  将全量迁移宣布完成；运行时与图形证明必须随真正的实现批次提供。

## Evidence

- `R/visualization.R:1440-1545,1604-1738` 与 `R/modeling.R:1680-1767`：
  LR 实时聚合、pathway 缓存、阈值等号及 cell/group 分支差异。
- `R/analysis.R:913-946`：缓存 LR 选择与本次阈值不一致；
  `R/analysis.R:213-223`：中心性实际写回；
  `R/SpatialCellChat_class.R:1349-1356,1477-1493,1526-1533`：访问器边界。
- `R/analysis.R:296-308,565-594,715-735,1249-1273,1345-1571,3081-3089`：
  contribution 归一化、比较检验、NMF 输入、相似性/嵌入与相反方向变化过滤。
- `R/visualization.R:4248-4280`：pattern dot 的最大乘积分数；
  `R/analysis.R:3467-3557,4016-4168`：cell topic 与方向场。
- `tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.Rmd:298-480`：
  表达背景、群体网络、角色、热点、共现和 topic 的用户分析链。
- 2026-10-20 主代理使用 R 4.4.3 与项目库运行标准输入数值探针，退出码 0。
  实际调用 `computeCentralityLocal()`、`.sc_communication_field_layer()`、
  `computeCommunProbPathway()` 和旧式 plain-list 输入的 `extractEnrichedLR_internal()`；
  后者不是最终 schema 的公开入口验收。另对 contribution、差异过滤、旧绘图门控执行
  源码公式对照，不将其表述为运行了相应绘图函数。数值与复现输入见分析语义卡的验证节。
  尚未运行名称重排/缺失层、全部公开入口或图形交互验收。
- 独立 reviewer 对本提案及父审计与源码的关系给出 Approve：只同意语义门禁方向，
  不代表维护者已裁定具体统计更正，也不代表语义卡或全量迁移已验收。
  审阅指出第一实施批次可围绕两个表达入口和共享 LR 选择链，不应只验证 `features=`
  成功便宣布整个入口完成；空 LR 交集仍返回基因的行为必须另行裁定。
- 既有 implemented Notes 中的历史测试不计入本次验证；本次不创建永久回归测试，
  以后按已裁定的消费者可见含义设计测试，不锁定已定位的错误输出。
- 2026-10-21 的 PP1 真实输入、最终 schema 子链、实际图形与 Seurat 量纲探针详见
  [调用方方案 §7](../../../docs/PP1_DOWNSTREAM_REFACTORING_PLAN.md#7-本轮运行证据与限制)。
  已有评分/两指标角色图的局部成功，与 Gi 形状、旧存储、Seurat 桥接的实际失败分开记录；
  不将源码映射当历史 niche API 等价证明，不将缺依赖的完整统计/社区路径计为通过。

## Acceptance criteria

1. 两份语义卡覆盖本提案 Consumer impact 的主要分析链，每个入口或有明确差异说明的
   函数组均包含八项语义字段；源码事实、静态风险、运行证据和推荐变化可区分。
2. 每个争议项列出证据、受影响消费者、建议处理方向和未决边界；在维护者裁定之前，
   不将其写成已批准的新默认行为，不将已知缺陷作为应保持的黄金输出。
3. 语义反例能区分：同边缘和而不同路由；正通讯强度而零净向量；`p == thresh`；
   缓存阈值与本次阈值不同；原始量与显示缩放；名称顺序置换与缺失层。
   每批实现只要求其受影响语义的反例，不能用未执行场景冒充通过。
4. 真正迁移时必须有最终 11-slot 输入、数值/名称对齐对照和实际返回/图形观察；
   plotly 对象构造成功不能替代浏览器交互，recordedplot 必须可重放，空间图必须核对
   点位/坐标朝向。不需要图形的分析路径运行实际函数并观察数据即可。
5. 公开 API 与所有受影响调用方、roxygen/Rd、教程及状态文档随实施批次同步；
   若需要改变已实施 centrality/producer 契约，先关联后继决策，不回写旧决策。
6. 本提案在关联迁移单元及父审计要求的消费者均完成语义核定、实现、证据关联和文档
   同步前保持 `proposed`；不能以本次文档交付关闭迁移。
7. 创建或转换 Note 后运行
   `python .agents/skills/agent-notes-governance/scripts/validate_agent_notes.py --root .agents/notes`；
   对新增或更改的 Markdown 相对链接检查目标存在。验证既有 Note 集合时不使用 `--strict`。
8. PP1 范围只有在 Visium 输入/生产者、Gi→metadata→Lee、topic→community→Seurat、
   图/CSV 及保存重载均有链级证据后才可关闭；[B0–B7 与 T01–T22](../../../docs/PP1_DOWNSTREAM_REFACTORING_PLAN.md#5-细化迁移批次与完成门禁)
   是实施与验收设计，不是本次通过清单，也不缩减父审计的其它消费者范围。
