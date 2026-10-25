# 下游可视化逐函数语义卡

## 阅读边界与证据级别

本页是迁移前的科学语义审阅材料，不是已完成迁移、已修复缺陷或已验证可运行的 API 手册。决策边界见 [下游通讯分析与可视化语义契约（proposed）](../.agents/notes/proposed/2026-10-20-downstream-communication-semantics-contract.md)；父迁移审计仍为 proposed。分析量的定义、NMF/topic 分解、中心性、相似性和生产者前置条件见 [分析语义卡](ANALYSIS_SEMANTICS.md)，本页只补充图形消费时发生的选择、变换和解释。

- **源码事实**：具名函数和行范围支持的实际读取、计算、绘图调用和返回分支。记录旧字段不表示最终 11-slot 对象已经支持该入口。
- **[INFERENCE]**：静态推导的科学解释、故障风险或预期差异，不是运行结果。
- **建议/未决**：第 9 节 V-* 问题需维护者裁定；当前缺陷不是必须保留的兼容性基线。本页不批准修改默认值、统计方法、签名或缓存布局。
- 每卡固定八项。**第八项都是尚未执行的入口级验收设计**，不是测试结果。本页作者没有运行 R、build、tests、lint、formatter 或浏览器检查。
- 主代理本次完成的孤立数值/helper/公式探针集中在 [分析文档第 8 节](ANALYSIS_SEMANTICS.md#8-本次运行证据与未验证范围)。它们支持局部机制审阅，不构成任何图形入口、整条分析链或浏览器交互的通过证据；既有审计中的历史运行结果也不计为本次验收。
- **2026-10-21 PP1 调用方增补**：[PP1 下游重构与测试方案](PP1_DOWNSTREAM_REFACTORING_PLAN.md) 明确 Gi metadata→热点类别→Lee、group/cell 层选择、niche/Seurat、场图与文件导出的实际调用契约；补充最终 matrix 坐标下 Gi wrapper 的已观察失败，以及评分图/角色热图的实际观察。本页既有判别反例仍是待验收设计。

## 共同坐标、量纲与阅读规则

令 $P_{ij\ell}$ 为某 dataset、某节点分辨率下的通讯权重，行是 sender，列是 receiver，层是 LR 或 pathway。以下轴必须分开：

| 轴/量 | 含义 | 不允许的替代解释 |
|---|---|---|
| `net` / `netP` | LR / pathway 层 | 不是 cell / group 的另一种写法 |
| cell / group | 单细胞或 spot节点 / 群节点 | `individual` 图名不自动表示 cell 分辨率 |
| dataset | 样本、条件或保存的对象实例 | joint 网络不是按 dataset 嵌套的旧列表 |
| 表达矩阵 | gene × cell/spot，或 metadata 数值 | 表达、共表达不等于通讯权重或显著通讯 |
| group 网络空间位置 | 每群坐标分量的中位数 | 群间弧线不是真实单细胞路径 |
| 中心性、topic loading、Gi、Lee L | 不同来源、不同统计单位的量 | 不能统一标为“通讯强度” |
| 显示坐标与缩放 | 坐标交换、y 翻转、相对色域、图形大小 | 不产生新空间测量、物理流速或跨图可比绝对量 |

**共同生产者边界**：cell→group 的 `avg/sum` 有支持度门控和各自分母，不是绘图时把细胞边任意相加；详见分析卡的“轴、量纲与共同前提”。group pathway 生产者清除 `p >= thresh`，cell pathway 分支不施加同一 group p 值门控；`netP` 固定了生产时选择。`netVisual_aggregate()` 则从 LR 在绘图时清除 `p > thresh` 后求和，不能仅因标题是 pathway 就机械改读 `netP`。中心性当前由生产者写入 `slot$centr` / `slot$centr.cell`，不是 `group$centr`（[生产者](../R/analysis.R#L213-L223)）。mode 和参数也应分开：`misc$.mode` 与 `params()` 访问的 `misc$.param` 不是同一项。

**坐标边界**：下列多数旧二维入口交换坐标前两列并反转 y；field/flow 同时变换向量分量。它们常按行位置绑定值与坐标，不能把“行数相同”当作名称已对齐。人工 3D 展示在 `plotFeatures()`、`plotStatistics_Gi()` 中把 `z` 设为 `three.dim.z`，没有测得的 z 轴；不能据此宣称三维组织、第三空间轴的热点或三维迁移轨迹。

## 1. 表达背景，而非通讯网络

### 1.1 `spatialFeaturePlot`

依据：[选择与连续表达](../R/visualization.R#L437-L579)、[二值表达](../R/visualization.R#L580-L694)、[复合物表达](../R/modeling.R#L2531-L2543)、[缓存 LR 选择](../R/analysis.R#L913-L946)。

1. **科学问题**：选定基因、LR 对或 pathway 相关基因在组织哪里表达；binary 模式检查同一 spot 的 ligand/receptor 表达组合，不检查两个细胞间的通讯边。
2. **输入单位与轴**：`features/signaling/pairLR.use` 三选一。连续模式从旧 `data.signaling/data` 取 gene × cell 值，或整体改从 metadata 列取值；点始终是 cell/spot。`do.group` 只控制 LR 富集筛选来源，不把表达图点变成群。坐标要求两列，交换后反转 y。
3. **计算定义**：连续模式 pathway 经 `extractEnrichedLR()` 取基因；显式 LR 的 enriched group 分支清除 `p>thresh` 后要求该层存在正边，cell 分支只检查已有正边。连续模式展开 complex 为亚基；binary 模式用 `computeExpr_LR()`，complex 用亚基几何均值，再按严格 `> cutoff` 分为 ligand、receptor、Both、None。连续模式 `<= cutoff` 置 NA。
4. **视觉编码**：连续图每 feature 一面板、颜色为该 feature 表达/metadata 值，NA 灰色；binary 颜色是上述四类，点大小不是通讯量。每面板自有连续色域；`show.legend.combined` 不是共享数值范围。Both 表示同一位置两者过表达阈值，不表示检验通过。
5. **信息损失**：丢失 sender→receiver 配对、网络权重和 p 值；binary 再丢失表达幅度。展开亚基图与复合物几何均值图不是同一测量。按行绑定也不修复名称错位。
6. **统计边界**：不做表达差异或通讯新检验。pathway 选择最终依赖缓存 `LR.sig`，本次 `thresh` 不一定改变返回 LR；显式 LR 分支才按本次 prob 检查。无 LR 交集但仍返回基因的 helper 行为列为缺陷，不能解读为显著信号（V-003）。
7. **前置计算及返回**：需要表达、坐标、注释；LR/pathway 模式还需 DB/LR，enriched 路径需相应网络或缓存。返回 ggplot 或 patchwork；源码 `ncol<=1` 时只取 `gg[[1]]`，即使有多个 feature 也不是完整多图列表（V-004）。旧表达槽读取是迁移门禁。
8. **判别反例（未执行）**：保持所有细胞表达不变，只把远距离通讯边置零，表达图不应被解释为通讯仍存在；令某复合物一个亚基为零、另一个为正，单亚基连续图可明亮而复合物表达为零，binary 分类不应将该复合物判为阳性。缓存阈值 0.05、本次 0.01 与显式 LR 分支应分别核对。

### 1.2 `plotly_spatialLRpairPlot` 与显示 helper `plotFeatures`

依据：[LR/feature 选择及真实返回](../R/visualization.R#L7289-L7423)、[普通/3D helper](../R/visualization.R#L7457-L7776)、[2D trace](../R/visualization.R#L7779-L7897)。

1. **科学问题**：交互查看最多两个 feature 或一个 LR 的表达位置、数值与 cell id；不是画 LR 通讯边。`plotFeatures()` 本身只消费三个数值列，不知道输入是表达还是其他分数。
2. **输入单位与轴**：wrapper 读旧 `object@data` 或 metadata，行是 cell/spot、每 feature 一个图；坐标两列交换。`do.group` 控制 enriched 筛选，不控制节点分辨率。helper 的 `normal/2d/3d` 是渲染方式；wrapper 固定调用 `method="2d"`。
3. **计算定义**：显式 LR enriched 的 group 分支清除 `p>thresh`，cell 分支检查正边。随后直接匹配 ligand/receptor 名称到表达行；complex 展开语句被注释，**没有**与 `spatialFeaturePlot(binary)` 相同的 complex 几何均值。`cells.highlight` 只改变显示透明度。helper 3D 的 `z=three.dim.z` 是人工平面。
4. **视觉编码**：颜色为 feature 值，hover 包括坐标、cell id、score；非选中点较透明。2D y 反向且等比例；3D 是固定 z 平面的 scatter3d，不是 score 高度，也不是测得的三维坐标。各 feature 的色域独立；计算出的 `max.feature.data` 没有传为共享色域。
5. **信息损失**：没有通讯方向、边权和统计证据；复杂 LR 可只留下可匹配的单个 partner，不能把一张图当作完整 LR 展示。人工 z 与 camera 只是布局。
6. **统计边界**：无新检验；enriched 不检验用户圈选细胞。highlight 后色条/trace 是否一致需实际浏览器检查，不能由对象类推断。helper 按位置重命名输入列并按固定 trace 下标设置色条，边界类别需验证（V-004、V-018）。
7. **前置计算及返回**：wrapper `do.plot=TRUE` 返回 **plotly 对象的 list**，不是单个 plotly；FALSE 返回可展示 feature 的**数目**，不是表达数据。helper normal 返回 ggplot/patchwork，2d/3d 返回 plotly。需可用表达/metadata、坐标，LR 分支需 DB/LR，enriched 时另需网络。
8. **判别反例（未执行）**：LR 的 receptor 是未作为表达行保存的 complex，核对是否只出现 ligand；选中所有点、零点及一个点，检查 hover id、透明度和色域。改变 `three.dim.z` 只能移动显示平面，不能改变数值或产生测得 z。

关联消费者：`plotly_spatialLRpairPlot_shiny()` 固定 `enriched.only=FALSE`，用上述数目创建面板、用 `customdata` 关联选择，并返回 `shiny::shinyApp()` 对象（[选择与面板](../R/visualization.R#L7936-L8044)、[返回](../R/visualization.R#L8111-L8114)）。它不是通讯重新推断入口；本次不改 observer 生命周期。

### 1.3 `plotGeneExpression`：violin / dot / bar

依据：[wrapper](../R/visualization.R#L6632-L6669)、[dot](../R/visualization.R#L6703-L6732)、[violin](../R/visualization.R#L6760-L6799)、[bar](../R/visualization.R#L6885-L6942)。

1. **科学问题**：LR/pathway 相关基因在各细胞群中的表达分布、表达比例或汇总表达；它为通讯推断提供表达背景，不给出通讯网络验证。
2. **输入单位与轴**：旧 `data.signaling` 被交给 `Seurat::CreateSeuratObject(counts=...)`；meta 与 idents 定义分组。gene × group 是 dot/bar 的显示轴，violin 保留组内单细胞分布。若 features 与 signaling 同给，源码警告并优先 features。
3. **计算定义**：pathway 基因来自 `extractEnrichedLR()` 缓存选择。violin 经 Seurat VlnPlot helper；dot 经 Seurat DotPlot，默认 `scale=TRUE`、截断范围 `[-2.5,2.5]`；bar 按 `truncatedMean/triMean/median` 汇总 Seurat 的 data 层，默认截尾均值 `trim=0.1`。wrapper 没有执行新的表达归一化，实际取哪一 assay 层及版本行为须随迁移确认。
4. **视觉编码**：violin 宽度/形状展示分布，颜色为组/条件；dot 点大小为表达细胞比例、颜色为 Seurat 返回的平均/缩放表达，默认可翻转轴；bar 条高是选定表达汇总。它们的颜色/大小均不是通讯 p 值或边权。
5. **信息损失**：dot/bar 丢失组内异质性及空间位置；dot 的按基因缩放不保留跨基因绝对量，bar 不能区分少数极高表达与更多低表达的全部差异。
6. **统计边界**：默认没有独立生物学重复的差异检验；`barplot_internal` 的可选统计不能把已汇总 group 数值变成重复样本。缓存 enrichment 仍受 V-003 约束。
7. **前置计算及返回**：需可构建的 Seurat 对象和相应 assay 层、meta、features；返回 dot 的 ggplot，或 violin/bar 的 patchwork。meta 名称不一致时当前只覆盖行名而不重排；bar 无 split 的 matrix `$gene` 赋值及最后 feature 复用 `df.use` 是静态风险（V-005），不保证该分支能完成。
8. **判别反例（未执行）**：一组 10% 细胞高表达、另一组 100% 低表达，令平均相等，dot 大小应不同、violin 分布不同；独立置换 meta 行顺序可暴露“改名代替对齐”。至少两基因检验 bar 最后一面板是否真正来自最后基因。

## 2. LR 网络与布局

### 2.1 `netVisual_aggregate`

依据：[完整入口](../R/visualization.R#L1440-L1545)、[pathway 生产者](../R/modeling.R#L1680-L1767)。

1. **科学问题**：一个 pathway 所属、当前绘图筛选后的 LR 在群间汇总成什么有向网络；不是从 pathway 缓存机械取一层。
2. **输入单位与轴**：`LR$LRsig` 按 pathway 查找 LR，与旧 `net$prob` 第三轴取交集；矩阵行 sender 群、列 receiver 群。没有 cell 分辨率参数。
3. **计算定义**：从 LR 概率清除 `pval > thresh`，去掉总和为零的 LR，再按 LR 轴求和 $A_{ij}=\sum_{\ell\in\mathrm{selected}}P_{ij\ell}$。之后交给布局 renderer 做 source/target、top 等显示筛选。**不是**默认读 `netP`；后者的生产时点、阈值等号都不同。
4. **视觉编码**：circle/hierarchy 用有向边宽；chord 用带宽；spatial 用群中位数位置间弧线。节点默认等大；`vertex.weight=NULL` 使用群细胞数；spatial 的 incoming/outgoing 字符选项另取 `netP$centr`，不保证与本次 LR 门控同一时点。
5. **信息损失**：求和隐藏 LR 组成、相反变化与各 LR 证据差异；group 早已丢失单细胞路由。图上群间大边可能由一个 LR 或多个较弱 LR 贡献。
6. **统计边界**：阈值只过滤既有 LR p 值；求和不产生 pathway p 值、因果效应或独立重复。`p==thresh` 在本入口保留，与 pathway producer 和 bubble 可能不同（V-002）。
7. **前置计算及返回**：需群 LR 网络、p 值和注释；spatial 需坐标/idents，centrality 权重选项另需相应缓存。hierarchy/circle/chord 返回 recordedplot；spatial 返回 ggplot。当前旧字段及 vector `vertex.weight` 条件判断需迁移审阅（V-001、V-006）。
8. **判别反例（未执行）**：两条 LR 的 p 分别恰等于阈值、严格小于阈值，比较实时 LR 求和与先前 `netP` 缓存；再保持汇总矩阵不变交换 LR 贡献，aggregate 图应相同而 individual 图不同。

### 2.2 `netVisual_individual`

依据：[完整入口](../R/visualization.R#L1604-L1738)。

1. **科学问题**：同一 pathway 中各个 LR 的群间网络有何差异；**individual 指 LR 对，不指单个细胞**。
2. **输入单位与轴**：与 aggregate 相同的旧 group LR 张量，可再由 `pairLR.use` 取子集；每面板是一层 LR，节点仍是群。
3. **计算定义**：清除 `p>thresh`、保留总和非零的 LR，逐层交给 renderer，不跨 LR 求和；默认 `edge.weight.max=max(prob)` 来自所有选中 LR，使同一调用的面板使用共同边宽分母。
4. **视觉编码**：标题为 LR 展示名，边/带宽为该 LR 群间权重；四种布局差异见 2.3。circle/chord/空间面板列表与 hierarchy 的双 receiver 分区不是同一种组合图对象。
5. **信息损失**：仍无细胞间路由，p 值仅作门控而不直接编码；弱 LR 可因共享宽度范围显得不明显，但这不等于不显著。
6. **统计边界**：无 LR 间差异检验；逐图比较边宽需保持共同上限与同样筛选。独立调用各自自动取最大值时，最粗边不能直接比较。
7. **前置计算及返回**：需同 aggregate；hierarchy 返回一个 recordedplot，circle/chord 返回 recordedplot 的 list；spatial 分支尝试返回 ggplot 的 list，但引用本函数未定义的 `labels`、`idents.use`，不能承诺实际运行成功（V-006）。roxygen 的单一 recordedplot 返回描述不覆盖全部分支。
8. **判别反例（未执行）**：两个 LR 使用不同 sender→receiver 路由却有相同总权重，逐 LR 图应区分；扩大每群的细胞数但保持群网络固定，不应把新增细胞画成 individual 网络节点。空间分支须单独运行并检查坐标与标注。

### 2.3 Renderer：`netVisual_circle` / `netVisual_hierarchy1/2` / `netVisual_chord_cell_internal` / `netVisual_spatial`

依据：[hierarchy1](../R/visualization.R#L1795-L1928)、[hierarchy2](../R/visualization.R#L1985-L2121)、[circle](../R/visualization.R#L2176-L2299)、[spatial](../R/visualization.R#L2380-L2509)、[chord](../R/visualization.R#L3567-L3698)。

1. **科学问题**：对已经确定的有向加权网络选择可读的排布；renderer 不负责重新推断通讯。四种布局分别突出整体拓扑、receiver 分区、群间带状联系、群级地理摘要。
2. **输入单位与轴**：circle/hierarchy/spatial 输入 sender × receiver 矩阵；chord 也接受 source/target/prob 表。节点含义取决于调用方；aggregate/individual 传群矩阵。spatial 另用每个 cell/spot 的坐标和 labels 计算群节点位置。
3. **计算定义**：circle/hierarchy/spatial 先用全矩阵 `quantile(net,1-top)` 清除严格较小值，再筛 source/target；分位数包含零，不等于精确保留 top 比例正边。hierarchy1/2 只展开各自 `vertex.receiver` 列并复制 target 节点作布局。spatial 每群分别取 x/y 中位数，默认去掉自环。chord 保留正边，按 sector 顺序、scale、reduce 渲染；保持孤立 sector 时插入不可见极小 fake link，并可能令 `scale=TRUE`。aggregate/individual 的 chord 调用未传 `top`，不应假设四布局过滤相同。
4. **视觉编码**：circle 的圆周位置无空间单位，边色为 sender；hierarchy 的位置是 receiver 分区结构，不是时间先后；chord 的带宽与 sector 宽度受 `scale` 影响；spatial 的弧线端点是群中位数，背景散点才是 cell/spot。前三类非空间位置都不能解释成距离。circle/hierarchy/spatial 边宽通常为 `0.3 + weight/max * edge.width.max`（scale 开启）或 `0.3 + weight * edge.width.max`；节点大小还有布局各自偏移项。
5. **信息损失**：top、去孤立点、去自环、chord reduce 会隐藏结构；布局和图形缩放会隐藏绝对大小。群空间中位数可落在没有任何细胞的空隙中，弧线不是单细胞路由或真实分子运输轨迹。
6. **统计边界**：不计算 p 值；`top` 是权重显示分位数，不是显著性阈值。chord 的 fake link 仅为布局保留，不是观测到的弱通讯；不能参与科学计数。
7. **前置计算及返回**：输入名称、颜色、顶点权重必须对齐；hierarchy 需有效 receiver 分区，spatial 需对应 labels/coordinates。circle/hierarchy/chord 在图形设备绘制并返回 recordedplot，spatial 返回 ggplot。hierarchy2 的筛选分支读取不在参数中的 `object@idents`，是待修问题而非契约（V-006）。
8. **判别反例（未执行）**：同一群分布在组织两端，中位数位于空白区，空间弧线必须标成群摘要；给矩阵大量零及少量相同正边，top 不会恰保留指定数量。加入纯自环，circle 与默认 spatial 应明确呈现不同过滤结果；chord 与 circle 的 top 差异须实图核对。

## 3. 通讯表、矩阵与条件差异

### 3.1 `netVisual_bubble`

依据：[表选择及变换](../R/visualization.R#L3079-L3300)、[图形及返回](../R/visualization.R#L3302-L3406)、[实际表生产门控](../R/analysis.R#L2405-L2448)。

1. **科学问题**：哪些 LR 在指定群对存在推断通讯，以及不同 dataset 的这些条目如何并列；不是比较数据集差异的统计检验。
2. **输入单位与轴**：每点是 source→target 群对 × LR（comparison 时另加 dataset）。wrapper 固定 `slot.name="net"` 调用 `subsetCommunication()`；本函数没有 cell 分辨率参数。
3. **计算定义**：当前表 producer 对 LR 网络输入清除 `p>=thresh` 并只展开正边，p 缺失则输出 NA。bubble 再把 p 分为 `>0.05`、`(0.01,0.05]`、`<=0.01` 三档；零 prob 置 NA，保留 `prob.original`，颜色量变成 $-1/\log(prob)$，并按显示分位数截断。`max.dataset/min.dataset` 筛选显示值的相对极值，非差异显著性；缺群对用 NA 占位。
4. **视觉编码**：横轴为有向群对（及 dataset），纵轴为 LR 展示名；颜色是变换、截断后的 prob，大小是 p 档而非效应大小。默认色条常只标 min/max，不能据色差还原原始倍数。p 档标签的 `<` 文案没有准确反映实际包含等号。
5. **信息损失**：低于门控的边消失，p 档丢失精确值，分位截断丢失尾部大小；空白混合缺层、无正边、被过滤或 NA。无原始 sender/receiver 的 cell 路由。
6. **统计边界**：图上 p 来自各自网络推断，不是跨 dataset 差异 p。无 p 的边不能被自动视作显著。`thresh` 与固定 p 档阈值是两套设置；比较方向与显著性不可混同（V-002、V-009）。
7. **前置计算及返回**：需 `subsetCommunication` 可消费的网络与 LR 注释；comparison 仍沿用旧 dataset 布局。返回 ggplot，或 `list(communication=df, gg.obj=g)`；表内 `prob` 已是显示量，原始值在 `prob.original`。单数据集非标准权重分支引用未限定的 `prob.original`，末尾无条件用 `pairLR.use` 重设 LR levels，需核定默认调用是否丢标签（V-009）。
8. **判别反例（未执行）**：对同一边设 `p==thresh`，与 aggregate 比较；给 `prob=0.1,0.2` 检查颜色不保持 1:2。令一个 dataset 缺层、另一个真实零，返回表应可追溯而非宣称“显著下降”；省略 `pairLR.use` 仍须看到正确 LR 轴。

### 3.2 `netVisual_heatmap`

依据：[输入分支](../R/visualization.R#L2734-L2807)、[过滤、注释及 draw](../R/visualization.R#L2810-L2911)。

1. **科学问题**：显示群间通讯矩阵，或第二 dataset 相对第一的 count/weight 差值；它保留 sender–receiver 路由，比只看总量更细。
2. **输入单位与轴**：行 sender 群、列 receiver 群。single 可选 signaling 层或聚合 count/weight；merged 分支读旧 dataset 嵌套聚合矩阵，忽略 signaling 所暗示的逐层选择。`slot.name` 决定 single signaling 是 LR 还是 pathway。
3. **计算定义**：single LR 清除 `p>thresh`，pathway 直接读缓存；single measure 直接读预聚合矩阵。merged 逐位置 `obj2-obj1`，没有在这里按名称重排。可筛 source/target、行列、去对角；侧栏为 `rowSums(abs(mat))` 和 `colSums(abs(mat))`，**差值情况下是绝对变化和**，不是净增加。通常零单元置 NA 显示白色。
4. **视觉编码**：矩阵颜色为原量或有符号差值；负/正值分色，中间零色。顶部/右侧条长是对应 receiver/sender 的绝对值和；不是平均强度，也不是显著性。
5. **信息损失**：聚合 measure 隐藏 LR 身份；零转 NA 混合未显示零和缺失；自动色域使不同调用不能凭相同色深比较绝对量。侧栏丢失增减方向。
6. **统计边界**：差值无差异检验和置信区间。源码 `sum(mat)==0` 会拒绝正负完全抵消但非空的差异网络；以有符号边缘和判孤立也可能误判（V-007）。
7. **前置计算及返回**：需相应 group prob、aggregateNet count/weight 或可比较 dataset 矩阵。末表达式是 `draw(ht1)`，不是显式 `return(ht1)`；返回的是 draw 的结果（通常为已绘制 HeatmapList，[INFERENCE]，须按安装版本实测）。`cluster_columns` 实际用了 `cluster.rows`；remove.isolate 重复删除、行注释复用列名等列为 V-008。
8. **判别反例（未执行）**：差值矩阵含 `+2,-2` 总和零却有明显重路由，不能验收成“无通讯”；把两个 dataset 的群名顺序反转，差值应按名称核对；分别设置行聚类与列聚类，并只显示非方形子矩阵，检查真实轴与注释。

### 3.3 `netVisual_diffInteraction`

依据：[完整入口](../R/visualization.R#L2572-L2692)。

1. **科学问题**：第二 dataset 相对第一在群间总体边数/强度哪里增加或减少；不识别具体变化 LR，也不建立条件因果关系。
2. **输入单位与轴**：旧 `net[[comparison]][[measure]]`，measure 为 count/weight 或相应 merged 聚合；节点是群，边是有向群对。两矩阵必须具有同名同序节点，代码没有自动对齐步骤。
3. **计算定义**：$D=A_2-A_1$，source/target 筛选后可去孤立点；按 `quantile(abs(D),1-top)` 去小差值，再以有符号权重建图，设颜色后把边权取绝对值计算宽度。
4. **视觉编码**：默认红为第二条件增加、蓝为减少；宽度为 `abs(D)` 的原量/相对缩放，箭头仍表示 sender→receiver，不表示变化时间；圆周布局非空间。开启 edge label 时在取绝对值前保存有符号标签。
5. **信息损失**：没有原始两个网络大小、LR 组成、cell 路由；同样差值可能来自完全不同基线。top 隐藏较小变化，独立最大值缩放妨碍跨图比较。
6. **统计边界**：没有差异 p 值；count 差是有效边/层数量的上游定义，不能说成细胞数差。`remove.isolate` 用有符号行列和，可把正负抵消误作无边（V-007）。
7. **前置计算及返回**：需已聚合、可比且名称对齐的两个 dataset 网络；当前旧 nested 布局需迁移。绘制 base 图并返回 recordedplot，不返回差值矩阵。
8. **判别反例（未执行）**：反转 `comparison=c(1,2)` 为 `(2,1)`，红蓝应反转而绝对边宽保持；同一节点两边 `+1,-1` 不应在科学上被认作孤立。节点顺序置换必须在数值和图中都被发现。

## 4. 中心性空间分数与空间统计

### 4.1 `spatialVisual_scoring`

依据：[缓存读取](../R/visualization.R#L724-L770)、[单量和合并编码](../R/visualization.R#L772-L933)、[中心性写回](../R/analysis.R#L213-L223)。

1. **科学问题**：某信号的 incoming/outgoing 或其他已计算中心性在组织中如何分布；不是直接画边，也不是识别统计热点。
2. **输入单位与轴**：一次仅一个 signaling；缓存轴是 measure × node × layer。`do.group=FALSE` 读 `centr.cell`，TRUE 读 `centr` 后把每群分数复制到该群所有 cell/spot；群模式不是每点独立测得分数。
3. **计算定义**：只实现 `method=="centrality"` 的值来源。普通模式使用缓存值；binary 用 `value>0`。merge 要求恰两个 measure，分为均零、仅第一正、仅第二正、两者均正；连续 merged 的 Both 点颜色值为两个分数的算术均值。
4. **视觉编码**：位置为 cell/spot 坐标，颜色为分数或四分类；merge 的三个非零类别各自有色阶，不是一条共同量尺。outdeg/indeg 标成 outgoing/incoming，但其他中心性单位不同；灰色不等于没有所有形式通讯。
5. **信息损失**：out/in strength 丢失具体路由；group 复制抹去群内异质性；binary 抹去强弱，merge 均值不能分辨 `(1,9)` 与 `(5,5)`。
6. **统计边界**：不新增 p 值或热点检验。正中心性不是显著性；若缓存来自不同阈值，图形不会重新校正。非负分类逻辑不能无条件用于可为负的自定义分数。
7. **前置计算及返回**：先运行 `netAnalysis_computeCentrality()`；实际生产落点仍是 `centr/centr.cell`。未 merge 时 combined=TRUE 返回 patchwork，FALSE 返回 ggplot list；merge 返回单 ggplot。group 按 `levels(idents)` 与 centr 的位置对应，需名称核对（V-010）。
8. **判别反例（未执行）**：同样每节点 out/in 总量但不同路由的两个网络，此图可相同而网络图不同；同群两个点分数悬殊时，group 模式应显示同值而 cell 模式不同。反转 centrality 节点顺序测试对齐，merge 两个不等量纲 measure 不能被默认当作可比较平均。

### 4.2 `spatialGiPlot`

依据：[分数聚合与两次调用](../R/visualization.R#L8841-L9006)、[返回分支](../R/visualization.R#L9008-L9064)。

1. **科学问题**：由单细胞/spot 中心性形成的输入量，哪些局部邻域显示高值/低值聚集；不是直接给某条通讯边做空间检验。
2. **输入单位与轴**：`centr.cell[measure,cell,signaling]`，可选多个信号层；空间统计的观察单位是 cell/spot，邻域为二维 kNN。group.highlight 控制图形高亮，未限定计算样本。
3. **计算定义**：默认 `do.binary=TRUE` 先把每个 measure×cell×layer 的正值改为 1，再跨 signaling 求和；因此默认输入是该细胞**正分数信号层数**，不是总通讯强度，也不是对最终 Gi 二值化。FALSE 才跨层求连续分数和；每 measure 交给 `plotStatistics_Gi()`，并另以 computation 模式再算一次取数据。
4. **视觉编码**：底层 normal/3d 颜色为 Gi/Gi* 统计量，p.value 可影响点形；scoring 则画输入分数，density/contour 为其后处理。多个 measure 组合；图名的 incoming/outgoing 由 measure 推出。
5. **信息损失**：跨信号求和失去是哪条 LR/pathway；默认二值化使一条强信号与一条弱正信号同贡献 1。热点显示也不保留原始边路由。
6. **统计边界**：Gi 的局部 p 值不等于中心性显著性、某 LR 的通讯 p 值或热点间差异 p；本 wrapper 没有多重比较校正。改变 k、binary、信号集合改变科学问题，不能作为纯渲染改动。
7. **前置计算及返回**：需 `centr.cell`、坐标与 idents。默认 `return.object=TRUE` 将 Gi 和 p 写入返回对象的 meta，`show.plot` 仅控制显示；FALSE 返回 patchwork 或 plotly，meta 修改不随该返回暴露。保存 p 固定读 `Pr(z != E(Gi))`，而 Gi* helper 使用星号键，属于 V-011。超过两个非 ggplot 面板的组合也只取前两个，须裁定支持范围。
8. **判别反例（未执行）**：一个 cell 在一层分数 100、另一个在三层各 0.1，binary 与连续输入的排序应相反；分别运行 Gi/Gi* 核对返回 meta 的统计量和 p 列，不能只观察图是否返回。非高亮组仍应参与邻域计算。

### 4.3 `plotStatistics_Gi`

依据：[Gi/Gi* 计算与普通图](../R/visualization.R#L9087-L9208)、[人工 3D](../R/visualization.R#L9210-L9237)、[density/contour/computation](../R/visualization.R#L9321-L9415)。

1. **科学问题**：给定二维位置上的一个数值变量，局部邻域是否呈现高值或低值聚集；该函数不要求输入来自通讯，也不证明空间因果。
2. **输入单位与轴**：三列按位置改名为 x、y、score；每行一个空间观察点。`n` 是每点 kNN 邻居数；`group.by/group.highlight` 是显示注释，不是统计分组。
3. **计算定义**：`knearneigh(k=n,longlat=FALSE)`→`knn2nb()`；gstar=TRUE 另 `include.self()`。`localG(score,nb2listw(...,style="B"))` 给标准化局部统计量，保存为 `statistic_G`，internals 拼入数据。density 对 Gi 做 `calculate_density(...,method="wkde")`；contour 对 Gi 做 kriging 图层；scoring 画原 score，但仍先执行 Gi 计算。
4. **视觉编码**：normal 的颜色为 Gi，设 p.value 后以严格 `p<p.value` 控制实/空心；3d 同类颜色/点形，z 恒为 `three.dim.z`。density 色为平滑后的量，contour 为插值面；它们不等同原始每点 Gi。非高亮透明度主要在 3d trace 分支使用。
5. **信息损失**：邻域压缩原始分数和邻居身份，平滑/插值再丢失局部细节；固定 z 并不提供三维样本位置。点形只保留阈值二分类。
6. **统计边界**：Gi/Gi* 使用不同的 self 包含定义；p 来自 localG internals 的局部检验，当前函数没有跨点/跨信号校正。高/低色值不自动等于显著；p=NULL 时没有显著性区分。3D hover 用 Accept/Reject 表述阈值分支，不能将其解读成接受“真实通讯存在”的假设（V-011）。
7. **前置计算及返回**：需有限且适合 kNN 的二维位置、合法 n、非退化 score 及 spdep；density/contour 另需相应 helper/绘图依赖。normal/scoring/density/contour 返回 ggplot，3d 返回 plotly，computation 返回含输入与 internals 的 data.frame。退化输入及依赖版本的 p 列名称尚须实测。
8. **判别反例（未执行）**：孤立高值点在“邻居不高”与“周围也高”两种场景应有不同 Gi；Gi* 纳入自身后可再变化。统一平移人工 z 不得改变 Gi/p。相同 Gi 图色但一组 p 未过阈值，必须从点形/图例辨别，不能把颜色视为显著性。

### 4.4 `spatialLeePlot`

依据：[权重选择及委托](../R/visualization.R#L8368-L8548)。

1. **科学问题**：两套细胞类别或连续属性是否在同一位置/邻域共现；不同 `lw.type` 问的是不同尺度的关联，不是有向通讯。
2. **输入单位与轴**：group.by1 默认 idents，group.by2 为 metadata 列或直接变量；每元素/矩阵行对应一个 cell/spot。图上轴是两组变量/类别，不是 sender 和 receiver。
3. **计算定义**：用户 lw 直接传给 `plotStatistics_Lee()`；Itself 传 coords，由 helper 建恒等权重；contact.range 用空间距离缓存的 adj.contact 并以 style B 转 listw；interaction.range 由 d.spatial 生成 P.spatial，置对角为 1 后 style W。缓存参数不匹配会重算距离。交换 x/y 坐标不把自权重变成邻近权重。
4. **视觉编码**：由 helper 选择 data/dot/heatmap；颜色/点大小为 Lee L。`cutoff.Lee` 是 L 的显示阈值，不是 p 值。行列转置关系见 4.5。
5. **信息损失**：全局类别对统计丢失哪些具体细胞贡献、局部异质性及通讯路由；不同邻域半径可能产生同样 L，但生物解释不同。
6. **统计边界**：wrapper/helper 都没有生成 Lee 显著性检验。Itself 在非退化输入下相当于同点相关的特殊情形（[INFERENCE]），不能叫“发现空间邻近效应”；contact/interaction 的权重含义须随结果披露。
7. **前置计算及返回**：需坐标、idents/meta，距离分支需有效空间尺度参数；不要求先做通讯推断。返回 helper 的 data.frame、ggplot 或 Heatmap 对象。给定 lw 时必须保证其顺序与变量一致，当前 wrapper 不代做名称对齐。
8. **判别反例（未执行）**：两个互斥细胞类各自成团，Itself 可负相关而接触邻域可能出现边界共现；固定标签、扩大邻域范围后检查 L 的变化，不得把范围参数只当作视觉平滑。

### 4.5 `plotStatistics_Lee`

依据：[输入转换、Lee 计算、阈值](../R/visualization.R#L8580-L8708)、[数据/点图](../R/visualization.R#L8711-L8743)、[heatmap 返回](../R/visualization.R#L8745-L8816)。

1. **科学问题**：以给定空间权重比较两套变量的全局双变量空间关联；不是 Gi 那种逐点热点，也不是 LR 显著性。
2. **输入单位与轴**：factor 转 one-hot（cell × category），numeric 转单列，matrix/Matrix 的行是 cell、列是变量。每一 x 变量与每一 y 变量组合计算一个 L；coords 与 lw 只能提供一个。
3. **计算定义**：coords 模式只按行数构造恒等矩阵，并不从坐标距离寻找邻居；每对变量调用 `spdep::lee(...,zero.policy=TRUE)`，只取 `$L`。默认 `cutoff=TRUE` 将 `L<=0` 置 NA。numeric 输入会覆盖请求的 type 为 dot。内部 `dfLee` 为 y×x，随后转置为 `mat`。
4. **视觉编码**：data 分支返回 y×x；heatmap 用 `mat`，故行 x 变量、列 y 变量；dot 使用转置后的表，横轴对应 y 输入、纵轴对应 x 输入。颜色与大小都是 L，尽管局部变量/Heatmap name 写着 proportion，不是组成比例。默认白色包含被 cutoff 隐藏的负关联和零。
5. **信息损失**：每类别对一个数，不能定位具体热点；默认隐藏负/零 L 丢失排斥信息。numeric 强制图形输出会丢失调用者期待的数据接口。
6. **统计边界**：**没有 p 值、置换检验或显著性校正**，`cutoff.Lee` 是效应量阈值；不能把“只剩正值”描述为显著正共现。与 Gi 的局部统计、p 点形必须严格区分。
7. **前置计算及返回**：需已按名称/顺序对齐的变量和 lw/coords；非 numeric 且 type=data 返回 data.frame，dot 返回 ggplot，heatmap 返回未由此函数 draw 的 Heatmap。`x.levels/y.levels` 在转置后对应的实际轴与参数名存在错位风险（V-012）。
8. **判别反例（未执行）**：用两个互斥类别让 L 为负，cutoff 开/关应暴露“隐藏”而不是“零关联”；用 2 个 x 类、3 个 y 类追踪每个轴名；请求 numeric 输入 type=data，应记录源码实际强制 dot，而不把它当成期望的新契约。

## 5. 方向场、流线与浏览器粒子

### 5.1 `netVisual_CommunField` / `netVisual_CommunFieldZoomIn`

依据：[基础场图](../R/visualization.R#L5477-L5740)、[ZoomIn 委托](../R/visualization.R#L6460-L6537)、[场生产者](../R/analysis.R#L4016-L4168)。

1. **科学问题**：每个 cell/spot 的某信号在局部主要向哪个方向输出/输入；ZoomIn 是同一场的总览与局部裁切，不是另一尺度重新推断。
2. **输入单位与轴**：场是 cell × 2 分量 × signaling；outgoing 向量位于 sender，incoming 位于 receiver，但两者指向都沿 sender→receiver 的方向。旧基础入口读 `slot$field`；当前生产者实际写 `slot$cell$field`，两者尚未切换一致（V-013）。
3. **计算定义**：生产者对每节点的有限正边按权重降序排列，以累计权重占比选边：`sparse=TRUE` 保留累计量严格 `< top × 总量` 的前缀，空选且 top>0 时取第一条；FALSE 在正边不超过两条时全取，否则用 `<=`。这不是 renderer 的 top 分位数。对选中边把“通讯权重 × 单位位移向量”求和，零长度位移不贡献向量。绘图只选 signaling，交换场与坐标两列，计算模长，去掉零模箭头；ZoomIn 调两次基础图，分别高亮网格窗口、裁切该窗口，不重算场。
4. **视觉编码**：箭头表示净方向与显示缩放后的模长，颜色通常为 sender/receiver 节点类别；背景为 cell/spot 标签。零向量点可有背景但无箭头。y 与 dy 的绘图变换必须一起观察；zoom anchor 是显示窗口参数。
5. **信息损失**：方向求和丢失每条边、邻居数与路径；等权相反方向会抵消，**有正通讯强度也可无箭头**。场不是传输总量，也不是组织中真实迁移路线。
6. **统计边界**：不产生新 p 值；top 与 min.mag 是选边/显示参数，非显著性检验。方向模长不能代替 row/column sum；缺箭头不证明无通讯。
7. **前置计算及返回**：先算 cell 网络及 `computeCommunField()`；需坐标、idents 同序。基础返回 ggplot，ZoomIn 返回 cowplot 组合的 ggplot 类对象；旧字段及 `my_as_sparse3Darray()` 路径仍须迁移，不能保证最终 schema 可直接运行。
8. **判别反例（未执行）**：中心点等权连向左右两点，保留全部边时总 outgoing 为正而净向量为零；只留右边后应向右。缩放/裁切不能把零净向量解释成无通讯，也不能把 window 外的贡献悄悄从场中重算删除。

### 5.2 `netVisual_CommunFieldGrid` / `netVisual_CommunFlow`

依据：[Grid 聚合与插值](../R/visualization.R#L5810-L6010)、[Grid 渲染](../R/visualization.R#L6027-L6090)、[Flow 聚合与显示参数](../R/visualization.R#L6153-L6370)、[Flow 渲染](../R/visualization.R#L6387-L6448)、[网格尺度和归属](../R/spatial.R#L187-L283)。

1. **科学问题**：局部方向场在较粗空间尺度上的主要方向如何组织；Grid 用采样箭头，Flow 用插值场上的流线，不是两种通讯推断模型。
2. **输入单位与轴**：读取 `slot$cell$field[[pattern]]` 的 SparseChatArray；一个 signaling、cell 的二维向量。网格单位来自显示坐标，若有 ratio 可报告物理 cellsize；仅把 sf CRS 设为 3857 不会赋予组织坐标真实地图单位。
3. **计算定义**：网格 cellsize 为用户值或最近邻距离来源，再乘 grid.resolution；按 sf 相交关系归格，可有边界点多格归属。`grid.field.mean` 名称虽含 mean，实际只支持**各分量 median 或 sum**。对 u/v 分别作 Barnes 插值，再在有效区域外置零。Grid 按 skip/min.mag 采样；Flow 用 metR streamline，默认 L、min.L、种子数来自网格尺度/占用格数。
4. **视觉编码**：两者背景颜色为细胞类别；Grid 的箭头宽度映射 `after_stat(norm_mag)`；Flow 的线宽映射插值向量模长并缩放到指定范围。流线弯曲是插值后积分显示，L/res/n/jitter 改变可见线段，不代表观测时间或真实细胞轨迹。
5. **信息损失**：cell→grid 分量汇总、插值、采样逐层损失局部路由和绝对尺度；同格相反方向抵消。median 与 sum 对细胞密度敏感程度不同；分量 median 也不等于选出一条真实 cell 向量。
6. **统计边界**：无新检验或物理速度标定；箭头粗细、流线密度不能当细胞数或显著性。当前 `sum` 对空格给零、随后 `na.omit()` 仍可保留这些格，可能与 median 的有效域不同（[INFERENCE]，V-013）；不能以“插值后连续”推断组织空洞真实有流。
7. **前置计算及返回**：需当前 cell-level field、匹配坐标/idents、网格与 oce/metR/sf 依赖；两者返回 ggplot，不修改场缓存。占用格、边界、退化坐标与全零场需要实图验证，不能把已存在新存储读取当作视觉验收。
8. **判别反例（未执行）**：同格放入两条反向等权向量并改变 cell 密度，比较 median/sum；在组织中央留空洞，确认插值与流线不被误读为跨空洞通讯。对已知水平场，检查 y/向量变换未转成垂直或反向，改变 L 只改显示段长。

### 5.3 `netVisual_CommunFlowWidget`

依据：[R 插值与有效域](../R/visualization_widget.R#L94-L182)、[payload 准备](../R/visualization_widget.R#L185-L376)、[widget 返回](../R/visualization_widget.R#L415-L487)、[浏览器采样](../inst/htmlwidgets/spatialcellchat-commun-flow.js#L86-L141)、[粒子步进](../inst/htmlwidgets/spatialcellchat-commun-flow.js#L203-L262)。

1. **科学问题**：交互观察一个**固定**通讯方向场的空间走向、局部净流与总强度；动画是方向编码，不是时间序列数据，也不是模拟出的生物动力学。
2. **输入单位与轴**：cell-level SparseChatArray 场、一个 LR/pathway、incoming/outgoing；坐标至少两列但只用前两列。严格要求场/坐标/idents 的 cell 名称与顺序相同，二维坐标唯一；粒子不是输入 cell。
3. **计算定义**：交换 x/y、向量分量并翻转显示 y/dy；网格分量 median/sum、仅占用格参与插值，Barnes 后保留 occupied mask、按 field.grid.max 抽样。全零场有专门零网格；缺 oce 的非零场明确报错。浏览器在有效四角网格内双线性采样，按 visual speed、帧间隔及最大模长归一化移动粒子。payload 同时保存 current=向量模长，mass=原 cell prob 层 row/column sum（若可用），两者不相等。
4. **视觉编码**：背景点色为细胞类；移动方向表示净向量场方向，粒子亮度/线宽随相对模长变化，粒子数、速度、trail 是显示参数。hover 区分 current 与可用的 mass；没有以秒为单位的生物过程速率。Pause/Play、速度、粒子数及缩放是交互控件。
5. **信息损失**：选边成场、网格、插值、抽样丢失路由；相反边抵消后 mass 可正而 current 为零。动画轨迹可经过非 cell 位置，不能解释为某个真实细胞移动或 ligand 沿该曲线运输。
6. **统计边界**：没有显著性检验；particle.speed 不改变通讯估计。mask 和无效值处理是显示约束；不能把静止区与未计算区合并解释。没有 measured z，且只显示二维投影。
7. **前置计算及返回**：需当前 `computeCommunField()` 结果、坐标/名字对齐与打包 JS/CSS；返回 `htmlwidgets::createWidget()` 创建的 `spatialcellchat-commun-flow` htmlwidget，不是 plotly、视频或时间步结果。R payload 成功不证明 browser 正确播放。
8. **判别反例（未执行）**：用有正 mass 但零 current 的对称场，hover 应保留两量区别、无伪造运动；改变粒子速度不改变 payload 通讯量。浏览器实际检查 Pause/Play、圈定有效域、缩放后点位、hover id、重绘和窗口 resize，不能以 widget class 代替验收。

## 6. Topic 消费图

以下只记录消费者的选择、汇总与显示；`identifyCellTopics()` 的中心性输入、NMF 分解与缩放见 [分析语义卡](ANALYSIS_SEMANTICS.md)。loading 不是 LR 通讯权重，也不是基因表达。

### 6.1 `spatialTopicPlot`

依据：[完整入口](../R/visualization.R#L7123-L7267)。

1. **科学问题**：某个 incoming/outgoing 通讯 topic 在哪些 cell/spot 上有较高 loading；不是该 topic 的具体通讯路径。
2. **输入单位与轴**：`slot$topic[[pattern]]$cell` 为 cell × topic，绘图转成 topic × cell；一个面板一个 topic，位置仍是 cell/spot。topic 标签由 `topics.key` 加序号重建，不是稳定的生物类别名。
3. **计算定义**：直接读取既有 cell loading；可将 `<=cutoff` 置 NA，不在此处重拟合、重归一化或重分配 topic。topics 子集仅用于显示。
4. **视觉编码**：颜色为该 topic loading，各面板独立色阶，点大小固定；灰色为 NA/阈值以下。共享图例外观不等于共享色域。
5. **信息损失**：不显示该 topic 的 LR 成分或 sender→receiver；每 topic 单独图不展示联合成员关系。相同色深跨 topic 不能当相同绝对贡献。
6. **统计边界**：loading/cutoff 没有 p 值意义；高 loading 是分解坐标，不证明此位置通讯显著。独立拟合的 Topic_1 不保证是同一因子。
7. **前置计算及返回**：先算对应 topic 与坐标；combine=TRUE 返回 patchwork，FALSE 返回 ggplot list。数据与坐标按位置绑定，需检查 cell 名称一致；不回写缓存。
8. **判别反例（未执行）**：把一个 topic 的全体 loading 乘 100，独立色域可令图外观近似不变，但数值图例必须变化；交换 topic 次序或 cell 行序，必须区分“重命名/错配”与真正空间热点改变。

### 6.2 `netVisual_TopicComposition`

依据：[输入和聚合](../R/visualization.R#L4766-L4918)、[归一化、删减与返回](../R/visualization.R#L4920-L5038)、[thresholdedMean](../R/modeling.R#L2896-L2903)。

1. **科学问题**：topic/其他属性在细胞群中的相对构成或关联，取决于输入分支和归一化轴；不是群间通讯强度。
2. **输入单位与轴**：默认 x 为 idents；SCC 的 y=`topicIn/topicOut` 从 **net** topic 取 cell loading，Seurat 从 reduction 取值；也接受 feature/category × cell 矩阵。常见结果 prop 行是 topic/y 类、列是 x 群；dot 把矩阵第一维映射为横轴、第二维为纵轴，与 heatmap 行/列方向不同。
3. **计算定义**：按离散群汇总时用 `thresholdedMean`：非零比例 `<trim` 则置零，否则取普通均值，**不是 trim 掉两端数值**。两矩阵分支用 `x.group %*% t(y.group)` 后转置，是加权重叠，不是群均值。`scale.rows=TRUE` 除行和；FALSE **除列和**，不是关闭缩放。NA 归一化结果置零。
4. **视觉编码**：heatmap 颜色、dot 颜色与大小为归一化后的比例；不同分母对应“某 topic 分布在哪些群”或“某群由哪些 topic 组成”。行缩放后两行相同色深只表示各自行内份额，不表示相同总 loading/通讯量。
5. **信息损失**：群汇总丢失 cell 异质性和位置，归一化丢失绝对量；均值再行归一化也不等于按群细胞数计算的全组织 topic 总量份额。矩阵乘积受 cell 数及输入缩放影响。
6. **统计边界**：没有组成差异检验、置信区间或 enrichment p；cutoff.prop 是显示选择阈值。当前 remove.isolate 由行最大值得到 idx，却删除列，不能把这个错轴行为写成科学规则（V-014）。
7. **前置计算及返回**：需 topic/矩阵与 cell 分组对齐；metadata 分支仍含 `@meta.data`，并不都适用于最终 SCC。dot 返回 ggplot，heatmap 返回 Heatmap；不返回未归一化 prop。矩阵、离散分组等分支应分别核对，不能仅用一个正方阵证明接口正确。
8. **判别反例（未执行）**：topic 原量两行分别 `(1,3)` 与 `(100,300)`，行归一化均为 `(0.25,0.75)`，图不应宣称两 topic 总强度相同；用 3 topics×2 groups 暴露错轴删除。令非零比例恰等于 trim，核对保留的是均值而不是截尾均值。

### 6.3 `netVisual_TopicSignaling`

依据：[输入读取](../R/visualization.R#L5072-L5148)、[截取、聚合、排序及返回](../R/visualization.R#L5151-L5294)。

1. **科学问题**：每个 topic 中哪些 LR、ligand 或 pathway 的 loading 居前，用来解释 topic 内容；不是 LR 显著性排序，也不是预测靶基因。
2. **输入单位与轴**：SCC 取 `slot$topic[[pattern]]$signaling`（signaling × topic），也尝试支持 Seurat/matrix。每面板一个 topic，横轴为显示特征的 rank，纵轴为绝对 loading。
3. **计算定义**：先每 topic 选 nfeatures（默认按绝对值，balanced 可取两端），**再**按 feature_cat 做展示聚合；ligand/pathway 是在选中子集内求 loading 均值，不是全体 LR 求和。再按绝对值降序，标注 nTop；LR 名字通过 DB 映射展示名。
4. **视觉编码**：曲线/点高度是 absolute loading；颜色仅区分符号，标签为 top feature。标准非负 topic 下没有正负效应的对照意义；外部有符号矩阵的“负”也不等于抑制通讯。
5. **信息损失**：先截 top 再聚合会漏掉同 ligand/pathway 的低 loading LR；均值忽略条目数，绝对值去掉符号大小的方向解释；没有空间分布。
6. **统计边界**：无 p、FDR 或稳定性选择；“marker”在这里是 loading 排名，不是正式 differential/enrichment 检验。独立 topic 拟合间不可直接按序号比较。
7. **前置计算及返回**：需 topic signaling loadings 与匹配 DB；combine=TRUE 返回 patchwork，FALSE 为 ggplot list。`species` 的 match.arg choices 未在签名列明，Seurat 分支引用未定义的 reduction/projected，外部 db 接受形状与后续 `$interaction` 使用存在不一致，列入 V-015；这些不能被“支持多个对象类”的文案掩盖。
8. **判别反例（未执行）**：同 pathway 有许多中等 LR、另一 pathway 只有一个高 LR，改变 nfeatures 可改变聚合均值排名；应明确这是选中子集的 loading，而不是 pathway 总通讯强度。使用 `slot.name="netP"` 的 pathway loadings 时也要检查 DB 映射是否仍错误按 LR 名处理。

## 7. Pattern 与网络嵌入消费图

### 7.1 `netAnalysis_river`

依据：[cell–pattern 部分](../R/visualization.R#L3978-L4101)、[pattern–signaling 部分](../R/visualization.R#L4103-L4146)、[返回组合](../R/visualization.R#L4197-L4203)。

1. **科学问题**：哪些群和信号共享某个 outgoing/incoming 通讯 pattern；这是因子结构摘要，不是观察到的物质流或 source→target 网络。
2. **输入单位与轴**：读 `slot$pattern[[pattern]]$pattern$cell/signaling` 的 Contribution 长表；轴为 CellGroup→Pattern 和 Pattern→Signaling。群 pattern 输入如何聚合和归一化见分析卡，不在本图重算 NMF。
3. **计算定义**：两张表各自将 `Contribution<cutoff` 置零；可筛群/信号，do.order 以贡献矩阵距离层次聚类排序。无 signaling 表时只画群–pattern；有表时分别作两张流图后并排，不计算群→pattern→signal 的联合三元权重。
4. **视觉编码**：流带宽为各自的 Contribution，颜色为群/pattern/信号类别；左右相接的 pattern 标签是因子对应，不保证两边宽度遵循同一质量守恒分母。
5. **信息损失**：阈值隐藏混合模式；每侧边缘 loading 不包含具体 sender–receiver 路由。外观看似贯通不能证明某群向某靶群经该 pathway 输送多少信号。
6. **统计边界**：cutoff 不对应 p 值；因子及其标签无显著性证明。图形“flow”与 rankNet 的总权重、空间场流线均不同。
7. **前置计算及返回**：需 `identifyCommunicationPatterns()` 的对应结果、ggalluvial；返回 ggplot 或 cowplot 组合图，不返回联合重构矩阵。pattern 参数选择和两侧类别配色需实际绘图核对。
8. **判别反例（未执行）**：两张网络具有相同 out/in 边缘和但不同路由，可产生相同 pattern 输入与 river；另令 cutoff 两侧改变某一小 loading，观察带消失而原网络未变化，避免把显示门控当新统计结论。

### 7.2 `netAnalysis_dot`

依据：[完整计算与编码](../R/visualization.R#L4229-L4306)。

1. **科学问题**：每群–信号组合最强地通过哪个保留 pattern 联系；不是完整 NMF 矩阵重构，也不是直接通讯概率点图。
2. **输入单位与轴**：群 pattern 的 cell 表、signaling 表及其原输入 data；横轴 signaling，纵轴群，无具体 receiver/sender 对。
3. **计算定义**：默认 cutoff=1/模式数；先分别门控两个 Contribution，再按 Pattern join、逐项相乘。原 data 中为零的群–信号组合被置零；按拼接 id 取 `top_n(1,Contribution)`，即**保留最大 pattern 乘积**，并列最大可保留多行。不是对 pattern 求和的 $WH$；零乘积改 NA。
4. **视觉编码**：点大小为最大乘积分数，颜色/填充为群类别，不是 p 值；坐标轴是群–信号而非 sender–receiver。
5. **信息损失**：舍弃次强模式，无法还原总贡献或混合度；输入零值屏蔽使图也不等同未经约束的因子乘积。并列点可能重叠。
6. **统计边界**：无 p 值、置信区间；小点不等于不显著。拼接 id 无分隔符会造成名称碰撞的风险（V-016），不能把错配得分当科学定义。
7. **前置计算及返回**：需群 pattern 的两侧表和 data，返回 ggplot。pathway.show/group.show 是显示子集；单个缺失 pathway 的补零条件与多个缺失时不同，需做名称/缺层验收。
8. **判别反例（未执行）**：同一组合两个 pattern 的乘积为 0.3、0.2，大小应反映当前最大 0.3 而不是 0.5；两个乘积同为 0.3 时观察重复点。用群/信号名 `AB+C` 与 `A+BC` 发现无分隔键碰撞，不能只比对输出行数。

### 7.3 `netVisual_embedding` / `netVisual_embeddingZoomIn` / `netVisual_embeddingPairwise` / `netVisual_embeddingPairwiseZoomIn`

依据：[single](../R/visualization.R#L4337-L4404)、[single zoom](../R/visualization.R#L4431-L4486)、[pairwise](../R/visualization.R#L4518-L4604)、[pairwise zoom](../R/visualization.R#L4633-L4727)。

1. **科学问题**：哪些信号网络在 functional/structural 相似性空间中邻近、属于同一聚类；pairwise 展示条件间的网络实例，而不是细胞的表达 UMAP。
2. **输入单位与轴**：一个点=一个 LR/pathway 网络；联合图一个点=layer–dataset 实例。x/y 是已计算 dr 的两维，无空间长度或细胞位置意义。functional/structural 的定义与嵌入拟合见分析卡。
3. **计算定义**：读取 similarity 的 dr、group、matrix，不在此处重拟合；默认从 similarity 列和等于 1 识别待去掉孤立层。另从 prob 各层求总权重，除所有显示层最大值作点大小。pairwise 拼接各 dataset 总量，并以 `layer--dataset` 名称拆分标签；zoom 只是按 cluster 分面，仍用同一全局相对总量。
4. **视觉编码**：位置是相似性嵌入、颜色是网络 cluster、大小是归一化层总权重；pairwise 形状是 dataset。ZoomIn 每 cluster 一面板，重新适配坐标范围，不能跨面板目测距离作定量比较。
5. **信息损失**：二维嵌入压缩网络差异，丢失边及 LR 组成；相对大小抹去绝对总量单位。相邻点不必总量相同，颜色相同不等于信号功能已被生物学证明一致。
6. **统计边界**：没有 cluster 显著性、距离置信区间或跨条件检验；top.label 是标注控制，不是显著性选择。不同独立嵌入坐标的旋转/布局不能当生物变化。
7. **前置计算及返回**：先算对应 similarity、embedding、clustering，且 dr/group/prob 的 layer 名称必须对齐。single/pairwise 返回 ggplot，zoom 返回 cowplot 组合图；pairwise 仍读旧 dataset 网络列表。源码按位置拼接数据，single top.label<1 的 data.label 分支、single zoom 未 match.arg(type) 及缺层删除需单独审阅（V-017）。
8. **判别反例（未执行）**：置换 prob 层序但保持带名 dr 不变，正确对应的点大小/标签应不变；同总权重不同路由的网络可点同大但位置不同。缺一层、孤立层及两个 dataset 重名 pathway 必须分别定位，不得用补零或行位置假装已对齐。

## 8. 后续视觉验收门禁（未执行）

本页八字段卡和数值探针只定义后续迁移的可证伪标准，不宣称当前全链可运行。每批只验证其受影响路径，但不能用一种布局替代其他返回/渲染分支。

1. **数据与轴**：从最终 11-slot 对象运行实际入口；核对 layer、cell/group、dataset 名称与顺序。包含名称置换、缺层、零/NA、p 等于阈值、缓存阈值不同的场景。`SparseChatArray` 必须逐层/稀疏归约，不能物化 cell×cell×layer 稠密张量；`marginSums(...,3)` 是跨层矩阵和，不是逐层总量。
2. **返回契约**：实际观察 list、ggplot、patchwork、recordedplot、Heatmap/draw 返回、SpatialCellChat 更新对象、plotly/htmlwidget 的差别；检查数据字段语义，不能只断言继承某个 class。新数据接口或返回变化需另行决策，不在本卡自动批准。
3. **静态图**：实际打开输出或渲染到图形设备；recordedplot 必须重放，ComplexHeatmap 必须 draw。核对点位、坐标朝向、发送/接收方向、自环、群中位数、边宽/色域、图例和裁切。正负抵消网络、全零场、空白组织区域、仅一个面板/多个面板都要观察。
4. **交互图**：实际浏览器核对 hover 的 cell id/数值、选中与反选、透明度、全部/无选中、色条、缩放及 resize。plotly 对象构造成功不能代替浏览器；widget 还需核对 Pause/Play、显示速度与粒子数量变化不修改科学量、零净向量与正 mass 的区别。
5. **科学对照**：表达不变而通讯改变；同边缘和不同路由；同总量不同 LR 组成；正强度零净向量；行缩放相同但原量不同；Gi 的局部 p 与 Lee 无 p 的区别；固定 z 的人工 3D。验收描述须同时指出“应该保持的科学量”和“只改变的显示量”。
6. **缺陷分流**：已知故障/错轴/误导统计不能作为黄金输出冻结。先分类为存储迁移、实现/文档不一致、科学/API 语义更正，得到裁定后再建立正确断言；本页和关联 Note 保持 proposed 边界。

## 9. 稳定问题表：证据、消费者与未决边界

V-* 是本页稳定引用标识，不表示已修复，也不是已批准的新默认值。证据为源码事实；影响中的风险若未运行验证则标 `[INFERENCE]`。实施时不得把任一条旧错误当作兼容要求。

| ID | 证据与当前事实 | 受影响消费者 | 建议方向（未批准） | 开放决策/验收边界 |
|---|---|---|---|---|
| V-001 | [表达旧槽](../R/visualization.R#L452-L454)、[aggregate 旧 prob](../R/visualization.R#L1459-L1472)、[heatmap 旧 options/mode 与 dataset](../R/visualization.R#L2764-L2793)、[pairwise 旧列表](../R/visualization.R#L4527-L4534)。读取路径混合新旧结构，不代表 final schema 可直接运行。 | 表达、aggregate/individual、bubble、heatmap、diff、embedding、部分 topic。 | 每批按三个独立轴迁移所有读取与名称逻辑；复用当前访问器/稀疏结果，不加旧槽 shim。 | dataset 存档与 joint 入口、可选 p 值、稀疏输出逐项核定；中心性实际落点不能擅自改成 group$centr。 |
| V-002 | [aggregate `p>`](../R/visualization.R#L1469-L1476)、[pathway `p>=`](../R/modeling.R#L1691-L1698)、[bubble table `p>=`](../R/analysis.R#L2431-L2445)。实时 LR、生产缓存和表过滤不同。 | aggregate/individual、heatmap、bubble、所有 pathway 图。 | 明示筛选时点及阈值等号，做同输入数值/图形对照。 | 统一阈值属于科学语义更正；本轮不预先选择保留 `>` 还是 `>=`，也不批准 aggregate 改读 netP。 |
| V-003 | [extractEnrichedLR_internal](../R/analysis.R#L913-L946) 的 gated prob 不决定最终 LR.sig 交集；空交集也未清空后续取基因的 pairLR。 | spatialFeaturePlot 的 signaling、plotGeneExpression 及 LR 选择消费链。 | 分开“缓存富集选择”和“本次实时筛选”，裁定空 LR 与 geneLR 返回一致性。 | 修改缓存/阈值语义与修复空结果要分开归因；不把空 LR 非空基因当作正确富集证据。 |
| V-004 | [连续/二值 complex 差异](../R/visualization.R#L518-L523)、[binary complex](../R/visualization.R#L621-L627)、[plotly 直接匹配](../R/visualization.R#L7367-L7387)、[ncol=1 仅首图](../R/visualization.R#L574-L577)、[plotly 固定 trace 下标](../R/visualization.R#L7878-L7893)。 | spatialFeaturePlot、plotly_spatialLRpairPlot、plotFeatures、Shiny adapter。 | 明确“亚基表达/complex 量/缺 partner”三者，确保 feature 面板完整、trace 与 cell id 一致。 | 缺 complex 应报错、展开还是计算复合物需选择；修复不能悄悄把表达换成通讯量。实际浏览器核对全/空高亮。 |
| V-005 | [meta 改名不重排](../R/visualization.R#L6640-L6647)、[bar matrix 与最后 df.use](../R/visualization.R#L6918-L6936)。[INFERENCE] 可能错配分组或展示前一个基因/失败。 | plotGeneExpression、barPlot 及相关教程。 | 以名称对齐为不变量，明确 assay 层及多基因返回内容。 | 先确认来源矩阵是否归一化和 Seurat 版本层行为；不是借文档更换表达统计量。 |
| V-006 | [individual spatial 未定义 labels/idents.use](../R/visualization.R#L1714-L1719)、[hierarchy2 外部 object](../R/visualization.R#L2004-L2024)、[spatial weight 条件](../R/visualization.R#L1516-L1528)；chord 与其他 renderer 不共享 top 过滤。 | aggregate/individual 四布局及直接 renderer。 | 逐布局确认局部参数来源、过滤时点、返回类及向量权重处理，不靠调用者全局变量。 | 是否统一 top 支持、返回 list 与组合图是 API 决定；群中位数路线解释必须保留其摘要性质。 |
| V-007 | [diff 孤立判定](../R/visualization.R#L2616-L2624)、[heatmap 零总和拒绝](../R/visualization.R#L2891-L2897)。有符号边缘和/全局和可抵消。 | diffInteraction、比较 heatmap。 | 区分“没有非零边”和“净变化为零”，名称对齐后再比较。 | 修复空网络判定不能变成新的差异显著性检验；正确断言需包含正负相抵但非空矩阵。 |
| V-008 | [heatmap 重复删点](../R/visualization.R#L2834-L2843)、[按列构建行注释](../R/visualization.R#L2880-L2889)、[列聚类用错参数](../R/visualization.R#L2898-L2901)。 | netVisual_heatmap 子集、remove.isolate、cluster.cols。 | 独立行列名称和索引，修复重复删除与参数失效。 | 非方形子集、单行/单列与全零场景都需实际 draw；侧栏的绝对变化和文案需核定。 |
| V-009 | [bubble p 档与变换](../R/visualization.R#L3139-L3150)、[固定 pairLR.use 重排](../R/visualization.R#L3298-L3300)、[p 图例等号](../R/visualization.R#L3345-L3346)。非标准 prob 分支还读未限定 `prob.original`。 | bubble single/comparison、return.data。 | 原量、显示量、p 档显式区分，保留无 p 语义，修正名称/边界显示。 | prob 非概率分布值是否允许、如何显示 `>=1`/NA、默认 LR 排序需裁定；不靠硬编码替代科学输入。 |
| V-010 | [scoring group 位置复制](../R/visualization.R#L751-L768)、[merge 均值](../R/visualization.R#L875-L891)。 | spatialVisual_scoring。 | 按群名对齐并明确 group 值复制，merge 标注量纲与分母。 | 是否允许任意 measure 合并、不同尺度 measure 是否可均值化需明确；不能只让图成功返回。 |
| V-011 | [Gi 默认先二值化](../R/visualization.R#L8900-L8907)、[固定保存 Gi p 键](../R/visualization.R#L9000-L9003)、[Gi/Gi* 不同键](../R/visualization.R#L9114-L9139)、[hover Accept/Reject](../R/visualization.R#L9215-L9223)。 | spatialGiPlot、plotStatistics_Gi 与 meta 下游。 | 明示 binary 输入含义，按实际 Gi/Gi* 名称保存；图例/hover 用明确的零假设和阈值措辞。 | binary 默认、局部 p 多重校正及超过两个面板支持属不同决策；未获批准前不改统计方法。 |
| V-012 | [Lee 仅取 L](../R/visualization.R#L8660-L8676)、[cutoff 与转置](../R/visualization.R#L8688-L8703)、[numeric 强制 dot](../R/visualization.R#L8610-L8632)。 | spatialLeePlot、plotStatistics_Lee、type=data 消费方。 | 清楚区分效应量与显著性，核对 x/y、矩阵行列、levels 和返回契约。 | 是否保留默认隐藏负值、是否允许 numeric data 返回、是否新增检验均须单独裁定；不能把 cutoff 叫显著性。 |
| V-013 | [旧 field 读取](../R/visualization.R#L5515-L5516)、[新生产者](../R/analysis.R#L4164-L4168)、[Grid 空格和 mask](../R/visualization.R#L5927-L6010)、[widget 只用 occupied](../R/visualization_widget.R#L290-L311)、[多格归属](../R/spatial.R#L271-L283)。 | CommunField/ZoomIn、FieldGrid、Flow、FlowWidget。 | 完成基础场旧路径切换；分别核对分量 median/sum、空格有效域与边界点重复贡献。 | 是否统一静态图与 widget mask、边界归属策略、零向量显示是独立选择；不能将净流为零改成无通讯。 |
| V-014 | [composition thresholdedMean 汇总](../R/visualization.R#L4863-L4873)、[归一化与错轴删除](../R/visualization.R#L4922-L4932)。 | TopicComposition 的群与矩阵分支。 | 说明行/列比例条件、非零支持门控，修复行索引删列。 | `scale.rows=FALSE` 当前仍归一化，若需原量视图须 API 决策；均值比例不能更名成组织总量。 |
| V-015 | [TopicSignaling 输入分支](../R/visualization.R#L5094-L5143)、[先 top 再均值聚合](../R/visualization.R#L5205-L5239)。 | SCC/Seurat/matrix topic signaling、netP loadings。 | 明确受支持对象/DB 形状和 LR/pathway 名称域，去掉未定义变量路径。 | 改成全量 pathway 总 loading 会改变问题，不能当机械修复；species/default 和 netP 支持需裁定。 |
| V-016 | [dot 阈值乘积与 top_n](../R/visualization.R#L4248-L4282)。最大值而非 WH 求和、并列多行、无分隔键。 | netAnalysis_dot 及 pattern 解读。 | 文档明确最大模式分数；键按两列组合且显式处理并列/缺层。 | 是否保留 max 或改 sum 是科学/API 决策；仅修键碰撞不能悄悄重定义分数。 |
| V-017 | [single 按位置拼接](../R/visualization.R#L4344-L4358)、[single top.label 分支](../R/visualization.R#L4373-L4388)、[zoom type](../R/visualization.R#L4431-L4449)、[pairwise 拼接](../R/visualization.R#L4557-L4571)。 | 四种 embedding 图。 | 使用显式 layer–dataset 名称对齐 dr/group/total，处理缺层/孤立层/单层；核对标注和 type 默认分支。 | 缺层是缺失还是零、独立拟合还是 joint 拟合必须披露；不从图形位置推断数据集效应。 |
| V-018 | [人工 z](../R/visualization.R#L7672-L7684)、[Gi 固定 z](../R/visualization.R#L9230-L9237)、[独立 feature 色阶](../R/visualization.R#L553-L565)、[topic 色阶](../R/visualization.R#L7207-L7228)、[粒子步长](../inst/htmlwidgets/spatialcellchat-commun-flow.js#L222-L233)。 | 表达、topic、Gi、空间场及浏览器展示。 | 图例/说明区分原量、相对显示量、人工坐标和视觉速度；用真实图和 browser 检查。 | 是否引入共享色域、真实 3D 或物理时间模型均需另立范围；本次不批准，也不把构造对象成功当验收。 |

## 覆盖与本次交付边界

本页覆盖 23 张八字段卡：表达三入口及 helper；aggregate/individual 和四类 renderer；bubble/heatmap/diffInteraction；空间 centrality scoring；Gi wrapper/helper 与 Lee wrapper/helper；基础场/ZoomIn、Grid/Flow、htmlwidget；topic spatial/composition/signaling；river/dot 与四种 embedding 消费者。Shiny LR 选择适配器列作关联消费者，不扩展其生命周期改造。

本次只有语义文档，不修改任何 R 实现、签名、默认值、tests、NAMESPACE、man、依赖或已实施 Note。当前最重要风险是旧存储读取、阈值/缓存不等价、名称与轴错配、实际返回分支差异及未完成的图形/浏览器验证。数值机制证据见分析文档第 8 节；上述判别反例和视觉门禁仍必须在真正迁移批次执行。
