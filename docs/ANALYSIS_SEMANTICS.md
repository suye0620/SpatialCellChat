# 下游分析逐函数语义卡

## 阅读边界与证据级别

本页是迁移前的科学语义审阅材料，不是已完成迁移或已验证可运行的 API 手册。对应决策为 [下游通讯分析与可视化语义契约（proposed）](../.agents/notes/proposed/2026-10-20-downstream-communication-semantics-contract.md)；图形端的完整编码、空间解释与交互边界见 [可视化语义卡](VISUALIZATION_SEMANTICS.md)。父迁移审计仍为 proposed。本页不批准更改统计方法、签名、默认值或缓存位置。

- **文档意图**：源码 roxygen 声明的用途；与执行语句不同时分别列出。
- **源码事实**：下列具名函数及源码行范围直接支持的读取、计算和返回行为。旧字段路径仍按源码记录，不能据此推定最终 11-slot 对象可以直接运行这些入口。
- **[INFERENCE]**：静态推导的解释、风险或预期区分，不是运行结果。
- **建议/未决**：问题表中的方向需维护者裁定；已知缺陷不应成为必须保持的黄金输出。
- **2026-10-21 调用方增补**：[PP1 下游重构与测试方案](PP1_DOWNSTREAM_REFACTORING_PLAN.md) 将本页语义落实到 `visium_PP1.R` 的实际 IO、默认参数、metadata/CSV/Seurat 消费链；包含 Visium 混合权重与硬标签聚合的区别、现有 topic/community 生产者与 Seurat 桥接契约，以及新的真实输入/最终 schema 探针。该增补不改变 A-* 的未决状态。
- 每张卡固定八项。第八项均为**拟议的迁移验收场景，尚未作为相应入口的验收检查执行**；不是测试结果，也不提供看似可直接运行的占位代码。2026-10-20 的孤立数值/公式探针单列于第 8 节；2026-10-21 的 PP1 子链与 Seurat 探针见上述调用方方案的运行证据节，均不代表整张卡或函数链已经验收。

## 轴、量纲与共同前提

令 $P^{d,r}_{ij\ell}$ 为数据集 $d$、节点分辨率 $r$ 下的非负通讯权重：行 $i$ 是 sender，列 $j$ 是 receiver，层 $\ell$ 是 LR 或 pathway。`net/netP`、`cell/group`、dataset 是三个独立选择轴；$P$ 不是表达矩阵，也不能仅因变量名 `prob` 就解释成互斥事件的概率分布。

| 量 | 轴与定义 | 不能替代的量 |
|---|---|---|
| 层总强度 $F_\ell=\sum_{ij}P_{ij\ell}$ | 每个 LR/pathway 一个值；`rankNet(weight)` 的 raw information flow | 不是显著性、熵或独立生物学重复数 |
| 正边数 $C_\ell=\sum_{ij}1[P_{ij\ell}>0]$ | 每层节点对的数量 | 不等于稀疏存储条目数；pathway 层计数也不是其 LR 数 |
| outgoing $O_{i\ell}=\sum_jP_{ij\ell}$ | sender × layer | 丢失 receiver 路由 |
| incoming $I_{j\ell}=\sum_iP_{ij\ell}$ | receiver × layer | 丢失 sender 路由 |
| 中心性 $C_{mi\ell}$ | measure × node × layer | out/in strength、PageRank、betweenness 等不是同一单位 |
| 群模式输入 $X_{g\ell}$ | 群 × layer，边缘和后按每层最大值缩放 | 不是 cell topic 的原始输入，也不是整个邻接张量 |
| cell topic 输入 $A_{\ell c}$ | layer × cell，由 `centr.cell` 的 outdeg/indeg 转置取得 | 不是表达 topic，也不是群平均模式 |
| 相似性 $S_{ab}$ 与嵌入 $Y_{a q}$ | $a,b$ 是**信号网络**；联合分析中是 layer–dataset 实例 | UMAP 点不是细胞、细胞群或空间位置 |

**稀疏归约**：[marginSums.SparseChatArray](../R/SpatialCellChat_class.R#L138-L163) 的 `margin=1/2` 分别给节点 × 层的 outgoing/incoming，`margin=3` 给跨层求和后的节点 × 节点矩阵；后者不等价于 `apply(P, 3, sum)` 的逐层总量。迁移不能为复刻旧 `apply` 而物化 cell × cell × layer 稠密张量。

**cell→group 不是任意求和**：[computeAvgCommunProb_LR_Sum / Avg](../R/modeling.R#L1378-L1500) 先按群聚合，再乘表达比例及 sender/receiver 支持门控；`Avg` 的分母来自将稀疏层已有存储项置 1 后聚合的计数，不是所有可能细胞对，也不是群大小简单乘积。门控中的 `min.cells.sr` 实际用关联条目计数聚合，不能仅按参数文案解释成去重细胞数。这里记录生产者定义，不在分析/绘图层另造分母；显式存储零与真实正边是否等价须随生产者审阅。

**Visium 混合比例是另一生产契约**：[Visium helpers](../R/modeling.R#L957-L1007) 用 spot×type 比例 D 计算 `t(D) %*% P %*% D`；avg 的分母为 `t(B) %*% E %*% B`，其中 B=1[D>0]，E 将 P 的已存储项置 1。它不是把每个 spot 指派给一个硬标签再调用普通 avg。真实 PP1 的群/比例不一致、截断后的质量与群列置换定义见 [PP1 §2.2](PP1_DOWNSTREAM_REFACTORING_PLAN.md#22-群定义混合比例与-visium-专属生产者)；下游迁移不能静默更换该科学问题（A-18）。

**阈值与缓存**：[computeCommunProbPathway](../R/modeling.R#L1680-L1767) 的 group 分支在生产时清除 `pval >= thresh`，cell 分支聚合已有 cell LR 层而不应用同一 group p 值门控。`netP` 不是任意新阈值下实时重算的结果。LR 实时聚合图及其 `p > thresh` 边界见 [可视化语义卡](VISUALIZATION_SEMANTICS.md)。`params()` 实际读取 `misc$.param`，**不是** `misc$.mode`（[访问器源码](../R/SpatialCellChat_class.R#L1526-L1533)）；迁移中的 mode、样本存档和 joint 网络必须分别处理。

## 1. 筛选、贡献与信息流

### 1.1 `netAnalysis_contribution`

依据：[R/analysis.R:237–397](../R/analysis.R#L237-L397)。

1. **科学问题**：文档意图是某 pathway 的各 LR 对总体信号贡献多少；实际是选中 LR 张量经全局 min–max 变换后的相对份额，不总等于原始通讯强度占比。
2. **输入单位与轴**：从 `LR$LRsig` 按 pathway 找 LR，再与网络层取交集；`do.group=TRUE` 读旧 `net$prob/pval`，FALSE 读旧 `prob.cell`。`vertex.receiver` 选择 receiver 轴，不是 sender。
3. **计算定义**：group 清除 `p>thresh`，cell 无该检验；剔除总和为零的候选 LR。对保留张量计算 $Q=(P-\min P)/(\max P-\min P)$，再取 $c_\ell=\sum_{ij}Q_{ij\ell}/\sum_{ij\ell}Q_{ij\ell}$。receiver 分区的两个分子仍共用**全体**分母，不在每区重归一化。全常数时分母为零，后续 NA 份额置零。
4. **视觉编码**：条长表示 $c_\ell$；无 receiver 分区时按贡献降序，补齐用的零条只为布局。分区时显示 All、Hierarchy1、Hierarchy2 三面板，不是三次独立模型。
5. **信息损失**：丢失 LR 内的 source–target 路由和绝对总量；减去最小值还可能改变 LR 间份额。分区图保留“该部分占全体”的意义，而非“该部分内部组成”。
6. **统计边界**：只消费既有 p 值，没有 LR 贡献的置信区间或差异检验；零贡献可来自门控、min–max 或退化归一化，不能一律解释为无生物学作用。
7. **前置计算及返回**：需 LR 注释与对应分辨率的 LR 权重，group 还使用 p 值。通常返回图；`return.data=TRUE` 返回 `LR.contribution` 和 `gg.obj`。源码在分区时反复覆盖 `df`，最终表是 **Hierarchy2** 的正贡献，不是 All/两区的完整数据；未分区才是总贡献表（A-03）。
8. **判别反例（未执行）**：两 LR 层分别有两个权重 `(1,2)` 与 `(3,4)`：raw 份额为 `3:7`，源码 min–max 后为 `1:5`。另设 receiver 两区各贡献一半，区内条长之和应为 `0.5` 而非 `1`；据此区分分母和返回表指向。

### 1.2 `rankNet`：single/comparison × weight/count × group/cell

依据：[参数文案与 single 分支](../R/analysis.R#L401-L554)、[comparison、统计及返回](../R/analysis.R#L555-L831)。

1. **科学问题**：排序一个或多个数据集的 LR/pathway 总强度或正边数；“information flow”在这里是权重和，不是信息论流量。比较分支还可请求边配对的 Wilcoxon 检验。
2. **输入单位与轴**：`net` 层是 LR，`netP` 层是 pathway；group 读旧 `prob/pval` 并清除 `p>thresh`，cell 读旧 `prob.cell`。single 读一个网络；comparison 把所选 slot 当 dataset 列表。source/target 筛选先把非选中边置零。
3. **计算定义**：`weight` 的 raw 值为 $F_\ell$；`count` 先转为 $1[P>0]$ 再求和。group weight 的 `contribution.scaled` 是 $-1/\log F_\ell$，无穷/负值分支再按 raw 排序分配显示值；cell weight 和 count 的 scaled 等于 raw。comparison 按层名并集补零，relative 是后续条件相对第一个条件的 **raw 比值**，经 `format(...,digits=1)` 转为数值；不是 scaled 比值。
4. **视觉编码**：single 固定画 scaled，源码不按 `show.raw` 切换；comparison 非堆叠由 `show.raw` 切换 raw/scaled，堆叠用 raw 经 `position="fill"` 变成条件份额。颜色区分 dataset；标签颜色依据 relative、`tol` 及可选 p 值。堆叠 50% 线只对两个条件的等份有直观意义。
5. **信息损失**：同一总强度/边数可对应不同路由；补零混合“缺层”和“确实零值”；scaled 条高不保留 raw 倍数关系。source/target 子集不能自动保证检验也只比较同一子集。
6. **统计边界**：仅两条件且 `do.stat=TRUE` 时检验；单位是展平的**网络边**，不是独立生物学样本。`prob.list` 在 source/target 筛选前保存；删除两列总和为零的行后做 paired Wilcoxon，不按节点名配对，也未在该函数校正多重比较。少于等于三条可用边、有缺层 NA 或检验给出 NA 时源码写 `p=0`，这是待纠正问题，不是显著性证据。`j` 还在本段赋值前被用于尺寸检查，不能承诺统计分支可运行（A-04）。
7. **前置计算及返回**：需要一致的权重生产定义和已核定 dataset 入口；返回 ggplot，或 `signaling.contribution` 表与 `gg.obj`；表区分 `contribution`、`contribution.scaled`、comparison relative 及请求的 `pvalues`。返回前 raw/scaled 取绝对值。当前 dataset 嵌套读取不能直接替代最终对象的样本存档（A-01）。
8. **判别反例（未执行）**：层 A 一条权重 `0.8`、层 B 四条各 `0.1`，weight 排 A 前而 count 排 B 前；同时核对 raw `0.8:0.4` 与 scaled 不保持 `2:1`。比较只改未选中的 receiver 边，可区分显示筛选与检验筛选；三条有效边或缺一层必须触发“无法检验”的裁定，而非把旧 `p=0` 当预期显著结果。

### 1.3 `extractEnrichedLR` / `extractEnrichedLR_internal`

依据：[R/analysis.R:836–946](../R/analysis.R#L836-L946)。

1. **科学问题**：文档意图是提取 pathway 内显著 LR 及组成基因，供表达背景检查；不是对选中细胞重新检验通讯。
2. **输入单位与轴**：候选从 `LR$LRsig` 按 pathway 查找，再按 group/cell 网络 LR 层及 `LR.sig`/`LR.sig.cell` 缓存筛选。基因通过 DB complex/geneInfo 展开；即使 `enriched.only=FALSE`，候选仍来自 `LR$LRsig`，不是无条件列出数据库全部 LR。
3. **计算定义**：`enriched.only=TRUE` 的 group 分支虽计算 `prob[pval>thresh] <- 0`，最终 LR 名称取 `intersect(pairLR.name, LR.sig)`，**没有用本次门控后 prob 重选**。cell 同样用缓存。wrapper 对每 pathway 调 helper；旧 merged 分支对 dataset 取并集，多个 pathway 结果再串接。
4. **视觉编码**：自身不画图；`plotGeneExpression` 等消费返回基因画表达，颜色/点大小属于表达量或表达比例，不能改称本次通讯显著性。
5. **信息损失**：丢失具体 sender–receiver 边、权重、p 值及 dataset 来源；并集掩盖条件特异性。复合物拆基因也不证明所有亚基同时形成有效复合物。
6. **统计边界**：`thresh` 参数不保证改变返回 LR；缓存生成时点决定 enrichment 选择。空交集分支只发消息而未清空 `pairLR`，随后仍可提取候选基因；[INFERENCE] 空 LR 与非空基因不应被解释为已证实的信号（A-02）。
7. **前置计算及返回**：enriched 路径要求 `aggregateNet` 生成缓存；helper 返回 `list(geneLR, pairLR.name.use)`，wrapper 文档约定返回 LR data.frame，或 `list(pairLR, geneLR)`。`enriched.only=FALSE` 的 helper 初始 `pairLR.name.use` 是 data.frame，而 TRUE 分支变为向量；wrapper 的字符转换须单独核定，不能仅凭文案保证内容一致。mode 仍读旧 `options$mode`。
8. **判别反例（未执行）**：缓存在 `0.05` 生成，某 LR 的唯一边 `p=0.04`；本次 `thresh=0.01`。源码缓存选择仍可能返回该 LR，实时重筛应排除它。再构造缓存交集为空，核对 LR 与基因是否同为空，而不只检查返回对象形状。

### 1.4 `subsetCommunication` / `subsetCommunication_internal`

依据：[wrapper](../R/analysis.R#L2296-L2370)、[稀疏层/表输入、门控及过滤](../R/analysis.R#L2405-L2545)、[pathway 汇总及返回](../R/analysis.R#L2549-L2604)。

1. **科学问题**：把指定群间、LR/pathway 或 DEG 条件的推断通讯整理为可追溯表；不是重算通讯模型。
2. **输入单位与轴**：内部非表输入读 `net$group$prob`，逐层展开 source × target；没有 cell 分辨率选择参数。也可输入已有通讯 data.frame。`pairLR.use` 为 pathway 列时 wrapper 改选 `netP`；它与 `signaling` 互斥。数值群索引按 `levels(idents)` 解释。
3. **计算定义**：LR 层有 p 值时清除 `p>=thresh`，只输出正权重；无 p 值列时输出 `NA`。netP 读已生产的 pathway 层，不再门控。**已有表输入绕过这一门控和正边展开**；随后按注释、dataset、DEG 指标及 source/target 过滤。DEG p 值用 `<=`；正/零 logFC 阈值取 `>=`，负阈值取 `<=`，pct 用 `>=`。netP 最终按 source–target–pathway 求权重和，p 值取算术均值。
4. **视觉编码**：自身不画图；输出供 bubble、chord、差异筛选等消费。LR 的 `prob/pval` 与 DEG 的 `ligand.* / receptor.*` 是不同统计来源；netP `pval=NA` 不能画成显著零 p 值。
5. **信息损失**：零边和被门控边不在表中；pathway 汇总丢失 LR 组成及其证据差异；已有表与网络输入不具有天然相同筛选语义。pathway 行用首个匹配 LR 提供注释的中间步骤不赋予 pathway 单一 ligand/receptor 身份。
6. **统计边界**：没有新检验。直接 netP 输入的 NA 表示无独立 pathway p 值；若从已有 LR 表聚合得到 p 值均值，它不是合法的联合 pathway 检验（A-16）。缺失 p 值、缺失 DEG 与零值不应混同。
7. **前置计算及返回**：LR group 需 `computeAvgCommunProb`，pathway 需 `computeCommunProbPathway`；DEG 过滤需相应映射列。single 返回 data.frame；文档意图 merged 返回逐 dataset 列表，但源码 merged 自动分支固定遍历 `object@net`、仍按旧样本嵌套解释，非最终样本存档契约。wrapper 即使传 `net` 仍访问 object 的 mode/LR/idents，不能把 `object=NULL` 视为已支持的独立表接口。
8. **判别反例（未执行）**：同一边 `p==0.05`、权重正，LR 网络输入在 `thresh=0.05` 被排除，而预制表输入不因这个参数自动删除。pathway 缓存生产阈值 `0.05` 后再请求 `0.01`，不应假称已重算；两 LR 的 `p=0.001,0.9` 求均值也不能变成 pathway 显著性的证据。

## 2. 中心性与角色

### 2.1 `netAnalysis_computeCentrality` / `computeCentralityLocal`

依据：[local 指标](../R/analysis.R#L1-L58)、[生产者读取、门控、返回位置](../R/analysis.R#L63-L232)。

1. **科学问题**：描述每个信号网络中节点的发送、接收、连接与拓扑重要性；“Mediator/Influencer”是网络角色标签，不是已验证的因果调控身份。
2. **输入单位与轴**：单层是有向加权 node × node 邻接矩阵，结果为 measure × node × layer。对象可选 group/cell、net/netP；也可直接输入二维或三维 `net`，显式 net 分支不应用对象路径的 p 值门控。
3. **计算定义**：group LR 先把零权重的 p 值设 1，再清除 `p>=thresh`；group netP 消费已门控层，cell 消费 `cell$prob`。local 的 unweighted out/in 是正边计数，weighted out/in 是 igraph strength；另有 PageRank。全量分支增加 hub、authority、eigen、betweenness、flowbet、info；betweenness 前把图边权换成 `1/weight` 距离，flowbet/info 则用原矩阵调用 sna。
4. **视觉编码**：自身不画图；角色 network/heatmap/scatter 与 cell topic 消费各指标。颜色或位置应标明具体 measure，不把所有中心性都称为“通讯强度”。
5. **信息损失**：out/in 仅保留边缘和；相同 strength 不保证相同 betweenness 或路径。把各层非线性中心性求和通常也不等于在跨层合并邻接图上重算中心性。
6. **统计边界**：这些是条件于推断网络的描述量，不产生新的显著性。local 捕获 flowbet/info 错误后返回零，会把失败与真实零混合；这是 A-06 待裁定行为，不是应固定的科学约束。
7. **前置计算及返回**：默认 producer `do.group=FALSE, degree.only=TRUE` 只生成五项；local 自身默认 `degree.only=FALSE`。对象路径返回更新对象，**group 写 slot 顶层 `centr`，cell 写 `centr.cell`，不是 `group$centr`**；无 object 返回三轴数组。默认角色 network 请求 flowbet/info，须全量 group 中心性，不能直接消费 producer 默认结果。
8. **判别反例（未执行）**：三节点顺时针与逆时针等权环具有相同 out/in strength，却有不同具名路由；证明角色边缘和不还原网络。另一场景仅一个 LR 边 `p==thresh`，group net 应清除而显式 net 输入不自动清除；另把无法计算与真实零中心性分开审阅。

### 2.2 `netAnalysis_signalingRole_network`

依据：[R/analysis.R:2608–2673](../R/analysis.R#L2608-L2673)。

1. **科学问题**：在每个给定信号网络内，哪些群在某种中心性上相对突出？
2. **输入单位与轴**：读 `slot$centr[measure,node,signaling]`；逐层形成 measure × group 矩阵，只用 group 结果，不读 `centr.cell`。
3. **计算定义**：每个信号、每种 measure 单独除以该行最大值；默认 outdeg/indeg/flowbet/info 改名 Sender/Receiver/Mediator/Influencer。各行是相对最大值缩放，不是行和为一。
4. **视觉编码**：热图行是角色、列是群，颜色是行内相对 importance；每个 signaling 单独 draw。白色还可来自零行除零后的缺失。
5. **信息损失**：不同 measure、不同信号的绝对量级不可用同色比较；`(1,2)` 与 `(100,200)` 会有相同色块。
6. **统计边界**：不做角色显著性检验。全零或计算失败的行不能因显示白色就被解释为已证明无作用。
7. **前置计算及返回**：先生成含所需 measure 的 group `centr`，默认 measure 需要 `degree.only=FALSE`。函数循环 `draw`，没有显式返回热图列表；应按绘图副作用审阅，不承诺 ggplot 返回。
8. **判别反例（未执行）**：两个 pathway 的 outdeg 分别为 `(1,2)`、`(100,200)`，热图应同色而 raw 中心性相差百倍；同加一个全零 info 行，检查是否与合法最小值区分。

### 2.3 `netAnalysis_signalingRole_scatter`

依据：[R/analysis.R:2676–2799](../R/analysis.R#L2676-L2799)。

1. **科学问题**：各群在选中信号集合内总体更像发送端还是接收端？
2. **输入单位与轴**：一个点是一个群；从 `centr` 选 x/y measure，再沿 layer 求和。指定 signaling 先与旧 `object@netP$pathways` 相交，即使 `slot.name="net"` 也使用 pathway 名称集合（A-01）。
3. **计算定义**：$x_g=\sum_\ell C_{x,g,\ell}$、$y_g=\sum_\ell C_{y,g,\ell}$。默认是总 out/in strength；点大小另调用 `aggregateNet(...,return.object=FALSE,remove.isolate=FALSE)$group$count`，取行和+列和−对角线，自环在大小中只计一次。
4. **视觉编码**：x/y 是所选中心性之和，色是群，shape 是用户分类，大小是聚合 LR 连接计数；大小不必与中心性使用相同生产阈值，且不是细胞数量。
5. **信息损失**：去掉各信号/各 partner 的分布；两个路由不同的群可落在同一点。非线性 measure 的跨层和也不是合并网络的中心性。
6. **统计边界**：无比较检验或误差条。点大小来自另一次 aggregateNet，不能由大小推导 x/y 的显著性；必须核对门控来源一致性。
7. **前置计算及返回**：需 group `centr`，还需 aggregateNet 所消费的 LR 网络；返回 ggplot。源码把中心性向量按位置赋给 `levels(idents)`，未显示按 node 名重排，名称次序属于迁移门禁（A-07）。
8. **判别反例（未执行）**：把单群自环计数设为 2，应按 `2+2−2=2` 编码大小，而默认 out/in 坐标各含该自环。保持中心性不变，仅改变 aggregateNet 门控下的弱 LR，大小可改变但坐标不变，不能解释成强度同步变化。

### 2.4 `netAnalysis_signalingRole_heatmap`

依据：[文档意图](../R/analysis.R#L3164-L3173)、[实现](../R/analysis.R#L3192-L3278)。

1. **科学问题**：哪些群相对主导某信号的 outgoing/incoming，总体强度分布又如何？
2. **输入单位与轴**：从 group `centr` 取 outdeg/indeg，形成 signaling × group 矩阵。`pattern="all"` 是二者相加；用户指定但不存在的 signaling 先补零行。
3. **计算定义**：outgoing/incoming 用相应 measure；all 用 $O+I$，自环在该和中计两次。保留 `mat.ori`，色值按每行最大值缩放，零改 NA。顶部条是 raw 列和；右侧条却是 raw 行和经过 $-1/\log$ 及异常值重排后的显示量，而非文案所称的原始总强度。
4. **视觉编码**：色块是行内 relative strength，顶条是每群 raw 总量，右条是 transformed 每信号总量；三者尺度不同。白色合并了零及缺失/补零后的不可定义缩放。
5. **信息损失**：色块不能比较信号绝对强弱；out/in 合并丢失方向；补零使缺失层与真零难区分。
6. **统计边界**：无显著性检验。色块亮不等于更显著，右条显示变换也不产生统计标准化（A-05）。
7. **前置计算及返回**：需 group `centr` 的 outdeg/indeg；返回 `ComplexHeatmap::Heatmap` 对象，函数不主动 draw。
8. **判别反例（未执行）**：两信号 `(0.1,0.2)` 与 `(1,2)` 色块相同，但 raw 顶条不同，右条不可按 `1:10` 解读；再加入不存在的层与确实零层，核对它们是否应使用同一种白色。

### 2.5 `netAnalysis_diff_signalingRole_scatter`

依据：[R/analysis.R:2803–2964](../R/analysis.R#L2803-L2964)。

1. **科学问题**：两个条件间，哪些**细胞群**的总体发送/接收角色增强或减弱？
2. **输入单位与轴**：一个点是群；文档接受 merged 对象或对象列表（先 merge），旧源码按 dataset 嵌套读取两份 group `centr`。默认 signaling 是两条件 pathway 并集。
3. **计算定义**：每条件对 signaling 名称对齐，缺失层补零；先作条件 2−条件 1，再沿 layer 求和得 $\Delta O_g,\Delta I_g$。`signaling.exclude` 先排层，`idents.exclude` 在点表生成后排群。node 轴赋值仍按位置，不能把统一群标签误当已完成名称对齐。
4. **视觉编码**：x/y 是中心性差值，零线标示无净变化，正号指第二条件更大；颜色是群、shape 是用户分类，点大小固定。
5. **信息损失**：各 pathway 的相反变化会在总和里抵消；缺失层补零把未测/未产出与无通讯合并。
6. **统计边界**：没有差异检验、重复层级或置信区间；位于原点只说明选定指标净变化为零，不是所有信号不变。非 degree measure 的差值仍保持该 measure 单位，不应沿用 strength 标签。
7. **前置计算及返回**：两条件分别计算相同定义的 group 中心性；当前旧 merged 路径需迁移至真实 dataset 来源。返回 ggplot，不返回检验表。
8. **判别反例（未执行）**：某群 pathway A outgoing 增 `2`，B 减 `2`，点的 x 应为零，但逐层变化非零。把第二条件中心性 node 顺序置换而保留名称，数学差值应不变；源码位置赋值风险须被识别，不能用只看维度的检查替代。

### 2.6 `netAnalysis_signalingChanges_scatter`

依据：[R/analysis.R:2968–3161](../R/analysis.R#L2968-L3161)。

1. **科学问题**：固定一个群，两个条件间哪些**信号层**改变了其 outgoing/incoming？与上一卡的“一个点一个群”不同。
2. **输入单位与轴**：从两条件 `centr` 取 `idents.use`；一个点是 pathway/所选信号，坐标是条件 2−条件 1 的两个 measure。文案/分支出现 single 对象提示，但后续仍使用两条件索引，不能据此承诺 single 模式可运行（A-01）。
3. **计算定义**：层并集补零后取差。specificity 用两条件某方向乘积为零且和非零来识别“仅一条件非零”；形状再分 incoming/outgoing/both，颜色结合差值符号分 Shared、条件1/2 specific。源码用 `rowSums(mat.diff)!=0` 留点，因此 $(+a,-a)$ 会被删除。
4. **视觉编码**：x/y 是方向差值，色是零/非零与差值符号派生的条件类别，shape 是方向特异类别，点大小固定。`top.label` 依据绝对坐标分位数选择标签，不是按 p 值排序。
5. **信息损失**：不能看到该群与哪个 partner 发生改变；缺层补零可能制造“specific”；相反等幅的真正角色交换被现有过滤丢失（A-13）。
6. **统计边界**：“specific”只表示阈值后支持/零值状态，不是统计学条件特异性。差异方向与显著性必须分开。
7. **前置计算及返回**：需两条件同指标 group `centr`、可比群身份及信号映射；返回 ggplot。缺失层与真实零、single 分支和比较条件数均待核定。
8. **判别反例（未执行）**：固定群的某信号从 `(out=1,in=3)` 变为 `(3,1)`，差为 `(2,−2)`；科学上是角色交换，现有过滤却删除它。此场景应作为更正裁定，而不能要求迁移永远删除该点。

## 3. 群通信模式与单细胞 topic

### 3.1 `selectK`

依据：[R/analysis.R:950–1009](../R/analysis.R#L950-L1009)。

1. **科学问题**：文档意图是以 NMF consensus 的 Cophenetic/Silhouette 稳定性辅助选择**群模式**数；不是自动证明真实生物程序数。
2. **输入单位与轴**：旧 `slot$prob`，outgoing 对 receiver 求和得 group × layer，incoming 对 sender 求和；每列除以该层最大值，再剔除行和为零的群。
3. **计算定义**：意图对 `k.range` 做多次 NMF rank estimation；实现准备了 `data`，调用却是 `nmfEstimateRank(data_sr, range=20:30,nrun=30L,seed=seed.use)`。`data_sr` 未在该函数定义，且 `k.range/nrun` 没有传入该调用。不能把文档意图当作当前已实现选择算法（A-08）。
4. **视觉编码**：若估计结果可得，横轴 rank，纵轴 Cophenetic 或 Silhouette，颜色/分面区分指标；曲线不是检验 p 值。
5. **信息损失**：输入只保留群的边缘强度相对构成，路由和绝对层强度已丢失；稳定性评价不提供完整模型适配信息。
6. **统计边界**：同一矩阵多随机起点的共识稳定性不是跨生物学样本复现性；当前错误调用不能作为有效选择证据。
7. **前置计算及返回**：需 group 层权重；文档/末尾代码意图返回 ggplot，不写入自动选出的 k。`identifyCommunicationPatterns(k=NULL)` 只报错要求人工选择，不自动消费此图。
8. **判别反例（未执行）**：用户只请求 ranks `2,3`、`nrun=4`，应能区分“尊重参数的 rank 评估”与当前硬编码 `20:30/30`；两个路由不同但 outgoing 边缘矩阵相同的网络应给相同的输入评估问题，而非声称检测到路由差异。

### 3.2 `identifyCommunicationPatterns`

依据：[R/analysis.R:1223–1327](../R/analysis.R#L1223-L1327)、[scaleMat](../R/utilities.R#L998-L1024)。

1. **科学问题**：寻找群与信号层共同参与的非负相对模式；不是直接聚类细胞，也不是重建全部 sender–receiver 路由。
2. **输入单位与轴**：$X_{g\ell}=O_{g\ell}/\max_gO_{g\ell}$ 或对应 incoming；群 × layer。剔除零行后以指定 k 分解。没有新的 p 值门控，使用已有网络状态；全零列除零边界需要核定。
3. **计算定义**：`NMF::nmf(...,method="lee",seed="nndsvd")` 得 $X\approx W H$，$W$ 为 group × pattern，$H$ 为 pattern × layer。存储前 `scaleMat(W,"r1")` 令每群的 pattern 份额和为 1；`scaleMat(H,"c1")` 令每层的 pattern 份额和为 1。分别归一化后的因子不再承诺原始 $WH$ 重构。
4. **视觉编码**：可选两热图分别是群–pattern 与信号–pattern 的归一化 loading；其消费者 `netAnalysis_river/dot` 见下两卡及 [可视化语义卡](VISUALIZATION_SEMANTICS.md)。不是各 pathway raw 强度热图。
5. **信息损失**：先边缘化丢掉对端，逐层 max 缩放丢掉层量级，低秩分解及 loading 归一化进一步压缩；pattern 编号也没有跨独立拟合的天然对应关系。
6. **统计边界**：loading 不是通讯 p 值，k 个 pattern 不是独立统计显著程序。零行/零列和因子非唯一性影响解释（A-09）。
7. **前置计算及返回**：需用户指定 k 和 group 层权重；返回更新对象，`slot$pattern[[outgoing/incoming]]` 包含 `data=data0` 与长表 `pattern$cell/signaling`。`data0` 是 max-scaled 输入，不是原始通讯张量；`k.range` 不在拟合中自动选 k。
8. **判别反例（未执行）**：把一个 pathway 所有边同时乘 100，其列 max 归一化后的 $X$ 不变，raw information flow 却增百倍。三节点顺/逆环边缘和相同，群模式不能识别路由反转；应另用边支持相似性或网络图。

### 3.3 `netAnalysis_river`

依据：[R/visualization.R:3978–4104](../R/visualization.R#L3978-L4104)、[后半与返回](../R/visualization.R#L4103-L4203)。

1. **科学问题**：哪些群和信号在某 NMF pattern 上有较高归一化 loading？
2. **输入单位与轴**：读群模式的 `pattern$cell`（group–pattern）及 `pattern$signaling`（pattern–layer）长表；不是 cell topic 的两矩阵。outgoing/incoming 选择已经拟合好的不同模型。
3. **计算定义**：分别把两表 `Contribution<cutoff` 置零，默认 `0.5`；两表独立画两段 river 并拼接，没有计算完整的 group–pattern–layer 联合流。无 signaling 表时只画群–pattern。sources/targets 在有 signaling 的分支都筛同一 CellGroup 轴，不新增 sender–receiver 轴。
4. **视觉编码**：带宽是各段的归一化 loading；颜色标记群、pattern 或 signaling，顺序可由 loading 距离聚类决定。两段同名 pattern 的宽度不构成守恒的物质/通讯流。
5. **信息损失**：低于 cutoff 的混合参与被隐藏；独立归一化的两侧不能相乘解读成真实边权，筛选也不重新拟合模式。
6. **统计边界**：cutoff 是展示规则，不是显著性或 posterior probability 阈值。
7. **前置计算及返回**：需 `identifyCommunicationPatterns` 的相应 pattern 表；返回组合图对象。输出不更新 NMF 模型；细粒度布局语义见可视化文档。
8. **判别反例（未执行）**：某群在两模式均为 `0.5`，默认保留两条；改为三模式均 `1/3`，默认全部隐藏，尽管群有通讯。某 pattern 的群侧总 loading 与信号侧总 loading 不相等，不能要求河流跨两图守恒。

### 3.4 `netAnalysis_dot`

依据：[R/visualization.R:4229–4306](../R/visualization.R#L4229-L4306)。

1. **科学问题**：以最强的共同 pattern loading 概括某群与某信号的模式关联；不是完整 NMF 重建矩阵。
2. **输入单位与轴**：读归一化 group × pattern 的 $\widetilde W$、pattern × layer 的 $\widetilde H$ 长表及群模式输入 `data0`；图的行是群，列是 signaling。
3. **计算定义**：默认 cutoff 为 `1/k`，两因子分别低于 cutoff 置零；按 pattern join 后乘积。原 `data0` 对应群–信号为零则把乘积置零，再按群–信号取 `top_n(1,Contribution)`，即 $D_{g\ell}=\max_a(\widetilde W_{ga}\widetilde H_{a\ell})$ 的显示分数，**不是** $\sum_a\widetilde W_{ga}\widetilde H_{a\ell}$；并列最大可保留多行。
4. **视觉编码**：点大小是最大乘积分数，颜色是群；零分数转 NA 因而不显示。pathway/group 筛选改变显示集合，不重新训练。
5. **信息损失**：非最大 pattern 的贡献被丢弃；点缺失可能来自原输入零，也可能来自 cutoff，不等于未推断到通讯。源码用无分隔符拼接群/信号作键，[INFERENCE] 名称碰撞还可能误合并不同组合（A-10）。
6. **统计边界**：点大小不是原始权重、重建期望或 p 值。max 与 sum 是不同科学摘要，不能在存储迁移时悄然互换。
7. **前置计算及返回**：需相应 group pattern 与 `data0`；返回 ggplot，不写回重建矩阵。若未来选择 sum，应作为独立语义变更，而不是声称修复存储读取。
8. **判别反例（未执行）**：同一群–信号在两个 pattern 上乘积均为 `0.25`，max 应为 `0.25`，sum 为 `0.5`，且当前 `top_n` 会保留并列行。另将原输入设零但因子乘积非零，当前 masking 应仍隐藏该点。

### 3.5 `identifyCellTopics` 与 rank 辅助 `suggestK1`

依据：[R/analysis.R:3447–3557](../R/analysis.R#L3447-L3557)、[suggestK1](../R/analysis.R#L3812-L3923)、[CellChat2Seurat](../R/utilities.R#L712-L745)。

1. **科学问题**：按每个细胞的发送/接收信号强度分解出 topic，描述细胞间通讯程序异质性；不是群模式的重命名，也不是基因表达 topic。
2. **输入单位与轴**：从 `centr.cell` 取 outdeg 或 indeg，转为 layer × cell 矩阵 $A$；默认 `slot.name="net"` 因而层是 LR，选择 netP 才是 pathway。默认不缩放；`do.scale=TRUE` 转 Seurat 后 `ScaleData(do.center=FALSE)`，没有在此调用 NormalizeData。
3. **计算定义**：一个 `topic.k` 时 `RcppML::nmf` 分解 $A\approx w\,diag(d)\,h$，存 $W=w\,diag(\sqrt d)$（layer × topic）、$H=(diag(\sqrt d)h)^T$（cell × topic），因此 $A\approx WH^T$；不使用群模式的 r1/c1 归一化。多个 k 调 `suggestK1`：`p>0` 时各抽约 p 比例的行与列，仍在**完整 A** 上拟合，再在抽出的子矩阵计算 Frobenius 残差/元素数；不是留出子矩阵训练。
4. **视觉编码**：单 k 本函数不画 topic 分布，后续 topic 空间图/热图消费 H/W；多 k 显示 rank 评价图。topic 空间图点代表 cell，但其颜色是 topic loading，而不是网络 UMAP 的点语义。
5. **信息损失**：out/in 边缘强度先丢失对端路由；低秩 topic 再丢失残差。H/W 不是概率，topic 编号在独立拟合间无天然对应关系。可选缩放会改变层之间的相对权重。
6. **统计边界**：文档称 rank 分支为 cross validation，但完整训练后抽样算残差不构成 out-of-sample 验证。抽行和抽列意味着约 $p^2$ 的元素覆盖，而不是 p 比例独立留出元素；helper 接受 L1 却未传给内部 nmf，单 k 拟合才传 L1（A-15）。不能用该误差宣称泛化能力。
7. **前置计算及返回**：需 `do.group=FALSE` 中心性。单 k 写 `slot$topic[[pattern]]$cell/signaling`；多 k 只写 `slot$tmp$suggestK[[pattern]]` 并返回对象，不自动选择/生成 topic。`do.scale` 还依赖旧 `CellChat2Seurat` 对 data/坐标的读取，最终 schema 兼容尚需核定。
8. **判别反例（未执行）**：两群均值相同但各群内细胞分成强/弱两类，群模式可相同，cell topic 输入仍保留细胞差异。把一个信号层乘百倍，默认 topic 输入改变，但群 max-scaled 模式不变；rank 评估中给抽中子矩阵加噪声，需区分“噪声已参与训练”与真正留出验证。

### 3.6 `identifyTopicCommunities` 与 Seurat 桥接

依据：[社区生产者](../R/analysis.R#L3595-L3790)、[CellChat2Seurat](../R/utilities.R#L712-L761)、[PP1 topic/community 消费链](PP1_DOWNSTREAM_REFACTORING_PLAN.md#33-niche--当前-topiccommunityseurat-api)。当前存在公开社区算法；旧 niche API 的历史等价性尚未证明，不能用网络 `netClustering` 替代。

1. **科学问题**：按细胞/spot 的发送或接收 topic loading 划分通讯社区；可选空间分支进一步纳入位置。它不是对 LR/pathway 网络聚类，也不等同于既有细胞类型。
2. **输入单位与轴**：同一 slot/pattern 的 H 为 cell×topic，W 为 LR/pathway×topic；`centr.cell` 的 indeg/outdeg 转成通讯特征×cell assay。H、W、中心性、metadata、坐标必须按身份一致，而不是只检查维度相等。
3. **计算定义**：H 的前 nPrograms 个 topic 经 FindNeighbors 构造 SNN，再 FindClusters；默认 seurat 分支不使用 spatialWeight。显式 resolution 直接聚类，否则扫描分辨率并估计社区数；autoStardust、weightStardust、fusion 的空间组合各异。`topic.k` 与 resolution 控制不同阶段；seed.use 当前未下传。
4. **视觉编码**：默认 do.plot 输出社区空间图与 marker 热图；空间点颜色是 ICR/OCR 类别，不是 H 的连续 loading。后续 Seurat UMAP 的点仍是 cell/spot，不是信号网络。
5. **信息损失**：低秩 H 再经邻居图和离散社区标签，不能由社区恢复完整通讯路由。nPrograms 会丢弃未选 topic；社区编号不具有跨拟合天然身份。空间分支与纯 H 分支也不能互称等价。
6. **统计边界**：当前将原始通讯强度放入 Seurat data 层再计算 markers；该层的 log-data 假设影响 fold change 和 AverageExpression，原量均值相同并不保证这些输出相同。需核定 assay 层、变换及效应量（A-17），不能直接套 RNA 归一化。空间分支距离矩阵的邻居解释和自动 resolution 下界另需核定；社区整链未运行。
7. **前置计算及返回**：需 cell 中心性与固定 k 的 H/W。返回更新 Chat；写 meta$ICR/OCR，保持 Chat@idents 不变；return.seurat=TRUE 时保存 topic[[pattern]]$seurat.obj，含 label.ICR/OCR、reduction 和 markers。桥接器即使显式传 counts 仍读取旧 @data。topic 重拟合、net/netP 共用 meta 键、return.seurat=FALSE 保留旧 srt 的状态风险须按 PP1-10 验收，不能只改槽名。
8. **判别反例（未执行）**：固定 H 改坐标，seurat 分支不应被误认为应用了 spatialWeight；原量均值相同而分布不同的两组区分线性均值与 expm1 后均值；依次重拟合 topic、切换 net/netP、关闭 Seurat 保存并重载，检验 H/W、labels、markers 是否同源。完整场景见 PP1 T15–T18；既有量纲探针不替代社区验收。

## 4. 网络相似性、嵌入与分类

### 4.1 `computeNetSimilarity` / `computeNetSimilarityPairwise`

依据：[single](../R/analysis.R#L1331-L1403)、[pairwise](../R/analysis.R#L1408-L1508)、[buildSNN](../R/analysis.R#L1671-L1686)、[结构距离及辅助](../R/analysis.R#L1996-L2165)。

1. **科学问题**：哪些信号网络具有相似的通讯支持或拓扑形态？`functional` 在这里衡量**同名节点对上正边支持的重叠**，不直接衡量基因功能、通路作用或权重强弱。
2. **输入单位与轴**：single 从旧 `slot$prob` 取每层网络，输出 layer × layer；pairwise 拼接所选 dataset 的所有层，节点标签为 `layer--dataset`，不是先把同名 pathway 平均。functional 要求节点身份、顺序及方向可比较；structural 使用各图的拓扑分布，可处理不同节点数但不代表群身份匹配。
3. **计算定义**：可选 `thresh` 是非零权重分位数裁剪比例，**不是 p 值阈值**；single 在全张量计算分位数，pairwise 每 dataset 各算。二者均先二值化 $B_{ij\ell}=1[P_{ij\ell}>0]$。functional 原始相似性为 $J(a,b)=|E_a\cap E_b|/|E_a\cup E_b|$；structural 为 `1-computeNetD_structure`，后者组合节点距离分布、网络节点分散度、排序 alpha-centrality 及补图的分布差异，默认权重 `0.45/0.45/0.1`。最终保存的是原始相似性与 SNN **逐元素相乘**的结果，不是仅 Jaccard 或 `1-D`。
4. **视觉编码**：自身不画图；供 netEmbedding/netClustering 与 embedding 图消费。点的空间接近来自后续相似性剖面降维，不是通讯距离或组织距离。
5. **信息损失**：权重二值化丢掉幅度；functional 保留具名边支持但不保留其强度，structural 更不承诺区分具名路由；SNN 再抑制非近邻关系。pairwise 没有先按 node 名重排矩阵（A-07）。
6. **统计边界**：不做相似性显著性检验。functional 的 `0/0` 被改零、对角设一；structural 的 Inf/NA 距离被改零，等效相似性一，可能把不可计算当完全相似（A-11）。这些退化规则不是经核定的生物学含义。
7. **前置计算及返回**：需对应粒度的网络、可比名称与权重口径；返回更新对象，写 `slot$similarity[[type]]$matrix[["single" 或 comparison.key]]`。这里 type 两分支、single/pairwise 四条路径均需迁移；joint 顶层网络不能假作旧 dataset 列表。
8. **判别反例（未执行）**：同一支持的两图，全部正边从 `0.01` 放大到 `10`，无裁剪时 functional 原始 Jaccard 仍为 1；改成三节点反向环则边缘和相同而支持不重叠，Jaccard 为 0。仅置换矩阵行列顺序且保留名称，数学相似性应不变；空图与计算失败也不能默认同义。

### 4.2 `netEmbedding` / `runUMAP`

依据：[netEmbedding](../R/analysis.R#L1512-L1571)、[runUMAP](../R/utilities.R#L2022-L2069)、[single 图消费](../R/visualization.R#L4337-L4370)、[pairwise 图消费](../R/visualization.R#L4518-L4571)。

1. **科学问题**：把信号网络之间的相似性剖面压缩到低维，发现网络层之间的局部组织结构；不是给细胞做转录组 UMAP。
2. **输入单位与轴**：输入先前保存的相似性方阵；每个点是一条信号网络，pairwise 是一个 layer–dataset 实例。根据 mode/comparison 找 matrix；默认删除列和恰为 1 的点，即在既有对角一、非负相似性的定义下，删除只有自身支持的孤立信号。
3. **计算定义**：对过滤后的 Similarity 运行 UMAP，而不是显式传预计算距离 `1-S`。`umap-learn` 经 `runUMAP` 转置输入，默认 `metric="correlation"`、seed 42；`uwot` 直接使用 Similarity，源码没有显式指定 metric/seed。对称矩阵转置不改数值，但两个后端的参数语义不能据此视作一致（A-12）。
4. **视觉编码**：自身不画；`netVisual_embedding`/Pairwise 以 Y 的二维坐标定位，颜色来自网络 cluster，点大小另用该层总权重/最大总权重；pairwise 另标 dataset。因此位置、颜色、大小分别来自嵌入、分类和强度三个过程。
5. **信息损失**：降维不保留全部距离，默认移除孤立信号；UMAP 坐标旋转、镜像或局部挤压不具有直接生物单位。点大不必更远，也不必更相似。
6. **统计边界**：嵌入距离不是效应量或 p 值；独立拟合坐标不能直接相减比较。后端/seed/邻居数会影响显示；这是模型选择而非检验。
7. **前置计算及返回**：先算指定 type/comparison 的 similarity；返回更新对象，写 `$similarity[[type]]$dr[[comparison.key]]`。后续图必须把 Y、group 和层总量按网络实例名称对齐；只共享行数不够。
8. **判别反例（未执行）**：只放大同一支持网络的权重而不改二值支持，similarity 输入不变，但 embedding 图的点大小会变。删除一个孤立信号后另将 prob 的层顺序打乱，位置标签和点大小必须仍按名称归位，不能由位置拼接产生错误解释。

### 4.3 `netClustering`

依据：[R/analysis.R:1575–1653](../R/analysis.R#L1575-L1653)、[eigengap/Laplacian 辅助入口](../R/analysis.R#L1699-L1762)。

1. **科学问题**：给嵌入空间中的信号网络分组；不是细胞类型聚类，也不是检验同组通路功能等价。
2. **输入单位与轴**：读取 `similarity[[type]]$dr[[comparison.key]]` 的 Y，行是信号网络，列是嵌入维度；不是读取原始相似性方阵做所有分支的聚类。
3. **计算定义**：kmeans 指定 k 时在 Y 上 `nstart=10`；未指定 k 时遍历 `2:min(N−1,10)`，平均不同 k 的同簇矩阵形成 consensus，用 computeEigengap 的 upper_bound 选最终 k。spectral 分支同样把 **Y** 当 A 构造 `diag(rowSums(A))−A` 并做特征分解；普通 N×2 嵌入不是方阵，[INFERENCE] 该分支的邻接/嵌入语义不成立，不能当作已可用谱聚类（A-12）。
4. **视觉编码**：写出的 group 供 embedding 图着色；未指定 k 的 kmeans 可保存 eigengap PDF，spectral 分支直接画 eigenvalue 序列。颜色是算法分组标签，不是条件或强度。
5. **信息损失**：用离散 cluster 压缩连续嵌入；分类还继承前面二值化、SNN 与 UMAP 的损失。cluster 编号本身无大小顺序或跨独立拟合对应关系。
6. **统计边界**：自动 k 的共识矩阵来自不同 k 的划分，不是生物学重复的统计置信度；谱分支的输入定义需要先裁定，不能仅为让矩阵运算通过而改形状。
7. **前置计算及返回**：需同 type/comparison 的 embedding；返回更新对象，写 `$similarity[[type]]$group[[comparison.key]]`。与 Y 的网络实例身份必须保持对应。
8. **判别反例（未执行）**：相同信号网络换一个 UMAP 邻居参数后 kmeans 分组可能变化，不能解读成生物变化。给 spectral 常见的五个网络二维 Y，需裁定它应消费五阶相似性图而非把 5×2 坐标伪装成邻接矩阵。

### 4.4 `rankSimilarity`

依据：[R/analysis.R:1785–1847](../R/analysis.R#L1785-L1847)。

1. **科学问题**：在同一次联合嵌入中，哪些同名信号在两条件的坐标位移较大？名称虽为 similarity，排序量实为低维距离。
2. **输入单位与轴**：`comparison1` 选已联合嵌入的 dataset 集，`comparison2` 是其中两条件的索引；从 `layer--dataset` 的 Y 行名拆出身份，仅取两条件共同出现的 layer。
3. **计算定义**：按共同 pathway 名显式重排两份坐标，计算 $d_\ell=\|Y_{\ell,d_2}-Y_{\ell,d_1}\|_2$，按距离排序；不是原始 Jaccard、结构距离、raw information flow 差或 fold-change。
4. **视觉编码**：横向条长为嵌入空间欧氏距离，标签是共同 pathway；没有方向正负，不能从条长读出哪个条件增强。
5. **信息损失**：只保留共享且未被嵌入移除的信号，条件特有/缺失信号不参与排名；把向量位移压成长度，丢掉方向和原始网络改变细节。
6. **统计边界**：没有检验。距离依赖联合嵌入的网络集合、后端和参数，不能当作跨实验可比的绝对效应量。
7. **前置计算及返回**：需相应联合 similarity→embedding；返回 ggplot，数据表在图数据中，而非单独返回统计对象。dataset 名仍从旧 slot 列表取，行名解析与新存档键需要统一（A-01/A-07）。
8. **判别反例（未执行）**：某 pathway 仅条件2存在，可能是重要差异，却不会出现在 shared-only 排名；共同 pathway 在两个独立 UMAP 中仅整体旋转产生的位移也不能交给此函数解释，必须使用同一次联合嵌入。

## 5. 表达差异注释

### 5.1 `netMappingDEG`

依据：[R/analysis.R:3283–3387](../R/analysis.R#L3283-L3387)、[DEG 对照定义](../R/utilities.R#L1503-L1533)、[统计列及筛选](../R/utilities.R#L1534-L1587)、[当前存储位置](../R/utilities.R#L1624-L1636)。

1. **科学问题**：把 ligand 所在 source 群和 receptor 所在 target 群的表达差异证据注释到已推断通讯行；不是推断新的通讯差异或因果效应。
2. **输入单位与轴**：通讯行是 source–target–LR，DEG 行是 cluster–gene；来源是 `features.name + ".info"`。上游 Wilcoxon 可比较群内指定 dataset 对其余 dataset，或某群对其它细胞；必须保留实际对照定义。Moran's I 等空间特征统计不能假装含同样的 logFC/pct DEG 列。
3. **计算定义**：先 `subsetCommunication(thresh)`，按 `source.ligand` / `target.receptor` 与 `clusters.features` 匹配。单基因直接抄 pvalues/logFC/pct.1/pct.2；复合物展开亚基，对**匹配到的亚基**各统计列取算术均值，缺失留 NA。映射不重新计算表达或通讯检验。
4. **视觉编码**：自身不画图，返回表供 subsetCommunication 的 DEG 过滤及差异展示消费；通讯 `pval` 与 ligand/receptor 的 `pvalues` 必须区分。
5. **信息损失**：复合物均值掩盖亚基方向不一致、缺失和极端限制；一条 DEG 统计可能复制到很多通讯边，不增加独立证据数。仅 cluster–gene 键不包含 dataset，不能把合并后的行标签当作已解决对照身份。
6. **统计边界**：上游同时产生 raw `pvalues` 和 Bonferroni `pvalues.adj`，本函数映射的是 raw 列；复合物 p 值均值不是复合物联合显著性检验，平均 logFC 也不证明复合物功能同步变化（A-14）。DEG 的不显著/未保留基因产生缺失，不等于零变化。
7. **前置计算及返回**：需合适的 DEG 统计表、LR 注释、group 通讯。文档意图返回附八个 ligand/receptor 统计列的 data.frame；当前读取旧 `object@var.features`，而上游已写 `misc$.var.features`，还需核定 single data.frame 与 merged list 的分支区分，不能仅迁移字段就宣布端到端可用（A-01/A-14）。
8. **判别反例（未执行）**：复合物两个亚基 logFC 为 `+2/−2`，平均为零，却不是“两个都不变”；p 值为 `0.001/0.9` 的均值也不是该复合物显著性。source 与 target 互换但基因不变，应映射到不同群的 DEG，防止只按基因名连接。

## 6. 待裁定问题与迁移前置门禁

以下 ID 用于后续 Note/迁移批次引用。建议方向均**未批准实施**；只阻断受影响路径，不要求先改完所有问题才能迁移独立消费者。证据指向当前源码，不把旧缺陷当兼容性要求。

| ID | 证据与可观察问题 | 受影响函数/消费者 | 建议方向（未决建议） | 必须裁定的边界 |
|---|---|---|---|---|
| **A-01** | [旧字段/merged 读取](../R/analysis.R#L856-L898)、[subset merged](../R/analysis.R#L2319-L2365)、[旧 DEG 位置](../R/analysis.R#L3296-L3314) 对照 [新 DEG 生产](../R/utilities.R#L1624-L1636)；部分分析已新读 group/cell，其他仍读 prob/prob.cell/options/var.features | contribution、rankNet、extractEnrichedLR、subsetCommunication merged、模式/相似性、两差异图、netMappingDEG、CellChat2Seurat | 按实际生产者/消费者成对切换，使用真实样本存档而非把 joint 网络按 dataset 下标解释；保留 `centr/centr.cell` 当前落点 | 每入口支持的 single/merged/cell/group 组合；哪些路径单纯存储迁移，哪些需改变 API 承诺。不得借 `params()` 读取 mode |
| **A-02** | [提取 helper](../R/analysis.R#L913-L946) 计算 gated prob 却只按缓存 LR.sig 选择；空交集可能仍提基因；FALSE 分支名称类型与 TRUE 不同 | extractEnrichedLR*、plotGeneExpression 等基因消费者 | 明确“缓存 enrichment”或“本次重筛”的唯一含义，统一 LR 与基因来源并追溯缓存阈值 | 是否改变 thresh 的作用、空集返回与 enriched.only=FALSE 语义；与 aggregate/可视化阈值等号差异分别裁定 |
| **A-03** | [min–max](../R/analysis.R#L296-L308)、[分区与返回](../R/analysis.R#L354-L393) 改变 raw 份额，返回最后分区数据 | netAnalysis_contribution | 区分 raw 相对贡献、显示归一化与 receiver 分区，表与图应指向同一命名内容 | 是否移除 min–max；全常数输出；返回 All/分区表的 API，不能在存储迁移暗改 |
| **A-04** | [统计前缓存](../R/analysis.R#L565-L595)、[检验及失败 p=0](../R/analysis.R#L715-L735)；展平边配对未按名字、失败被显著化 | rankNet(do.stat)、依赖其 pvalues/标签颜色者 | 无法检验必须与显著区分；让检验对象与显示 source/target 子集相同，先名称配对 | 边是否可作有效重复、缺层如何处理、多重比较、失败/NA 表达方式；不能以旧 p=0 作为通过条件 |
| **A-05** | [rankNet scaled](../R/analysis.R#L497-L548)、[role heatmap 三尺度](../R/analysis.R#L3232-L3268)；raw、scaled、relative 与右条文案不一致 | rankNet、signalingRole_heatmap、下游导出表 | 显示变换与原始统计量显式分列/标注，数值推断用 raw | 是否保留 $-1/\log$、大于等于 1/常数/全零的显示策略、single show.raw 的契约 |
| **A-06** | [local 异常归零](../R/analysis.R#L29-L54)、[producer 默认](../R/analysis.R#L83-L91)、[network 默认指标](../R/analysis.R#L2633-L2648) | centrality、角色 network/heatmap/scatter、cell topic | 明确 measure 集与失败状态；默认调用链能满足请求指标，不把失败零当科学零 | 是否调整默认值或要求调用者显式选择；失败值 NA/状态与消费者呈现；另迁存储落点须独立决策 |
| **A-07** | [role 节点按位置赋值](../R/analysis.R#L2725-L2733)、[比较层补零](../R/analysis.R#L2882-L2913)、[functional 直接矩阵相乘](../R/analysis.R#L1463-L1476)、[embedding 显示位置拼接](../R/visualization.R#L4357-L4358) | 角色图、rankNet 配对、computeNetSimilarityPairwise、embedding 图、rankSimilarity | node/layer/dataset 三轴按名称对齐；真实零与缺失层分开；用 permutation 反例审阅 | 缺层补零是否科学合理、不同细胞组成的可比性、联合实例命名与名称分隔符 |
| **A-08** | [selectK 实参](../R/analysis.R#L971-L1009)：局部准备 data，实际用未定义 data_sr 和硬编码 rank/nrun | selectK → identifyCommunicationPatterns 的 k 选择 | 按核定的输入矩阵/用户参数评估 rank，并明确人工选择而非自动最优 | 稳定性方法是否保留、空/低秩输入处理；不能把未运行的图意图视为有效证据 |
| **A-09** | [群 max 缩放与因子归一化](../R/analysis.R#L1253-L1273)、[scaleMat r1/c1](../R/utilities.R#L1021-L1022) | selectK、identifyCommunicationPatterns、river/dot | 将输入归一化、拟合因子、展示份额分清；明确全零列与零行规则 | 是否保留原因子供重建、零/NA 处理与跨拟合 pattern 对齐；不要误写成 cell topic 契约 |
| **A-10** | [dot 乘积、mask、top_n](../R/visualization.R#L4248-L4282)：max 而非 sum，拼接键可能碰撞、并列可多行 | netAnalysis_dot 及对“contribution”的说明 | 先裁定“最强共同 pattern”还是“加总重建”；用成对身份区分群–层，不依赖歧义字符串 | max/sum 是否改义、并列如何展示、mask 与 cutoff 的科学解释；旧形状/图像不等于正确基线 |
| **A-11** | [functional/structural 退化值](../R/analysis.R#L1356-L1395)、[pairwise](../R/analysis.R#L1453-L1501)；结构失败距离归零成相似性一 | computeNetSimilarity*、embedding/clustering | 明确支持相似性与权重无关；把不可计算与完全相似区分，记录裁剪与 SNN 两阶段 | 空图比较、结构失败策略、支持/权重是否另立量；本次不添加新相似性模型 |
| **A-12** | [embedding 后端](../R/analysis.R#L1546-L1564)、[runUMAP](../R/utilities.R#L2022-L2068)、[spectral 读取 Y](../R/analysis.R#L1607-L1647) | netEmbedding、netClustering、rankSimilarity、embedding 图 | 为“相似性剖面嵌入”和“基于图的谱聚类”分别核定输入；后端参数可追溯 | correlation/euclidean/预计算距离是否变更、随机性、谱聚类究竟消费 S 还是 Y，不能只做维度补丁 |
| **A-13** | [rowSums 差过滤](../R/analysis.R#L3081-L3089) 删除方向相反等幅改变 | netAnalysis_signalingChanges_scatter | 根据每方向是否有变化判断保留，而不是总和是否零；另核定 specificity 命名 | 是否展示零变化层、缺层/零值与“特异”定义；$(a,-a)$ 不应固化为必须隐藏 |
| **A-14** | [DEG 单基因/复合物映射](../R/analysis.R#L3315-L3384)、[上游 raw/adjusted](../R/utilities.R#L1555-L1587)；p 值平均、亚基缺失、对照身份与表/list 区分 | netMappingDEG、subsetCommunication 的 DEG 过滤 | 把注释摘要与统计检验分开；保留亚基证据、对照键和缺失，不把均值称联合显著性 | raw/adjusted 列选择、复合物汇总规则、dataset 对照语义与 single/merged 返回处理 |
| **A-15** | [topic rank 分支](../R/analysis.R#L3524-L3554)、[suggestK1](../R/analysis.R#L3838-L3910)：全矩阵训练后抽子矩阵算误差，L1 未下传 | identifyCellTopics 多 k、suggestK1 | 如保留算法，应称训练重建误差抽样评价；若要交叉验证，须单独核定真正留出训练协议 | prop 的元素/行列单位、是否改算法或文案、L1 与单 k 拟合一致性、误差归一化 |
| **A-16** | [网络门控/表 bypass](../R/analysis.R#L2412-L2449)、[pathway p 均值](../R/analysis.R#L2549-L2560) | subsetCommunication*、消费 pathway 表显著性的图 | 保留 pathway 无独立 p 值的事实，表输入是否再筛选须明示；不能把 LR p 均值充当 pathway p | `thresh` 在网络/表/netP 三入口的作用、NA/缺失与既有表字段、signaling/pairLR 筛选的一致性 |
| **A-17** | [社区 assay/marker](../R/analysis.R#L3766-L3786) 将原始通讯 counts 放入 data；[PP1 量纲探针](PP1_DOWNSTREAM_REFACTORING_PLAN.md#7-本轮运行证据与限制) 已观察线性均值与 Seurat 均值/效应量不同 | identifyTopicCommunities、CellChat2Seurat、AverageExpression、markers 与热图 | 先核定原量层、变换、效应量和特征名映射；不只替换 slot，也不默认套 RNA normalization | 线性均值或变换后均值、marker 检验/效应量、Seurat 版本；社区全链仍未验收 |
| **A-18** | [Visium 入口](../R/modeling.R#L738-L798) 仍需旧存储；[混合比例公式](../R/modeling.R#L957-L1007) 与硬标签 avg 不同，PP1 六列 D 对九群 | computeAvgCommunProb_Visium 及 PP1 全部 group/pathway 下游 | 保持混合比例问题，先核定群映射、质量和 null model；迁移真实生产者，不静默替换成普通 avg | 六类/九类定义、截断后是否归一化、分母中存储零、群列置换的科学含义；见 PP1-02 |

## 7. 迁移批次应补充的证据

以下是**未执行的验收设计**，不是本次通过清单。具体决策理由、负责人及批准状态仍留在 [proposed Note](../.agents/notes/proposed/2026-10-20-downstream-communication-semantics-contract.md)。

1. **单位与归约**：同边缘和、不同路由网络同时比较 contribution/rank/role/pattern 与 functional similarity；证明哪些结果应相同、哪些不应相同。cell→group 必须遵循生产者分母与支持门控，不用全张量总和冒充群统计。
2. **阈值与缓存**：覆盖 `p==thresh`、缓存阈值不同于本次阈值、LR 与 netP、group 与 cell、已有表与网络入口；明确哪些调用本来就不能重新检验。无法计算与 NA 不得转换为显著零 p 值。
3. **原始与显示量**：检查 raw/scaled/relative 字段，角色热图三种尺度、贡献的分区分母、dot 的 max 与 sum、群模式与 cell topic 的不同归一化。不能只检查图对象构造成功。
4. **身份对齐与缺失**：置换同名群/层/dataset 的存储顺序，观察数值及标签应不变；分别引入“真实零层”“未生产层”“条件特有层”，裁定补零与省略是否符合该函数问题。DEG 对照身份不能由表行数替代。
5. **统计及模型边界**：网络边不是生物学重复；DEG p 值不是通讯差异 p 值；训练残差不是留出验证；UMAP 距离不是网络原始距离。统计更正、算法更正和纯存储迁移分别记录，失败基线不得作为兼容成功。
6. **运行与图形证据**：真正实施时由集成负责人用最终 11-slot 对象运行受影响入口，观察实际数值、名称、返回表及图形；全程保持稀疏层处理。第 8 节的局部探针不覆盖整包/可选包或图形设备，不得引用本页作为全函数可运行证明。

## 8. 本次运行证据与未验证范围

本节记录 **2026-10-20** 主代理回传的运行结果：**stdin 数值探针退出码 0**；环境为 R 4.4.3、项目 renv 库，Matrix 1.7.4、igraph 2.3.2、dplyr 1.2.1；加载 `R/SpatialCellChat_class.R`、`R/utilities.R`、`R/modeling.R`、`R/analysis.R`、`R/database.R`。以下区分实际函数调用与仅复现源码公式；没有新增 R 实现、测试文件或持久探针脚本。后续 PP1 子链的实际图形观察与新失败证据另见 [PP1 方案的运行证据](PP1_DOWNSTREAM_REFACTORING_PLAN.md#7-本轮运行证据与限制)，不扩大本节探针的覆盖范围。

| 编号与证据类型 | 输入/调用范围 | 主代理报告的实际观察 | 可支持与不能支持的结论 |
|---|---|---|---|
| P-01：实际函数 `computeCentralityLocal` | 三节点正向/反向有向环，正边权重均为 2；[源码](../R/analysis.R#L9-L58) | 两网络不同，各节点 outdeg/indeg 均为 2 | 支持“边缘强度相同不等于路由相同”；未运行所有角色图、群/细胞 producer 路径或相似性入口 |
| P-02：实际函数 `.sc_communication_field_layer` | 中心分别向坐标 `(−1,0)`、`(1,0)` 的邻居发权重 1，`top=1,sparse=FALSE`；[源码](../R/analysis.R#L4016-L4077) | outgoing 总强度 2，净向量 `(0,0)` | 支持“正通讯强度可以产生零净向量”；不是空间流图、widget 或浏览器验收。解释见 [可视化语义卡](VISUALIZATION_SEMANTICS.md) |
| P-03：实际函数 `computeCommunProbPathway` | 单 LR、同一 P1 pathway；两边权重 `0.4/0.6`、p 值 `0.05/0.01`；[源码](../R/modeling.R#L1680-L1767) | `thresh=0.05` 时 group 总量 `0.6`、cell 总量 `1`；`thresh=0.10` 时 group 总量 `1`；`pval=NULL` 时 group 总量 `1` | 支持 group 的 `>=` 边界、无 p 值门控和 cell 分支不同；未调用 netVisual_aggregate。旧可视化 `p>0.05` 公式另算得总量 `1`，仅是公式对照 |
| P-04：实际 helper `extractEnrichedLR_internal` | 旧式 plain-list 隔离输入；唯一 LR1 的 p 值 `0.9`，缓存 LR.sig 含 LR1；[源码](../R/analysis.R#L913-L946) | `thresh=0.01` 与 `0.95` 均保留 LR1；清空 LR.sig 后 selected 数为 0，但 geneLR 仍为 L、R | 支持 A-02 的缓存选择和空集合不一致；不是 wrapper 或最终 11-slot 对象的端到端运行证明 |
| P-05：仅源码归一化公式 | 全正 2×2×2 数组，两层分别全为 1 和全为 2；[对应公式](../R/analysis.R#L296-L308) | raw 份额为 `1/3,2/3`，min–max 后为 `0,1` | 支持 min–max 可改变 LR 份额；**未运行** netAnalysis_contribution 或其图形/返回分区分支 |
| P-06：仅源码保留条件公式 | 差值向量 `(2,−2)`；[对应过滤](../R/analysis.R#L3087-L3089) | “任一方向改变”为 TRUE，但 `rowSums!=0` 保留条件为 FALSE | 支持 A-13 的等幅相反变化丢失；**未运行** netAnalysis_signalingChanges_scatter |

上述退出码只说明该批探针完成，不等于全部 A-* 问题已解决。尚未提供各分析入口的完整迁移回归、NMF/UMAP 实际拟合、rankNet 统计分支、真实返回图形、recordedplot 重放、空间点位方向或浏览器交互证据；这些仍须随受影响实现批次补充。所有卡片的判别场景仍是待实施的验收设计，不能因相近局部机制已被探针观察就宣称整张卡通过。
