# 可视化迁移问题记录：P3 场与分解

日期：2026-11-10

## P3-001：单细胞通讯场图读取旧 field 路径

位置：`netVisual_CommunField()` 从 `slot$field[[pattern]]` 取值，新版 `computeCommunField()` 将 `SparseChatArray` 写入 `slot$cell$field[[pattern]]`。原代码还尝试把所有层转换为三维对象，而绘图实际只需要一个 signaling 层。

影响：基础场图和调用它的 `netVisual_CommunFieldZoomIn()`；现有的 grid 和 flow 使用新路径。修复只读取选定的细胞×方向矩阵，保持箭头方向、筛选与绘制逻辑，检查细胞名与坐标对应；必须使用真实生产者计算后的字段进行实际绘图测试。

首次测试未运行：项目 renv 库缺少 `metR` 与 `oce`，已有 grid 测试报告 SKIP。两者均已在 DESCRIPTION Suggests 中声明；需安装后重新执行，不把 SKIP 算作通过。

补充：安装依赖后，9 细胞专项 fixture 的基础场图与放大图均由 `computeCommunField()` 生成真实输入并成功保存 PDF；但产生 `No shared levels found` 颜色尺度警告，须检查有无细胞群被错误着色。项目 lockfile 尚未记录新装依赖，环境可复现性仍待处理。

## P3-002：pattern 和 similarity 上游仍使用旧 group 概率与 mode

`identifyCommunicationPatterns()`、`computeNetSimilarity()` 从已删除的 `slot$prob` 读三维数组；`netEmbedding()`、`netClustering()` 访问已删除的 `options` slot。结果是 river/dot 和 embedding 即使只修绘图入口，也无法通过现有生产链生成输入。应将 group 的各 SparseChatArray 层按原轴顺序还原为小规模 group×group×signal 数组，保留 NMF、相似度、SNN、UMAP 和聚类的原算法；mode 对齐 `misc$.mode`。这些迁移不应将 cell×cell 稀疏矩阵稠密化。comparison producer 与消费者现已在 P2 回归中覆盖；空间图仍受 dataset-local 坐标约束。

## P3-003：单对象 embedding 绘图读取旧概率，且可能按位置错配

`netVisual_embedding()` 和 `netVisual_embeddingZoomIn()` 用旧的 `slot$prob` 汇总每条 pathway 的点大小，并假定概率层顺序、嵌入坐标和聚类结果完全相同。新版应从 `group$prob` 求每层总强度，按嵌入坐标中的 pathway 名称对齐；点大小仍为各层总强度除以最大值。缺失名称或无法对齐时应报明确错误，不能静默错配。

## P3-004：源码直载测试未加载 ggalluvial

模拟对象的 NMF producer 已计算出 pattern，但 river 测试在 `to_lodes_form()` 处报未找到函数。该函数通过包 NAMESPACE 导入，源码直载测试须显式附加 `ggalluvial` 后重试；此现象本身尚不足以认定包安装后的入口有同样错误。

## P3-005：模拟对象的 embedding 聚类数边界

实际 producer 链在模拟对象上运行到 k-means 时，测试固定使用 `k = 2` 报 `number of cluster centres must lie between 1 and nrow(x)`。需检查实际 embedding 的 pathway 数量，按数据规模使用合法 k；若只能使用一个 cluster，要明确该测试不能覆盖多 cluster 分面。

重测：模拟对象只有 2 条 pathway，使用 k=1 后两个 embedding 图可绘制；该用例不覆盖多 cluster。

## P3-006：通讯场放大区域边界和颜色尺度

`netVisual_CommunField()` 的 `cut()` 默认不包含最小坐标，使位于坐标边界的细胞在 zoom 区域被误排除；颜色映射使用仅含有非零箭头所属群的 palette，但背景点可以属于其他群，导致其颜色缺失。保持原区域分割和箭头算法，仅让最小边界纳入并保留全部点的群颜色；需实际渲染确认。

边界修正后专项测试暴露另一问题：选中区域有细胞但没有非零向量时，向 `metR::geom_arrow()` 传入零行数据会在绘图阶段报 `invalid 'arrow' argument`。该区域应保留背景点而不添加空箭头层，不能人工制造向量。

## P3-007：TopicComposition 元数据 slot 和孤立项轴

`netVisual_TopicComposition()` 直接访问 `object@meta.data`，新 SpatialCellChat 使用 `@meta`；`remove.isolate` 查找低比例行却删除列。迁移需按对象类型读取元数据并删除对应行，不能更改均值、归一化或 topic 权重计算。

## P3-008：Seurat TopicSignaling 引用未定义变量

`netVisual_TopicSignaling()` 已取得指定 reduction 的 `feature.loadings`，但紧接着又用未定义的 `reduction` 与 `projected` 重新取值。保留前者并去除错误引用，不更改载荷排序。

进一步检查发现显式 `db` 分支按 interaction data.frame 检查列，但之后作为 `CellChatDB$interaction` 使用；应先补足旧行为要求的名称列，再按 interaction name 建立索引并包装成数据库列表。测试须使用真实 Seurat reduction 与实际 PDF，而不是仅断言代码不再引用未定义变量。

验证：使用合成 Seurat 对象的实际 reduction 与 interaction data.frame 生成 PDF，输出在 `tests_dev/figures/P3/simulated/netVisual_TopicSignaling-Seurat.pdf`；此测试不代表已覆盖三套真实 Seurat 数据。

## P3-009：pattern 和 embedding 上游永久关闭全局警告

`selectK()`、`identifyCommunicationPatterns()`、`netEmbedding()` 与 `netAnalysis_dot()` 调用 `options(warn = -1)` 且不恢复，使后续绘图和统计过程的警告在整个 R 会话中消失。删除这些全局更改，不改变计算；还需检查其他绘图入口的同类问题。

## P3-010：真实稀疏场无法转换，且位置未校验

`netVisual_CommunFieldGrid()` 与 `netVisual_CommunFlow()` 直接将选中场层转换为 data.frame。模拟对象的 Grid 和 Visium 的 Flow 实测在 `dgCMatrix` 上报错；先前小型测试的数据恰好没有暴露此问题。两函数也未验证场向量与坐标的细胞顺序，相同长度但顺序不同时会将方向对应到错误细胞。修复方向是只将选中的细胞×2 方向矩阵转为普通矩阵，并与基础场图一致地检查行名顺序；不改变网格聚合、向量方向或流线算法。三对象真实绘图仍需重测。

## P3-011：全量 Slide-seq 网格场图耗时

全量 29,536 细胞、LAMININ outgoing、默认网格参数生成 122,138 个网格（28,821 个非空），实际 PDF 成功，但本机本次绘制耗时约 273 秒，且出现 `st_centroid` 属性警告。该结果不是性能达标证明；需要进一步剖析空网格与插值/几何计算的成本后再决定如何优化，不能在迁移中未经验证地改变网格或场的统计含义。

## P3-012：全量场图可读性与流线积分

全量 Slide-seq 默认参数的 Grid 和 Flow 已保存 PDF，但预览显示高密度背景点遮挡场向量；Flow 额外警告仅执行一步积分，预览主要呈现直线而非可辨识的流线。两项状态修正为 `VISUAL-REVIEW`，不计为验收通过。应检查尺度与插值参数，并在不更改场计算意义的前提下调整呈现；不得通过丢弃细胞或伪造更多积分步数掩盖问题。对应预览保存在 `tests_dev/figures/P3/slideseq/`。

对照：保留全部 29,536 个细胞，设置 `point.size = 0.22`、`image.alpha = 0.4`，Flow 另设 `res = 8`。Flow 对照不再出现仅一步积分警告，预览中曲线较默认图清楚，但仍需检查流线密度和解释；Grid 对照仍只有少量明显箭头，需核对场值与网格显示密度。两张对照图仅证明对应 PDF 可以写出，单项汇总均标为 `VISUAL-REVIEW`；不改变默认图的状态，也不据此认定问题已解决。


## P3-013：通讯流浏览器测试与网格抽样诊断

浏览器 CDP 实测模拟、Visium、完整 Slide-seq 三套通讯流 widget 的初始化、细胞数量、动画移动、暂停、恢复及粒子数、透明度、缩放控件，结果见 `tests_dev/figures/P3/commun-flow-widget-browser-summary-2026-11-10.csv`。截图是在修改控件后捕获；hover、窗口尺寸变化及方向的科学正确性不在本次覆盖范围。

Slide-seq LAMININ outgoing 的诊断中，29,536 个细胞有 19,474 个非零场向量；插值网格 2,401 个点中有 540 个非零。默认 skip=4 时统计层保留 100 行，不跳格保留 2,401 行；保留行可能含零向量，不能直接当作可见箭头数量。无跳格对照见 `tests_dev/figures/P3/slideseq/netVisual_CommunFieldGrid-no-skip.pdf`。对照预览仍有大量细小点且箭头不普遍清晰，不能把跳格视为唯一原因。这是显示抽样和场稀疏性诊断，不证明默认图已可读，也不能证明矢量的生物学正确性。


## P3-014：网格向量聚合重复提取 sf 对象

位置：`netVisual_CommunFieldGrid()` 和 `netVisual_CommunFlow()` 中按网格聚合 dx/dy 的环节。原实现为每个网格、每个分量分别对完整 `sf` 对象取行；完整 Slide-seq 生成 122,138 个网格，包含 93,317 个空网格，产生大量重复的空间对象子集开销。现在仅在此步骤预先提取 dx/dy 数值列，再按原始成员索引调用同一个 `median` 或 `sum`；网格、成员关系、空网格、缺失值、插值和向量方向均未改动。

验证：`tests_dev/test-grid-visualization.R` 对空集合、零值、缺失值、无穷值和重复成员逐项比较原 sf 子集法，median 与 sum 均完全一致，并实绘两个入口。`tests_dev/benchmark-grid-field-aggregation.R` 的同进程单次合成对照（2,000 点、6,000 网格）中，median 聚合旧/新分别约 3.58/0.12 秒，sum 约 3.22/小于计时分辨率的 0.01 秒；仅为此聚合步骤，不是整个绘图的加速倍数。

完整 Slide-seq 的 LAMININ outgoing 本次 Grid 和 Flow 对象构造分别约 23.97 和 23.79 秒，未包含字段计算和 PDF 输出。对照 PDF 以相同 `pdftoppm` 参数转为 PNG 后，优化前后 Grid 与 Flow 的 PNG 哈希均相同；对应原始 Grid 对照和 Flow 对照分别见 `tests_dev/figures/P3/slideseq/netVisual_CommunFieldGrid-no-skip.pdf` 和 `tests_dev/figures/P3/slideseq/netVisual_CommunFlow-adjusted.pdf`。旧的约 273 秒记录覆盖范围不同，不能计算整体加速倍数。显示问题和 `st_centroid` 警告仍待处理。
