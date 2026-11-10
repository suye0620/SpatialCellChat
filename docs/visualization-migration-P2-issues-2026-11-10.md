# 可视化迁移问题记录：P2 网络图

日期：2026-11-10

## P2-001：aggregate / individual 网络图读取旧概率字段

位置：`netVisual_aggregate()` 和 `netVisual_individual()` 在查询 pathway 对应 LR 后读取 `object@net$prob`、`object@net$pval`；新结构为 `object@net$group$prob` 与 `object@net$group$pval`，每层是稀疏矩阵。

原因和本质：只是替换字段名不够，因为 `SparseChatArray` 不支持旧代码的三维逻辑赋值。旧算法先按调用时 `pval > thresh` 将对应边置零，再保留至少有一条非零边的 LR；aggregate 将保留的 LR 按细胞群对求和，individual 分 LR 绘制。这个边界的 `>` 不能擅自变成 `>=`，也不能直接使用按此前阈值计算的 `netP$group$prob`。

影响：两个入口的 circle、hierarchy、spatial、chord 布局；`netVisual_chord_cell()` 和旧总入口 `netVisual()` 存在相似的读取方式，后续按顺序检查。

修复方向：逐 LR 读取 group 概率和 p 值，保留旧阈值与 LR 排列，再生成小规模 group × group × LR 数组供现有 renderer 使用；不稠密化 cell-level 矩阵。为单层、阈值相等、无 p 值、无匹配和空间布局补测试。

## P2-002：individual 的空间布局引用未定义 labels

位置：`netVisual_individual(layout = "spatial")` 从对象读取坐标后把未定义的 `labels` 传给 `netVisual_spatial()`。

修复方向：与 aggregate 空间布局一致，使用当前对象 `@idents`，不改变绘制算法。

## P2-003：individual 空间布局还引用未定义的 idents.use

`netVisual_individual()` 没有 `idents.use` 形参，但 spatial 分支把该变量传给 `netVisual_spatial()`，调用时报 `object 'idents.use' not found`。沿用 renderer 的默认 `NULL`，不新增过滤语义。

## P2-004：源码直载 runner 未加载 circlize

`netVisual_chord_cell_internal()` 使用通过包 `NAMESPACE` 导入的 `circos.clear` 等符号；直接 source 源码的 runner 未附加 `circlize`。修复测试环境后重新出图，不修改 chord 算法。

## P2-005：Visium individual 多面板超过固定 PDF 画布

位置：`tests_dev/render_viz_figures.R` 为所有基础图固定使用 9×7 英寸的 PDF 设备；`netVisual_individual()` 的 hierarchy 布局按 LR 数量设置行数，circle 和 chord 默认每行两个 LR。

实际现象：Visium 首个 pathway 的 circle、hierarchy 和 chord 报 `figure margins too large`；同一数据的 aggregate 和 individual spatial 通过。绘图设备太小可能导致此错误，但仍需在根据完整 LR 面板数量调整设备尺寸后重新测试，不能据此认定函数已通过。

处理方向：不缩减 LR，不更改统计或绘图布局；测试设备按完整面板行数确定高度。若超过 PDF 设备上限，明确记录，不能把部分面板视为完整验收。

重测：首个 pathway 有 88 个 LR，circle 和 chord 使用 9×154 英寸页面后通过；spatial 实际生成 88 页。hierarchy 按 3.5 英寸/行需要 308 英寸，超过 PDF 页面限制，测试脚本主动报错；继续检查较紧凑的 2.2 英寸/行是否可以完整绘制。

结果：hierarchy 在 9×193.6 英寸页面上完成 88 行绘制。Slide-seq 完整 29,536 细胞对象的 aggregate/individual 四种布局均运行成功，共 8/8；其首个 pathway 有 16 个 LR。页面很大，仍需人工查看图像可读性，运行成功不等于全面视觉验收。

## P2-006：单对象 heatmap 和 barplot 读取旧网络结构

位置：`netVisual_heatmap()`、`netVisual_barplot()` 使用旧的 `net$count/weight`、`net$prob`、`netP$prob` 和旧的 merged 判定；新聚合矩阵位于 `net$group$count/weight`，单层概率位于 `net$group$prob`、`netP$group$prob`。

本质和影响：新对象的 group 本身是 list，旧的 `is.list(object@net[[1]])` 误判为 merged；`options$mode` 新对象通常不存在。heatmap 的 LR 模式仍须保留 `pval > thresh` 剔除的旧边界，barplot 原先没有这个门控，不应顺手添加。comparison 的生产者缺失是独立问题，不能声称已验收。

额外缺陷：heatmap 的 isolate 索引删除两次，`cluster.cols` 未传给列聚类；barplot 在未指定来源和目标时未创建 `df.net`。修复后必须验证单对象默认调用、各概率层、阈值边界及绘图。

## P2-007：bubble 单对象判断和强度变换引用错误

位置：`netVisual_bubble()` 用 `is.list(object@net[[1]])` 判断 merged，会将新对象的 group 分支误判；强度变换遇到无穷或负值时引用未定义的 `prob.original`，实际数据在 `df.net$prob.original`。

影响：误报对象模式；概率达到或超过 1 时该分支报错。先保留原变换和排序规则，修正数据来源；comparison 另行核验。

首轮迁移测试补充：新类完全没有 `options` slot，不只是缺少其中的 `mode`。最初沿用 `object@options$mode` 的判断导致三个入口均报 `no slot of name "options"`；现改为读取 `object@misc$mode`，须重新验证。

更正：已有 `subsetCommunication()` 使用 `object@misc$.mode`，不是 `misc$mode`；入口已对齐这一既有约定。

## P2-008：heatmap 单一正值时把颜色函数拼接成字符串

模拟对象默认 count 网络实际绘制时报 `invalid color name 'function ...'`。只有一种正值时，代码对 `circlize::colorRamp3()` 返回的函数执行 `c("white", color.heatmap.use)`，将函数转成文本。连续颜色函数本来就能映射网络中的零值，应仅对颜色向量执行原有的拼接。

## P2-009：bubble pathway 筛选依赖 searchPair 的错误列名

模拟对象的 `LR$LRsig` 前四列为 ligand、receptor 等；`searchPair(pair.only = TRUE)` 选出了 interaction_name、pathway_name、ligand、receptor 四列，随后却用输入表的前四列名重命名。`subsetCommunication_internal()` 按结果的 interaction_name 筛选时可能取到错误列或错误值，bubble 因此报告无通讯。应在该调用点直接依据 `LR$pathway_name` 选择 interaction_name，保持 pathway 和阈值原语义；其他 `searchPair` 调用方需单独检查。

模拟 fixture 缺少展示名称 `interaction_name_2`，该列此前被补为 NA；改为缺列时使用已有的 `interaction_name` 作为显示名，不改变通讯值和筛选。

## P2-010：Visium bubble 返回对象但 PDF 渲染失败

真实 Visium 首个 pathway 的 `netVisual_bubble()` 调用返回 ggplot，但 `ggsave()` 报 `missing value where TRUE/FALSE needed`。错误发生在实际绘图阶段，尚未定位；该项维持 FAIL，不能根据返回对象类型判定成功。

定位：Visium fixture 未进行 permutation，全部 2344 个通讯记录的 pval 都是 NA；旧代码仍用 `range = c(min(pval), max(pval))` 创建点大小比例尺，导致渲染报错。另一个独立问题是未传 `pairLR.use` 时纵轴顺序由空名单构造，全部 LR 标签被转换成 NA。修复方向：仅当完全没有 p 值时使用固定点大小、不显示 p 值图例，不将 NA 伪装成显著；保留实际 LR 标签。混合缺失 p 值的绘图语义仍需专项验证。

首次保存成功后检查 PDF 文本，发现 88 个 LR 行挤在 5.5 英寸高度，标签严重重叠。测试输出需按实际行列数扩展画布，再检查可读性；保存成功本身不能作为视觉验收。

## P2-011：Slide-seq 整批 runner 运行期间脚本被修改

全量 Slide-seq 的 11 项均在输出中显示 PASS，但测试脚本仍在运行时被编辑，进程末尾出现 `unexpected ','` 并以非零码退出。该批记录不作为最终通过依据；应在脚本稳定后重新执行。单个图像的存在也不代替完整运行结果。

## P2-012：chord cell 入口读取旧概率结构

位置：`netVisual_chord_cell(signaling = ...)` 使用已删除的 `net$prob/pval`。其 `netP` 模式的旧行为是按本次阈值过滤 LR 后实时求和，而不是直接使用预计算 pathway 概率；`net` 模式则逐 LR 画图。迁移应复用已验证的 group LR 筛选，不更换上述行为。

## P2-013：chord gene 入口读取旧三维数组

`netVisual_chord_gene()` 用 `net$prob/pval` 生成完整稠密三维表，新对象使用 group 的 SparseChatArray。必须逐层处理 group 矩阵并保持原 `pval > thresh` 门控，按 interaction_name 关联 LR 信息。不能直接调用使用 `pval >= thresh` 的 `subsetCommunication()`，否则阈值相等时结果不同。原先经 `searchPair()` 的 pathway 筛选也存在 P2-009 的列名问题。

测试注册过程中曾遗漏两个右括号，导致 runner 解析失败，尚未进入绘图函数；已修正，需重跑。

稳定重跑：Visium 与全量 Slide-seq 的当前注册项均为 13/13，进程退出码 0。人工预览发现 Visium 的 88 LR pathway 在 9×7 英寸 gene chord 页面上，sector 标签显著重叠；运行成功不表示其视觉可读性通过。结果保留在 `tests_dev/figures/P2/visium/` 供检查，后续需要合理处理高密度标签，不缩减数据冒充完整绘图。

对照：保留完整 pathway 与权重，使用 20×20 英寸页面、`lab.cex = 0.4`、`small.gap = 0.25` 另生成 `netVisual_chord_gene-adjusted.pdf`。预览显示部分标签仍重叠，因此该参数调整不能视为完全修复，默认图的状态仍为 `VISUAL-REVIEW`。

## P2-014：旧总入口 netVisual 的导出分支仍依赖旧字段

`netVisual()` 使用 `net$prob/pval` 以及稠密三维数组赋值；spatial 分支还引用未定义的 labels、idents.use 和 edge.weight.max，且调用返回 ggplot 的 renderer 后没有显式打印。该函数直接管理 SVG、PNG 和 PDF 设备，不能仅验证返回值；修复时须核查每种格式和所有布局的实际输出与设备清理。

处理：数据选择复用逐 LR 的既有阈值逻辑；spatial 分支按全部 LR 构造 ggplot 列表并组合，再导出 individual/aggregate 文件，分别使用对应边权上限。对超多 LR 的大尺寸位图和 PDF 设备限制仍需专项验证。

## P2-015：merged comparison producer 与空间坐标语义

状态：comparison producer 已实现，不再是“缺少生产者”。`mergeSpatialCellChat()` 将各输入的通讯结果保留在 `misc$.datasets` 快照，并同步带前缀的 cell 轴；真实 producer 的模拟对象回归覆盖 `rankNet`、`subsetCommunication`、差异图、heatmap/barplot/bubble、pairwise embedding、signaling-role analysis，以及 subset/lift/update 生命周期。

仍须注意两个边界。第一，`count.merged`/`weight.merged` 是按 `mergeInteractions()` 提供的 group 映射聚合出的独立矩阵，不是普通 `count`/`weight` 的别名；输入对象应先分别调用 `mergeInteractions()`，缺少聚合矩阵时 comparison 会报错。第二，空间数据集的坐标通常是各自局部坐标。merged 顶层不再拼接这些坐标，而以 `images$coordinate.scope = "dataset-local"` 标记并将坐标保留在快照；空间图入口明确拒绝把 merged 容器当作统一组织平面。

当前 producer-specific 回归使用两个独立生成的 simulated objects。它验证新 producer 和 comparison 消费者的契约，不等同于三个原始测试 object（simulated、Visium、Slide-seq）的 merged biological comparison。三套 single object 的既定绘图迁移产物仍按各自测试结果检查；需要共同空间坐标的 merged spatial plot 不在当前支持范围。

验收记录位于 `tests_dev/test-merge-spatialcellchat.R`。重跑该测试、原 analysis/P2/P3/P4 迁移脚本、`R CMD INSTALL` 和 `git diff --check`；若环境缺失 R 或某项失败，须保留为显式未验收项。

## P2-016：绘图函数关闭整个 R 会话的警告

hierarchy、circle、spatial、diffInteraction 等 renderer 调用 `options(warn = -1)` 后不恢复，导致后续操作的警告也被隐藏。迁移期间移除这些全局选项更改，保留绘图和统计计算；测试需重新运行以观察此前可能被掩盖的警告。

## P2-017：contribution 和 rankNet 的稀疏 slot 迁移

位置：`R/analysis.R` 的 `netAnalysis_contribution()` 和 `rankNet()`。此前这两个入口分别读取旧的 `net$prob` / `net$pval`、`net$prob.cell` 以及旧 comparison 容器；新对象的概率现在从 `net` 或 `netP` 下的 `group$prob` / `group$pval` 和 `cell$prob` 逐层读取。

本质：单对象路径是存储迁移，不能因为函数返回 ggplot 就自动视为完成；还必须确认 LR 选择、群级 p 值门控、归一化、weight/count 语义、receiver 分面和 cell-level 稀疏性没有被改变。comparison 路径则需要真实的 merged object producer 来保证 dataset、group、LR/pathway 和 cell 轴的名称对齐。迁移总表曾将 `rankNet` 描述为 bubble 的依赖，但源码中没有发现 `netVisual_bubble()` 直接调用 `rankNet()`，不应据此推断 bubble 单对象结果失效。

当前结果：`netAnalysis_contribution()` 的 group、cell 和 receiver 分面，以及 `rankNet()` 的 `netP` weight/count、`net` cell count 已切换到 `SparseChatArray` 新 slot；三套真实 fixture（simulated、Visium 905 spots、完整 Slide-seq 29,536 cells）均已成功出图。回归脚本为 `tests_dev/test-analysis-network-migration.R`，产物和摘要分别位于 `tests_dev/figures/P2/<dataset>/` 与 `tests_dev/figures/summary-all-P2-*.csv`。

当前结果：single-object 数值路径保留原有 LR 选择、p-value 门控、归一化及 count/weight 语义；comparison 路径现由 `mergeSpatialCellChat()` 提供真实快照，并在 `tests_dev/test-merge-spatialcellchat.R` 中通过 producer 驱动回归。该回归是两个 simulated object 的合成对比，不代表真实样本间的生物学比较，也不替代三套 single-object 图像验收。

剩余边界：对真实 Visium/Slide-seq 数据是否存在可解释的共同空间坐标系，不能由相同单位或像素比例推断。merged 空间对象因此保留 dataset-local 快照并拒绝顶层空间绘图；非空间 comparison 可正常消费 merged 快照。


## P2-018：merged 差异热图的正负变化相互抵消

位置：`R/visualization.R::netVisual_heatmap()` 的空结果判断及 `remove.isolate` 分支。merged 差异是第二个数据集减第一个数据集，既可正也可负。旧代码使用 `sum(mat) == 0` 判断整张图无变化；例如 A→A 增加 1、A→B 减少 1 时总和为零，尽管两条变化都存在仍报错。同理，`remove.isolate = TRUE` 用行列的带符号和判断群是否孤立，可能删除收发变化互相抵消的群。

修复：只有每个元素都为零时才按空矩阵处理；只有对应行和列都没有非零元素时才删去孤立群。差值、方向、阈值和色阶均未重算。全零检查提前至色阶构造之前，使全零输入返回原有的无通讯提示，而不是先因色阶断点重复报错。

验证：`tests_dev/test-merged-heatmap-edgecases.R` 使用真实 `mergeSpatialCellChat()` 生产的模拟对象，覆盖 +1/−1 总和抵消、纯负变化、四条变化的行列抵消及全零；在热图矩阵中核对非零数值。对模拟、Visium、完整 Slide-seq 的单对象热图重新出图并核对非零群矩阵，汇总见 `tests_dev/figures/P2/heatmap-cancellation-regression-2026-11-10.csv`。比较边界图保存在 `tests_dev/figures/merged/heatmap-edgecases/`。

## P2-019：merged 比较图实绘与图例待核对

此前 producer 回归对一些比较入口仅断言返回类型，未保存实际图。本次基于模拟对象及其经过明确修改的副本，使用真实 `mergeSpatialCellChat()` 生产快照，输出 10 张比较 PDF 至 `tests_dev/figures/merged/producer-comparison/`；这不是独立生物样本的比较，也不是用户的视觉验收。

另观察到 merged `netVisual_heatmap()` 的差值直接来自第二数据集减第一数据集，图例仍固定写 `Relative values`；当所有差值非负时，默认两色映射将零映射为第一种颜色。该说明和色阶可能令读者误解，需要结合已有绘图约定单独评估，本次没有改变配色或差值定义。
