`computeCommunProb()` 前的通用准备链是：

```text
原始 counts（或已标准化表达）
  → normalizeData()                         # 仅 raw counts
  → createSpatialCellChat()
  → 设置/筛选 DB（可选 subsetDB()）
  → subsetData()                            # 必需
  → preProcessing()                         # 可选：ALRA 补全
  → identifyOverExpressedGenes()            # 自动特征选择时必需
  → identifyOverExpressedInteractions()     # 生成 LRsig；通常必需
  → projectData()                           # 可选：仅 raw.use = FALSE
  → computeCommunProb()
```

### 必需主链

1. **标准化输入**
   - raw counts：`normalizeData()`；新构造器也可通过 `input.assay = "raw"` 自动完成。
   - 已标准化表达：使用 `input.assay = "norm"`，不要重复标准化。

2. **构造对象**：`createSpatialCellChat()`
   - 需要 expression、`meta`、`group.by`。
   - 空间数据还需 coordinates；有物理尺度时需要 `spatial.factors = list(ratio, tol)`。
   - 当前新 schema 的对应位置：`assay$raw/norm`、`meta`、`idents`、`images`

3. **设置并筛选数据库**

将物种对应的 CellChatDB 写入对象；空间数据通常保留以下 interaction 类型：

```r
chat@DB <- subsetDB(
  CellChatDB.mouse,
  search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
  non_protein = FALSE
)
```

新 object 中数据库仍位于 `DB`，核心字段为 `DB$interaction`、`DB$complex`、`DB$cofactor` 和 `DB$geneInfo`。

4. **提取 signaling genes：`subsetData()`**

这是通常不可跳过的准备步骤。它根据 `DB$interaction`、complex subunits 和 cofactors，从完整表达矩阵提取后续计算所需的基因：

```text
assay$norm → assay$signaling
```

旧版对应关系：

```text
@data → @data.signaling
```

即使使用完整数据库，也建议执行 `subsetData()`，因为后续 LR 表达计算应基于 signaling layer，而不是整个转录组。

空间对象的数据库如果缺少 `annotation`，旧版会将 interaction 默认标记为 `Secreted Signaling`；但正式流程最好在调用前提供正确 annotation，以免把接触依赖信号错误当作扩散信号。

5. **可选表达补全：`preProcessing()`**

该步骤主要用于 ALRA dropout 补全，不是所有数据集的硬性前置条件：

```r
chat <- preProcessing(chat)
```

旧版默认处理 `data.signaling`；新 object 应处理：

```text
assay$signaling
```

如果后续要严格复现某个已有对象，不能无条件启用 ALRA。ALRA 会改变表达值和后续概率；应根据该对象是否已有补全结果、测序深度和复现要求决定。

6. **识别过表达 signaling genes：`identifyOverExpressedGenes()`**

这是自动选择候选 signaling features 的步骤：

```r
chat <- identifyOverExpressedGenes(
  chat,
  selection.method = "wilcox"
)
```

当前基线函数支持的 `selection.method` 为：

```text
"wilcox"    按细胞群做差异表达
"moransi"   基于空间坐标的 Moran's I
"meringue"  基于空间表达模式的 MERINGUE
```

结果不是供 `computeCommunProb()` 直接读取的表达矩阵，而是候选 feature 元数据。新 object 应存放在：

```r
misc(chat, ".var.features")
```

或等价的 `chat@misc$.var.features` 中。空间方法还需要：

```r
images$coordinates
```

如果用户已经有明确的 feature 集，也可以不运行自动统计筛选，直接将该 feature 集传给后续 LR 选择函数；但这属于手动特征边界，不再是完整的自动 over-expression 流程。

7. **识别过表达 LR interactions：`identifyOverExpressedInteractions()`**

该步骤把 signaling features 与数据库 interaction/complex 进行匹配，并生成最终参与概率计算的 LR 表：

```r
chat <- identifyOverExpressedInteractions(
  chat,
  variable.both = TRUE
)
```

输出位置：

```text
LR$LRsig
```

`variable.both = TRUE`：ligand 和 receptor 两侧都必须满足 feature 条件；通常更严格。

`variable.both = FALSE`：允许仅一侧满足过表达条件，适用于希望保留更多候选 interaction 的场景。

complex 不能只检查 complex 名称本身，还必须确认其 subunits 能在 signaling 表达层中解析。最终 `LR$LRsig` 的顺序应保持稳定，尤其是混合了 secreted、ECM 和 cell-cell contact interaction 时。

8. **可选 PPI 平滑：`projectData()`**

只有计划以投影后的表达执行通信概率计算时才需要：

```r
chat <- projectData(
  chat,
  adjMatrix = ppi,
  alpha = 0.5,
  normalizeAdjMatrix = "rows"
)
```

旧版路径：

```text
@data.signaling → @data.project
```

新 object 路径：

```text
assay$signaling → assay$smooth
```

对应关系：

```r
computeCommunProb(chat, raw.use = TRUE)   # 使用 signaling
computeCommunProb(chat, raw.use = FALSE)  # 使用 smooth/projected expression
```

不要为了满足“完整流程”而无条件运行 `projectData()`；`raw.use = TRUE` 时它不是必要步骤。

9. **确认空间推断参数**

`computeCommunProb()` 会根据对象中的空间信息计算距离约束。通常不需要用户手动先调用 `computeCellDistance()`：

```text
images$coordinates
images$spatial.factors$ratio
images$spatial.factors$tol
```

这些信息用于内部计算：

- `interaction.range`：扩散信号允许的最大距离；
- `contact.range`：接触依赖信号的距离阈值；
- `scale.distance`：距离缩放；
- `distance.use`：是否将空间距离作为概率约束；
- `contact.dependent`：是否按 interaction annotation 区分接触信号和扩散信号；
- `contact.dependent.forced`：是否强制所有 LR 按接触依赖模式处理。

手动运行 `computeCellDistance()` 主要用于检查距离矩阵或提前发现坐标、ratio、tol 问题，不是 `computeCommunProb()` 前的必需步骤。

## 最小标准流程

```r
chat <- createSpatialCellChat(
  object = counts,
  meta = meta,
  group.by = "cell_type",
  input.assay = "raw",
  datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = ratio, tol = tol)
)

chat@DB <- subsetDB(
  CellChatDB.mouse,
  search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
  non_protein = FALSE
)

chat <- subsetData(chat)
chat <- identifyOverExpressedGenes(
  chat,
  selection.method = "wilcox"
)
chat <- identifyOverExpressedInteractions(
  chat,
  variable.both = TRUE
)

chat <- computeCommunProb(
  chat,
  distance.use = TRUE,
  contact.dependent = TRUE
)
```

## 含 ALRA/PPI 的可选流程

```r
chat <- subsetData(chat)
chat <- preProcessing(chat)                 # 可选：ALRA
chat <- identifyOverExpressedGenes(chat, selection.method = "meringue")
chat <- identifyOverExpressedInteractions(chat, variable.both = FALSE)
chat <- projectData(chat, adjMatrix = ppi)  # 可选：仅 raw.use = FALSE 时需要
chat <- computeCommunProb(chat, raw.use = FALSE)
```

## 不属于 `computeCommunProb()` 前置步骤的函数

以下步骤应在概率计算之后执行：

```text
filterProbability()
filterCommunication()
computeCommunProbPathway()
aggregateNet()
netAnalysis_computeCentrality()
```

最终依赖关系可以压缩为：

```text
assay$norm
  + DB
  + idents
  + images$coordinates/spatial.factors
      ↓
assay$signaling
      ↓
misc$.var.features
      ↓
LR$LRsig
      ↓
computeCommunProb()
```

因此，最核心的直接前置条件是：`assay$signaling` 已存在、`idents` 有效、`LR$LRsig` 已生成；空间对象还必须有合法的 `images$coordinates` 和空间尺度信息。
