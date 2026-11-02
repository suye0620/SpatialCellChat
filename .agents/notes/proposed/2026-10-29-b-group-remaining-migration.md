# Agent Note: B 组剩余 6 个函数的迁移方案

Status: proposed
Governance: v1
Date: 2026-10-29
Decision type: migration
Scope: R/utilities.R (`setIdent`, `updateClusterLabels`, `identifyOverExpressedLigandReceptor`), R/modeling.R (`computeAveExpr`), R/spatial.R (`communicationDistPlot`, `communicationDistPlot2`), R/SpatialCellChat_class.R (重复的 `setIdent`), NAMESPACE, man/, tests_dev
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## 1. 现状

B 组 15 个函数中 9 个已可用（`projectData`、`asSeurat` 及其链路）。剩余 6 个全部读已删字段，因此**当前调用即报错**，只是因为零调用者所以没暴露。

| 函数 | 位置 | 读的旧字段 | 包内/测试调用者 | 导出 |
|---|---|---|---|---|
| `setIdent` | `R/utilities.R:994` | `net$prob` / `net$pval`（L1005–1021） | `updateClusterLabels`（2 处） | 是 |
| `setIdent` | `R/SpatialCellChat_class.R:1326` | — | 无（被上一份覆盖） | 同名 |
| `updateClusterLabels` | `R/utilities.R:1044` | `@meta` / `@idents`（L1046–1062） | 无 | 是 |
| `computeAveExpr` | `R/modeling.R:2269` | `slot(object, "data.signaling" / "data")`（L2270、L2279）、`@idents` / `@meta`（L2290、L2296） | 无 | 是 |
| `identifyOverExpressedLigandReceptor` | `R/utilities.R:1531` | `@var.features`（L1541–1544、L1583） | 无 | 是 |
| `communicationDistPlot2` | `R/spatial.R:11` | `@options$parameter`、`@images$.distance`、`net$prob.cell`（L18、L23、L32） | 无 | 是 |
| `communicationDistPlot` | `R/spatial.R:95` | 同上（L100、L110、L128、L131、L136） | 无 | 是 |

**字段映射基础**（已核实新容器的能力）：

| 旧 | 新 | 容器 |
|---|---|---|
| `net$prob` / `net$pval` | `net$group$prob` / `net$group$pval` | `SparseChatArray`（group × group × layer） |
| `net$tmp$*` | `net$group$count` / `net$group$weight` | group × group 矩阵 |
| `net$centr`（组级） | 同名 | measure × group × layer 数组 |
| `@net$prob.cell`（组级） | `net$cell$prob` | `SparseChatArray`（cell × cell × layer） |
| `netP` 下同理 | `netP$group$*`、`netP$centr` | 同上 |
| `@meta` / `@idents` | `meta()` / `idents()` | data.frame / factor |
| `slot(object, "data.signaling")` | `assay(object, "signaling")` | dgCMatrix |
| `@var.features` | `misc$.var.features` | list |
| `@options$parameter` | `misc$.param$communication` | list（含 `interaction.range`、`all.diffusible`、`all.contact.dependent`、`nLR1`、`layer.names`） |

`SparseChatArray` 已支持 `[`（含 `drop`）、`[[`、`dimnames<-`、`dim<-`、`names<-`、`t`、`marginSums`，所以下标重排和改名都不需要新增容器方法。

## 2. 逐函数方案

### 2.1 `setIdent`（迁移）——本批唯一的重活

**为什么不能只换个字段**：改 `idents` 等于改分组轴，对象里有 5 类东西挂在这个轴上，漏掉任何一类都会让 `validObject()` 失败（validator 会拿 `net$group$prob`、`net$group$count` 的 dimnames 与 `levels(idents)` 比对），或者更糟——落进 validator 不管的 `centr`，静默错位。

要在 `net` 和 `netP` 两个槽里同步：

| 层 | 形状 | 重排 | 改名 |
|---|---|---|---|
| `group$prob` / `group$pval` | `SparseChatArray` | `arr[idx, idx, ]` | `dimnames(arr) <- list(new, new, old[[3]])` |
| `group$count` / `group$weight` | 矩阵 / dgCMatrix | `mat[idx, idx, drop = FALSE]` | `dimnames(mat) <- list(new, new)` |
| `centr` | measure × group × layer | `centr[, idx, , drop = FALSE]` | `dimnames(centr)[[2]] <- new` |

抽成一个内部函数 `.sc_reindex_groups(object, idx = NULL)`：`idx` 为整数向量时重排，为 `NULL` 时改名。

判定逻辑：

```
if (identical(old.levels, new.levels))            什么都不做
else if (setequal(old.levels, new.levels))        idx <- match(new, old); 重排
else                                              改名（要求组数相同，否则先 stop）
```

结尾统一 `validObject(object)`。这样"半更新"的对象出不了这个函数。

同时删除 `R/SpatialCellChat_class.R:1326` 的重复定义：C locale 下 `SpatialCellChat_class.R` 先于 `utilities.R` 加载，它自类重写起就被覆盖，签名还不同（`ident.use = "ident"`），留着只会误导。

**风险**：`setequal` 判定为改名时，若新旧组数不同，`dimnames<-` 会报一个不直观的长度错误。方案是在动手前显式检查并给出"组数必须不变"的错误。

### 2.2 `updateClusterLabels`（迁移，随 2.1）

改动点：

- `object@idents` → `idents(object)`，`object@meta[[...]]` → `meta(object)[[...]]`；
- `new.cluster.metaname` 写入后，`as.factor(meta[[col]])` 会丢掉 `names`，而 validator 要求 `names(idents) == rownames(meta)` —— 必须手动补 `names(value) <- rownames(md)`；
- 尾部 `computeCommunProbPathway` → `aggregateNet` → `netAnalysis_computeCentrality(slot.name = "netP")` 三个调用保持不动（这三个函数已迁移）；
- 保持 `setIdent(display.warning = FALSE)`，避免连锁警告。

### 2.3 `computeAveExpr`（迁移）

- `slot.name` 的取值从 `c("data.signaling", "data")` 改为 `c("signaling", "norm")`；
- 旧值仍接受，但走一次 `warning()` 映射（`data.signaling → signaling`、`data → norm`），避免已有脚本直接断；
- `slot(object, slot.name)` → `assay(object, slot.name)`，空层报错时点名 `subsetData()`；
- `group.by` 分支的 `object@meta[[group.by]]` → `meta(object)[[group.by]]`。

### 2.4 `identifyOverExpressedLigandReceptor`（迁移）

- 3 处 `object@var.features` → `object@misc$.var.features`（读 `features.info`、写 `features.LR`）；
- 顺带修两个既有 bug（与迁移无关，但会一起动这段代码）：
  1. `features` 非 NULL 时，`markers.all` 还没赋值就被 `subset()` 使用；
  2. `for (i in 1:nrow(...))` / `1:length(...)` 在空输入时迭代两次（`1:0`），改成 `seq_len()` / `seq_along()`；
- `DB$interaction` 缺 `ligand`/`receptor` 时给明确错误；
- 结尾 `validObject()`。

### 2.5 / 2.6 `communicationDistPlot`、`communicationDistPlot2`（迁移）

> 2026-10-31 起本节由 `2026-10-31-spatial-c2c-dist-plot.md` 接管：两个函数连同 `spatialCCCDistPlot` 一并三合一为 `spatialC2CDistPlot()`，下面"只迁移不合并"的方案作废。

**先更正上一版的判断：这两个不需要新增功能，是我把两个事实记错了。**

| 上一版说的 | 实际 |
|---|---|
| `interaction.range` 只剩在 `images(object, ".distance")$.parameters` | `misc$.param$communication$interaction.range` 有存（`R/modeling.R:460`） |
| `all.diffusible` / `all.contact.dependent` 从未存进对象 | 两个 flag 都在 `misc$.param$communication`（`R/modeling.R:463–464`），`nLR1`、`layer.names` 也在 |

`object@images$.distance` 这一行**根本不用改**：`images` 槽还在，`.distance` 仍是它的键，只是可以顺手写成 `images(object, ".distance")`。

真正要改的只有三类：

1. **参数来源** —— `object@options$parameter[[x]]` → `misc$.param$communication[[x]]`；缺失时报错并点名 `computeCommunProb()`。
2. **概率数组** —— `object@net$prob.cell` → `object@net$cell$prob`；spill 模式下后者是 NULL，单独报错。
3. **三元组** —— `prob.cell$i` / `$j` / `$k` 已不存在，用一个 helper 从 `SparseChatArray` 每层重建：

```r
.sc_cell_prob_triplets <- function(prob) {
  layers <- unclass(prob)
  do.call(rbind, lapply(seq_along(layers), function(k) {
    s <- summary(layers[[k]])          # dgCMatrix -> data.frame(i, j, x)
    data.frame(i = s$i, j = s$j, k = k)
  }))
}
```

`communicationDistPlot` 里"可扩散 / 接触依赖"的两分桶，旧代码用 `max(which(annotation %in% 可扩散))` 自己算阈值；迁移时优先读 `misc$.param$communication$nLR1`（正是 `computeCommunProb` 当时用的那个数），缺失时再退回按 `LR$LRsig$annotation` 逐层判定。

这两个函数**保留导出**，不改 NAMESPACE，不动 `man/`。

## 3. 需要你拍板的五件事

| 编号 | 问题 | 我的建议 | 影响 |
|---|---|---|---|
| **D1** | `communicationDistPlot` / `2` 归档还是迁移 | **迁移**（推翻上一版的归档建议） | 上一版误以为数据形态变了；实际参数都在 `misc$.param$communication`，只剩三元组要重建，约 60–80 行 |
| **D2** | `setIdent` 遇到组数变化（如 4 组→3 组）怎么办 | **直接报错拒绝** | 现有对象结构无法表达"合并两个组"；拒绝比静默错位安全 |
| **D3** | `setIdent` 是否同步 `net$centr` | **同步** | 不同步则 `centr` 的 dim2 与 idents 脱节，下游按名取会错 |
| **D4** | `computeAveExpr(slot.name=)` 旧值 `"data.signaling"`/`"data"` | **保留但告警** | 直接删更干净，但会让旧脚本硬断 |
| **D5** | 距离图的横轴上限 / 分桶阈值来源 | `nLR1` + `interaction.range` 都读 `misc$.param$communication` | 读不到就报错点名 `computeCommunProb()`，不要静默给个默认值 |

## 4. 测试计划

新增 `tests_dev/test-labelling.R`，用真实 pipeline 构造对象（20 细胞 4 组、3 对 LR）：

```
createSpatialCellChat → 挂 DB/LRsig → computeCommunProb → computeAvgCommunProb
                     → computeCommunProbPathway → aggregateNet
```

断言覆盖：

- **setIdent 重排**：`idents` 层级顺序、`net$group$prob` / `netP$group$prob` 行列名、`count` 矩阵行列名、`centr` 的 dim2、按组名取到的值不变、`validObject` 通过；
- **setIdent 改名**：新标签写入、层名（dim3）不变、数值不变、`netP` 同步；
- **setIdent 守卫**：组数变化报错、`ident.use` 不存在报错、标签含 `0` 报错、无网络的裸对象也能用；
- **updateClusterLabels**：新列进 `meta`、`idents` 换成新名、三层结果被重建、`validObject` 通过；
- **computeAveExpr**：默认读 `signaling`、`norm` 层可读、与手工分组中位数一致、旧值告警、空层报错点名 `subsetData`；
- **identifyOverExpressedLigandReceptor**：写进 `misc$.var.features$features.LR`、复合体亚基归到复合体名下、DB 外基因被丢弃、`features` 过滤生效、缺 `misc$.var.features` 时报错；
- **communicationDistPlot / 2**：在同一个真实对象上，`enriched.only = TRUE` 和 `= FALSE` 各画一次，返回 `ggplot` / `patchwork` 且不报错；`signaling.type` 取 `All` / `Secreted` / `Contact` 都能出图；`net$cell$prob` 为 NULL（spill 模式）时报错；缺 `misc$.param$communication` 时报错点名 `computeCommunProb()`；
- **导出状态**：四个迁移函数的 export 仍在，`man/` 下的页面没动。

回归另跑：seurat-bridge 41、projectData 19、computeCommunProb 61、computeAvgCommunProb 162、filterCommunication 22、normalizeData 9、pre-inference 15、class 22、accessors 22、cpp-prob-layer 20、helpers 16、SparseChatArray 41。

## 5. 提交拆分

| 顺序 | 内容 | 说明 |
|---|---|---|
| 1 | `setIdent` + `.sc_reindex_groups` + 删 class.R 重复定义 + `updateClusterLabels` | 一个逻辑单元：分组轴一致性；附 test-labelling 的对应断言 |
| 2 | `computeAveExpr` + `identifyOverExpressedLigandReceptor` | 两个小改，可独立回滚 |
| 3 | 迁移 `communicationDistPlot` / `2` + 新增 `.sc_cell_prob_triplets` | 纯 `R/spatial.R` 内改动，不摘 export；和前面两步无耦合，可单独回滚 |
| 4 | 决策笔记转 implemented + `docs/FUNCTION_INVENTORY.md` 更新 | 收尾 |

## 6. 验收标准

- B 组 15 个函数在 `R/` 下**不再有任何一处**读 `net$prob` / `net$pval` / `@var.features` / `@options$parameter` / `net$prob.cell` / `@data.signaling`；
- 六个函数的导出状态不变（`NAMESPACE` 不动、`man/` 不动），不留 `R/legacy/` 归档；
- `tests_dev/test-labelling.R` 全绿，且其余回归无新增失败；
- `R CMD INSTALL` 的失败信息**仍然只有**那 4 个 undefined exports（不新增、不减少）；
- 所有 R 文件 `parse` 通过。
