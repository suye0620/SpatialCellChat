# Agent Note: B 组最后两个函数（`computeAveExpr` / `identifyOverExpressedLigandReceptor`）

Status: implemented
Governance: v1
Date: 2026-11-03（当日实施）
Decision type: migration
Scope: R/modeling.R（`computeAveExpr`）、R/utilities.R（`identifyOverExpressedLigandReceptor`）、man/、tests_dev/test-b-group-tail.R
Owner: maintainer
Impact: low（两个函数零调用者，无导出变化）
Supersedes: §2.3 / §2.4 of `2026-10-29-b-group-remaining-migration.md`（口径一致，仅细化）

## 1. `computeAveExpr`（R/modeling.R:2269）

### 1.1 改动

| 现状 | 改为 |
|---|---|
| `slot.name = c("data.signaling", "data")` | `slot.name = c("signaling", "norm")` |
| `slot(object, slot.name)` | `assay(object, slot.name)` |
| `labels <- object@idents` | `labels <- idents(object)` |
| `labels <- object@meta[[group.by]]` | `labels <- meta(object)[[group.by]]` |

空层直接报错并点名 `subsetData()`：

```r
if (is.null(data.use) || sum(dim(data.use)) == 0) {
  hint <- if (identical(slot.name, "signaling")) "; run `subsetData()` first" else ""
  stop("assay$", slot.name, " is empty", hint, call. = FALSE)
}
```

### 1.2 顺带清掉一段死代码

```r
labels <- object@idents
if (!is.factor(labels)) {
  message("Use the joint cell labels from the merged CellChat object")
  labels <- object@idents      # 和上面一模一样
}
```

两个分支赋同一个值，`labels` 若还不是 factor，后面 `factor(labels)` 统一处理。整段删掉，只留 `labels <- idents(object)`。

### 1.3 待定：是否把 `"smooth"` 也列进来（D-B1）

`assay` 层的合法值是 `raw / norm / scale / smooth / signaling`，其中 `signaling` 与 `smooth` 都是 `subsetData()` 产生的基因子集（`smooth` 是 `raw.use = FALSE` 时 `computeCommunProb` 消费的那层）。

旧函数只有 `data.signaling` / `data` 两个选项，迁移到 `signaling` / `norm` 是等价映射。但用户若想对比"平滑后"的表达，`smooth` 是唯一入口。加进去是 `c("signaling", "smooth", "norm")`，一行的事。

### 1.4 待定：旧取值是否保留兼容（D-B2）

按"重构不留过渡层"，**直接删除**，不设 `warning()` 映射。旧脚本传 `"data.signaling"` 会得到 `match.arg` 的标准报错。

## 2. `identifyOverExpressedLigandReceptor`（R/utilities.R:1676）

### 2.1 改动

| 现状 | 改为 |
|---|---|
| `object@var.features[[features.name]]`（读） | `object@misc$.var.features[[features.name]]` |
| `object@var.features[[features.name.LR]] <- ...`（写） | `object@misc$.var.features[[...]] <- ...` |
| `if (is.list(object@var.features))` 检查 | 直接检查 `misc$.var.features[[...]]` 是否为 NULL |

### 2.2 三个实际 bug（与迁移无关，但在这段代码里）

1. **`features` 分支先用后定义**
   ```r
   } else {
     features.use <- features
     rm(features)
     markers.all <- subset(markers.all, subset = features %in% features.use)
   }
   ```
   `markers.all` 此时还不存在；且 `rm(features)` 之后又用 `features`。传 `features` 必崩。
   改法：先无条件取 `markers.all <- misc$.var.features[[...info]]`（NULL 就报错），再按 `features` 过滤。

2. **`1:nrow()` / `1:length()` 在空输入上迭代两次**
   `markers.all` 为 0 行、或 `complex_input` 没有 subunit 列时，`1:0` 会跑两轮。改成 `seq_len()` / `seq_along()`。

3. **`DB$interaction` 缺列时报错难懂**
   `select(interaction_input, ligand, receptor)` 在列缺失时抛 dplyr 的管道错误。前面加一句显式检查。

### 2.3 语义确认（不需要决策，写进文档）

`features` 参数的含义是"只用这些标记基因"，即对 `markers.all` 做行过滤，不是"直接把这些基因当作结果"。复合体亚基仍会被映射成复合体名（`SU1` → `CPLX`）。

## 3. 测试（新增 `tests_dev/test-b-group-tail.R`）

`computeAveExpr`：
- 默认读 `signaling`，行名 = 该层基因；
- `slot.name = "norm"` 读到全量基因；
- 与手工分组中位数一致；
- 空 `signaling` 报错且消息含 `subsetData`；
- `group.by` 指向 meta 列时按该列分组；
- 传旧值 `"data.signaling"` 报 `match.arg` 错误（确认没有兼容层）。

`identifyOverExpressedLigandReceptor`：
- 结果写进 `misc$.var.features$features.LR`；
- 直接 LR 基因保留、复合体亚基映射成复合体名、DB 外基因丢弃；
- `features` 参数过滤生效（这条覆盖 bug 1）；
- `misc$.var.features` 为空时报错点名 `identifyOverExpressedGenes()`；
- `DB$interaction` 缺列时报错可读。

## 4. 提交

一个提交（两个函数 + 测试 + 文档），信息：`refactor: migrate the last two B-group helpers off the deleted slots`。

## 5. 实施结果

两条决策都按建议落地：**D-B1 加入 `"smooth"`**（`c("signaling", "smooth", "norm")`）、**D-B2 直接删除旧取值**，不留兼容层。

`computeAveExpr` 落地时多做了两件事：

1. 空层判定从 `sum(dim(...)) == 0` 改成 `length(dim(...)) != 2L || any(dim(...) == 0L)`——前者在以基因为空、细胞非空时和为 0 会漏判；
2. `group.by` 指向不存在的列时显式报错（原实现会得到 `NULL` 然后在 `factor()` 里抛一个难懂的错误）。

`identifyOverExpressedLigandReceptor` 的三个 bug 全部修掉，另外补了两条前置校验：`misc$.var.features[[...]]` 必须是带 `features` 列的 data.frame；`DB$complex` 为 NULL 时按空表处理（否则 `complex_input[, grepl(...)]` 会抛下标错误）。

### 验证

`tests_dev/test-b-group-tail.R` **18 项**：`computeAveExpr` 的三层读取、与手工分组中位数一致、`group.by` 生效、旧取值报 `match.arg` 错误、空 `signaling` 点名 `subsetData()`、空 `smooth` 点名 `projectData()`；`identifyOverExpressedLigandReceptor` 的写回位置、复合体映射、DB 外基因丢弃、`features` 过滤（覆盖原 bug 1）、缺标记表、DB 缺列、complex 表缺失。

回归：seurat-bridge 41、projectData 19、computeCommunProb 61、computeAvgCommunProb 162、filterCommunication 22、normalizeData 9、comm-distance 26、setident 31、class 22、accessors 22、spatial-image 15、pre-inference 15、computeCellDistance 8、cpp-prob-layer 20、helpers 16 全绿；两个 Rd `checkRd` 通过；14 个 R 文件 parse 通过。

**B 组到此清空（15/15）。**

## 6. 待确认（已批准）

| | 问题 | 建议 |
|---|---|---|
| **D-B1** | `computeAveExpr` 是否接受 `"smooth"` | **加**（`smooth` 是 `raw.use = FALSE` 的唯一入口，不加就没法对平滑表达求均值） |
| **D-B2** | 旧取值 `"data.signaling"` / `"data"` 是否保留兼容 | **不保留**，直接删（符合"不留过渡层"） |
