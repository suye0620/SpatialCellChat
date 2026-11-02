# Agent Note: `setIdent()` 按"分组轴一致性"重设计（方案 A）

Status: implemented
Governance: v1
Date: 2026-11-02（当日实施）
Decision type: api
Scope: R/utilities.R（`setIdent`、`.sc_same_partition`、`.sc_reindex_group_layers`、`updateClusterLabels`）、R/SpatialCellChat_class.R（删重复定义）、man/、tests_dev/test-setident.R
Owner: maintainer
Impact: high（语义变更：不再静默容忍"换分组"；移除一个参数）
Supersedes: none
Superseded by: none

## 1. 问题回顾

`setIdent` 改的是分组轴，不是普通属性。挂在轴上的五类结构：

| 结构 | 形状 | validator 检查 |
|---|---|---|
| `net$group$prob` / `pval` | group × group × layer（`SparseChatArray`） | ✅ 与 `levels(idents)` 比对 |
| `net$group$count` / `weight` | group × group 矩阵 | ✅ 同上 |
| `net$centr`、`netP$centr` | measure × group × layer | ❌ **完全不查** |
| `netP$group$prob` / `pval` | 同 net | ✅ |

旧实现把**重排、改名、换分组**三件事塞进一个 API：前两件事结果可保留，第三件事在数学上无法迁移（group 级数组必须重新聚合）。判据是"旧组名是否包含于新层级"，既读错字段名，也区分不了"改名"和"换了分组但名字恰好重合"。

## 2. 方案 A：先判"是不是同一个划分"

一句话：**同一划分 → 同步重排/改名；不同划分 → 丢弃分组级结果、保留细胞级结果。**

### 2.1 判据（可严格判定）

两个标注是同一个划分（partition）当且仅当：每个旧标签只对应一个新标签，且每个新标签只来自一个旧标签。

```r
.sc_same_partition <- function(old, new) {
  tab <- table(old, new) > 0
  all(rowSums(tab) == 1L) && all(colSums(tab) == 1L)
}
```

额外前提：`length(levels(old)) == length(levels(new))`，否则 dimnames 无法一一映射。

这个判据同时覆盖：
- 纯置换（顺序变、名字不变）
- 纯改名（双射换标）
- 改名 + 换序的混合（旧实现的两个分支会漏判这种）

明确排除：合并（多旧 → 一新）、拆分（一旧 → 多新）、重新分组。

### 2.2 分支一：同一划分 → 重排 + 改名，一条代码路径

不再区分"重排"和"改名"两个分支，统一：

```r
idx <- match(new_levels, old_levels)   # 新第 k 位取旧的第 idx[k] 位

arr <- arr[idx, idx, ]                        # prob / pval
dimnames(arr) <- list(new, new, dimnames(arr)[[3]])

mat <- mat[idx, idx, drop = FALSE]            # count / weight
dimnames(mat) <- list(new, new)

centr <- centr[, idx, , drop = FALSE]         # measure × group × layer
dimnames(centr)[[2]] <- new
```

两个槽都要走：`net` 与 `netP`。

**不动的东西**：`net$cell$prob`、`netP$cell$prob`（细胞级，与分组无关）、`net$tmp`、`misc`、`dr`、`LR`、`DB`。

为什么统一成一条路径：混合情形（A→X 且换序）下，"先判重排还是改名"的旧逻辑必然二选一，两种都不对；而"先按位置重排、再覆盖成新标签"对三种情形都正确。

### 2.3 分支二：不同划分 → 精确丢弃分组级结果

丢弃：`net$group`（整个 list）、`netP$group`、`net$centr`、`netP$centr`。
保留：`net$cell$prob`（细胞级概率完全可复用）。

随后 `message()` 给出确切的重跑配方：

```
computeAvgCommunProb() → computeCommunProbPathway() → aggregateNet()
→ netAnalysis_computeCentrality(slot.name = "netP")
```

**这是对原 A 方案的一处修正，见 §10 的 D-A1。**

### 2.4 空层级

组数必须相等。没有细胞落在其中的层级，其 dimnames 位置无法从数据推断，按"新旧各自剩余的层级升序配对"。这条写进文档，避免用户以为是按名字匹配。

## 3. API

```r
setIdent(object, ident.use = NULL, levels = NULL)
```

- `ident.use`：取 `meta(object)` 的某一列作为新 idents；
- `levels`：设置层级顺序（可只给部分，`factor()` 语义）；
- **移除 `display.warning`**：它承载的那句"`computeCommunProb` 之后的所有计算都要重跑"在原语义下恒真、在新语义下恒假（同一划分时结果完全有效），留着只会误导；
- 赋值 `object@idents` 时必须补 `names(value) <- rownames(meta)`——`as.factor(meta[[col]])` 会丢 names，validator 要求 `names(idents) == rownames(meta)`；
- 结尾 `validObject()`。

## 4. `updateClusterLabels` 跟着简化

| 旧行为 | 新行为 |
|---|---|
| 改标签后重跑 `computeCommunProbPathway` → `aggregateNet` → `netAnalysis_computeCentrality` | **不再重跑**。同一划分下这些数值完全不变，重跑是纯浪费，还多一个失败点（要求 `net$group$prob` 存在） |
| 允许 `mapvalues` 的 many-to-one（合并标签） | **要求双射**：`new.cluster.name` 长度必须等于旧组数且无重复；否则报错并指出这是换划分，需先 `computeAvgCommunProb()` 重跑 |

其余保留：`old.cluster.name`、`new.order`、`new.cluster.metaname` 写入 meta。

## 5. 删除重复定义

`R/SpatialCellChat_class.R:1326` 的 `setIdent(chat, ident.use = "ident")`。C locale 下该文件先于 `utilities.R` 加载，它自类重写起就被覆盖，且完全不含分组轴逻辑。删掉后 `R/` 下只应存在一份 `setIdent`。

## 6. 测试计划（新增 `tests_dev/test-setident.R`）

在真实 pipeline 对象上（`computeCommunProb` → `computeAvgCommunProb` → `computeCommunProbPathway` → `aggregateNet` → 中心性）：

1. **置换**：`levels` 换序 → 五类结构 dimnames 全部按新序、按组名取到的值不变、`validObject` 通过；
2. **改名**：双射换标 → 数值不变、层名（dim3）不变、`centr` 的 dim2 同步；
3. **混合**：改名 + 换序同时 → 结果等于"先改名再换序"；
4. **合并**（`levels` 少一个）→ `net$group` / `netP$group` / `centr` 全部消失，`net$cell$prob` 仍在，message 含 `computeAvgCommunProb`；
5. **拆分**（一个旧组拆成两个）→ 同 4；
6. **空层级**：删除某组所有细胞后改名，dimnames 仍对齐；
7. **names 约束**：`ident.use` 路径后 `names(idents) == rownames(meta)`；
8. **无网络对象**：只有 idents 变化，不报错；
9. `updateClusterLabels`：双射 → 成功且不触发重跑；非双射 → 报错；
10. `R/` 下只有一份 `setIdent` 定义（`grep` 计数为 1）。

## 7. 文档

- `man/setIdent.Rd` 重写：放"同划分 / 不同划分"两条分支的语义表 + 丢弃清单；
- `man/updateClusterLabels.Rd`：删掉"会重跑"的描述，加上"必须双射"。

## 8. 与其他 B 组任务的关系

本方案只覆盖 `setIdent` / `updateClusterLabels`。剩余两个（`computeAveExpr`、`identifyOverExpressedLigandReceptor`）是单纯字段替换，另行处理。

## 9. 验收标准

- `R/` 下 `setIdent` 定义唯一；
- `net` / `netP` 的 `group$prob` / `pval` / `count` / `weight` / `centr` 在置换与改名后与"重新按新组计算"的 dimnames 完全一致；
- 不同划分时对象仍 `validObject`（分组级被清空、细胞级保留）；
- `tests_dev/test-setident.R` 全绿；回归无新增失败；`R CMD INSTALL` 仍只报 4 个 undefined exports。

## 10. 实施结果（与方案一致）

三条决策全部按建议落地：D-A1 丢弃分组级结果、D-A2 移除 `display.warning`、D-A3 `updateClusterLabels` 停止重跑。落地时补了三个内部函数：

| 函数 | 职责 |
|---|---|
| `.sc_same_partition(old, new)` | 交叉表判定是否双射 |
| `.sc_group_level_map(old, new, old.levels, new.levels)` | 旧层级 → 新层级的完整映射；观测层级按内容配，空层级按剩余顺序配 |
| `.sc_reindex_group_layers` | 单一路径重排五类结构 |
| `.sc_has_group_layers` / `.sc_drop_group_layers` | 判定与清空分组级结果 |

### 实施中踩到的一个陷阱（已写进测试注释）

`net$centr` 与 `net$centr.cell` 在同一个 list 里，而 `$` 会**部分匹配**：当只有 `centr.cell` 时，`object@net$centr` 会静默返回 `centr.cell`（measure × cell × layer），看起来"有值"，实际轴完全不同。此外 `netAnalysis_computeCentrality()` 默认 `do.group = FALSE`，写的是 `centr.cell`；要拿到组级 `centr` 必须显式 `do.group = TRUE`。

因此实现里一律用 `[[ "centr" ]]` 精确取，测试 fixture 也用 `do.group = TRUE` + 精确索引。

### 验证

新增 `tests_dev/test-setident.R` **31 项**：置换（五类结构 + 值不变 + validObject）、改名（层名不变 + names 对齐）、改名+换序混合与分步等价、合并/拆分（分组级清空、细胞级保留、message 点名重跑）、空层级映射（直接测 helper）、无网络对象、`updateClusterLabels` 双射与非双射报错、`R/` 下 `setIdent` 定义唯一。

回归：seurat-bridge 41、projectData 19、computeCommunProb 61、computeAvgCommunProb 162、filterCommunication 22、normalizeData 9、comm-distance 26、class 22、accessors 22、spatial-image 15、pre-inference 15、computeCellDistance 8、cpp-prob-layer 20、helpers 16 全绿；`man/setIdent.Rd`、`man/updateClusterLabels.Rd` `checkRd` 通过；14 个 R 文件 parse 通过。

## 11. 需要你确认的三件事（已于 2026-11-02 批准）

| | 问题 | 我的建议 | 另一面 |
|---|---|---|---|
| **D-A1** | 不同划分时：丢弃分组级结果，还是直接报错 | **丢弃并提示重跑**（保留 `net$cell$prob`） | 报错更"硬"，但用户就没有合法路径改分组了——除非重建对象，而细胞级概率明明可以复用 |
| **D-A2** | 是否移除 `display.warning` | **移除** | 保留则要给它编一个新含义，实际没人用 |
| **D-A3** | `updateClusterLabels` 是否停止重跑那三个函数 | **停止** | 保留重跑 = 每次改名都白算一遍通路聚合与中心性 |
