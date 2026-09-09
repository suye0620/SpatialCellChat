# SpatialCellChat 数据结构重构计划（历史讨论稿）

> 2026-06-25 · 早期计划文档
>
> **状态：已过时，不作为当前实现依据。** 本文档记录的是早期 v3 讨论版，使用 `expression`、`communication`、`spatial` 等命名；当前开发应以 `DATA_STRUCTURE.md` 的最终 11-slot schema 和 `REFACTORING_PLAN.md` 的实施顺序为准。当前最终设计中，`.mode`、`.datatype` 位于 `misc` 内部，通信结果位于 `net` / `netP`，`SparseChatArray` 作为其中 `prob` / `pval` 的底层类型。

---

## Context

当前 SpatialCellChat 的 S4 class 设计继承自 CellChat v1/v2，从组级计算演进到单细胞分辨率后，原有数据结构已不适合新的计算规模。核心问题：`sparse3Darray` 格式对 GPU 不友好、`net` 是裸 list 无类型约束、`tmp` 字段泛滥、多份数据冗余。在 Phase 1 重构之前，必须先打好数据结构的基础。

---

## SparseChatArray — 新的 3D 稀疏容器

### 定义（~30 行，放在 `R/SpatialCellChat_class.R` 顶部）

```r
#' @title SparseChatArray
#' @description 一个命名稀疏矩阵集合，对外表现得像三维数组，支持 [,,k] 索引。
#' 底层是 list_of_dgCMatrix，可直接供 cuSPARSE 等 GPU 库使用。
SparseChatArray <- function(x, ..., dimnames = NULL) {
  if (!is.list(x)) stop("x must be a list")
  if (length(x) == 0) stop("x must be non-empty")
  if (!all(vapply(x, inherits, logical(1), "dgCMatrix")))
    stop("all elements must be dgCMatrix")
  dims <- unique(lapply(x, dim))
  if (length(dims) > 1) stop("all matrices must have the same dimensions")
  if (!is.null(dimnames)) names(x) <- dimnames[[3]]
  for (i in seq_along(x)) {
    rownames(x[[i]]) <- dimnames[[1]]
    colnames(x[[i]]) <- dimnames[[2]]
  }
  structure(x, class = "SparseChatArray")
}

`[.SparseChatArray` <- function(x, i, j, k, drop = TRUE) { ... }
`dim.SparseChatArray` <- function(x) c(nrow(x[[1]]), ncol(x[[1]]), length(x))
`dimnames.SparseChatArray` <- function(x) list(rownames(x[[1]]), colnames(x[[1]]), names(x))
length.SparseChatArray <- function(x) length(unclass(x))
names.SparseChatArray <- function(x) names(unclass(x))
print.SparseChatArray <- function(x, ...) { ... }
```

存储层面 = list_of_dgCMatrix（GPU 友好），接口层面 = 3D 数组（`x[,,k]` 兼容现有下游代码）。

**关于重命名说明** — "Sparse" 表示稀疏存储，"Chat" 表明第三维度（LR 通路层）与通讯结果绑定，"Array" 是对外表现为三维数组。与 SpatialCellChat 命名一脉相承。

---

## 辅助函数：操作记录

在 `R/utilities.R` 中新增一个轻量级辅助函数，供所有主要函数调用：

```r
.log_operation <- function(object, funcname, params = list()) {
  entry <- list(
    function = funcname,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    params = params,
    version = as.character(packageVersion("SpatialCellChat"))
  )
  object@log <- c(object@log, list(entry))
  return(object)
}
```

使用方式（在每个主要导出函数的 return 之前）：

```r
computeCommunProb <- function(object, ...) {
  # ... 计算逻辑 ...
  object <- .log_operation(object, "computeCommunProb",
    params = list(interaction.range = interaction.range, ...))
  return(object)
}
```

设计原则：

- 只记"操作了什么"，不记"数据本身"（数据在 slot 里）
- 参数只记关键参数，不记巨大矩阵
- 纯追加，不删除。用户可以通过 `object@log` 查看全流程
- 类比 Seurat 的 `seu@commands`

---

## 新 SpatialCellChat S4 class

### Slot 设计

按照之前的计划，我们要开始对spatial CellChat动手了，它的核心数据结构是什么

### `@log` 结构示例

```r
object@log[[1]]
# $function: "preProcessing"
# $time    : "2026-06-25 14:32:00"
# $params  : list(slot.name = "data.signaling", ...)
# $version : "0.1.0"

object@log[[2]]
# $function: "computeCommunProb"
# $time    : "2026-06-25 14:35:12"
# $params  : list(interaction.range = 250, ...)
# $version : "0.1.0"
```

### 与旧结构对应关系

| 旧 slot                           | 新 slot                                     | 备注                           |
| -------------------------------- | ------------------------------------------ | ---------------------------- |
| `@data`                          | `@expression$norm`                         | 标准化后数据                       |
| `@data.raw`                      | `@expression$raw`                          | 可选，用完即删                      |
| `@data.signaling`                | `@expression$signaling`                    |                              |
| `@data.scale`                    | `@expression$scale`                        |                              |
| `@data.project`                  | `@expression$project`                      |                              |
| `@net$prob.cell` (sparse3Darray) | `@communication$cell` (SparseChatArray)    | **格式变更，核心改动**                |
| `@net$prob` (array)              | `@communication$group` (SparseChatArray)   | 3D array → 稀疏                |
| `@net$pval` (array)              | `@communication$pval` (SparseChatArray)    | 同上                           |
| `@netP$prob` (array)             | `@communication$pathway` (SparseChatArray) | 同上                           |
| `@net$count/sum/weight`          | `@communication$network$count/weight`      | sum/weight 合并为 weight        |
| `@net$centr`                     | `@communication$network$centr`             |                              |
| `@net$LR.sig`                    | `@communication$LRsig`                     |                              |
| `@options$misc`                  | `@log` + `@options`                        | misc 拆分，元数据归 log，实参归 options |
| `@net$tmp`                       | **删除**                                     | 各字段各归其位                      |
| `@net$tmp$prob.cell` (list)      | **删除**                                     | 去掉双存储                        |
| `@net$tmp$Lavg/Ravg/...`         | 直接参数传递或 `@options$tmp.LR`                  | 计算中间结果不存入 object             |
| `@images`                        | `@spatial`                                 | 重命名，结构明确                     |
| `@net` list                      | **删除**                                     | 拆入 communication             |
| `@netP` list                     | **删除**                                     | 拆入 communication             |

---

## 波及范围

### 必须改动的文件

| 文件                          | 改动范围                                                                                                                      |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `R/SpatialCellChat_class.R` | class 定义重写，SparseChatArray 定义，create/merge/update/subset/lift 全部改                                                         |
| `R/modeling.R`              | computeCommunProb：去掉 `my_as_sparse3Darray`，直接 `SparseChatArray(prob.cell_)`；aggregateNet：写新格式；computeCommunProbPathway：同上 |
| `R/analysis.R`              | netAnalysis_computeCentrality：改 `centr` 索引；identifyCommunicationPatterns：改 `netP` 索引                                      |
| `R/visualization.R`         | 所有 `prob.cell[,,k]` → 通过 SparseChatArray 的 `[` 方法兼容                                                                       |
| `R/utilities.R`             | `my_as_sparse3Darray` 保留（其他地方还在用），`preProcessing` 改 slot 路径                                                               |
| `R/spatial.R`               | `object@images` → `object@spatial`                                                                                        |

### 不改，保持不动

- `R/data.R` — 内置数据库，无关
- `R/database.R` — 数据库操作，无关
- `src/` — C++ 代码，无关
- `R/RcppExports.R` — 自动生成，无关

### 兼容模式：SparseChatArray `[` 方法

所有下游代码写 `prob.cell[,,k]` 的，换了 SparseChatArray 后照常工作。唯一需要改的是用 `slot()`、`object@net`、`object@images` 直接访问的结构层级路径——文本替换可在 IDE 里批量完成。

---

## 实施顺序

```
 1. 写 SparseChatArray 类 + 方法    — R/SpatialCellChat_class.R
 2. 写新 SpatialCellChat class 定义   — R/SpatialCellChat_class.R
 3. 更新 createSpatialCellChat        — 写新 slot
 4. 更新 modeling.R 核心管线          — 输入用新格式
 5. 更新 utilities.R（preProcessing） — 改 slot 路径
 6. 更新 aggregateNet                 — 写新 communication 结构
 7. 更新 analysis.R                   — 改索引
 8. 更新 spatial.R                    — @images → @spatial
 9. 更新 visualization.R              — 适配新索引
10. 更新 merge/update/subset/lift    — 适配新结构
```

每个步骤后 `devtools::check()` — 不改动函数签名。

---

## 撤销路径

保留旧 class 定义注释在文件的 `# @deprecated` 区块中。如果新结构暴雷，可以恢复旧定义 + 写一个 `as.spatialcellchat.v2()` 转换器，无需 git reset。
