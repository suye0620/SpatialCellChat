# SpatialCellChat 重构计划

> 2026-08-03 · 当前实施基线
> 基于合作者确认的最终数据结构
>
> **权威说明**：当前目标 schema 以 `docs/DATA_STRUCTURE.md` 的“最终数据结构”节为准，共 11 个顶层 slot。`.mode` 和 `.datatype` 归入 `misc`，不是顶层 slot。本文件中的实施顺序和技术约束均以此为前提。

---

## 最终数据结构

```r
SpatialCellChat <- methods::setClass("SpatialCellChat",
  slots = c(
    assay    = "list",       # raw, norm, scale, smooth, signaling
    images   = "list",       # coordinates, coordinate.system, spatial.factors, rasters, fov, caches
    meta     = "data.frame", # 细胞元数据
    idents   = "factor",     # cell group labels for all cells
    features = "data.frame", # 基因元数据
    LR       = "list",       # LR 通讯条目信息
    dr       = "list",       # pca, umap, ...
    net      = "list",       # LR-level: cell$prob/pval/count/weight/centr/field, group$prob/pval/count/weight/centr
    netP     = "list",       # pathway-level: 结构同 net
    DB       = "list",       # interaction, complex, cofactor, geneInfo
    misc     = "list"        # .mode, .datatype, .param, .log, .var.features, .datasets
  )
)
```

共 **11 个 slot**（对比旧 CellChat v2 的 14 个、v3 讨论版的 9 个）。

关键变更：
- `assay$project` → `assay$smooth`
- `colData`/`rowData` → `meta`/`features`（还原旧命名，降低 co-authors 迁移成本）
- `signalingData` 撤销，归入 `LR`
- `netData` 撤销，归入 `net`/`netP`
- `images` 恢复独立 slot，`.distance` 移入 `images`
- `.mode`、`.datatype` 按合作者确认保留在 `misc` 内部，不作为顶层 slot

---

## SparseChatArray

通讯概率矩阵（nC × nC × N）的传统 3D 稠密表示在单细胞规模下不可接受。设计一种基于 `list_of_dgCMatrix` 的轻量 S3 类，对外暴露 3D 数组索引语义。

在最终 schema 中，`SparseChatArray` 不是独立 S4 slot，而是通信结果字段的底层类型，主要用于：

```text
net$cell$prob / net$cell$pval
net$group$prob / net$group$pval
netP$cell$prob / netP$cell$pval
netP$group$prob / netP$group$pval
```

`net` 保存 LR-level 结果，`netP` 保存 pathway-level 结果；两者都采用 `$cell` / `$group` 两层结构。`count`、`weight`、`centr`、`field` 等二维或派生结果仍放在对应的 `net` / `netP` 子结构中。

```r
SparseChatArray <- function(x) {
  structure(x, class = "SparseChatArray")
}
```

当前已实现的方法：`[`、`[[`、`dim`、`dimnames`、`length`、`names`、`print`、`t`、`marginSums`、`as.data.frame`。

当前待扩展：
- `Ops.SparseChatArray` — 逐层算术运算
- 进一步补充非方阵、dimnames、二次转置和异常输入测试

`marginSums` 的跨第三维求和已接入 Rcpp `cpp_sum_layers`；性能测试使用 `microbenchmark` 比较 C++ 实现和 R `Reduce(+)`，具体结果见 `tests_dev/test-SparseChatArray.md`。

### Rcpp 跨层求和：逐列稠密累加器

输入 N 个 n×n dgCMatrix，输出一个 n×n dgCMatrix。

```
对每列 j:
  1. val[n] = 0, mark[n] = false
  2. 遍历所有层在该列的非零行:
     val[row] += x, mark[row] = true
  3. sort(rows), 写入结果，清理 mark
```

- O(总非零元)，无中间膨胀
- 每列 val[n] (8n bytes) + mark[n] (n bytes) 暂存，n=50000 时 ~450KB
- 可加 OpenMP 按列并行

---

## 并行计算与 CLI 提示

### setEnvironment — 入口初始化

```r
setEnvironment <- function(workers = 4,
                           future.globals.maxSize = 10000 * 1024^2) {
  if (workers == 1) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = workers, gc = TRUE)
  }
  options(future.globals.maxSize = future.globals.maxSize)
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  
  cli::cli_h1("SpatialCellChat")
  if (workers == 1) {
    cli::cli_alert_info("Sequential mode")
  } else {
    cli::cli_alert_success("Parallel with {workers} workers")
  }
}
```

### my_future_lapply — 无感复用

```r
my_future_lapply <- function(X, FUN, ...,
                             future.seed = TRUE,
                             .progress = TRUE) {
  n <- length(X)
  if (.progress && n > 1 && future::nbrOfWorkers() > 1) {
    p <- progressr::progressor(along = X)
    f <- function(x) { res <- FUN(x); p(); res }
  } else {
    f <- FUN
  }
  future.apply::future_lapply(X, f, ..., future.seed = future.seed)
}
```

`progressr::handlers("cli")` 走的条件（condition）机制——workers 发信号 → future 带回主进程 → cli 渲染进度条。与 `future::multisession` 完全兼容，已验证通过。

### 设计原则

1. **`setEnvironment` 做一次完整初始化**（策略、进度条样式、全局选项），运行时无副作用
2. **`my_future_lapply` 不重复打印策略信息、不重复设 handlers**
3. **cli 输出统一** — 进度条用 `progressr + handler_cli`，状态消息用 `cli::cli_alert_*`

---

## Phase 0/1/2 已冻结并实现

### Phase 0：数据不变量

- `idents` 顶层 slot 始终为覆盖全部细胞的 factor，保持 `object@idents` 的直接访问兼容；样本级信息归档在 `misc$.datasets`。
- `misc$.mode` 只取 `single` 或 `merged`；merged 模式的样本级结果归档于 `misc$.datasets`，顶层 `net` / `netP` 只表示 joint 结果。
- 空间坐标、空间因子、可选的 Visium raster 图像和空间缓存统一存于 `images`；距离与接触矩阵位于 `images$.distance`。
- `images$rasters` 保存命名 raster 记录；每条记录包含 `image`、`scale.factors`、`spot.coordinates`、可选 `spot.radius`/`source`，并规范化 `coordinate.system` 与 `transform`。
- `assay$norm` 是 genes × cells 的主表达矩阵，`meta`、`features`、坐标和所有通讯矩阵必须与其行列名称对齐。
- 任意对象写入后都应通过 `validObject()` 或 `validateSpatialCellChat()` 检查；输入数据改变时，依赖的通讯或空间缓存必须清除。

### Phase 1：SpatialCellChat S4 基础设施

已完成：

- 11-slot `SpatialCellChat` 类定义。
- `setValidity("SpatialCellChat", ...)` 结构校验。
- `validateSpatialCellChat()` 显式校验入口。
- 面向最终 schema 的 `show()` 方法。
- `createSpatialCellChat()` 的 matrix / Seurat / SingleCellExperiment 输入分流，以及 RNA / spatial 对象初始化。
- 开发测试：`tests_dev/test-SpatialCellChat-class.R`。

Phase 2 已建立统一访问器；下一步进入 Phase 3，先实现 subset、merge、updateObject，保证对象变更时所有 slot 同步。`modeling.R`、`analysis.R`、`spatial.R`、`visualization.R` 的业务逻辑迁移仍待后续阶段完成。

## Phase 2：访问器已完成

已建立并测试以下统一访问器：

```r
assay(object, layer = NULL)
meta(object, columns = NULL)
idents(object)
images(object, key = NULL)
LR(object, key = NULL)
DB(object, key = NULL)
communication(object, slot.name, resolution, measure)
params(object, key = NULL)
```

对应的替换写入接口使用同名 `<-` 形式。所有写入操作都会调用 `validObject()`；表达矩阵、元数据、空间信息、通信矩阵和参数均在 accessor 层保留结构校验。`object@idents` 始终是直接 factor，不存在 `joint` 或 `samples` 子字段。

开发测试：`tests_dev/test-SpatialCellChat-accessors.R`。

## Visium 图像与 spatialDimPlot

已参考 Seurat 4.4.0 的 `Read10X_Image()`、`Load10X_Spatial()`、`VisiumV1` 和 `GetTissueCoordinates()` 契约实现：

- `readSpatialImage()` 从 Visium `spatial/` 目录读取 raster、scale factors 和 tissue positions，返回规范化 raster 记录。
- `createSpatialCellChat()` 接受 VisiumV1、raster list；Seurat spatial 对象默认导入第一张图像到 `images$rasters$main`。
- `images$rasters` 使用命名普通 list 保存可选图像记录；分析坐标始终独立保存在 `images$coordinates`。
- `spatialDimPlot()` 在存在可对齐 raster 时默认叠加背景；`image = FALSE` 保持坐标散点图；`image.alpha` 控制透明度，`crop` 控制绘图区范围。
- 开发测试：`tests_dev/test-SpatialCellChat-spatial-image.R`。

## 重构策略

### 基本原则

1. **按模块逐个重构，不一把抓**
2. **叶节点优先** — 先改没有依赖的底层（SparseChatArray 方法），再改上层（计算函数→下游分析→绘图）
3. **行为不变性** — 每重写一个函数，对照原输出确保结果一致
4. **先测试再整合** — 重写完的函数在 `tests_dev/` 下用 testthat 验证后，才合并到 `R/`
5. **roxygen2 完整手写** — 不跨文件引用 `@template`，每个函数独立完整，降低 co-authors 维护负担

### 实施顺序

```
Phase 1: 数据结构基础设施
  1.1 完善 SparseChatArray 方法（marginSums、as.data.frame）
  1.2 写 Rcpp 跨层求和函数 cpp_sum_layers
  1.3 写并测试 my_future_lapply / my_future_sapply / setEnvironment

Phase 2: 表达式层
  2.1 更新 createSpatialCellChat（新 slot 路径）
  2.2 更新 preProcessing（assay$project → assay$smooth）

Phase 3: 计算管线
  3.1 重构 computeExpr_* 辅助函数
  3.2 重构 computeCommunProb（去掉 my_as_sparse3Darray，去掉 tmp）
  3.3 重构 computeCommunProbPathway
  3.4 重构 aggregateNet（将 count/weight/centr/field 写入对应的 net/netP 子结构）

Phase 4: 下游分析
  4.1 重构 netAnalysis_computeCentrality
  4.2 重构 rankNet、identifyCommunicationPatterns 等

Phase 5: 空间与可视化
  5.1 重构 spatial.R（@images 路径还原）
  5.2 重构 visualization.R（新索引适配）

Phase 6: 收尾
  6.1 更新 mergeSpatialCellChat / subsetSpatialCellChat
  6.2 清理废弃函数（my_as_sparse3Darray、spatstat.sparse 引用）
  6.3 全面 devtools::check()
```

### 测试策略

```
tests_dev/
├── test-SparseChatArray.R
├── test-utils-parallel.R
├── test-computeCommunProb.R
└── ...
```

每个测试脚本独立 `source()`，不注册到包。模块稳定后移入 `tests/testthat/`。

### 撤销路径

旧 class 定义的关键代码在 git 历史中可回溯（commit `9804d1d`、`1c72d15`）。不在文件中保留废弃注释块，git 本身即安全网。

---

## 参考

- 数据结构设计详细说明：[DATA_STRUCTURE.md](DATA_STRUCTURE.md)
- 旧版重构计划（2026-06-25，v3 讨论版）：[refactoring-plan.html](refactoring-plan.html)
