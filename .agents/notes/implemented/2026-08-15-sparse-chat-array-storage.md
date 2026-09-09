# Agent Note: SparseChatArray 的按层稀疏存储

Status: implemented
Date: 2026-08-15
Decision type: architecture

## Problem

通讯概率的自然形状是 `cell × cell × LR/pathway`。在单细胞规模下构造稠密 3D 数组会产生不可接受的内存占用；全局三元组压缩又会增加按层访问和下游聚合的复杂度。

## Decision

使用轻量 S3 类 `SparseChatArray` 表示三维通讯结果。底层是等维、带名称的 `dgCMatrix` 列表，每个元素代表一个 LR 或 pathway 层。

对外保留 3D 数组式语义：

```r
arr[i, j, k]
dim(arr)
dimnames(arr)
```

当前已建立的基础方法包括 `[`, `[[`, `dim`, `dimnames`, `length`, `names`, `print`, `t`, `marginSums` 和 `as.data.frame`。跨第三维聚合使用 Rcpp `cpp_sum_layers`，避免先把所有层展开成稠密对象。

`net` / `netP` 中的 cell-level 和 group-level `prob` / `pval` 使用该类型；二维 `count`、`weight` 和派生中心性结果继续存放在对应子结构中。

## Alternatives considered

- 原生稠密 3D array：在 `nC × nC × N` 场景下不可接受。
- `spatstat.sparse::sparse3Darray`：引入额外生态依赖，全局索引使单层切片和按需迭代不够直接。
- `DelayedArray` / `HDF5Array`：对当前小中型数据集引入过重的 Bioconductor 依赖。
- 3D COO / `dgTMatrix` 平铺：按层访问和聚合不如独立 `dgCMatrix` 直接。

## Consequences

单层访问为直接列表索引，适合按 LR/pathway 层并行处理；跨层计算必须使用 `marginSums` 或等价的稀疏实现，禁止重新引入 `apply(prob, ..., 3, ...)` 造成稠密膨胀。

`Ops.SparseChatArray`、更多非方阵边界和异常输入测试属于后续扩展，不改变当前底层选择。

## Evidence

- `docs/DATA_STRUCTURE.md`
- `docs/REFACTORING_PLAN.md`
- `R/SpatialCellChat_class.R` 中的稀疏数组校验与访问路径
- `tests_dev/test-SparseChatArray*`
- `src/SpatialChat_Rcpp.cpp`
