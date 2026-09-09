# Agent Note: 不使用稠密或重型第三方 3D 通讯存储

Status: rejected
Date: 2026-08-15
Decision type: architecture

## Problem

通讯结果具有 `cell × cell × LR/pathway` 的三维形状。稠密数组和全局三元组结构在大规模空间转录组数据上会制造内存峰值或增加不必要的索引转换。

## Decision

拒绝以下持久化存储方案：

- 原生稠密 3D array；
- 继续以 `spatstat.sparse::sparse3Darray` 作为核心新存储；
- 为当前对象默认引入 `DelayedArray` / `HDF5Array` 等重型外部层。

当前目标是 `SparseChatArray`：命名 `dgCMatrix` 列表 + 3D 索引语义 + 稀疏跨层聚合。

## Why rejected

- 稠密表示的空间复杂度是 `O(nC^2 × N)`，在单细胞规模下不可接受。
- 全局稀疏三元组不利于直接按 LR/pathway 层访问。
- 重型外部依赖会增加安装、版本和部署成本，且不是当前主要瓶颈的必要解法。

## Consequences

新增计算或绘图代码不得通过 `as.array()`、`apply(..., 3, ...)` 等方式无条件物化三维数据。需要跨层求和时使用 `marginSums()` 或专用稀疏实现。
