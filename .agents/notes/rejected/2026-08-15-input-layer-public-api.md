# Agent Note: 不保留 input.layer 公开号名

Status: rejected
Date: 2026-08-15
Decision type: simplification

## Problem

`input.layer` 容易与 Seurat 的内部 layer / slot 机制混淆，也不能清楚表达用户输入的是 raw counts 还是 normalized expression。

## Decision

公开 API 不使用 `input.layer`，不保留同名兼容别名。使用：

```r
input.assay = c("norm", "raw")
```

Seurat 的 `counts` / `data` 和 SCE 的 `counts` / `logcounts` 只作为内部对象读取映射，由 `input.assay` 决定。

## Why rejected

- `layer` 是外部框架的实现术语，raw/norm 才是本包的输入语义。
- 同一参数名在 Seurat v4/v5 读取路径中含义不同。
- 同时保留两个参数会形成长期 API 歧义和文档分叉。

## Consequences

构造器、日志、参数记录和回归测试统一使用 `input.assay`。内部 `layer.name` 可以继续存在，但仅表示读取 Seurat 数据所需的 `counts` 或 `data` 名称，不对外暴露。
