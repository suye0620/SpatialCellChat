# Agent Note: 构造器 input.assay 与归一化契约

Status: implemented
Date: 2026-08-15
Decision type: architecture

## Problem

构造 `SpatialCellChat` 时，输入矩阵可能是 raw counts，也可能已经归一化。仅根据数值范围或对象类型猜测输入语义，会把 normalized 数据误当 raw，或让 raw 输入绕过必须的归一化流程。

## Decision

`createSpatialCellChat()` 使用显式参数：

```r
input.assay = c("norm", "raw")
normalize = NULL
scale.factor = 10000
do.log = TRUE
```

契约如下：

- `input.assay = "raw"` 默认 `normalize = TRUE`；raw 输入必须是有限、非负数据。
- `input.assay = "norm"` 默认 `normalize = FALSE`；normalized 输入不自动写入 raw 层。
- `raw + normalize = FALSE` 被拒绝。
- `norm + normalize = TRUE` 被拒绝。
- raw normalization 使用按细胞 library size 缩放，并可执行 `log1p`；结果写入 `assay$norm`。
- 归一化后清除依赖旧表达层的 `assay$scale`、`assay$smooth` 和 `assay$signaling` 派生结果。
- Seurat 映射为 `input.assay = "raw"` → `counts`，`"norm"` → `data`；SCE 映射为 `counts` 与 `logcounts`。
- 归一化参数写入 `misc$.param`，操作记录写入 `misc$.log`。
- 内部变量 `layer.name` 仅表示外部 Seurat 的 `counts` / `data` 层，不作为公开 API 参数。

`normalizeData()` 同时支持矩阵输入和 `SpatialCellChat` 输入。对对象输入，它从 `assay$raw` 生成 `assay$norm`，并清理派生层。

## Alternatives considered

- 用 `input.layer` 作为公开参数：名称容易与 Seurat 的内部 `layer` 概念混淆，且不能表达 raw/norm 的生物学语义。
- 根据矩阵是否包含负值自动判断 raw/norm：对中心化或特殊平台数据不可靠。
- raw 输入允许 `normalize = FALSE`：会产生缺少 `assay$norm` 的半初始化对象。

## Consequences

构造器调用方必须明确输入语义；不再保留 `input.layer` 兼容别名。归一化行为、参数和派生层清理都可以被回归测试固定。

## Evidence

- `R/SpatialCellChat_class.R`
- `R/utilities.R`
- `tests_dev/test-normalizeData.R`
