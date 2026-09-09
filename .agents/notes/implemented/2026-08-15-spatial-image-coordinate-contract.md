# Agent Note: 空间图像、坐标和距离缓存契约

Status: implemented
Date: 2026-08-15
Decision type: architecture

## Problem

空间分析坐标、图像像素坐标、空间比例因子和距离缓存曾分散在多个旧路径中。绘图函数在切换 raster 或距离阈值后可能复用错误缓存，坐标轴方向也容易依赖隐含的 x/y 交换。

## Decision

统一使用 `images` 保存空间相关数据：

- `images$coordinates` 是分析坐标，独立于 raster 像素坐标。
- `images$coordinate.system` 描述坐标单位、校准状态和轴方向。
- `images$spatial.factors` 保存坐标到工作单位的转换参数。
- `images$rasters` 是命名 raster 记录；每条记录包含图像、scale factors、spot coordinates，并可包含 spot radius、来源和坐标变换。
- `images$.distance` 保存空间距离与接触邻接缓存，并携带计算参数。

`spatialDimPlot()` 在存在可对齐 raster 时叠加背景；`image = FALSE` 保持纯坐标散点图；`image.alpha` 和 `crop` 控制显示。

距离缓存只有在当前 interaction/contact range、ratio、tol 和坐标与缓存参数一致时才复用；参数不一致或旧缓存缺少必要元数据时重新计算。`computeCellDistance()` 使用精确范围查询构造稀疏距离矩阵，避免先物化 `n × n` 稠密矩阵。

## Alternatives considered

- 把分析坐标直接存入 raster 记录：会把分析坐标和图像像素坐标耦合。
- 继续从旧的 `images$result.computeCellDistance` 读取：无法可靠判断缓存参数是否匹配。
- 每次绘图无条件重算距离：避免缓存错误但浪费大型空间数据的计算成本。

## Consequences

新增图像来源必须通过统一的 raster 记录和坐标系统契约进入对象。所有绘图函数应从 accessor 或统一预处理层获取坐标，不应各自实现隐含轴变换。

## Evidence

- `docs/REFACTORING_PLAN.md`
- `docs/VISUALIZATION_AUDIT.md`
- `R/spatial.R`
- `tests_dev/test-SpatialCellChat-spatial-image.R`
- `tests_dev/test-computeCellDistance.R`
