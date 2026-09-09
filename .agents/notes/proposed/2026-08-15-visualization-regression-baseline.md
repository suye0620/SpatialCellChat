# Agent Note: 最终 schema 可视化回归基线

Status: proposed
Date: 2026-08-15
Decision type: testing

## Problem

目前部分可视化审计已经有 source-level smoke test，但 `spatialFeaturePlot()` 和 `plotly_spatialLRpairPlot()` 在最终 11-slot 对象上的迁移仍未完成。若只测试旧对象或内存 helper，可能把旧 slot 读取误判为兼容成功。

## Proposal

建立一个最小的最终 schema fixture，覆盖：

- raw/norm assay；
- meta、idents、features；
- images coordinates 与可选 raster；
- LR、net、netP 的最小 cell/group 结果；
- 1 个全零 layer、1 个含 NA 的 layer、多个 feature/pathway；
- 空 highlight 和单 cell highlight。

测试分三层：

1. **数据契约**：所有绘图函数只从 accessor 或最终 slot 读取。
2. **静态 ggplot smoke**：返回对象可打印、坐标比例和缺失值策略稳定。
3. **Plotly smoke**：2d、3d、3D stack 的 trace 数量、scene、colorbar、hover 和 selection 状态稳定。

旧 slot 访问应作为明确失败或迁移阻塞被捕获，不能加入隐式兼容路径。

## Acceptance criteria

- fixture 不依赖旧 `@data*` 或旧通信 slot。
- 每个已迁移函数至少有一个最终 schema 成功案例和一个边界案例。
- 测试输出能区分“未迁移错误”和“真实绘图回归”。
- 完整 `devtools::check()` 所需的可选依赖缺失时，targeted smoke test 仍能独立报告结果。

## Evidence

- `docs/VISUALIZATION_AUDIT.md` 第 7、8 节
- `tests_dev/audit_visualization_findings.R`
- `tests_dev/test-SpatialCellChat-spatial-image.R`
