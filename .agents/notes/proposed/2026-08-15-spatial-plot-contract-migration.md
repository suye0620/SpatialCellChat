# Agent Note: spatialXXXPlot 统一数据契约与迁移

Status: proposed
Date: 2026-08-15
Decision type: architecture

## Problem

`R/visualization.R` 中的 `spatialXXXPlot`、feature/LR 绘图和 plotly 路径重复实现坐标变换、颜色尺度、主题和数据整理。部分函数仍直接读取旧的 `@data`、`@data.signaling` 或旧通信结构，无法在最终 11-slot schema 上稳定运行。

## Proposal

按以下顺序重构：

1. 先通过 `assay()`、`communication()`、`LR()` 和 `images()` 等访问层迁移数据来源。
2. 所有空间绘图先生成统一的 long-format `plot_data`，至少包含 cell id、分析坐标、feature/layer、value、group 和 highlight 状态。
3. 把分析坐标到绘图坐标的轴方向变换集中到 helper；raster pixel 坐标通过显式 transform 转换，不在每个函数中重复 x/y 交换和 `scale_y_reverse()`。
4. 收敛 feature/topic/scoring 的主题、颜色尺度和缺失值策略；避免重复内联 theme 和重复构造 colorbar。
5. ggplot2 使用 `after_stat()`、`aes()` / `.data` 和 `inherits()` 等稳定语法，减少 `aes_string` 与 ggplot 内部 panel API 的依赖。
6. plotly 与 ggplot2 共享同一份标准化数据和颜色域，交互层只负责 trace、hover、scene 和 selection。

## Acceptance criteria

- `spatialFeaturePlot()` 与 `plotly_spatialLRpairPlot()` 不再访问旧 slot。
- 普通空间图、raster 背景图和 plotly 图使用同一坐标语义。
- 单 feature、多 feature、全零值、NA、highlight 为空和大 spot 数都有稳定行为。
- 旧 slot 迁移完成后删除兼容分支，而不是继续叠加 fallback。
- 迁移后的函数有针对最终 schema 的回归测试。

## Alternatives considered

- 只在每个旧函数内部添加 `if (new schema)` 分支：会复制数据访问和坐标语义，延长旧结构生命周期。
- 先统一视觉主题、最后迁移数据结构：视觉重构会反复修改，且无法验证真正的 schema 契约。

## Risks

`visualization.R` 规模较大；应按数据访问、坐标预处理、主题/色图、单函数迁移分批完成，每批保留独立 smoke test。

## Evidence

- `docs/VISUALIZATION_AUDIT.md` 第 4、6、8 节
- `docs/REFACTORING_PLAN.md` Phase 5
