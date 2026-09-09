# Agent Note: Gi/Lee 默认语义与 Plotly 生命周期加固

Status: proposed
Date: 2026-08-15
Decision type: testing

## Problem

空间统计函数的数值实现大体正确，但默认值和交互生命周期会改变用户对结果的解释：Gi 是否二值化、Lee 是否隐藏负关联、highlight 变量命名、plotly trace 数量和 Shiny observer 累积都需要固定契约。

## Proposal

在迁移函数前锁定以下 API 和测试行为：

### Gi

- 明确 `do.binary` 的含义：默认继续二值化，或改为连续中心性；无论选择哪一个，都必须在文档和测试中明确，不能让用户误以为是连续 Gi。
- 修复 `spot_highlight` 语义反转，使用 `!(spot_labels %in% group.highlight)` 表达非高亮状态。
- 中心性跨 pathway 聚合所得的是 incoming/outgoing 总分，文档中明确它不是单 pathway 统计。
- 用 `inherits(x, "ggplot")` 代替脆弱的 `class(x)[[1]]` 判断。

### Lee

- 评估并固定 `cutoff.Lee`：负空间关联是合法输出，默认隐藏所有 `L <= 0` 可能误导解释；推荐默认不截断，或明确保留当前行为的理由。
- 对 `lw.type = "Itself"` 明确说明其退化为 Pearson 相关。
- 数值型 label 输入触发 dot 模式的行为写入文档。

### Plotly / Shiny

- 不通过 `pl$x$data[[N]]` 硬编码 trace 索引设置 colorbar；构建 trace 时显式传入，或按稳定 trace name 定位。
- highlight 为空时不创建空 trace。
- 将 `legendgroup` 机器 key 与 `legendgrouptitle` 显示文本拆开。
- 将 HTML hover 文本统一转义。
- `plotly_spatialLRpairPlot_shiny()` 使用 `reactiveVal` 保存选中 cell 集合，避免在每次 selection 事件中嵌套注册 observer。

## Acceptance criteria

每个默认值和生命周期修复至少有一个边界测试，覆盖：连续/二值 score、负 Lee 值、空 highlight、trace 数量变化、重复 selection 事件和恶意 HTML 元数据。

## Evidence

- `docs/VISUALIZATION_AUDIT.md` 第 3、5、6、8 节
- `tests_dev/audit_visualization_findings.R`
