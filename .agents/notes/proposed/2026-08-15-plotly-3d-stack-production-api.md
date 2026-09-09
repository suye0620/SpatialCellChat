# Agent Note: Plotly 生产级 3D stack

Status: proposed
Date: 2026-08-15
Decision type: feature

## Problem

当前 plotly 的 `3d` 模式主要是把点放在固定 z 平面，属于空间 3D 外壳，不是真正的多 feature/pathway z-stack。多层分别建 trace 时还会产生多个独立 colorbar 和不一致的层间色域。

## Proposal

新增生产级 3D stack 绘图 API，语义为：每个 feature/pathway 是一个水平切片，沿 z 轴堆叠，x/y 使用同一套空间坐标。

核心契约：

- 输入为坐标表和 `cell × layer` 数值表；layer 名称必须稳定且唯一。
- 将所有层合并为一张 long table，使用单个 `scatter3d` trace，共享一根 colorbar。
- 颜色范围在所有层上全局计算；存在非负数据时以 0 作为下界，保证层间可比。
- z 轴使用 `tickmode = "array"`，`ticktext` 显示 feature/pathway 名称。
- `z.axis.space`、点大小、颜色方案和 hover 字段是显式参数。
- hover 显示 layer、cell id 和 value；不直接拼接未经转义的用户元数据。
- 默认压低 z 轴 aspect ratio，避免层数增加后空间坐标被拉成长条。
- 大 spot 数据支持 `scattergl` / WebGL 路径或明确的点数上限。

## Acceptance criteria

- 至少 3 层和 300 个点可构建 plotly 对象，只有一根颜色标尺。
- z 轴刻度与 layer 名称一一对应。
- 所有层共享同一色域，层间不发生独立归一化。
- 空层、全 NA 层、单层和全零数据有明确错误或显示行为。
- 生产 API 加入 NAMESPACE 前，先完成最终 schema 输入、hover 转义、大数据渲染和交互选择测试。

## Current evidence

`tests_dev/audit_3dstack_prototype.html` 已验证 3 层 × 300 点的原型。该原型证明构建路径可行，但尚未定义生产 API，也不是当前导出的包函数。

## Alternatives considered

- 每个 layer 建一个独立 trace：会生成多根 colorbar，且默认色域可能各自归一化。
- 把 z 轴设为数值 value：会混淆 layer 位置和观测值，无法表达切片语义。
- 继续把所有 feature 分别输出为二维图：失去层间空间比较和统一交互。
