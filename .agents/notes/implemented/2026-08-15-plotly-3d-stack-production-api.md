# Agent Note: Plotly 生产级 3D stack

Status: implemented
Governance: v1
Date: 2026-11-10
Decision type: api
Scope: `plotly_spatialStackPlot()` in R/visualization.R, NAMESPACE, man, and P4 visualization/browser tests
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

当前 plotly 的 `3d` 模式主要是把点放在固定 z 平面，属于空间 3D 外壳，不是真正的多 feature/pathway z-stack。多层分别建 trace 时还会产生多个独立 colorbar 和不一致的层间色域。

## Decision

新增 `plotly_spatialStackPlot(object, layers, ...)` 生产 API。每个输入 layer
是一张水平切片，沿 z 轴按列顺序堆叠；x/y 来自单个 `SpatialCellChat`
对象的共享空间坐标，z 只编码 layer 顺序，不表示物理深度。

核心契约：

- 输入为坐标表和 `cell × layer` 数值表；layer 名称必须稳定且唯一。
- 将所有层合并为一张 long table，使用单个 `scatter3d` trace，共享一根 colorbar。
- 颜色范围在所有层上全局计算；存在非负数据时以 0 作为下界，保证层间可比。
- z 轴使用 `tickmode = "array"`，`ticktext` 显示 feature/pathway 名称。
- `z.axis.space`、点大小、颜色方案和 hover 字段是显式参数。
- hover 显示 layer、cell id 和 value；不直接拼接未经转义的用户元数据。
- 默认压低 z 轴 aspect ratio，避免层数增加后空间坐标被拉成长条。
- 3D Plotly 的 `scatter3d` 使用 WebGL；函数不另行选择 `scattergl`，也不静默
  抽样或截断大数据。

坐标从
`object@images$coordinates` 的前两列读取；`layers` 必须是带唯一行名和唯一
非空列名的二维数值表，行名集合必须与对象细胞名完全一致，函数按细胞名重排
层数据。允许单层和部分 `NA`，但每层必须至少包含一个有限值；拒绝 `Inf`、
`-Inf` 和 `NaN`。所有 layer 使用一个全局 color domain；纯零输入用以 0 为
中心的有限窄区间避免退化色域，hover 仍显示原始零值。hover 文本对 layer 名
和 cell id 做 HTML 转义。

## Constraints and invariants

- The object must be a single-dataset `SpatialCellChat`; merged dataset-local
  coordinates are rejected by the shared spatial-coordinate guard.
- Coordinate row names and layer row names are unique, non-empty, and equal as
  sets; layer rows are reordered to coordinate order before plotting.
- Layer columns are unique, non-empty numeric names; every layer retains at
  least one finite value, and no `NaN` or infinite value is accepted.
- Layer order determines z tick order. Every finite layer-cell observation is
  represented once in the main trace; `NA` observations are omitted only at
  their own layer-cell positions.
- Color limits are computed once over all finite layers and are shared by the
  main trace. The API does not attach physical units to z or value.

## Acceptance criteria

- 至少 3 层和 300 个点可构建 Plotly 对象，只有一个数据 trace 和一根可见颜色标尺。
- z 轴刻度与 layer 名称一一对应。
- 所有层共享同一色域，层间不发生独立归一化。
- 空层、全 NA 层、单层和全零数据有明确错误或显示行为。
- 生产 API 已加入 NAMESPACE；help page、最终 schema 输入、hover 转义、大数据
  渲染和浏览器交互选择测试均已完成。
- `plotly_spatialStackPlot()` 已加入 NAMESPACE，帮助页和 P4 回归覆盖合法输入、按细胞名
  对齐、单 trace/单 colorbar、全局色域、z 轴层名、hover 转义以及非法/空层
  边界；3D trace 的浏览器 smoke 覆盖 Plotly 初始化和重布局。

## Current evidence

`tests_dev/test-plotly-spatial-stack.R` 验证 one/many layers、按名称对齐、NA/全零
和拒绝边界。`tests_dev/render_viz_p4.R` 在 simulated、905-cell Visium 和
29,536-cell Slide-seq fixtures 生成 3-layer HTML；Slide-seq 画出 88,608 个
layer-cell 点。`tests_dev/test-p4-browser-cdp.js` 的 15/15 Chrome CDP smoke
核验真实浏览器中的点覆盖、scene、layer ticks、单 colorbar 和 relayout。
`tests_dev/test-p4-shiny-cdp.js` 以真实鼠标拖拽检查 Shiny 两次选择替换及 Reset。
临时 library 的 `R CMD INSTALL` 和显式 namespace export 检查通过。

这份 Note 的接口只支持一个对象的坐标框架。它不定义跨切片配准、物理 z
单位、不同对象的细胞对应或组织层缺失策略；这些需要独立的数据契约。

## Consumer impact

- Producer/input: callers construct a numeric `cell × layer` matrix/data frame from
  expression, scores, or other named layer measurements; the function does not
  compute those measurements.
- Runtime consumer: `plotly_spatialStackPlot()` in `R/visualization.R`, using
  `object@images$coordinates` and caller-provided values.
- Public API: `NAMESPACE` exports the function; `man/plotly_spatialStackPlot.Rd`
  describes input, missing-value, and z-axis semantics.
- Tests/rendering: `tests_dev/test-plotly-spatial-stack.R`,
  `tests_dev/render_viz_p4.R`, `tests_dev/test-p4-browser-cdp.js`, and
  `tests_dev/test-p4-shiny-cdp.js`.
- Status records: `docs/VISUALIZATION_AUDIT.md`,
  `docs/visualization-migration-P4-issues-2026-11-10.md`, and
  `tests_dev/figures/VISUALIZATION_REVIEW.md`.

## Consequences

- All layers are materialized into a Plotly long data trace. Memory scales with
  finite `cells × layers`; callers should choose a deliberate set of layers.
- Missing numeric `NA` values omit only their corresponding layer-cell points;
  a layer with no finite values errors. Invalid coordinate/value alignment errors
  instead of silently changing cell identity.
- Layer z positions are ordinal display coordinates. No cross-object or physical
  tissue depth comparison is implied.
- WebGL avoids a point-count truncation policy, but does not guarantee responsive
  interaction for arbitrary layer counts. Existing three-object tests establish
  the verified 3-layer workload only.

## Alternatives considered

- 每个 layer 建一个独立 trace：会生成多根 colorbar，且默认色域可能各自归一化。
- 把 z 轴设为数值 value：会混淆 layer 位置和观测值，无法表达切片语义。
- 继续把所有 feature 分别输出为二维图：失去层间空间比较和统一交互。

## Validation

- `Rscript tests_dev/test-plotly-spatial-stack.R`: PASS.
- `Rscript tests_dev/test-viz-p4-migration.R`: PASS.
- `Rscript tests_dev/render_viz_p4.R simulated`, `... visium`, and `... slideseq`:
  all P4 widgets including stack HTML written successfully.
- `node tests_dev/test-p4-browser-cdp.js`: 15/15 PASS.
- `node tests_dev/test-p4-shiny-cdp.js` against
  `tests_dev/run-p4-shiny-browser.R`: mouse selection c1 -> c13 replaced rather
  than accumulated; Reset restored all 20 cells.
- All 13 top-level R files parse; source package installs into a temporary
  library and the installed namespace exports `plotly_spatialStackPlot()`.

## Evidence

- `R/visualization.R::plotly_spatialStackPlot()` contains named-coordinate
  alignment, finite-value validation, global color-domain construction, hover
  escaping, and single `scatter3d` construction.
- `tests_dev/test-plotly-spatial-stack.R` passes 18 focused checks for valid
  rendering and invalid/empty input boundaries.
- `tests_dev/figures/P4/{simulated,visium,slideseq}/summary.csv` records HTML
  generation for 20, 905, and 29,536 cells respectively; the Slide-seq stack
  contains 88,608 layer-cell points.
- `tests_dev/figures/P4/browser-summary-2026-11-10.csv` records 15/15 browser
  smoke cases as `PASS`; `shiny-browser-summary-2026-11-10.csv` records the
  simulated real mouse selection/replace/reset case.
