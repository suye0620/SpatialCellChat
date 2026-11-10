# 可视化迁移问题记录：P4 交互图

日期：2026-11-10

## P4-001：交互 LR 表达图读取旧 slot 且会稠密化全表达矩阵

`plotly_spatialLRpairPlot()` 读取已删除的 `object@data`、`net$prob`、`net$pval` 和 `net$prob.cell`，还会把完整表达矩阵转为稠密矩阵；全量 Slide-seq 风险很高。应读取 `assay$norm`，只提取所需基因，并复用静态表达图的逐 LR 稀疏筛选；group 模式保留本次阈值过滤，cell 模式只检查正的通信值。坐标、表达和 metadata 按细胞名对齐，不改变坐标转换和颜色映射。

Shiny 封装和 Plotly 输出需要分别实际验证；仅有 R 对象或静态 HTML 文件不代表浏览器交互已经通过。

## P4-002：Shiny 框选重复注册监听器

`plotly_spatialLRpairPlot_shiny()` 原来在每次 `plotly_selecting` 事件内新建一个 `input$update` 监听器；框选次数越多，Update 触发的监听器越多，旧选择也可能再次参与处理。另有一个未限定命名空间的 `observe()`，在没有附加 shiny 包的调用环境中可能找不到函数。调整为一个固定的 Update 监听器，在点击时读取当前选择，并显式使用 `shiny::observe()`；保持原有选择、反转和清除的用途。需要通过 Shiny 会话及浏览器交互分别验证。

## P4-003：交互维度图的坐标与标签可能错配

`plotly_spatialDimPlot()` 使用坐标行顺序，却按位置读取 `idents` 或 metadata 标签；当坐标行被重排而细胞名不变时，图不会报错，但细胞会显示为其他群。按坐标细胞名对齐标签，并拒绝缺失、重复或不匹配的名字，不改变绘图坐标与群配色。影响 2D、3D 和使用该图的 Shiny 页面。

## 当前证据与限制

三个对象各自生成了 2D/3D 维度图、3 层 spatial stack、单基因表达图与单指标 Gi 3D HTML；Slide-seq 每项数据 trace 均覆盖全部 29,536 个细胞，stack 覆盖 88,608 个 layer-cell 点。`tests_dev/test-p4-browser-cdp.js` 使用本地 Chrome CDP 对 2026-11-10 生成的 15 个页面执行了 Plotly 初始化、点数/细胞覆盖、3D scene、stack 单主 trace/单 colorbar/三层刻度和一次 `Plotly.relayout` 检查，汇总见 `tests_dev/figures/P4/browser-summary-2026-11-10.csv`，15/15 为 PASS。

`tests_dev/test-p4-shiny-cdp.js` 另行启动 simulated Shiny app，通过 Chrome CDP 的真实鼠标拖拽和按钮点击验证两次选择替换及 Reset：首次为 `c1`，第二次为 `c13`，两次不累加，Reset 恢复 20 个细胞；汇总见 `tests_dev/figures/P4/shiny-browser-summary-2026-11-10.csv`。这证明了 Shiny 选择器的浏览器交互契约，但没有替代在 Visium/Slide-seq 大对象上的性能验收。

## P4-004：三维堆叠生产入口与输入契约

生产入口为 `plotly_spatialStackPlot()`。它从单对象 `images$coordinates` 读取共享 x/y，接收带唯一细胞行名和层列名的 `cell × layer` 数值表，按细胞名对齐，沿 z 轴按层顺序堆叠；z 轴是层序的视觉编码，不解释为物理深度。所有层使用单个 `scatter3d` 主 trace 和全局色域/单 colorbar，不抽样或截断大对象；NaN、Inf、重复/缺失细胞名、重复/空层名和全 NA 层直接拒绝。三套 fixture 的 3 层 HTML 和 CDP 验收已通过。该入口不声称完成多切片组织配准；真实多切片配准仍属于不同的输入契约。
