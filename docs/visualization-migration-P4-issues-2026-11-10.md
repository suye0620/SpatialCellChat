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

三个对象各自生成了 2D/3D 维度图、单基因表达图与单指标 Gi 3D HTML；Slide-seq 每项数据 trace 均覆盖全部 29,536 个细胞。Chrome 静态截图已抽查模拟 2D/3D、Visium 表达图和 Gi 3D；未完成各图鼠标交互及 Shiny 浏览器端验收。Shiny 的连续两次框选、更新、反转和重置通过 `shiny::testServer()` 验证，不代表浏览器验收通过。所有 HTML 结果保留 `BROWSER-NOT-TESTED` 状态。

## P4-004：三维堆叠仅存在于审计原型

迁移计划要求 3D stack 生产化，但现有 `spatial3Dstack()` 只定义在 `tests_dev/audit_visualization_findings.R`，不是包中现有公开入口；计划未给出多层输入、坐标配准、细胞对应、z 轴单位与缺失切片的完整契约。把该原型复制到包内并用单样本自造多层，不能证明真实多切片数据的语义正确。此项尚未完成，需先确定真实输入及接口契约，再按真实数据验收。
