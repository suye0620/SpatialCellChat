# SpatialCellChat CLI 提示设计规范

## 已实现范围

当前约定覆盖 11-slot 推断链的九个函数：`computeCommunProb`、`filterProbability`、`computeAvgCommunProb`、`filterCommunication`、`computeCommunProbPathway`、`aggregateNet`、`netAnalysis_computeCentrality`、`subsetCommunication`、`relabelSpatialCellChat`，以及共享 `.cli()`、`my_future_lapply()` / `my_future_sapply()`。

这不是全包输出迁移完成声明。旧 Visium 分支及其他未迁移函数仍可能使用 `cat()`、`message()`、`cli.symbol()` 或 pbapply；其调用尚存，因此不删除 `cli.symbol()` 导出和 pbapply 依赖。

## 用户控制

全局选项 `SpatialCellChat.verbose` 默认是 `1L`：

| 级别 | 普通消息 | 进度权限 |
|---|---|---|
| 0 | 不显示 | 禁用 |
| 1 | 入口、关键输入、结果摘要 | 允许，仍服从 progressr |
| 2 | 额外显示分支、阈值、执行模式 | 同级别 1 |
| 3 | 另含显式 debug 消息 | 同级别 1 |

`computeCommunProb()` 和 `computeAvgCommunProb()` 的 `verbose` 保持 **logical**，不是整数分级参数。`FALSE` 禁用本次调用的普通消息和进度；`TRUE` 服从全局选项，不能覆盖全局静默。没有新增 `verbose=NULL` 或逐函数整数级别 API。没有局部参数的链路函数使用全局选项。

```r
# 一次性配置 future 策略、全局详细度和 CLI handler
setEnvironment(workers = 2L, verbose = 1L)

# 非交互 Rscript 默认不强制显示进度，需用户显式开启
options(progressr.enable = TRUE)

# 整条迁移链静默；真实 warning()/stop() 仍然可见
options(SpatialCellChat.verbose = 0L)

# 详细输出下，仅本次核心计算静默
options(SpatialCellChat.verbose = 2L)
chat <- computeCommunProb(chat, verbose = FALSE)
```

`setEnvironment()` 是显式配置入口，会设置 future plan、全局详细度和 progressr handlers。计算函数不更改这些全局配置。用户可以单独配置 future/progressr，不要求先调用 `setEnvironment()`；`options(progressr.enable = FALSE)` 始终关闭本链路的进度。

## 消息与 condition

- 普通状态统一调用 `.cli()`，在主进程发出；不要再拼接 `cli.symbol()` 或时间戳，也不要打印整张表或逐 LR 状态行。
- `header` / `subheader`、`info`、`success` 默认要求级别 1；`text` 要求 2；`debug` 要求 3。`.verbose` 仍是显式阈值覆盖，不改变现有语义。
- `.cli(.type="warning"/"danger")` 是默认阈值 0 的显示样式，不等价于 R warning/error condition。真实问题继续使用 `warning()` / `stop()`，不能为了统一外观而吞掉或替换 condition。
- 普通消息通过 cli 的 message/stderr 流输出。仅捕获 stdout 不能验证静默。
- Worker 只发出进度 condition 和真实 warning/error，不打印普通状态；不要把调试日志伪装成 warning。
- 局部 `verbose` 必须守住入口、参数摘要、嵌套调用和完成提示，不只包住主循环。

## 进度生命周期

共享 `.sc_progress_enabled(verbose=TRUE)` 合并局部 logical 权限与全局详细度；progressor 还检查 `progressr.enable` 和工作项数量。只有多于一个工作项时报告进度，不为零项或单项制造进度条。

`my_future_lapply()` / `my_future_sapply()` 在 sequential 和 multisession 下均使用作用域内的 progressr progressor。每项在 `FUN(...)` 返回后才更新；保留 `future.seed`、future 参数、结果命名及简化语义。封装不设置用户 handler，也不强制开启批处理进度。

串行计算保留原调度：使用 `with_progress()` 包住原循环；先计算/存储结果，再更新进度。`filterCommunication()` 仍为串行，不为显示进度而切换到 future。同步 C++ 内核只能在 LR 调用返回后更新，不能显示内部 permutation 的实时百分比。

| 函数/阶段 | 输出重点 | 进度单位 |
|---|---|---|
| `computeCommunProb` | 细胞/LR 规模；详细模式含距离、缓存、LR 分类、线程和落盘模式 | 已计算并写入内存或磁盘的 LR |
| `filterProbability` | LR 规模、抽样和分位数；详细模式含层来源 | 已过滤 LR |
| `computeAvgCommunProb` | 细胞群规模和聚合方法；详细模式含门限、置换配置 | 观测聚合、置换各自的有效 LR |
| `filterCommunication` | 过滤条件及排除计数 | 已处理的细胞级 LR，含被排除项 |
| `netAnalysis_computeCentrality` | 节点、信号层和指标数量 | 已完成的信号层 |
| `computeCommunProbPathway` / `aggregateNet` | 输入及输出网络规模 | 无额外进度条 |
| `relabelSpatialCellChat` | 详细模式显示编排边界；子函数负责自己的输出 | 无外层进度条 |
| `subsetCommunication` | 保持 getter 风格；保留必要提示和真实 warning | 无进度条 |

`computeCellDistance()` 与 `createCellCellContactMatrixFrom_dspatial()` 是静默计算助手；由 `computeCommunProb()` 按自己的详细度报告工作距离单位、阈值及缓存条目，避免嵌套输出穿透局部静默。

## 不变量与边界

- 详细度不改变概率、p 值、随机数序列、过滤阈值、调度、对象 schema 或落盘所有权。
- 内存安全检查不依赖 `verbose`。串行和并行计算均执行现有输出内存保护；并行模式在 worker 结果汇总后检查完整结果。
- `computeCommunProb(parallel=TRUE)` 的 LR 结果现在正确组装进 `net$cell$prob`；设置 `spill.dir` 时，调用方在 worker 计算完成后写入每个 LR 的 RDS，`filterProbability()` 可以正常消费并清理。并行落盘是“计算后写盘”，不是 worker 执行期间的流式写盘。

## 验证

`tests_dev/test-parallel.R` 覆盖消息流、完成事件、全局/局部静默、progressr opt-out、串行/多进程、参数/命名/简化、RNG 和 condition。两个核心计算套件保留局部静默回归。数值回归还包括 `test-filterCommunication.R`、`test-netpathway-aggregate.R`、`test-centrality-subset.R` 和 `test-computeCellDistance.R`。

真实小数据烟测须覆盖整条链路、普通/详细终端输出、串并行消费者及落盘消费，并与修改前的 `net` / `netP` / `idents` 结果比较；不比较时间和日志元数据。决策证据见 `.agents/notes/implemented/2026-10-19-compute-commun-prob-parallel-output.md`。
