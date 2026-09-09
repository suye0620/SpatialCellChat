# SpatialCellChat CLI 提示设计规范

> 2026-07-14 · 与 Claude Desktop 讨论后整理

---

## 设计目标

统一重构后所有函数的用户交互输出，消除当前代码中 `cat()` / `message()` / `cli.symbol()` / `cli::cli_alert_*` 四套并存的混乱局面。

---

## 核心原则

```
1. 一个入口初始化        → setEnvironment 统一配置策略/主题/进度条样式
2. 一种消息机制          → 全包走 cli::cli_alert_* 系列, 零 cat()/message()
3. 一条封装函数          → 只保留 .cli() 作为唯一调用点, 删除 cli.symbol()
4. 三级详细度            → verbose=0(静默) / 1(正常) / 2(详细), 全局可设
5. 进度条统一走 progressr → 用 my_future_lapply 内置, 不手工 cli_progress_bar
```

---

## 第 1 层：入口初始化 (`setEnvironment`)

```r
setEnvironment <- function(workers = 4,
                           verbose = getOption("SpatialCellChat.verbose", 1L),
                           future.globals.maxSize = 10000 * 1024^2) {
  if (workers == 1) future::plan("sequential")
  else future::plan("multisession", workers = workers, gc = TRUE)

  options(
    future.globals.maxSize = future.globals.maxSize,
    SpatialCellChat.verbose = verbose
  )

  progressr::handlers(global = TRUE)
  progressr::handlers("cli")

  if (verbose >= 1L) {
    cli::cli_h1("SpatialCellChat")
    if (workers == 1) cli::cli_alert_info("Sequential mode")
    else cli::cli_alert_success("Parallel with {workers} workers")
  }
}
```

**职责**：并行策略 / 全局选项 / 进度条样式 / 入口提示。一次性完成，运行时无副作用。

---

## 第 2 层：唯一消息封装 (`.cli`)

升级现有 `.cli()` 为唯一内部消息函数，增加 verbosity 门控：

```r
.cli <- function(text, .type = "info", ..., .verbose = NULL,
                 .env = parent.frame()) {
  lvl <- .verbose %||% switch(.type,
    info      = 1L,
    success   = 1L,
    warning   = 0L,
    danger    = 0L,
    header    = 1L,
    subheader = 1L,
    text      = 2L,
    progress  = 2L,
    debug     = 3L,
    1L)
  v <- getOption("SpatialCellChat.verbose", 1L)
  if (lvl > v) return(invisible(NULL))

  switch(.type,
    info      = cli::cli_alert_info(text, .envir = .env, ...),
    success   = cli::cli_alert_success(text, .envir = .env, ...),
    danger    = cli::cli_alert_danger(text, .envir = .env, ...),
    warning   = cli::cli_alert_warning(text, .envir = .env, ...),
    header    = cli::cli_h1(text, .envir = .env, ...),
    subheader = cli::cli_h2(text, .envir = .env, ...),
    text      = cli::cli_text(text, .envir = .env, ...),
    cli::cli_alert_info(text, .envir = .env, ...)
  )
}
```

### Verbosity 级别映射

| `.type` | 默认级别 | 语义 |
|---------|---------|------|
| `warning`, `danger` | 0 | 始终显示 |
| `header`, `subheader` | 1 | 函数入口 |
| `info`, `success` | 1 | 关键状态 |
| `text`, `progress` | 2 | 详细步骤 |
| `debug` | 3 | 调试信息 |

### 使用示例

```r
# 函数入口
.cli("Computing communication probabilities", .type = "subheader")        # lv 1

# 关键信息
.cli("Using {nrow(pairLR)} LR pairs", .type = "info")                     # lv 1
.cli("Done", .type = "success")                                           # lv 1

# 详细分支说明（仅 verbose >= 2 显示）
.cli("Contact mode for Cell-Cell Contact signaling", .type = "text")      # lv 2
.cli("Enforcing contact.dependent.forced = TRUE", .type = "text")         # lv 2

# 覆盖指定级别
.cli("Raw matrix dimension: {d}", .type = "text", .verbose = 3)           # lv 3
```

---

## 第 3 层：并行进度条 (`my_future_lapply`)

```r
my_future_lapply <- function(X, FUN, ..., future.seed = TRUE, .progress = TRUE) {
  n <- length(X)
  if (.progress && n > 1 && future::nbrOfWorkers() > 1) {
    p <- progressr::progressor(along = X)
    f <- function(x) { res <- FUN(x); p(); res }
  } else {
    f <- FUN
  }
  future.apply::future_lapply(X, f, ..., future.seed = future.seed)
}
```

**边界规则**：
- progressr 只负责**进度条渲染**，不负责状态消息
- 状态消息由 `.cli()` 在**主进程**中输出
- **workers 内部不调任何 `.cli()` 或 `cat()`** — future multisession 不保证捕获
- 调试信息在 worker 中用 `warning()` 走 condition 机制回主进程

---

## 第 4 层：函数内部模板

每个主要导出函数遵循固定结构：

```r
computeSomething <- function(object, ..., verbose = NULL) {
  # 入口
  .cli("Beginning operation", .type = "subheader")

  # 关键参数摘要 (lv 1)
  .cli("Parameters: range={range}, n_cells={n}", .type = "info")

  # 详细分支 (lv 2)
  .cli("Using distance-based constraints", .type = "text")

  # 并行计算（progressr 内置）
  result <- my_future_lapply(items, FUN = compute_one, .progress = TRUE)

  # 完成
  .cli("Operation complete", .type = "success")
  return(object)
}
```

---

## 第 5 层：Verbosity 级别

| 级别 | 行为 |
|------|------|
| 0 | 静默 — 仅 warning / danger |
| 1 | **正常（默认）** — header + info + success + 进度条 |
| 2 | 详细 — 额外显示分支决策和中间步骤 |
| 3 | 调试 — 所有 debug 信息 |

全局选项：`getOption("SpatialCellChat.verbose", 1L)`
函数级覆盖：`verbose = NULL`（走全局）或显式传入 0/1/2/3

---

## 第 6 层：show 方法

```r
setMethod("show", "SpatialCellChat", function(object) {
  .cli("SpatialCellChat", .type = "header")
  .cli("{n_cells} cells, {n_genes} genes", .type = "info")
  if (has_net) .cli("Communication network inferred", .type = "success")
  else .cli("Communication network not yet computed", .type = "info")
})
```

---

## 不做的事

| ❌ 当前反模式 | ✅ 目标做法 |
|-------------|-----------|
| `cat(cli.symbol(), "text\n")` | `.cli("text")` |
| `cat(paste0("... <<< [", Sys.time(), "]"))` | `.cli("Done", .type = "success")` |
| `cli.symbol("success")` + 拼接 | 直接 `.cli(.type = "success")` |
| `message("提示信息")` | `.cli("提示信息", .type = "info")` |
| workers 内 `cat()` | 删掉或改为 `warning()` 走 condition |
| 用 `cat()` 显示进度 | `my_future_lapply` 内置 progressr |
| `cat("Cell groups:", cli::col_red(...))` | `.cli("Cell groups: {.val {groups}}")` |

---

## 迁移策略

按 **函数依赖顺序** 从底层往上替换，不一把抓：

1. **先改基础设施** — `.cli()` 封装升级、删除 `cli.symbol()`
2. **再改工具函数** — `setEnvironment`、`my_future_lapply`
3. **再改核心计算** — `computeCommunProb` 等 modeling.R
4. **再改下游分析** — analysis.R
5. **最后改可视化** — visualization.R、spatial.R

每个函数替换时遵循：**删旧代码（cat/message/cli.symbol）→ 插入 `.cli()` → 验证**。不保留注释掉的旧输出行。
