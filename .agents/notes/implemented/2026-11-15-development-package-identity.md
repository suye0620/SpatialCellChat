# Agent Note: 开发包使用 SpatialCellChatDev 身份与上游基线共存

Status: implemented
Governance: v1
Date: 2026-11-15
Decision type: dependency
Scope: DESCRIPTION, NAMESPACE, compiled package identity, Rcpp DLL registration, tests, documentation, renv installation boundary
Owner: project maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

当前工作区的开发代码仍声明为 `Package: SpatialCellChat`，与上游基线包 `jinworks/SpatialCellChat@be88f300464cf970b36a6165ea7d2e41a47a7511` 使用相同的安装目录和 namespace 身份。两者无法在同一 R library 中安全安装；即使分别保存源码，测试进程也容易因同名 namespace 加载错误包。

当前开发包的包身份声明在 `DESCRIPTION:1`，动态库身份在 `NAMESPACE`、`R/SpatialCellChat-package.R` 和 `R/RcppExports.R`/`src/RcppExports.cpp` 中。已恢复的上游包是 15-slot `SpatialCellChat` 类，而当前开发代码是 11-slot schema；同名包会使基线对象迁移和开发包回归测试混用。

## Decision

将当前开发包的安装和 namespace 身份改为 `SpatialCellChatDev`，保留上游包的正式身份 `SpatialCellChat`。更新 DESCRIPTION、NAMESPACE、Rcpp 动态库注册、包资源解析、开发测试和开发文档中的包身份引用；保留 S4 类名 `SpatialCellChat` 及其对象兼容边界，不进行无证据的全局类名重命名。

开发包与上游包可以安装在同一个 R library，但涉及两个 `SpatialCellChat` 类定义的比较、旧对象迁移和基线对拍必须在独立 R 进程中执行。后续若需要同一 R 进程同时持有两种 S4 类，另行提出类名/继承方案，不由本次包名修改隐式解决。

## Constraints and invariants

- 上游包安装名和 namespace 仍为 `SpatialCellChat`，固定比较提交为 `be88f300464cf970b36a6165ea7d2e41a47a7511`。
- 开发包安装名、namespace 和已编译动态库注册名统一为 `SpatialCellChatDev`；不能只修改 DESCRIPTION 而遗留 `R_init_SpatialCellChat` 或 `PACKAGE = "SpatialCellChat"`。
- 当前 11-slot S4 类的槽位、对象字段和公共函数语义不因包重命名改变；S4 类名保留为 `SpatialCellChat`，以免破坏现有序列化对象和内部方法签名。
- 两个包在同一 library 中均可被发现和安装；开发包加载后不能从 `SpatialCellChat` namespace 解析开发代码。
- `SpatialCellChatDev` 的资源、htmlwidget 和 native routine 查找必须指向开发包自身。
- 本次变更不修改统计定义、通信结果、缓存语义或可视化算法。

## Alternatives considered

- 继续使用 `SpatialCellChat`：拒绝。会覆盖或遮蔽上游包，无法可靠进行基线比较和旧对象迁移。
- 只修改 DESCRIPTION，不修改 namespace、Rcpp 和资源路径：拒绝。安装可能成功，但 native routine、`system.file()` 和运行时 namespace 会指向错误身份。
- 将 S4 类名也改为 `SpatialCellChatDev`：本次拒绝。它会改变序列化对象的类身份、`methods::is()` 判断、方法签名和大量内部消费者，属于独立的对象 schema/API 迁移，应有单独决策和迁移方案。
- 在同一 R 进程同时加载两个包并依赖同名 S4 类自动共存：拒绝作为本次契约。S4 类注册是全局运行时边界，必须先证明类定义和方法分派不会互相覆盖。

## Consumer impact

- 生产和安装元数据：`DESCRIPTION`、`NAMESPACE`、`R/SpatialCellChat-package.R`、`R/RcppExports.R`、`src/RcppExports.cpp`、编译生成的 `src/symbols.rds`。
- 运行时资源：`R/visualization_widget.R` 中按包名查找 htmlwidget 资源的路径，以及任何 `system.file(package = ...)` 调用。
- 开发测试和脚本：所有加载当前源码包的 `tests_dev/*.R`、安装包 smoke test、基线比较脚本；上游对拍脚本必须继续显式指向 `SpatialCellChat`，开发脚本显式指向 `SpatialCellChatDev`。
- 文档和示例：README、开发文档、安装说明、测试说明中涉及当前开发包名称的内容。
- 配置与落盘边界：`renv` 中当前开发包的本地安装记录和独立上游基线 library；不把两个同名包重新写回同一安装目录。
- S4 类消费者：`methods::is(..., "SpatialCellChat")`、`setClass`、方法签名、对象读写和所有 11-slot 访问器；这些保留类名但必须在独立进程边界下使用。

## Consequences

开发包可与上游包并列安装，后续可以明确区分“上游旧算法/旧 schema”和“开发版 11-slot 实现”。安装和测试命令需要显式区分两个包名，现有直接写 `library(SpatialCellChat)` 的开发脚本必须更新，否则会继续加载基线包。

本次不解决同一 R 进程中的同名 S4 类共存；这会保留一个明确的后续边界，但避免用类名重命名破坏已有对象。Rcpp 动态库和 htmlwidget 资源身份必须同步修改，否则会出现安装成功但运行时加载错误的隐蔽故障。

## Evidence

- 当前开发包元数据：`DESCRIPTION:1` 声明 `Package: SpatialCellChat`。
- 当前开发包动态库声明：`NAMESPACE` 的 `useDynLib(SpatialCellChat, ...)`，`R/SpatialCellChat-package.R` 的 `@useDynLib SpatialCellChat`。
- 当前开发包已实现的 11-slot 类：[R/SpatialCellChat_class.R](/D:/AIProjects/SpatialCellChat/R/SpatialCellChat_class.R:323)。
- 当前开发包的占位迁移函数：[R/merging.R](/D:/AIProjects/SpatialCellChat/R/merging.R:453)。
- 上游基线记录：[.agents/notes/implemented/2026-08-23-spatialcellchat-upstream-baseline.md](/D:/AIProjects/SpatialCellChat/.agents/notes/implemented/2026-08-23-spatialcellchat-upstream-baseline.md)，固定提交为 `be88f300464cf970b36a6165ea7d2e41a47a7511`。
- 上游固定提交源码的 `R/SpatialCellChat_class.R`：15-slot 类和旧版 `updateSpatialCellChat()`。
- 已恢复的上游独立安装路径：`C:/Users/jinLab1/AppData/Local/Temp/SpatialCellChat-upstream-baseline-git-20261115212001/library/SpatialCellChat`。

## Acceptance criteria

- [x] 等价安装检查显示开发包名为 `SpatialCellChatDev`，且 `packageDescription("SpatialCellChatDev")$Package` 正确。
- [x] `SpatialCellChat` 与 `SpatialCellChatDev` 已安装在同一个临时 library，`find.package()` 返回两个不同目录。
- [x] 在独立 R 进程中分别加载两个包：上游类为 15-slot，开发类为当前 11-slot；开发包源码和已安装包均无 `SpatialCellChat::` 自调用。
- [x] 开发包的 Rcpp routine 能正常加载；`R_init_SpatialCellChatDev`、`_SpatialCellChatDev_*` 和 `NAMESPACE` 的动态库声明一致，`cpp_sum_layers()` 数值检查通过。
- [x] 开发包的 htmlwidget 静态资源从 `SpatialCellChatDev` 自身解析，`test-commun-flow-widget.R` 通过。
- [x] 需要使用当前实现的已安装 smoke test 改为 `SpatialCellChatDev`，`test-streamline-installed.R` 通过；上游基线仍保留为 `SpatialCellChat`。
- [x] 本 Note 在 Agent Notes validator 中报告为 `OK`；仓库级 validator 仍因既有 8 个旧 Note 的 20 个格式错误退出 1，这些错误与本次包身份变更无关，未被本次变更掩盖或修改。
