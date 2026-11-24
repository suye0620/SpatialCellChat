# Agent Note: 精简 Dev 包依赖声明

Status: proposed
Governance: v1
Date: 2026-11-23
Decision type: dependency
Scope: DESCRIPTION, NAMESPACE, R/legacy metadata, renv dependency tracking
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`SpatialCellChatDev` 仍沿用基线 `SpatialCellChat` 的部分依赖声明，但当前 11-slot 代码已经归档了旧表达和 `sparse3Darray` 工具链，并新增了 widget 与 streamline 实现。当前 `DESCRIPTION` 仍声明未被生产代码引用的 `Rfast`、`spatstat.sparse`、`scattermore`、`Rtsne`、`expm`、`progress`，还把仅用于开发或文档的 `rmarkdown`、`knitr`、`roxygen2` 放在 `Suggests`。`RcppEigen` 同时出现在 `LinkingTo` 和 `Suggests`。当前 `NAMESPACE` 还从已归档的 `R/legacy/` 收集了 `spatstat.sparse` 导入。

## Proposal

将 Dev 包的运行时声明按当前可构建的 `R/` 代码裁定。`Depends` 只保留 R 版本约束，将 `dplyr` 和 `igraph` 放入 `Imports` 并继续通过 `NAMESPACE` 提供代码所需符号。保留实际被生产函数调用的核心 `Imports`，删除 `Rfast` 和 `spatstat.sparse`。保留有可选生产入口的 `Suggests`，删除无生产引用的 `scattermore`、`Rtsne`、`expm`、`progress`、`RcppEigen`，并将开发和文档工具移至 `Config/Needs/renv` 或移除包级声明。移除归档文件中的过时 `spatstat.sparse` roxygen 导入，使生成的 `NAMESPACE` 与可构建源码一致。

## Constraints and invariants

- `tests_dev/` 的任何文件、数据和图像不删除，也不改变其测试依赖边界。
- `R/legacy/` 继续作为历史参考保留，但不产生构建包的运行时导入或导出。
- `SpatialCellChatDev` 的核心推断、11-slot 对象、`updateSpatialCellChat()`、三类测试对象和已验收的可视化入口继续可加载。
- `htmlwidgets`、`htmltools`、`memoise` 继续保留在 `Imports`，因为当前生产代码直接调用它们。
- `RcppEigen` 继续保留在 `LinkingTo`，C++ 源码仍必须可以编译。
- 可选入口缺少其 `Suggests` 包时，包加载和核心路径不能因此失败。

## Alternatives considered

- 继续完整复制基线的 `Imports`、`Depends`、`Suggests`。拒绝，因为静态调用证据显示旧工具链和若干旧扩展已经不再被生产源码使用，会增加安装和恢复负担。
- 把所有被生产函数调用的可选绘图与分析包都提升到 `Imports`。拒绝，因为这些入口已有可选功能边界，提升后会让基础安装强制带上 plotly、Shiny、field rendering、topic analysis 等扩展生态。
- 只修改 `DESCRIPTION`，保留旧 `NAMESPACE`。拒绝，因为生成的 namespace 仍会声明已归档 `spatstat.sparse`，导致包元数据与构建源码不一致。

## Consumer impact

生产消费者包括 `R/analysis.R`、`R/modeling.R`、`R/spatial.R`、`R/visualization.R`、`R/visualization_widget.R`、`R/SpatialCellChat_class.R`、`R/utilities.R`、`R/database.R`、`R/merging.R` 及生成的 `NAMESPACE`。开发消费者包括 `devtools`、roxygen 生成流程、`renv/settings.json`、`renv.lock`、README 的安装和构建说明。验证消费者包括三对象 smoke test、`tests_dev/test-updateSpatialCellChat.R`、已安装可视化测试、源码包构建和 `R CMD check`。`R/legacy/` 的归档文件只保留历史参考用途。

## Consequences

终端用户基础安装将少拉取未使用的旧依赖和开发工具。可选绘图、topic、SCE 和交互式功能继续按需安装对应 `Suggests`。开发者若要生成文档或渲染教程，应从 `Config/Needs/renv` 或 README 的开发环境步骤安装对应工具。依赖声明修改后需要重新生成 `NAMESPACE`，并核对 `renv.lock` 是否仍保持显式记录策略和开发环境可恢复性。

## Evidence

- 基线元数据：`renv/library/windows/R-4.4/x86_64-w64-mingw32/SpatialCellChat/DESCRIPTION`。
- Dev 元数据：`DESCRIPTION`。
- 当前 namespace：`NAMESPACE`，其中 `spatstat.sparse` 导入来自 `R/legacy/expression-toolchain.R`。
- 生产调用扫描：除 `R/legacy/` 和 `R/test_helpers.R` 后，`Rfast`、`spatstat.sparse`、`scattermore`、`Rtsne`、`expm`、`progress` 没有生产代码引用。
- 当前生产代码对 `memoise`、`htmlwidgets`、`htmltools` 有直接调用，对 `RcppEigen` 有 `src/SpatialChat_Rcpp.cpp` 的头文件依赖。
- `R CMD check` 当前在完整 Suggests 强制检查下报告多个未安装 Suggests，说明当前声明还混合了基础安装与未安装可选扩展。

## Acceptance criteria

- `DESCRIPTION` 的 `Depends`、`Imports`、`Suggests` 与裁定后的生产代码引用一致，且 `Rfast`、`spatstat.sparse`、无引用旧扩展、开发/文档工具不再错误地作为包运行时要求。
- 重新生成的 `NAMESPACE` 不再导入 `spatstat.sparse`，也不产生缺少声明的非基础包导入。
- 在项目 R 4.4.3 环境中完成源码构建、安装、包加载和 `R CMD check`，核心检查无 errors。
- `updateSpatialCellChat()` 的 3cts/4cts 迁移 smoke test 和至少一个 Visium 核心 workflow 继续通过。
- `git diff --check` 通过，`tests_dev/` 文件集合不减少。
- Agent Notes validator 对此 Note 通过。
