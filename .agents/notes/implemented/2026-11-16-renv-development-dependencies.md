# Agent Note: 保留开发工具和实际使用的可选依赖

Status: implemented
Governance: v1
Date: 2026-11-16
Decision type: dependency
Scope: DESCRIPTION, renv settings, renv.lock, development library
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

项目使用 explicit 快照，仅扫描强制依赖字段。devtools 是开发工具而非包运行时依赖，因此没有被记录；已安装并写入锁文件的 RcppML、metR、oce、spdep 及其传递依赖被 renv::status() 报告为未使用。直接按原规则快照会从锁文件移除这些可选功能依赖，影响可恢复性。

## Decision

保持 explicit 快照，在项目 DESCRIPTION 增加 `Config/Needs/renv`，声明 devtools 和当前已安装且列于 Suggests 的 RcppML、metR、oce、spdep，并将该字段加入 renv 的依赖字段；将它们的实际已安装版本和强制传递依赖记录到锁文件。保留既有 R 4.4.3、Bioconductor 3.22、GitHub 固定版本及其他已锁定依赖。

## Constraints and invariants

- devtools 仅作为开发环境依赖，不加入运行时 Imports/Depends。
- 不安装全部 Suggests；未声明的可选功能仍需按需安装。
- 保留既有锁定记录、GitHub SHA 和 R/Bioconductor 配置，不对整个环境 restore 或升级。
- 上游 SpatialCellChat 和 SpatialCellChatDev 仍可独立加载并保留原类定义。
- 同步状态必须通过实际 renv::status() 检查，不通过忽略依赖或删除有用包掩盖问题。

## Alternatives considered

- 启用所有 Suggests：会扩展到大量未安装的可选依赖，超出本次范围。
- 仅按原 explicit 规则快照：会将已使用的可选依赖从锁文件移除，且遗漏 devtools。
- 切换 all 快照：会将本地临时开发包和不相关安装包全部纳入锁文件，降低可复现性。

## Consumer impact

开发者的 devtools 命令、renv::status/snapshot/restore、项目 DESCRIPTION 与 renv/settings.json、锁文件消费者和依赖上述可选包的绘图/分析功能。Dev 包与基线包的安装、S4 类和算法契约不因本决策改变。

## Consequences

新增或移除需要稳定复现的开发/可选依赖时，须维护 `Config/Needs/renv`。其他 Suggests 不自动保证在干净环境可用。原有 Bioconductor 3.22 与 R 4.4 Windows 仓库可用性仍需单独处理，不能以状态同步替代干净恢复验证。

## Evidence

- `renv/settings.json` 的 explicit 快照和 Imports/Depends/LinkingTo 字段。
- 变更前 `renv::status()` 报告 11 个已安装、已记录但未判定使用的包；四个根包均列在 DESCRIPTION Suggests。
- 项目库安装 `devtools` 2.5.1 和 `pkgbuild` 1.4.8，后者识别 Rtools 4.4 并成功编译测试程序。
- `.agents/notes/implemented/2026-10-11-spatialcellchat-renv-r44-restore.md` 规定 R 4.4.3、Bioconductor 3.22 和 GitHub SHA 基线。

- 临时快照仅新增 38 条记录；全部原有锁定记录、仓库配置、R 和 Bioconductor 版本与两个 GitHub SHA 均保持不变。
- `renv::status()` 返回 `SYNCHRONIZED=TRUE`；项目库 devtools 从精简源码副本成功构建源码包。独立进程分别验证基线 15-slot 和开发版 11-slot。
- 普通快照因 Bioconductor 3.22 的 R 4.4 Windows 索引不可用而拒绝；在临时锁文件使用 `force=TRUE` 并核对全部差异后才更新项目锁文件。未验证全新环境恢复。

## Acceptance criteria

- [x] 项目 DESCRIPTION 和 renv 配置包含该开发依赖字段，且未改变包运行时依赖。
- [x] 锁文件保留原有包记录和 GitHub SHA，并记录 devtools 及其传递依赖；R 及 Bioconductor 版本不变。
- [x] 项目库中 devtools 与 Rtools 可用，独立进程检查上游 15-slot、开发版 11-slot 不变。
- [x] `renv::status()` 实际报告同步，或保留明确记录的可复现阻碍。
- [x] Agent Notes validator 对本 Note 无新增错误。
