# Agent Note: 轻量 Agent Notes 工作流

Status: implemented
Date: 2026-08-15
Decision type: process

## Problem

项目重构中会同时出现代码实现、阶段计划、算法取舍和被明确否决的方案。仅依赖聊天记录会丢失决策理由；把所有内容都写进 `AGENTS.md` 又会把硬性工作规则和项目历史混在一起。

## Decision

采用轻量三目录结构：

```text
.agents/notes/
├── proposed/
├── implemented/
└── rejected/
```

职责分别是：

- `proposed/`：已提出但尚未落地的方案。
- `implemented/`：当前代码已经采用、可作为长期依据的决策。
- `rejected/`：已评估并明确不采用的方案，保留拒绝理由以防重复讨论。

每个重大 Note 使用独立 Markdown 文件，至少包含：`Status`、`Date`、`Decision type`、`Problem`、`Decision/Proposal`、`Alternatives considered` 和 `Consequences`。

只为以下事项写 Note：数据结构、公开 API、算法原理、可视化架构、依赖选择、重大兼容性取舍、测试策略和明确拒绝的方案。机械改名、局部 typo 修复和一次性实验不强制写 Note。

`AGENTS.md` 继续承担每次会话必须遵守的硬规则；`docs/REFACTORING_PLAN.md` 继续承担阶段计划；Notes 只保存长期决策和其理由。

## Alternatives considered

- 复制 DeepSeek harness 的完整双语、sidecar、哈希冻结和自动归档系统：对当前 R 包规模过重。
- 把所有长期知识都放到 `AGENTS.md`：会让规则文件持续膨胀，并混入历史讨论。
- 只依赖 OMP 自动 memory：记忆是启发式背景，不能替代仓库内的可审计决策记录。

## Consequences

后续重大架构和 API 改动应同时更新对应 Note；Note 必须服从当前代码、用户最新指令和权威数据结构文档。若 Note 与代码冲突，应先修正 Note 或明确标为 proposed，不应把过时文字当作实现依据。
