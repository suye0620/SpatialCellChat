# Agent Note: 大型项目 Agent Notes 治理 Skill

Status: implemented
Governance: v1
Date: 2026-08-17
Decision type: process
Scope: .agents/skills/agent-notes-governance and .agents/notes
Owner: project maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

大型重构项目会同时产生架构、API、算法、依赖、测试和流程决策。只依赖聊天记录或任务列表会丢失选择理由、替代方案、消费者影响和验证证据；把这些内容写入通用规则文件又会让硬性规则与历史决策混杂。

## Decision

引入项目级 Skill `.agents/skills/agent-notes-governance`，把 Agent Notes 设计成可被其他 Agent 复用的决策治理流程。该 Skill 定义何时必须写 Note、Note 生命周期、v1 元数据契约、消费者映射、证据要求、审查门禁、实施门禁和渐进式 legacy 迁移策略。

## Constraints and invariants

- `implemented` Note 必须描述已经验证的当前事实，不能描述未来计划。
- `proposed` Note 不能作为已实施设计的证据。
- 一个事实只有一个权威来源：代码和测试定义行为，API 文档定义公开契约，Agent Note 定义决策理由。
- 删除或替换接口前必须列出生产代码、测试、示例、文档、配置、生成物和动态查找消费者。
- `Governance: v1` Note 必须通过 `.agents/skills/agent-notes-governance/scripts/validate_agent_notes.py`。

## Alternatives considered

- 继续使用现有轻量三目录 Notes，不增加 Skill：执行成本最低，但其他 Agent 无法稳定学习何时记录、如何审查和如何过渡状态。
- 直接照搬 DeepSeek Harness 的完整双语、sidecar、哈希冻结和自动归档机制：治理能力更强，但对当前项目过重，维护成本会压低记录率。
- 只写一份 docs 文档，不做 Skill：适合人读，不适合被 Agent 自动触发和执行。

## Consumer impact

- Agent 消费者：后续参与大型开发、架构审查、API 迁移、测试策略和依赖决策的 Agent 需要加载该 Skill。
- Note 消费者：`.agents/notes/{proposed,implemented,rejected}` 继续兼容 legacy Note；新改动优先使用 `Governance: v1`。
- 工具消费者：`new_agent_note.py` 生成 v1 Note；`validate_agent_notes.py` 校验 v1 Note 并对 legacy Note 输出 warning。
- 文档消费者：`references/note-contract.md` 是 Note 契约细节来源；`SKILL.md` 是执行入口。

## Consequences

新增 Skill 会让重大决策多一步记录和校验，但能降低跨会话遗忘、重复争论、错误状态标记和删除接口时漏掉消费者的风险。旧 Note 不强制批量迁移，避免一次性文档债务；严格模式留给未来项目完全采用 v1 后启用。

## Evidence

- `.agents/skills/agent-notes-governance/SKILL.md`
- `.agents/skills/agent-notes-governance/references/note-contract.md`
- `.agents/skills/agent-notes-governance/scripts/new_agent_note.py`
- `.agents/skills/agent-notes-governance/scripts/validate_agent_notes.py`
- `dist/agent-notes-governance.skill`
- `docs/DEEPSEEK_HARNESS_SPATIALCELLCHAT_ORIGINAL.md`
- `.agents/notes/implemented/2026-08-15-lightweight-agent-notes-workflow.md`

## Acceptance criteria

- Skill 元数据通过 `quick_validate.py`。
- `new_agent_note.py` 能生成受治理 v1 Note。
- `validate_agent_notes.py` 能拒绝未填写模板注释的 v1 Note，并接受填写完整的 v1 Note。
- 现有 legacy Notes 在非 strict 模式下只产生 warning，不阻塞渐进采用。
- `.skill` 包只包含 `SKILL.md`、`agents/`、`scripts/` 和 `references/`。
