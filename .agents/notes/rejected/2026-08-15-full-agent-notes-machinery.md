# Agent Note: 不照搬完整 DeepSeek Agent Notes 机制

Status: rejected
Date: 2026-08-15
Decision type: process

## Problem

DeepSeek harness 的 `.agents/notes` 包含双语副本、sidecar i18n 文件、生命周期检查、哈希冻结和自动归档等完整机制。SpatialCellChat 当前主要需要可审计的架构记忆，不需要同等复杂的文档基础设施。

## Decision

只采用三类目录：

```text
.agents/notes/proposed/
.agents/notes/implemented/
.agents/notes/rejected/
```

使用单份 Markdown Note；不在当前阶段引入双语同步、sidecar metadata、哈希 manifest 或强制每个小改动归档。

## Why rejected

- 维护机制本身会成为重构负担。
- 双语同步和冻结校验不是当前代码正确性的主要风险。
- 小规模项目使用重型流程会降低记录决策的实际执行率。

## Consequences

如果项目未来需要多人协作、外部贡献或自动审计，再单独提出 Note 讨论是否增加格式检查和归档机制；不能因为 DeepSeek harness 采用了这些机制，就默认迁移全部基础设施。
