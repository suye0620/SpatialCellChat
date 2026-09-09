# Agent Note: 不保留旧顶层 slot 与兼容别名

Status: rejected
Date: 2026-08-15
Decision type: simplification

## Problem

旧实现和历史讨论中出现过 `data`、`data.signaling`、`data.project`、`options`、`signalingData`、`netData` 等路径。为单个旧调用保留这些路径，会让最终 schema 长期同时维护两套对象模型。

## Decision

不恢复旧顶层 slot，不新增旧 slot 的 re-export、alias 或隐式 fallback。所有调用方迁移到最终 11-slot schema 和访问器后，删除旧路径。

旧 schema 在 git 历史中保留即可，不在源文件中维护大段废弃注释或影子结构。

## Why rejected

- 兼容层会掩盖迁移未完成，而不是暴露真实调用点。
- 同一对象可能出现两个来源不一致的表达矩阵或通信结果。
- 新函数继续复制旧索引逻辑，无法收敛访问契约。
- 测试矩阵和维护成本会随着每个旧路径叠加。

## Consequences

迁移期如果函数仍访问旧 slot，应报告为明确的迁移阻塞；不得把错误包裹成“兼容成功”。旧版本用户需要通过明确的对象迁移或版本边界进入新 schema。
