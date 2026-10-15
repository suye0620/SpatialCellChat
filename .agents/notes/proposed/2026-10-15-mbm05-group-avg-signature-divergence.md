# Agent Note: MBM05 全量规模 group 平均签名分歧

Status: proposed
Governance: v1
Date: 2026-10-15
Decision type: algorithm
Scope: `R/modeling.R` `computeAvgCommunProb`（Plan A kernel）、`src/SpatialChat_Rcpp.cpp` 群组聚合路径、benchmark 签名契约与下游 `net$group` 消费者
Owner: SpatialCellChat maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

全量 MBM05 基准（1563 LR × 29,536 cell，两个 worker 于 2026-10-15 首次全部成功）的签名对比显示：refactored 与 pinned-upstream 流水线的唯一分歧点在群组平均阶段。refactored `net$group$prob` 为 nnz 10,241 / sum 2.636712，upstream 为 nnz 10,352 / sum 2.640221——支撑条目差 111（−1.07%），数值和差 0.0035（−0.13%）。其余全部签名项逐位一致：cell 级（nnz 285,801,026、sum 14031.85）、pathway 级（nnz 89,174,034、sum 14031.85）、centrality 维度、pathway 列表。

Plan A（implemented Note 2026-09-28）的等价契约声明离散结果逐位一致、Prob ≤1 ulp，其证据为 30 层 fixture（max rel 1.43e-14、Pval 零翻转）。1,563 层量级下出现 111 个支撑条目差异与 0.13% 数值和差异，超出现有契约的已验证边界：要么门控/分母推理存在规模相关缺陷，要么分歧机制（如显式零保留、 CSC 求和序进入门控边界）尚未被识别。Rscript 兼容 shim 已被排除：两 worker 的 filter 阶段与 pathway 阶段 nnz/sum 完全一致，分歧严格始于 `computeAvgCommunProb`。

## Proposal

1. 冻结本次全量发现为新的等价性数据点；不再重跑 3.4 小时全管线复现它（artifacts 与 checksums 已保全）。
2. 为 benchmark worker 增加 mid-stage checkpoint（post-filter、post-avg），使 tensor 级逐条目对拍无需整轮重跑；现存 artifacts 只含签名，不含 group tensor。
3. 在同一 checkpoint 上做聚焦单阶段对拍：相同过滤输入下分别以 upstream `computeAvgCommunProb_LR_Avg` 与 Plan A kernel 重建群组输出，逐条目 diff，并对 111 个分歧条目分类（幅值量级、是否显式零、门控决策差异）。
4. 按分类结果裁定：若差异限于已文档化的数值噪声边界且门控决策一致，则把契约补全为"离散一致 + 全量规模数字界"；若门控或支撑语义存在计划外分歧，按 Plan A 缺陷处理，修复进入 Plan B 波次。

## Constraints and invariants

- 全量成功 artifact（`baseline_result.rds`，status "success"，2026-10-15 02:40:27）为不可覆写的原始证据；后续对拍不得破坏既有基准产物与 `historical_checksums.json` 既有条目。
- pinned upstream SHA `be88f300464cf970b36a6165ea7d2e41a47a7511` 不变；upstream oracle 对拍继续使用其自有函数。
- 任何契约修订必须同时陈述离散结果规则与全量规模实测数字界。
- Plan B 统计语义修正（空间约束 null、(b+1)/(N+1)、FDR）仍不在本 Note 范围内，除非分类证明存在缺陷需要连带处理。

## Alternatives considered

- 不测量、直接假设差异属已文档化 ulp 效应：拒绝——支撑条目（nnz）差异是离散差异，不是数值噪声。
- 立即重跑全管线并 dump tensor：拒绝作为第一步——3.4+ 小时成本可由 mid-stage checkpoint 后的单阶段聚焦对拍替代。
- 回退 Plan A kernel 为基线 crossprod 路径：拒绝——在分歧机制未识别前放弃 19× 加速不成比例。

## Consumer impact

- `net$group$prob` / `net$group$pval` 消费者：下游 filter 通信筛选、aggregateNet、centrality、组级可视化——数值传播严格始于该阶段（cell/pathway 签名证明上游边界未受影响）。
- Benchmark harness：`tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` 将获得 mid-stage checkpoint 能力。
- 治理契约：`.agents/notes/implemented/2026-09-28-compute-avg-commun-prob-plan-a.md` 的"组级 prob 的 1 ulp 契约"为被检对象；`.agents/notes/implemented/2026-10-13-mbm05-rscript-baseline-compatibility.md` 提供成功全量基线，交叉引用本 Note。
- 文档：DATA_STRUCTURE.md 的 11-slot 层级关系不受影响；PERMUTATION_TEST_AUDIT.md 的 Plan B 语义修正不因本 Note 自动生效。

## Consequences

- 若分类证明差异限于 ulp 契约内：全量基准签名解读口径更新（group 域允许已测数字界内的相对差），汇报与合作者沟通按新契约表述。
- 若发现计划外语义差：kernel 门控/den/显式零语义修复进入 Plan B，R 层测试契约（≤1e-12、identical 断言）同步修订，且需新增全量规模护栏。
- 全量签名对比从此成为 `computeAvgCommunProb` 改动后的常规验收证据：任何改动先以 mid-stage checkpoint 对拍支撑条目。

## Evidence

- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/signature_comparison.csv` 与 `signatures.csv`：12 项指标 10 项 TRUE；`group_nnz`（10241 vs 10352）与 `group_sum`（2.636712 vs 2.640221）FALSE；`comparison.txt` 记录 `All output signatures match: FALSE`。
- `baseline_result.rds`（3379 B，mtime 2026-10-15 02:40:27，status "success"，无 failure 字段；8 阶段 elapsed 3.44/7060.07/3069.20/2093.79/47.03/38.55/495.29/191.35 秒，workers 1/4/2/2/1/1/2/2）。
- 阶段级排除链：`filterProbability` 两侧 nnz 同为 285,801,026，`computeCommunProbPathway` 两侧 nnz 同为 89,174,034——分歧严格始于 averaging 阶段。
- Driver 于 2026-10-15 由双 artifact 重建全部对比 CSV/TXT/RDS，`Worker failures: none`；运行日志 `baseline_run4.log`、PID 记录 `baseline_run4.pid`、启动脚本 `launch_baseline_run4.ps1` 同目录留档。
- Plan A 契约原文与 30 层证据：`.agents/notes/implemented/2026-09-28-compute-avg-commun-prob-plan-a.md` "组级 prob 的 1 ulp 契约"（fixture ≤3.55e-16、MBM05 30 层 ≤1.43e-14、Pval 零翻转、契约断言 ≤1e-12）。
- 规模盲区确认：本 Note 之前所有等价性证据（160 项单元测试、2026-10-09 e2e、30 层冻结对拍）均为小规模 fixture；这是首次 1563 层双端成功对比。

## Acceptance criteria

- [ ] mid-stage checkpoint（post-filter、post-avg）加入 worker 且不改变各阶段计时语义；`historical_checksums.json` 既有条目保持匹配。
- [ ] 单阶段对拍在相同过滤输入上运行 upstream `computeAvgCommunProb_LR_Avg` 与 Plan A kernel，产出逐条目 diff 与 111 个分歧条目的分类（幅值、是否显式零、门控决策）。
- [ ] 分类结论写回本 Note 的 Decision/Evidence，并明确 Plan A 契约是否修订；若裁定为缺陷，另建修复 Note 并交叉链接。
- [ ] Agent Notes validator 报告 0 errors。
