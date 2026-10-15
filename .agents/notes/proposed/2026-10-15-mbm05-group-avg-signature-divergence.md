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

**串行复测（2026-10-16，约束"新旧代码两次串行运算"）将结论升级为逐连边证据**：cell 级 prob 张量在 `computeCommunProb` 与 `filterProbability` 两个断点均为**支持集与数值逐位一致**（407,065,015 与 285,801,026 条 max |Δx| = 0）；group 层 8×8×1,563 = 100,032 坐标全空间对比得到支撑差 148（仅旧）/37（仅新）、共同支撑 10,204 中 7,901 条 ≤3.6e-18（相对 ≤6.1e-14）、pval 差 1,033（**0↔1 完全翻转 80 + 1→0 37**，全部集中在 17 个 LR 层）、显著性(<0.05)集合 5,971→5,921（净 −50）。sum 差 0.003509 恰等于 only-old sum 0.007525 − only-new sum 0.004016，精确闭合。详见 `docs/SERIAL_COMPARISON_GROUP_AVG.md`。

Plan A（implemented Note 2026-09-28）的等价契约声明离散结果逐位一致、Prob ≤1 ulp，其证据为 30 层 fixture（max rel 1.43e-14、Pval 零翻转）。1,563 层量级下出现 111 个支撑条目差异、0.13% 数值和差异、以及 117 个 p 值完全翻转，超出现有契约的已验证边界：要幺门控/分母推理存在规模相关缺陷，要幺分歧机制（如显式零保留、 CSC 求和序进入门控边界）尚未被识别。Rscript 兼容 shim 已被排除：两 worker 的 filter 阶段与 pathway 阶段 nnz/sum 完全一致，分歧严格始于 `computeAvgCommunProb`；且新串行 vs 旧串行的 cell 级逐位一致证明差异不可能来自 filter/compute 端。

## Proposal

1. ✅ 冻结全量发现为新的等价性数据点：全套串行 artifacts 已落盘于 `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/`（含各阶段 tensor dump），不再重跑 3.4–7 小时全管线。
2. ✅ mid-stage dump 能力已加入 benchmark worker（`--dump-dir`，post-filter 对象亦落盘供单阶段重放）；并旁证既有 14 个历史条目 checksum 全部匹配。
3. ✅ 串行单阶段对拍完成并分类：支撑差集中在弱边（仅旧侧 x ∈ [2.2e-14, 5.45e-4]，仅新侧最大 7.64e-4），共同支撑差 ≤3.6e-18；pval 完全翻转与支撑差异在相同的 17 个 LR 层出现（HLA-E_KLRC2、IL18_IL18R1_IL18RAP、VTN_ITGA8_ITGB1 等），指向这些层存在大量**置换平局（T_perm ≡ T_obs 完全相等）**，平局的 `>=`/`>` 语义 + 1 ulp 偏差使 reject 计数 100↔0 全跳。
4. ⏳ 待维护者裁定（本 Note 保持 proposed）：是把 Plan A 契约修订为"全量规模下的量化界 + tie 平局语义显式化"，还是先在 Plan B 中做 per-layer 的 T 矩阵级归因与语义统一。若裁定为缺陷需修复，另建修复 Note 并交叉链接。

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
- **串行对拍证据集（2026-10-16）**：`tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/` 下 `{old,new}_prefilter_cell_prob`（7.8GB/4.8GB）、`{old,new}_postfilter_cell_prob`（5.5GB/3.4GB）、`{old,new}_postfilter_object`、`{old,new}_group_prob`、`{old,new}_group_pval`；对比脚本 `tests_dev/benchmarks/compare_serial_tensors.R`，逐连边分析 `analyze_group_diff1..6.R`（含 `analyze_group_edge_diff.R`）。
- cell 级对比结论：prefilter 与 postfilter 两个断点的支持集与数值**逐位一致**（common = 全部条目，only_old/only_new = 0，max |Δx| = 0），输出于 `{prefilter,postfilter}_cell_comparison.csv`。
- group 级结论：支撑差 148/37（弱边，仅旧侧 x ∈ [2.2e-14, 5.45e-4]、仅新侧最大 7.64e-4），sum 差 0.003509 与支撑 sum 精确闭合；共同支撑 10,204 条中 7,901 条 ≤3.6e-18（相对 ≤6.1e-14）；pval 全空间差 1,033（0→1 完全翻转 80、1→0 37，集中 17 个 LR 层；848 条双侧有效且位移 ≤0.5；significance 5,971→5,921），见 `group_summary.csv`、`group_edge_diff.csv` 与 `docs/SERIAL_COMPARISON_GROUP_AVG.md`。
- 串行执行剖面：old 26,204s vs new 1,385s（18.9×，旁证 Plan A 19× 加速声明）；old 串行签名与旧多进程 + shim 运行逐位一致（shim 全量等价旁证）。
- `baseline_result.rds`（3379 B，mtime 2026-10-15 02:40:27，status "success"，无 failure 字段；8 阶段 elapsed 3.44/7060.07/3069.20/2093.79/47.03/38.55/495.29/191.35 秒，workers 1/4/2/2/1/1/2/2）。
- 阶段级排除链：`filterProbability` 两侧 nnz 同为 285,801,026，`computeCommunProbPathway` 两侧 nnz 同为 89,174,034——分歧严格始于 averaging 阶段。
- Driver 于 2026-10-15 由双 artifact 重建全部对比 CSV/TXT/RDS，`Worker failures: none`；运行日志 `baseline_run4.log`、PID 记录 `baseline_run4.pid`、启动脚本 `launch_baseline_run4.ps1` 同目录留档。
- Plan A 契约原文与 30 层证据：`.agents/notes/implemented/2026-09-28-compute-avg-commun-prob-plan-a.md` "组级 prob 的 1 ulp 契约"（fixture ≤3.55e-16、MBM05 30 层 ≤1.43e-14、Pval 零翻转、契约断言 ≤1e-12）。
- 规模盲区确认：本 Note 之前所有等价性证据（160 项单元测试、2026-10-09 e2e、30 层冻结对拍）均为小规模 fixture；2026-10-16 串行全量对比是首次 1,563 层双端成功对比。

## Acceptance criteria

- [x] mid-stage dump（post-filter tensor/object、post-avg prob/pval）加入 worker，`--dump-dir` 可指定且不进 stage 计时；`historical_checksums.json` 既有 14 条全部匹配（2026-10-16 `verify_checksums.py` 复核）。
- [x] 串行全管线对拍完成（同 checkpoint、旧原样代码/新重构代码、无 shim），逐连边 diff 覆盖 filter 前后 cell prob 与 group prob/pval；111+148/37 分歧条目按幅值分类完毕（弱边 + tie 层重叠），见 docs/SERIAL_COMPARISON_GROUP_AVG.md。
- [x] 分类结论写回本 Note 的 Problem/Proposal/Evidence；**Plan A 契约修订与否尚未裁定**（裁定前保持 proposed）。
- [ ] Agent Notes validator 报告 0 errors（每次本 Note 变更后复核）。
