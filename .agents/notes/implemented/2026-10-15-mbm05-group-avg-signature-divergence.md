# Agent Note: MBM05 全量规模 group 平均签名分歧

Status: implemented
Governance: v1
Date: 2026-10-15
Decision type: algorithm
Scope: `R/modeling.R` `computeAvgCommunProb`（Plan A kernel）、`src/SpatialChat_Rcpp.cpp` 群组聚合路径、benchmark 签名契约与下游 `net$group` 消费者
Owner: SpatialCellChat maintainer
Impact: medium
Supersedes: .agents/notes/archived/2026-09-28-compute-avg-commun-prob-plan-a.md
Superseded by: none

## Problem

全量 MBM05 基准（1563 LR × 29,536 cell，两个 worker 于 2026-10-15 首次全部成功）的签名对比显示：refactored 与 pinned-upstream 流水线的唯一分歧点在群组平均阶段。refactored `net$group$prob` 为 nnz 10,241 / sum 2.636712，upstream 为 nnz 10,352 / sum 2.640221——支撑条目差 111（−1.07%），数值和差 0.0035（−0.13%）。其余全部签名项逐位一致：cell 级（nnz 285,801,026、sum 14031.85）、pathway 级（nnz 89,174,034、sum 14031.85）、centrality 维度、pathway 列表。

**串行复测（2026-10-16，约束"新旧代码两次串行运算"）将结论升级为逐连边证据**：cell 级 prob 张量在 `computeCommunProb` 与 `filterProbability` 两个断点均为**支持集与数值逐位一致**（407,065,015 与 285,801,026 条 max |Δx| = 0）；group 层 8×8×1,563 = 100,032 坐标全空间对比得到支撑差 148（仅旧）/37（仅新）、共同支撑 10,204 中 7,901 条 ≤3.6e-18（相对 ≤6.1e-14）、pval 差 1,033（**0↔1 完全翻转 80 + 1→0 37**，全部集中在 17 个 LR 层）、显著性(<0.05)集合 5,971→5,921（净 −50）。sum 差 0.003509 恰等于 only-old sum 0.007525 − only-new sum 0.004016，精确闭合。详见 `docs/SERIAL_COMPARISON_GROUP_AVG.md`。

Plan A（implemented Note 2026-09-28）的等价契约声明离散结果逐位一致、Prob ≤1 ulp，其证据为 30 层 fixture（max rel 1.43e-14、Pval 零翻转）。1,563 层量级下出现 111 个支撑条目差异、0.13% 数值和差异、以及 117 个 p 值完全翻转，超出现有契约的已验证边界。串行组件级复核已将 148/37 条支撑差异逐条匹配到 pinned upstream 的 `format(dataLR_temp[, -1], digits = 1) >= min.percent` 字符门控与新实现 `signif(mean, 1)` 数值门控之间的语义差异；共同支撑的约 1 ulp 尾差仍来自 CSC 单遍累加与 Matrix/BLAS 两段归约顺序不同。Rscript 兼容 shim 已被排除：两 worker 的 filter 阶段与 pathway 阶段 nnz/sum 完全一致，分歧严格始于 `computeAvgCommunProb`；且新串行 vs 旧串行的 cell 级逐位一致证明差异不可能来自 filter/compute 端。

## Decision

1. **采用 numeric `min.percent` gate**：对每个组的配体/受体阳性比例执行 `signif(mean, 1) >= min.percent`；R 侧按组预计算整数阈值，C++ observed/permutation kernel 使用 `cnt >= thr`。`signif` 的 `1` 是一位有效数字参数，不是比例阈值。
2. **将 pinned upstream 的 `format(..., digits = 1)` 字符比较定义为历史行为**，不再要求新实现复现其 group support 或 p-value 结果；全量 148/37 条支撑差异及相关 p 值传播作为已知、可解释的行为变更记录。
3. **保留 Plan A kernel 和置换检验方法**：RNG 生成、null 构造、单侧 `Pboot - Pnull > 0`、`nReject / nboot`、`Pval[Prob == 0] <- 1`、CSC 聚合和 OpenMP 输出隔离均不变。
4. **以 numeric reference 和数值界验收**：`tests_dev/reference-computeAvgCommunProb-numeric.R`、`.spec_boot` fixture 以及 MBM05 层子样本验证 group support/Pval 语义；共同支撑 `Prob` 允许已测 `≤1e-12` 归约差。全量 pinned-upstream artifact 保留为历史比较证据。

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
- 治理契约：本 Note 取代 `.agents/notes/archived/2026-09-28-compute-avg-commun-prob-plan-a.md` 中关于 pinned-upstream gate/Pval 严格逐位等价的过时边界；kernel 加速、RNG、null 构造、`Pboot - Pnull > 0` 和 `Pval[Prob==0] <- 1` 保持不变。
- 测试 oracle：`tests_dev/reference-computeAvgCommunProb-numeric.R` 是正式 numeric reference；`tests_dev/test-computeCommunProb.R` 的 format/signif fixture 仅记录历史差异。

## Consequences

- numeric gate 行为变更正式接受：148/37 条支撑差异及其 p 值传播不再作为兼容性失败；下游消费者继续读取相同的 `net$group$prob` / `net$group$pval` schema，但必须按 numeric gate 解释结果。
- 共同支撑的 `Prob` 仍采用 `≤1e-12` 数值界；严格 p 值比较的归约边界继续单独记录，不改 `>`、分母或 null。
- 全量签名对比仍是 `computeAvgCommunProb` 改动后的验收证据，但 pinned-upstream 差异按已裁定的 numeric contract 解释。

## Evidence

- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/signature_comparison.csv` 与 `signatures.csv`：12 项指标 10 项 TRUE；`group_nnz`（10241 vs 10352）与 `group_sum`（2.636712 vs 2.640221）FALSE；`comparison.txt` 记录 `All output signatures match: FALSE`。
- **串行对拍证据集（2026-10-16）**：`tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/` 下 `{old,new}_prefilter_cell_prob`（7.8GB/4.8GB）、`{old,new}_postfilter_cell_prob`（5.5GB/3.4GB）、`{old,new}_postfilter_object`、`{old,new}_group_prob`、`{old,new}_group_pval`；对比脚本 `tests_dev/benchmarks/compare_serial_tensors.R`，逐连边分析 `analyze_group_diff1..6.R`（含 `analyze_group_edge_diff.R`）。
- cell 级对比结论：prefilter 与 postfilter 两个断点的支持集与数值**逐位一致**（common = 全部条目，only_old/only_new = 0，max |Δx| = 0），输出于 `{prefilter,postfilter}_cell_comparison.csv`。
- group 级结论：支撑差 148/37（弱边，仅旧侧 x ∈ [2.2e-14, 5.45e-4]、仅新侧最大 7.64e-4），sum 差 0.003509 与支撑 sum 精确闭合；共同支撑 10,204 条中 7,901 条 ≤3.6e-18（相对 ≤6.1e-14）；pval 全空间差 1,033（0→1 完全翻转 80、1→0 37，集中 17 个 LR 层；848 条双侧有效且位移 ≤0.5；significance 5,971→5,921），见 `group_summary.csv`、`group_edge_diff.csv` 与 `docs/SERIAL_COMPARISON_GROUP_AVG.md`。
- 串行执行剖面：old 26,204s vs new 1,385s（18.9×，旁证 Plan A 19× 加速声明）；old 串行签名与旧多进程 + shim 运行逐位一致（shim 全量等价旁证）。
- `baseline_result.rds`（3379 B，mtime 2026-10-15 02:40:27，status "success"，无 failure 字段；8 阶段 elapsed 3.44/7060.07/3069.20/2093.79/47.03/38.55/495.29/191.35 秒，workers 1/4/2/2/1/1/2/2）。
- 阶段级排除链：`filterProbability` 两侧 nnz 同为 285,801,026，`computeCommunProbPathway` 两侧 nnz 同为 89,174,034——分歧严格始于 averaging 阶段。
- Driver 于 2026-10-15 由双 artifact 重建全部对比 CSV/TXT/RDS，`Worker failures: none`；运行日志 `baseline_run4.log`、PID 记录 `baseline_run4.pid`、启动脚本 `launch_baseline_run4.ps1` 同目录留档。
- Plan A 历史契约与 30 层证据：`.agents/notes/archived/2026-09-28-compute-avg-commun-prob-plan-a.md`（fixture ≤3.55e-16、MBM05 30 层 ≤1.43e-14、Pval 零翻转、契约断言 ≤1e-12）；当前 numeric contract 的直接验证见 `tests_dev/test-computeAvgCommunProb.R`。
- 规模盲区确认：本 Note 之前所有等价性证据（160 项单元测试、2026-10-09 e2e、30 层冻结对拍）均为小规模 fixture；2026-10-16 串行全量对比是首次 1,563 层双端成功对比。

## Acceptance criteria

- [x] mid-stage dump（post-filter tensor/object、post-avg prob/pval）加入 worker，`--dump-dir` 可指定且不进 stage 计时；`historical_checksums.json` 既有 14 条全部匹配（2026-10-16 `verify_checksums.py` 复核）。
- [x] 串行全管线对拍完成（同 checkpoint、旧原样代码/新重构代码、无 shim），逐连边 diff 覆盖 filter 前后 cell prob 与 group prob/pval；111+148/37 分歧条目按门控语义和幅值分类完毕（`format()`/`signif()` gate 差异 + 共同支撑 1 ulp 尾差），见 `docs/SERIAL_COMPARISON_GROUP_AVG.md` 与 `docs/COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md`。
- [x] 分类结论写回本 Note 的 Problem/Decision/Evidence；numeric `signif()` gate 已由维护者选择，旧 Plan A 的 pinned-upstream 离散等价边界由本 Note 取代。
- [x] `tests_dev/test-computeAvgCommunProb.R` 与 numeric reference 通过；R/C++ 注释、公开参数文档与优化文档同步记录 numeric contract。MBM05 层子样本的既有 2026-09-28 对拍证据保持有效，因本次仅重命名 reference 文件且未改变数值逻辑。
