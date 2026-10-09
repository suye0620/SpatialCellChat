# Agent Note: 模拟数据端到端一致性对拍（Wave 2b 验收）

Status: implemented
Governance: v1
Date: 2026-10-09
Decision type: testing
Scope: tests_dev/simulated_data; inference chain equivalence evidence (Wave 2b: filterCommunication / computeCommunProbPathway / aggregateNet / relabelSpatialCellChat)
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

Wave 2b 迁移（见 implemented note 2026-10-08-netpathway-aggregate-migration）的验收证据目前只有 fixture 级（手工小矩阵）测试与单元回归。缺少一条**真实规模、真实 spatial 参数**的端到端对拍：同一台机器、同一份模拟输入、同一组 seed，重构链路 vs 冻结基线包（be88f30）的完整推断产物逐槽位对拍。没有这个证据，"离散结果逐位一致、组级浮点 ≤1 ulp"的等价性契约只在合成 fixture 上成立过，不能宣称在真实数据规模（2000 细胞 × 35 LR × ~1960 万 nnz）上成立。

已有的 SpatialChat_1.rds 是 Linux 服务器产物，与 Windows 本机重跑结果无跨平台可比性（libm/BLAS 归约序差异不可移植）；它不能作为参照。

## Decision

采用"同机重建 + 同 seed 镜像 + 双重对拍"方案：

1. `tests_dev/simulated_data/run_s3_chain.R`：按 baseline_inference_compare.R 的构建序列从原始输入（countmat_1.csv / metadata_1.csv）重建 SpatialChat_3（5 项输入预检 bitwise），用重构后 R/ 源码跑 7 步链路（参数逐项镜像 _2@options$parameter），存 SpatialChat_3.rds。
2. `tests_dev/simulated_data/compare_s3_s2.R`：SpatialChat_3（11-slot）vs SpatialChat_2（legacy）逐层对拍，产出 e2e_consistency_layers.csv。
3. `tests_dev/simulated_data/rigor_check.R`：质疑驱动的补验脚本——层名锚定对齐、支持集断言、nnz 总账、按名匹配的 IEEE 精确值重验。

拒绝的替代方案（见 Alternatives）：跨平台对拍（_1）与"仅 fixture 测试不再做规模对拍"。

对拍判定标准（预注册）：
- 细胞级 prob：IEEE 逐位相等（`==` 断言，非 all.equal、无容差）；
- 组级 prob：≤1 ulp（≤1e-12 断言），因为 Rcpp 分层求和序 ≠ 基线 dgemm 归约序（Plan A 契约）；
- pval：整数计数（1/nboot 步进），|Δpval| ≤ 1 boot 且 pval<0.05 显著性翻转 = 0；
- 离散结果（count/LR.sig/pathways 成员与排序/稀疏支持集）：逐位一致。

## Constraints and invariants

- 输入同源性不变量：_3 重建后的 `assay$norm` ↔ `_2@data`、`assay$signaling` ↔ `_2@data.signaling`、coordinates、DB$interaction、LR$LRsig 必须逐位一致；否则对拍作废。
- 层对齐不变量：任何逐层比较必须经过名字锚定（_2 侧 `net$tmp$prob.cell` 的名字与顺序 == _3 侧 `dimnames(net$cell$prob)[[3]]` == `rownames(LRsig)`）；基线 sparse3Darray 的 z 维名为空，位置 k↔k 对齐不可作为证据，只能由 `net$tmp$prob.cell` 名单 + 逐 k nnz 佐证。
- seed 镜像不变量：filterProbability(nboot=100, seed=666)、computeAvgCommunProb(nboot=100, seed=1, avg.type="sum")——任何一侧 seed/参数不同，pval 比较即失效。
- pval 语义不变量：Plan B 统计修复实施前，新旧置换检验语义相同（`replicate(nboot, sample.int(nC,nC))`）；因此任何 pval 差异必须可归因于 1-ulp 浮点差穿越 `perm == obs` 边界的单 boot 翻转，出现 |Δpval|>1 或显著性单侧聚集即为实现回归。
- S4 嵌套槽不变量：对 S4 slot 内嵌套 list 的元素赋值（`x@misc$.var.features$features <- v`）不生效，必须整体写回（`vf <- x@misc$.var.features; vf$features <- v; x@misc$.var.features <- vf`）——与 relabelSpatialCellChat 的 idents 剥名同属一类陷阱，后续迁移脚本一律遵守。

## Alternatives considered

1. **与 SpatialChat_1.rds（Linux 服务器）对拍**：拒绝——跨平台 libm/BLAS 归约序差异不可移植，先验容差须放宽到 1e-8（2026-08-28 compare 脚本口径），无法区分"跨平台浮点"与"实现回归"两类信号；用户明确指示 _1 不参与。
2. **仅依赖 fixture 级测试，不做规模对拍**：拒绝——fixture 无法覆盖真实空间衰减核（exp/距离加权/AGAN）、2000×2000 稀疏模式和 1960 万 nnz 的累积效应；等价性契约宣称不完整。
3. **对拍 _2 加载后续跑**（而非重建 _3）：拒绝——_2 是旧类定义（无 misc 槽），链路函数已写新槽位；且"旧对象续跑"会把输入同源性问题与实现差异混在一起。同机重建 + 输入 bitwise 预检隔离了这两类信号。

## Consumer impact

- 测试资产：`tests_dev/simulated_data/run_s3_chain.R`、`compare_s3_s2.R`、`rigor_check.R`（新增，可复跑）；`baseline_inference_compare.R`（既有，作为 _2 的构建规格来源，未改动）。
- 产物：`tests_dev/simulated_data/SpatialChat_3.rds`（新链完整对象，validObject 通过）、`e2e_consistency_layers.csv`（逐层数值表）、`e2e_consistency_report.md`（报告）。
- 治理链接：implemented note 2026-10-08-netpathway-aggregate-migration 的 Equivalence Contract 由本 note 的 Evidence 补充规模级验证；docs/COMPUTATION_OPTIMIZATION.md §9.4 Wave 2b 行的 fixture 级证据声明保持有效（本 note 不改动该文档，规模级证据以本 note 为准）。
- 下游消费者：Plan B（置换检验统计修复）实施时以 SpatialChat_3.rds 为"新语义前"快照，可复用 compare 脚本量化语义变化的显式差值。
- 不影响任何运行时代码：本决策只新增测试资产，零生产代码变更。

## Consequences

- 维护成本：_2 由 renv 内冻结基线包（be88f30, Version 0.1.0）产生，重生成依赖 `baseline_inference_compare.R` + 该 pin；基线包升级时 _2 须重新生成并重跑对拍。
- 磁盘成本：SpatialChat_3.rds ~133 MB 量级，与 _2 同级；checkpoint 已清理。
- 复跑成本：全链 ~18 s（computeCommunProb 2.7 s + filterProbability 14.0 s + computeAvgCommunProb 0.7 s + 其余 <1 s），对拍 ~15 s，rigor_check ~12 s；总 <1 分钟，可纳入验收流程。
- 未来义务：Wave 后续（netAnalysis/子集调用方、Visium、subsetCommunication）迁移完成后，应扩展 compare 脚本覆盖对应新槽位，并保持本 note 的判定标准；Plan B 落地时本 note 的 pval 语义不变量将被 successor note 显式取代（届时需交叉链接）。

## Evidence

全链运行（run_s3_chain.R，2026-10-09）：
- 输入预检 5/5 bitwise PASS：count 矩阵、normalizeData→data、metadata 组 levels、subsetData→data.signaling、DB$interaction、LRsig interaction_name。
- 链路耗时：computeCommunProb 2.9 s（35 层 nnz=19,585,047）→ filterProbability 14.0 s → filterCommunication(links) 0.1 s → computeAvgCommunProb 0.8 s → filterCommunication(min.cells) 0.4 s → computeCommunProbPathway 0.1 s → aggregateNet 0.2 s；final validObject PASS。

逐槽位对拍（compare_s3_s2.R + rigor_check.R）：

| 槽位（新 vs 基线） | bitwise 层数 | max_abs | max_rel | 判定 |
|---|---|---|---|---|
| net$cell$prob vs net$prob.cell（35 层） | **35/35** | 0 | 0 | 逐位一致 |
| net$cell$count / weight | 2/2 | 0 | 0 | 逐位一致 |
| net$group$count | 1/1 | 0 | — | 逐位一致 |
| net$group$LR.sig / netP$pathways / netP$pathways.cell | — | — | — | 集合与排序逐位一致 |
| net$group$prob（35 层） | 0/35 | 4.83e-12 | **1.31e-14（≈0.06 ulp）** | 1-ulp 契约内 |
| netP$group$prob（35 层） | 0/35 | 4.83e-12 | 1.31e-14 | 同上 |
| net$group$weight | 0/1 | 9.10e-12 | 5.16e-15 | 同上 |
| net$group$pval（35 层） | 29/35 整层逐位 | Δ=±0.01 | — | 6 条目单 boot 翻转，显著性翻转 0 |

rigor_check.R（质疑驱动补验，全部 PASS）：
- 支持集：8 槽位 n_only_ref = n_only_new = 0（含 cell.prob 与 group.prob 全 35 层断言）。
- 层对齐：`_2@net$tmp$prob.cell` 名单与顺序 == `_3@net$cell$prob` dimnames[[3]] == rownames(LRsig)；逐 k nnz 证明 array k 序 == tmp list 序。
- 总账：两侧存储条目 3,022,361 == 3,022,361（含显式零），非零数一致。
- 按名匹配重验：3,022,361 个非零条目 IEEE 精确相等，max|delta| = 0。
- IEEE 边界：bitwise=TRUE 排除 NaN（NaN==NaN 为 FALSE）；`!=0` 过滤排除 -0 表示差异。

pval 差异归因：6 层（cL1-cR5、cL2-cR1、dL8-dR7、eL1-eR8、bL7-bR5、bL1-bR8）各 1-4 条目 Δ=±0.01（=1/100 boot）；pval<0.05 翻转 ref-only=0、new-only=0——组级 Prob 的 ≤1 ulp 差在 `perm == obs` 边界翻转了个别置换比较，属预期浮点边界效应，非实现回归。

## Acceptance criteria

- [x] 输入同源性：run_s3_chain.R 输入预检 5/5 bitwise PASS（同机重建可复现）。
- [x] 细胞级 prob IEEE 逐位一致：35/35 层，且经名字锚定重验（rigor_check.R，max|delta|=0 over 3,022,361 条目）。
- [x] 组级 prob ≤1 ulp：max_rel = 1.31e-14 ≤ 1e-12（35/35 层 + netP + weight）。
- [x] pval：29/35 层逐位；6 条目 |Δ|=1 boot；显著性翻转 0（预注册判定：|Δpval|≤1 且翻转 <0.5% → 无显著差异）。
- [x] 离散结果逐位一致：group.count、cell.count/weight、LR.sig、pathways 成员与排序、稀疏支持集（n_only=0）。
- [x] 复现命令：`Rscript tests_dev/simulated_data/run_s3_chain.R` → `Rscript tests_dev/simulated_data/compare_s3_s2.R` → `Rscript tests_dev/simulated_data/rigor_check.R`，三者均以 "OK/DONE" 结束、零 FAIL。
- [x] 治理：本 note 由 proposed 转 implemented（验证命令与证据齐备），validator 通过。
