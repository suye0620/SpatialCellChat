# Agent Note: preProcessing ALRA 填补适配最终 11-slot schema

Status: implemented
Governance: v1
Date: 2026-09-10
Decision type: architecture
Scope: `R/utilities.R` `preProcessing()`、`R/SpatialCellChat_class.R` assay 校验、`man/preProcessing.Rd`、`tests_dev/test-preProcessing.R`
Owner: SpatialCellChat maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`preProcessing()` 仍读写旧 14-slot 结构（`object@data.signaling`、`object@data`、`object@options`），在最终 11-slot 类（`assay`/`misc` 等）上必然报错，函数对当前对象完全不可用。此外基线实现存在多处与上游 ALRA 契约不符的缺陷：

- `k <- median(k.boot)` 对 10 个估计值取中位数可产生非整数 k，传入 `ALRA::alra()` 会导致列索引出错；
- `noise_start = min(round(0.8*K), 80)` 违反 ALRA `choose_k` 的 `noise_start <= K - 5` 硬约束，小矩阵必然 `stop`；
- `choose_k` 无奇异值间隙时返回 `k = -Inf`、k = 1 时 ALRA 内部 `diag()` 退化，均无守卫；
- rsvd 随机性未固定 seed，k 估计与填补结果不可复现。

同时，`.sc_validate_assay` 要求包括 `assay$signaling` 在内的所有非空层具有全基因维度，与 `subsetData()` 将 signaling 写成 DB 基因子集的语义直接矛盾：任何真实数据集（DB 信号基因 ⊂ 对象全基因）在 `subsetData → preProcessing` 标准流程中都会被 validator 拒绝。既有测试 fixture 恰好满足"DB 基因 = 全基因"，掩盖了该缺陷。

## Decision

`preProcessing()` 以 assay 层为目标重写，作为 ALRA 填补的稳定契约：

- 对象输入：`slot.name = c("signaling", "norm")`（默认 `signaling`），经 `assay()` 读取、`assay<-` 写回；矩阵输入返回 genes × cells 的 `dgCMatrix`。
- 参数显式化：`quantile.prob`（**ALRA 原生参数**，非本包定义；见下方语义澄清）、`seed.use = 1L`（固定 rsvd 随机性，k 估计与填补均可复现）。
- 秩估计：`K` 在 `min.dim < 100` 时取 `floor(0.9 * min.dim)`，否则沿用 ALRA 默认 100；`noise_start = max(1, min(round(0.8*K), 80, K-5))` 显式满足 ALRA 约束；k 取 10 次 `choose_k` 中位数并四舍五入为整数；守卫 `min.dim >= 7`、`2 <= k < min(dim)`。
- 结果取 ALRA 返回值的**命名**元素 `$A_norm_rank_k_cor_sc`，写出前恢复 `dimnames`。
- 对象写入路径：impute `norm` 时清除派生层 `scale`/`smooth`/`signaling`；`k` 记录到 `misc$.param$alra`（整数），操作写入 `misc$.log`。
- 输入语义：ALRA 要求 log-normalized 输入；对非负整数值输入发出 warning 提示先运行 `normalizeData()`。
- `.sc_validate_assay` 中 `assay$signaling` 改为基因子集语义：列必须精确等于对象细胞名，行必须是 `assay$norm` 基因顺序保持的子序列；其余层（`raw`/`norm`/`scale`/`smooth`）维持全维度校验。

### quantile.prob 语义澄清（ALRA 原生参数）

`quantile.prob` **是 ALRA 包 `alra()` 的原始参数**（KlugerLab/ALRA master：`alra(A_norm, k=0, q=10, quantile.prob=0.001, ...)`），不是 SpatialCellChat 定义的概念。ALRA 内部对低秩近似矩阵 `A_norm_rank_k` 的每个基因（列）计算 `quantile(x, quantile.prob)`，取绝对值作为该基因的自适应截断阈值，将低于阈值的元素置零，再对非零值做一、二阶矩匹配（均值/方差对齐到原始输入）。因此它控制**逐基因截断的激进程度**：分位越大 → 阈值越高 → 置零越多 → 输出越接近原始稀疏度；分位越小 → 保留越多低幅填补值。

ALRA 默认 `0.001`；本包基线（沿用自前身 SpatialChat）保持 `1e-5`——更宽松的阈值保留更多低幅填补值，输出更稠密。维护者确认保留 `1e-5` 默认（2026-09-10）：该超参数可能经过先前调参，无充分证据支持更改；已作为显式参数暴露，需要对齐 ALRA 默认的调用方可显式传 `quantile.prob = 0.001`。

## Constraints and invariants

- ALRA 输入必须是 log-normalized、dense、cells × genes 矩阵；本函数负责转置与 `as.matrix()` 具体化。
- `noise_start <= K - 5` 且 `K < min(dim)` 必须始终成立；`k` 必须是 `[2, min(dim))` 内的整数。
- 填补前后基因名与细胞名完全一致；输出非负（log 尺度）。
- 原始表达（非零）位点在填补后必须仍为非零（ALRA 恢复原始值保证）。
- `assay$signaling` 行序必须保持 `assay$norm` 基因顺序；任何写入路径经 `assay<-` 的 `validObject` 校验。
- 相同 `seed.use` 下重复运行产生逐元素一致的填补矩阵。

## Alternatives considered

1. 保留旧 slot 名（`data.signaling`/`data`）作兼容别名：违反 schema 治理决策"旧 slot 迁移后删除，不叠加兼容分支"，拒绝。
2. 将默认 `quantile.prob` 改为 ALRA 默认 0.001：改变填补语义且基线行为（1e-5）无证据表明错误；保留 1e-5 并参数化暴露，拒绝静默更改。
3. 在 validator 中为 signaling 豁免校验（返回空）：放弃真实不变量（细胞对齐、行序子序列），使 stale/错序层静默通过，拒绝。

## Consumer impact

- `R/utilities.R` `preProcessing()`：重写（唯一实现位置）。
- `R/SpatialCellChat_class.R` `.sc_validate_assay()`：signaling 子集语义。
- `man/preProcessing.Rd`：更新签名、参数与返回值文档。
- `tests_dev/test-preProcessing.R`：新增回归测试（19 项检查）。
- 全仓 grep 确认 `preProcessing()` 无内部 R 调用方；`man/preProcessing.Rd` 是唯一文档引用。
- 下游消费者：`identifyOverExpressedGenes()`、`identifyOverExpressedInteractions()`、`computeCommunProb()` 读取填补后的 `assay$signaling`，接口不变。
- 旧 `tests_dev` 调用（如有）引用 `slot.name = "data.signaling"` 的写法将报错——属 schema 治理预期的破坏性迁移。

## Consequences

- `preProcessing()` 在新 schema 上恢复可用且可复现；`subsetData → preProcessing` 标准流程在真实规模数据上验证通过。
- 旧 slot 用法（`slot.name = "data.signaling"`）不再接受——这是最终 schema 迁移的一部分，无兼容层。
- 全基因矩阵输入路径（`preProcessing(matrix)`，29,536 × 27,854 ≈ 6.6 GB dense）内存开销大；对象路径（subsetData 后，967 信号基因）dense 仅 ~230 MB。大矩阵成本是 ALRA dense 契约的固有限制，函数打印维度提示。
- `assay$signaling` 校验语义的收紧可能暴露其他写入路径的行序错误（当前未发现）。
- 未来若 ALRA 上游更改返回结构，命名访问 `$A_norm_rank_k_cor_sc` 比位置索引更早、更清晰地失败。

## Evidence

- ALRA 上游源码（KlugerLab/ALRA master，renv 锁定 SHA `f34d465`）：`choose_k(A_norm, K=100, thresh=6, noise_start=80, q=2)` 约束 `K > min(dim)` 即 stop、`noise_start > K-5` 即 stop；`alra(A_norm, k=0, q=10, quantile.prob=0.001)` 要求 dense matrix、返回命名 list。
- 本地安装探测（R 4.5.3，renv library）：`formals(ALRA::choose_k)`、`formals(ALRA::alra)` 与上游一致；`as.array(dgCMatrix)` 可行但语义含糊。
- `tests_dev/test-preProcessing.R`：19 项检查全部通过（层写入、维度/名保持、非负、原位保真、支持集不缩、k 整数记录、日志、validObject、同 seed 确定性、norm 路径清层、矩阵路径一致性、4 个拒绝路径）。
- 端到端验证（Slide-seq MBM05 rep1，27,854 基因 × 29,536 细胞）：`createSpatialCellChat(input.assay="raw") → subsetData → preProcessing(seed.use=1)` 全流程通过；ALRA 秩 k=49（2000 细胞子样本 k=47）；ALRA 本体 3.15 min；基因级 nonzero 均值相关（≥20 非零细胞基因，579/967）**0.9989**；非零占比 1.51% → 11.97%。
- 既有回归测试 `test-pre-inference-new-schema.R`、`test-SpatialCellChat-class.R`、`test-normalizeData.R` 在 validator 修改后全部通过。

## Acceptance criteria

- [x] `preProcessing()` 对象路径读写 `assay$signaling`/`assay$norm`，不再访问 `@data.signaling`/`@data`/`@options`。
- [x] `slot.name = c("signaling", "norm")`，非法值与空层（含 subsetData 提示）被明确拒绝。
- [x] `noise_start <= K - 5`、`K < min(dim)`、`k` 为 `[2, min(dim))` 整数的守卫全部生效（5×5 输入被拒绝）。
- [x] 相同 `seed.use` 重复运行产生逐元素一致的填补矩阵（测试断言通过）。
- [x] impute `norm` 后 `scale`/`smooth`/`signaling` 派生层被清除。
- [x] `misc$.param$alra$k` 为整数，`misc$.log` 含 `preProcessing` 条目，`validObject` 通过。
- [x] `.sc_validate_assay` 接受 `subsetData` 写出的基因子集 signaling 层，拒绝细胞不齐与行序错乱的层。
- [x] 全规模真实数据端到端流程通过，矩匹配相关系数 0.9989 佐证填补量级校准正确。
- [x] `man/preProcessing.Rd` 与实现签名一致；`tests_dev/test-preProcessing.R` 确定性回归通过。
