# Agent Note: computeAvgCommunProb 迁移至 11-slot 与置换检验 Plan A 加速

Status: implemented
Governance: v1
Date: 2026-09-28
Decision type: architecture + performance
Scope: `R/modeling.R` `computeAvgCommunProb`、`src/SpatialChat_Rcpp.cpp`（`cpp_group_avg_obs`/`cpp_group_avg_perm`）、`tests_dev/{test-computeAvgCommunProb.R, test-mbm05-avgperm.R, reference-computeAvgCommunProb-baseline.R}`
Owner: SpatialCellChat maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

1. `computeAvgCommunProb` 仍读写旧 14-slot 结构（`net$tmp$prob.cell`/`Lavg`/`Ravg`、`@options$parameter`），在已迁移的 11-slot 对象上入口即断（Wave 1 后 `net$tmp` 已删除）。
2. 置换段为推断管线剩余最大瓶颈：每 LR 每 boot 全套组平均（model.matrix + 2×aggregate + 4 次稀疏矩阵乘/dgemm）≈14.5 ms；500 有信号 LR × nboot=100 外推 ~25 min（COMPUTATION_OPTIMIZATION.md §9.4）。

## Proposal（Plan A，维护者 2026-09-27/28 裁定）

- **统计语义完全保持基线**：同 seed 同 permutation 生成调用（`set.seed` 位置逐字对齐 + 单次 `replicate(nboot, sample.int(nC, nC))`）、单侧计数 p = nReject/nboot、门控每个 boot 在置换标签上重算、`Pval[Prob==0] <- 1`。用户已确认：缺陷修复（(b+1)/(N+1)、空间约束 null、FDR 等）后置至 Plan B，与合作者确认后另行裁定（见 docs/PERMUTATION_TEST_AUDIT.md）。
- **kernel 化**：`cpp_group_avg_obs`（观测标签）/`cpp_group_avg_perm`（nboot 个置换一次调用），单遍 O(nnz·nboot) 分桶聚合取代 4 次矩阵乘 + dgemm；boots 维 OpenMP 并行（`nthreads`，输出列互斥，逐位与单线程一致）。
- **min.percent 门控整数阈值化**：决策 `signif(mean, 1) >= min.percent` 对 0/1 计数 cnt 单调非降，且组大小在标签置换下不变（分母 n_a 为每组建模常量）→ 每组一个整数阈值，R 侧用**基线同款 `mean()`+`signif()`** 二分探测求出（每组一次、~log₂n 次调用），kernel 内退化为 `cnt >= thr` 整数比较。避免移植 R `fprec`，门控逐位一致；单调性对 n=1..80 × τ∈{0.05,0.1,0.2,0.5} 穷举验证 + n=1000/29536 分层探针（test §2）。
- **supp 来源**：旧 `net$tmp$Lavg/Ravg` 仅经 `1*(dataLR>0)` 进门控；缓存 Ravg 含 cofactor 调整（`Rexpr*coA/coI`），但 coA/coI = 1+expr ≥ 1 恒正 → `>0` 模式为 cofactor 不变量 → 直接取 assay signaling 行 `>0`（复合物 = 全体亚基 `>0`，`geometricMean` 语义），门控输入与基线逐位一致。
- **tmp 层删除**（用户 2026-09-27 确认不再需要）：`LRsig.CCC/GGC.counts`、`LRsig.use.idx` 移入 `misc$.param$averaging`（grep 确认无下游消费者）。
- **输出落点**（DATA_STRUCTURE.md §层级关系）：`net$group$prob` / `net$group$pval`（SparseChatArray K×K×N）；参数 `misc$.param$averaging`；`.log_operation("computeAvgCommunProb")`。
- **新增参数**：`nthreads`（默认 1）、`verbose`；原签名参数序不变（`relabelSpatialCellChat` 调用兼容）。

### 等价性契约（关键裁定）

- **逐位一致**：RNG 流；全部门控决策（整数精确，任意累加序不变）；den = 块内存储条目数（含显式零，对应基线 `prob@x <- 1` 后 crossprod）；NaN→0；门控乘积；CCC/GGC 计数与两层活跃 LR 集合；Pval 数值（实测零翻转）。
- **顺序敏感唯一项：num 分块和**。基线为两段舍入（行内升序偏和 + `crossprod` 的 BLAS dgemm 归约序）；dgemm 归约序不可移植复现（且基线自身随 BLAS 线程数/指令集变化），kernel 采用文档化 CSC 遍历序（列升序、列内行升序）。实测 fixture 与 MBM05 30 层相对差 ≤3.6e-16（≈1 ulp），远低于基线自身运行间数值噪声；测试契约 ≤1e-12。
- 验收：`tests_dev/test-computeAvgCommunProb.R` 160 项全绿（kernel↔R 语义规范 `identical`、nthreads 一致、阈值穷举、冻结基线端到端对拍、拒绝路径、确定性）；`tests_dev/test-mbm05-avgperm.R` 端到端（数字见 Evidence）。

#### 组级 prob 的 1 ulp 契约（数学表述）

记号：每个 LR 层的细胞级概率为稀疏矩阵 $V \in \mathbb{R}^{n_C \times n_C}$（$v_{ij} \ge 0$），
分组 one-hot $D \in \{0,1\}^{n_C \times k}$，$D_{ia} = \mathbb{1}[g_i = a]$；
$\mathcal{J}_i$ = 第 $i$ 行的存储列索引（升序）。

**基线链路**（两段浮点求和，pre-migration `computeAvgCommunProb_LR_{Avg,Sum}`）：

$$
X = V D,\qquad X_{ib} \;=\; \sum_{\substack{j \in \mathcal{J}_i \cap b \\ \text{升序累加}}} v_{ij}
$$

$$
\mathrm{num}^{\mathrm{base}}_{ab} \;=\; (D^{\!\top} X)_{ab} \;=\; \sum_{i:\,g_i=a} X_{ib}
\quad\text{—— 由 BLAS dgemm 计算，累加序 = BLAS 内部序}
$$

分母（可精确复现）：

$$
\mathrm{den}_{ab} \;=\; \#\{\text{存储条目} (i,j) : g_i=a,\ g_j=b\}
$$

基线经 `prob@x <- 1` 后 `crossprod`，**显式零同样计数**；整数值、任意累加序不变 → 逐位一致。
`avg` 口径 $\mathrm{val}_{ab} = \mathrm{num}_{ab}/\mathrm{den}_{ab}$（$\mathrm{den}=0 \Rightarrow \mathrm{NaN} \Rightarrow 0$）；
`sum` 口径 $\mathrm{val}_{ab} = \mathrm{num}_{ab}$。最后乘 0/1 门控积 $G_{ab} = G^{\mathrm{percent}}_{ab} \cdot G^{\mathrm{sr}}_{ab}$。

**kernel 链路**（`sc_boot_block_avg`，CSC 单遍分桶，j 外层升序、列内 i 升序）：

$$
\mathrm{num}^{\mathrm{kernel}}_{ab} \;=\; \sum_{\substack{(i,j) \in \text{CSC 序} \\ g_i=a,\ g_j=b}} v_{ij}
\qquad \text{（同一实数多重集，不同求和顺序）}
$$

**等价类划分**：

1. **逐位一致**：permutation 矩阵、全部门控决策（整数阈值，任意累加序不变）、`den`、NaN→0、
   门控乘积、CCC/GGC 计数、两层活跃 LR 集合、**Pval 全部数值**（实测零翻转——kernel 内
   $T^{\mathrm{obs}}$ 与 $T^{(r)}$ 同舍入体系，离散比较自洽）。
2. **顺序敏感唯一项 `num`**：标准前向误差界
   $|\mathrm{fl}(\textstyle\sum_{t=1}^{m} v_t) - \sum_{t=1}^{m} v_t| \lesssim m\,\varepsilon\, \sum_t |v_t|$
   （$m$ = 块内条目数，$\varepsilon$ = 机器精度）。实测：fixture ≤3.55e-16、
   MBM05 30 层 ≤1.43e-14（相对差），测试契约断言 ≤1e-12。

**为什么 strict identical 不可达**：$\mathrm{num}^{\mathrm{base}}$ 的第二段归约序是 OpenBLAS
构建版本 / 线程数 / FMA 指令集的函数——**基线自身换环境重跑末位即漂移**，"与基线 strict 一致"
只对某一次具体基线执行有定义。复现该序的唯一途径是调用同一 BLAS 算子（即保留基线 crossprod
链路），与加速目标互斥。因此本契约以"离散结果逐位一致 + Prob ≤1 ulp"为**最终等价口径**，
与合作者沟通时按此表述。

## Constraints and invariants

- 相同 seed + 输入 → 输出逐位一致（R 层面重跑确定性；Prob 数值受 1 ulp 契约约束）。
- `computeAvgCommunProb_LR_{Avg,Sum}` 原样保留导出（Visium 变体与冻结基线对拍依赖）。
- 空层（nnz=0）被 CCC 计数剔除；GGC=0 层 pval 恒 1（基线语义）。
- `nboot < 1` 且 `do.permutation=TRUE` 显式报错（基线为隐晦崩溃；"非法输入明确报错"先例）。

## Alternatives considered

1. 移植 R `fprec` 实现 kernel 内 `signif`：拒绝——整数阈值化 + R 侧探针达到相同逐位保证，且无舍入逻辑移植风险。
2. bitwise 匹配 dgemm 归约序：不可行（BLAS 内核序随线程/FMA 配置变化，基线自身不稳定）。
3. 置换检验统计修正（空间约束 null、(b+1)/(N+1)、BH-FDR、门控条件化）：维护者裁定后置 Plan B，接口已在 PERMUTATION_TEST_AUDIT.md §5.3 预留。

## Consumer impact

- `R/modeling.R`：`computeAvgCommunProb` 重写；`LR_Avg`/`LR_Sum` 未动。
- 下游 `filterCommunication`/`computeCommunProbPathway`/`aggregateNet`/`netAnalysis_computeCentrality` 等仍读旧 `net$prob`/`net$pval` 槽——属下一波迁移（note 2026-09-10 Consumer impact 链）；组级结果新落点 `net$group$prob`/`net$group$pval`。
- `man/computeAvgCommunProb.Rd` 待 roxygen 重生成。

## Consequences

- MBM05 全量（1563 LR × nboot=100）：全量 125.66 s；基线外推 2372.6 s；加速比 19×（单线程；kernel 段 nthreads=8 再 ~2.5×）。
- 组级结果槽位变化：旧 `net$prob` → `net$group$prob`（SparseChatArray）；下游迁移波必须跟进，否则读旧槽得 NULL。

## Evidence

- `tests_dev/test-computeAvgCommunProb.R`：160 PASS（2026-09-28，R 4.5.3 / g++ 14.3.0 MinGW-w64）。
- 冻结基线对拍（be88f30 快照 plain-input 化，`computeAvgCommunProb_LR_Avg` 原函数驱动）：Prob 相对差 ≤3.55e-16；Pval 数值零差异（50 boots）。
- `tests_dev/test-mbm05-avgperm.R`（2026-09-28 实测）：computeCommunProb 400.0 s（4.07e8 nnz）→ filterProbability → computeAvgCommunProb 全量 125.66 s，validObject 通过；30 层冻结基线对拍 max rel = 1.43e-14、Pval 零翻转；基线 1518 ms/层；kernel 线程扩展 1/4/8 线程 = 50/30/20 ms（nnz=5e4 × nboot=100 × K=8 合成层）且逐位一致。
- `docs/PERMUTATION_TEST_AUDIT.md`（统计语义审计；Plan A 明确不修其中问题）。
- 修复记录：kernel 初版把 `perm` 条目误作组码直接索引（基线语义为 gather `group[permutation[,e]]`），nC>K 时越界写引发漂移性原生崩溃——由压力测试（60 随机 trial）暴露，已修复并加 perm 值域校验。

## Acceptance criteria

- [x] 11-slot 对象上全流程可运行（assay 读取、`net$group$prob`/`pval` 写入、params/log/validObject）
- [x] 门控、计数、活跃层集合、p 值与基线逐位一致
- [x] Prob ≤1 ulp 契约成立并有测试护栏（≤1e-12 断言）
- [x] 160 项测试全绿；kernel 与 R 语义规范 `identical`；nthreads 无关性
- [x] MBM05 全量耗时与加速比记录（125.66 s vs 外推 2372.6 s，19×；对拍与 validObject 全过）
