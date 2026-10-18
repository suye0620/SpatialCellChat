# MBM05 旧代码(串行) vs 新代码(串行) 逐层张量对比报告

日期:2026-10-16
决策范围: investigate the group-level (`computeAvgCommunProb`) signature divergence between the pinned-upstream pipeline and the refactored pipeline under fully-serial execution; document cell-prob equality boundaries (pre/post `filterProbability`) with value- and support-level evidence.

## 0. 结论(Key Findings)

1. **net 的 cell 级 prob(computeCommunProb 之后、filterProbability 之后)在旧代码串行 vs 新代码串行下逐位一致(bit-identical)**:支持集(support)与数值(values)零差异(common = 全部条目,only/both 差 0,max |Δ| = 0)。filterProbability 前后两个断点都成立。
2. **唯一差异锁定在 `computeAvgCommunProb` 群组平均层**:
   - 支撑条目差 111(旧侧独有 148、新侧独有 37);
   - 群组 prob 的 stored-sum 差 0.00351,**精确等于**两侧独有条目各自 sum 的差(0.007525 − 0.004016);
   - 共同支持的 10,204 个条目数值差 ≤ 3.6e-18(相对差 ≤ 6.1e-14,~1 ulp),7,901 个条目有该量级的尾数差,其余 2,303 个逐位相同;
   - p 值全空间差异 1,033(含 0↔1 完全翻转 80/37),完整翻转集中到 **17 个 LR 层**;显著性(<0.05)条目 5,971 → 5,921(净 −50;旧→不显著 90、不显著→显著 40)。
3. **因果归属:全部差异由重构策略(Plan A kernel 的 `computeAvgCommunProb`)导致**,与 `computeCommunProb`/`filterProbability`/`computeCommunProbPathway` 无关——cell 层两张 prob 张量在 filter 前后都逐位相同是铁证。
4. 数值后果在**生物量级上可忽略**(最大的旧侧独有条目值 5.5e-4,与典型强边 1e-3 同量级但数量少;数值差 ~1e-18);但在**算法契约文本上**,Plan A 的"离散结果逐位一致 + Pval 零翻转"声明(在 30 层 fixture 上验证)在全量 1,563 层规模被推翻——支撑条目离散差异与 p 值完全翻转都超出了该契约的已验证边界。

## 1. 实验设置

| 项 | 值 |
|---|---|
| 输入 | 同一 `tests_dev/test_data/MBM05Chat_pre_inference.rds` checkpoint + 同一 components |
| 旧代码 | pinned upstream package (SHA `be88f30...`),原样串行调用(future::sequential,无 shim) |
| 新代码 | refactor 后的 R/ 源码 + Rcpp kernel(`Rcpp::sourceCpp`),全部串行 `nthreads=1` |
| 执行方式 | 各自独立 Rscript 进程,隐藏控制台,同机同时段并行跑 |
| 参数 | benchmark 共同参数(Kh=0.5, n=1, distance.use, interaction.range=250, contact.range=10 等);filter nboot=100/seed 666;avg nboot=100/seed 1 |
| 张量 dump | 各阶段边界 `saveRDS(compress=FALSE)`:`{old,new}_prefilter_cell_prob`、`{old,new}_postfilter_cell_prob`、`{old,new}_postfilter_object`、`{old,new}_group_prob`、`{old,new}_group_pval` |

Serial 端执行脚本(含 dump 逻辑):`tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R`(--worker=new|baseline --baseline-workers=1 --dump-dir=...);对比脚本:`tests_dev/benchmarks/compare_serial_tensors.R`;分析脚本:`tests_dev/benchmarks/analyze_group_diff{1..6}.R`。

### 串行阶段耗时(serial vs serial)

| stage | old serial | new serial | 比值 |
|---|---:|---:|---:|
| create_normalize | 2.29s | 4.14s | — |
| computeCommunProb | 18,489s | 318s | 58× |
| filterProbability | 4,365s | 86s | 50× |
| computeAvgCommunProb | 2,519s | 73s | 34× |
| computeCommunProbPathway | 31s | 2.9s | 11× |
| aggregateNet | 38s | 4.8s | 8× |
| centrality_net | 664s | 785s | 0.85× |
| centrality_netP | 94s | 111s | 0.85× |
| **全程** | **26,204s (7.3h)** | **1,385s (23min)** | **18.9×** |

(附:Plan A Note 宣称的 19× 加速比,全量串行实测 18.9×,吻合。)

## 2. net 的 cell prob:数值与支撑逐位一致(含 filter 前后)

比较方式:两侧张量按统一 (i,j,k,x) 坐标三元组 → flat 索引 → 排序对齐 → 支持集交并 + 数值差。

### computeCommunProb 之后(即 filterProbability 之前)

| 指标 | old(串行) | new(串行) |
|---|---:|---:|
| stored entries | 407,065,015 | 407,065,015 |
| nonzero 支持 | 407,065,015 | 407,065,015 |
| stored sum | 1.65372033912306e4 | 1.65372033912306e4 |
| 共同支持条目 | 407,065,015 | — |
| 仅 old / 仅 new | 0 / 0 | — |
| 共同条目 max \|Δx\| | **0** | — |

→ **零位置逐位一致,数值逐位一致。**

### filterProbability 之后

| 指标 | old(串行) | new(串行) |
|---|---:|---:|
| stored entries | 285,801,026 | 285,801,026 |
| nonzero 支持 | 285,801,026 | 285,801,026 |
| stored sum | 1.40318482389692e4 | 1.40318482389692e4 |
| 共同支持条目 | 285,801,026 | — |
| 仅 old / 仅 new | 0 / 0 | — |
| 共同条目 max \|Δx\| | **0** | — |

→ filter 侧同样逐位一致(286M 条目,保存为 `*_postfilter_cell_prob.rds`,连同 filter 后完整对象 `*_postfilter_object.rds` 供日后单阶段重放)。

**推论**:compute 与 filter 两个阶段,重构将上游实现的输出**精确复现**(替换实现也复现了行为),支持/数值无一差;下游任何差异必然来自 averaging(Plan A kernel)。同时(附证)旧串行 vs 旧多进程(此前 stage_workers 运行,含 Rscript shim)签名逐位一致,cross-execute 的 shim 无数值影响。

## 3. group 层(computeAvgCommunProb)逐连边差异

### 3.1 支撑(连边)位置差异

| 项 | 值 |
|---|---:|
| total coords | 8×8×1,563 = 100,032 |
| 旧 nonzero prob 条目 | 10,352 |
| 新 nonzero prob 条目 | 10,241 |
| 共同支持 | 10,204 |
| 仅旧 | 148 |
| 仅新 | 37 |

148 个"仅旧"连边按值分布:x 范围 [2.2e-14, 5.45e-4];最大几条(全部 top-12):

| 层(LR) | (i,j) | x_old | pval_old |
|---|---|---:|---:|
| HLA-E_KLRC2 | 6→6 | 5.452e-4 | (pval==0,显著) |
| HLA-E_KLRC2 | 8→8 | 5.427e-4 | (pval==0) |
| TENM4_FLRT3_ADGRL1 | 8→8 | 4.406e-4 | 0.41 |
| HLA-E_KLRC2 | 8→6 | 3.531e-4 | (pval==0) |
| HLA-E_KLRC2 | 6→8 | 2.282e-4 | (pval==0) |
| INHABA_ACVR2A | 3→3 | 1.658e-4 | … |
| IL18_IL18R1_IL18RAP | 3→3 | 1.525e-4 | … |
| IL18_IL18R1_IL18RAP | 6→6 | 1.423e-4 | … |
| IL21_IL21R_IL2RG | 3→3 | 1.344e-4 | … |

37 个"仅新"连边:x 范围最大 7.64e-4(CD40LG_CD40 5→5),其余同量级或更小;pval<0.05 的 0 条。

**差异大不大**:只有弱连边在 1e-4 量级以下动;整体 sum 差 -0.0035(相对 0.13%),完全由支撑差本身解释(sum[only_old]=0.007525 − sum[only_new]=0.004016 = 0.003509 = 签名 sum 差,**精确闭合**)。

### 3.2 共同支持的数值差

| 项 | 值 |
|---|---|
| 共同条目 | 10,204 |
| 其中 \|Δx\|>0 | 7,901 |
| max \|Δx\| | 3.60e-18 |
| max 相对差 | 6.06e-14(最大值边上 ~2.4e-16) |
| 例 | x_old=1.891208e-3,x_new 同,dx=-8.7e-19 (rel 4.6e-16) |

→ 共同支持内的差异即 Plan A 已文档化的 **num 分块求和序 / 1-ulp 契约**,在此全量规模被复证(数值上无意义)。

### 3.3 p 值(pval)差异(全空间 100,032 个坐标)

存储语义:old `net$pval` 为 dense 8×8×1,563 数组(0 即显著);new 的 `net$group$pval` 为逐层 SparseChatArray,absent 条目 = 真实 p 值 0(存储约定掉显式零)。两侧全空间读取后比较:

| 项 | 值 |
|---|---:|
| pval 不同条目 | 1,033 |
| 完全翻转 0→1(old 最显著 → new 最不显著) | **80** |
| 完全翻转 1→0 | **37** |
| 0→(0,1) / 1→(0,1) | 1 / 0 |
| 双侧 valid 的量化位移(±k×0.01) | 848(最大 \|Δ\|≤0.5) |
| 显著性(<0.05 含 0)全空间 | old 5,971 → new 5,921(净 −50;旧→不显著 90,不显著→显著 40) |

**完整 0↔1 翻转集中于 17 个 LR 层**(多处同一层多组多边同时全翻转):IL18_IL18R1_IL18RAP、INHABA_ACVR2A、WNT2B_FZD9_LRP5、WNT6_FZD9_LRP5、IL21_IL21R_IL2RG、VTN_ITGA8_ITGB1、VTN_ITGAV_ITGB3、TNFSF13_TNFRSF17、TNFSF13_TNFRSF13B、VTN_ITGA2B_ITGB3、MIF_CD74_CXCR2、GDNF_GFRA1、TSLP_IL7R_CRLF2、**HLA-E_KLRC2**、PCDHGB6_PCDHGB6、IL24_IL20RA_IL20RB、IL12B_IL12RB1_IL12RB2。

**机制复核（基于组件级门控重算）**：支撑差异层与完整翻转层高度重叠，但不能据此直接归因为 `>=`/`>` 的置换平局。pinned upstream 的平均函数使用 `format(dataLR_temp[, -1], digits = 1) >= min.percent`；`format()` 返回字符型，当前 R 环境下会产生字符比较。新 kernel 使用数值 `signif(mean, 1) >= min.percent` 并预计算固定整数阈值。
以 MBM05 组件的表达支持和 group_edge_diff 逐条回算：36 个支撑差异 LR 层中，148/148 条 only-old 都是 `old gate=1,new gate=0`，37/37 条 only-new 都是 `old gate=0,new gate=1`。因此支撑差异的主因是门控语义不一致，而不是下溢。
门控差异还会传播到 `Pval[Prob == 0] <- 1`，解释 80 条 `0→1` 与 37 条 `1→0` 完全翻转。两侧 p 值比较符均为 `Pboot - Pnull > 0`；不存在已证实的“一侧 `>=`、另一侧 `>`”差异。
共同支撑上的另一类 p 值位移仍由求和顺序解释：Matrix/BLAS 两段归约与 C++ CSC 单遍累加可产生约 1 ulp 尾差，严格 `> 0` 比较在平局边界上可改变一个 boot 的 `0.01` 计数。完整复核见 [COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md](COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md)。

## 4. 机制判定与差异归因

1. 差异源头仍是重构后的 `computeAvgCommunProb`，而非 `computeCommunProb`、`filterProbability`、pathway 或 shim；cell-level pre/post-filter tensor 已逐位一致。
2. 共同支撑约 `1e-18` 的尾数差属于已文档化的求和顺序边界；148/37 条离散支撑差异对应正式采用的 numeric gate 与 pinned upstream 历史字符 gate 的有意行为边界。
3. p 值差异是两条链路叠加：门控变化导致观察边增删和 `Prob==0` 强制 p=1；共同边上的数值归约差异可在严格比较边界改变部分置换计数。
4. 维护者已裁定采用 `signif(mean, 1) >= min.percent` numeric semantics；Plan A 不再要求新实现复现 pinned upstream 的 `format()` 字符比较结果。
5. 生物量级评估保持不变：差异集中在弱边，group sum 相对差约 0.13%；后续验收以 numeric reference 和 `Prob ≤1e-12` 数值界为准。

## 5. 复现命令

```bash
# 串行重跑(含 dump)
"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R \
  --worker=new      --out=tests_dev/benchmarks/out/mbm05_serial_tensor_comparison \
  --baseline-workers=1 --dump-dir=tests_dev/benchmarks/out/mbm05_serial_tensor_comparison

# tensor 对比与 group 连边分析
"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" tests_dev/benchmarks/compare_serial_tensors.R tests_dev/benchmarks/out/mbm05_serial_tensor_comparison
"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" tests_dev/benchmarks/analyze_group_diff2..6.R
```

## 6. Artifacts 索引

| 文件 | 内容 |
|---|---|
| `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/{old,new}_*.rds` | 串行各阶段张量与对象 dump |
| `.../{prefilter,postfilter}_cell_comparison.csv` | cell 级支持+数值差异统计 |
| `.../group_summary.csv` / `.../group_edge_diff.csv` | group 级聚合 + 8,086 行差异连边详表 |
| `.agents/notes/implemented/2026-10-15-mbm05-group-avg-signature-divergence.md` | 该发现的治理 Note（implemented；numeric gate 已裁定） |
