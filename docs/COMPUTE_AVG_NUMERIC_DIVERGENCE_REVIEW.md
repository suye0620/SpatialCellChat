# `computeAvgCommunProb` 全量差异归因与新函数设计审视

> 状态：证据审查（不改变置换检验方法，也不改变 Plan A 的加速决策）
> 范围：MBM05 全量串行对拍中 `computeAvgCommunProb` 的 group-level `prob`/`pval` 差异，以及重构函数的实现假设、兼容性边界和潜在缺陷。
> 关联：[SERIAL_COMPARISON_GROUP_AVG.md](SERIAL_COMPARISON_GROUP_AVG.md)、[PERMUTATION_TEST_AUDIT.md](PERMUTATION_TEST_AUDIT.md)、[COMPUTATION_OPTIMIZATION.md](COMPUTATION_OPTIMIZATION.md)、`.agents/notes/implemented/2026-09-28-compute-avg-commun-prob-plan-a.md`

## 0. 结论摘要

这次差异不是单一原因，而是两个层次叠加：

1. **离散支撑差异的主因已经定位为 `min.percent` 门控语义不一致。** pinned upstream 在 `computeAvgCommunProb_LR_Avg` 中使用 `format(..., digits = 1) >= min.percent`；在当前 R 环境中，`format()` 返回字符型，比较会按字符语义进行。新实现则使用数值型 `signif(mean, 1) >= min.percent`，并把它预计算为固定整数阈值。两者不是等价变换。
2. **共同支撑内的约 `1 ulp` 尾差来自聚合求和顺序不同。** upstream 通过两段 `Matrix::crossprod`/BLAS 归约；新 kernel 按 CSC 顺序直接向 group bucket 累加。数学结果相同，但 IEEE-754 加法不满足结合律。
3. **p 值完全翻转主要由门控差异传播而来；部分 `0.01` 量化位移还受到浮点平局边界影响。** 两边的 p 值比较符都是 `Pboot - Pnull > 0`，没有证据表明一边使用 `>=`、另一边使用 `>`。此前把主要原因表述为“`>=`/`>` tie 语义差异”不够准确，现以本文件的归因为准。
4. **没有证据显示 cell-level 输入、`filterProbability`、随机置换流、worker 调度或 OpenMP race 是差异源。** prefilter/postfilter cell-level `prob` 在旧串行与新串行之间均支持集和数值逐位一致。
5. **新函数的 CSC kernel、显式零计数、置换 gather 和线程隔离设计本身没有被这次证据推翻。** 需要修正的是兼容性边界和测试 oracle，而不是推翻 Plan A 或重写假设检验方法。

---

## 1. 证据边界与不变决策

### 1.1 对拍边界

输入是同一 MBM05 pre-inference checkpoint；旧侧为 pinned upstream package SHA `be88f300464cf970b36a6165ea7d2e41a47a7511`，新侧为当前 R/Rcpp 实现。两侧均使用串行执行；新侧 kernel `nthreads = 1`。参数包括：

- `filterProbability`: `nboot = 100`, `seed.use = 666`；
- `computeAvgCommunProb`: `avg.type = "avg"`, `min.percent = 0.1`, `min.cells.sr = 10`, `nboot = 100`, `seed.use = 1`；
- 同一 cell-level `prob`、同一 8 组、同一 1,563 个 LR 层。

### 1.2 已确认的全量数字

| 层级 | 旧侧 | 新侧 | 结论 |
|---|---:|---:|---|
| prefilter cell `prob` | 407,065,015 stored/nonzero | 同值 | 支持集、数值 `max \|Δ\| = 0` |
| postfilter cell `prob` | 285,801,026 stored/nonzero | 同值 | 支持集、数值 `max \|Δ\| = 0` |
| group `prob` nonzero | 10,352 | 10,241 | 仅旧 148、仅新 37 |
| group 共同支撑 | 10,204 | — | 7,901 条有非零尾差 |
| 共同支撑最大 `\|Δprob\|` | — | `3.6049722235143e-18` | 相对差最大 `6.05822e-14` |
| group stored sum | `2.64022099022303` | `2.63671247907007` | 差异 `0.003509` |
| p 值差异 | — | 1,033 / 100,032 坐标 | `0↔1` 完全翻转 80/37 |
| 显著性 `< 0.05` | 5,971 | 5,921 | 净差 `-50` |

group sum 差异精确闭合：

```text
only_old sum - only_new sum
= 0.007525 - 0.004016
= 0.003509
= old group sum - new group sum
```

这说明 group 总和差异主要由离散门控导致的边增删产生，共同边的 `1e-18` 尾差没有形成可观的总量放大。

### 1.3 不在本次范围内改变的决策

以下内容保持既有 Plan A/Plan B 边界，不因本次数值审查而改变：

- 置换标签的生成方式和随机种子位置；
- 单侧统计量 `nReject / nboot`；
- 每次 permutation 重新计算 group gate；
- `Pval[Prob == 0] <- 1` 的旧行为；
- 空间约束 null、`(b+1)/(N+1)`、FDR 等统计改进；
- CSC 单遍聚合和 OpenMP 按 boot 隔离输出的加速原理。

本文件审查的是**实现是否兑现既有兼容性契约**，不是重新审议假设检验原理。

---

## 2. 差异源头已经锁定在 group average

`computeCommunProb` 和 `filterProbability` 的两个边界 dump 在旧串行/新串行之间逐位一致：

- prefilter cell `prob`：共同支撑 407,065,015 条，`max |Δ| = 0`；
- postfilter cell `prob`：共同支撑 285,801,026 条，`max |Δ| = 0`。

因此：

```text
cell input V_old == V_new bitwise
```

差异只能从 `computeAvgCommunProb` 开始。后续 pathway 和 centrality 的既有签名也没有显示出新的 cell/pathway 输入差异。

这条边界很重要：不能把本次 group-level 分歧归因于 filter 分位数、空间距离、RNG 或多进程导出。

---

## 3. 主因：`format()` 字符门控与数值 `signif()` 门控不一致

### 3.1 两套函数实际做的事情不同

pinned upstream 的平均函数使用：

```r
dataLR_temp <- 1 * (dataLR > 0)
dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
dataLR_percent <- 1 * (
  format(dataLR_temp[, -1], digits = 1) >= min.percent
)
```

当前重构代码的门控阈值位于 `R/modeling.R:87-102`：

```r
decision <- function(cnt)
  isTRUE(signif(mean(rep(c(1, 0), c(cnt, n - cnt))), 1) >= min.percent)
```

然后传入 C++，在 `src/SpatialChat_Rcpp.cpp:257-263` 中执行：

```cpp
gL[a] = (cntL[a] >= (double)thr[a]) ? 1.0 : 0.0;
gR[a] = (cntR[a] >= (double)thr[a]) ? 1.0 : 0.0;
```

新实现的数学假设是：

```text
min.percent gate = 一个只依赖 group size 和 count 的固定整数阈值
```

这个假设对数值型 `signif()` 门控成立，但对 upstream 的 `format()` 门控不成立。

### 3.2 `format()` 还会受整列和 permutation 影响

`format(data.frame, digits = 1)` 是按列格式化的。它是否使用固定小数或科学计数法，取决于该 LR 层各组比例的整体范围。因此旧门控不是单纯的：

```text
f(count, group_size)
```

而更接近：

```text
f(count, group_size, 同一 LR 列其它组的比例, R 的格式化规则)
```

在 permutation 中，组大小虽然不变，但每个组的表达比例会变，`format()` 选择的表示形式也可能改变。新 kernel 把阈值一次计算后复用于所有 LR 和所有 boot，这在数值 `signif()` 语义下合理，在 upstream 的实际 `format()` 语义下不等价。

### 3.3 已逐条核对支撑差异方向

使用 MBM05 组件中的表达支持、8 组标签和 `group_edge_diff.csv` 重新计算门控后：

| 支撑差异 | 条数 | 门控核对结果 |
|---|---:|---|
| `only_old` | 148 | 148/148 都是 `old gate = 1, new gate = 0` |
| `only_new` | 37 | 37/37 都是 `old gate = 0, new gate = 1` |

共有 36 个 LR 层产生支撑差异。这个结果把 185 条离散支撑差异从“可能由数值边界造成”提升为“与门控语义逐条吻合”的证据。

### 3.4 代表性例子

#### HLA-E_KLRC2：旧侧误保留低比例受体

`KLRC2` 在 8 个组中的阳性细胞数：

```text
0, 0, 3, 0, 2, 12, 0, 22
```

所有比例都低于 10%。旧列被格式化为类似：

```text
0e+00, 0e+00, 4e-03, 0e+00, 2e-03, 2e-02, 0e+00, 8e-04
```

在当前 R 环境下，字符比较使旧侧 8 个受体 gate 全部为 1；新侧数值比较使 8 个 gate 全部为 0。

结果：

```text
HLA-E_KLRC2, 6 -> 6
old prob = 5.45189615325912e-04
new prob = 0
old pval = 0
new pval = 1
```

#### MIF_CD74_CXCR2：成片 old-only

`CD74_CXCR2` 的阳性细胞数：

```text
0, 0, 3, 0, 3, 2, 0, 9
```

同样由于科学计数法字符串比较，旧侧受体 gate 全部通过，新侧全部不通过。该层产生 28 条 `only_old` 边。

例如：

```text
MIF_CD74_CXCR2, 1 -> 3
old prob = 2.11442257718820e-05
new prob = 0
old pval = 0.29
new pval = 1
```

#### CD40LG_CD40：新侧保留 9.7% 边界值

`CD40LG` 在 Myeloid 组（第 5 组）中的比例为：

```text
96 / 991 = 0.09687
```

旧侧固定格式是类似 `"0.097"`，与 `"0.1"` 比较后不通过；新侧执行：

```r
signif(0.09687, 1) == 0.1
```

因此新侧通过该 gate。

结果：

```text
CD40LG_CD40, 5 -> 5
old prob = 0
new prob = 7.63876840186550e-04
old pval = 1
new pval = 0
```

#### TENM4_FLRT3_ADGRL1：复合物低表达边界

`TENM4_FLRT3` 几乎只在 Tumor 组出现，支持数为：

```text
0, 0, 0, 0, 0, 0, 0, 21
```

比例约为 `0.00079`。旧科学计数法路径会保留旧侧边；新数值 gate 会删除这些边。代表性结果：

```text
TENM4_FLRT3_ADGRL1, 8 -> 8
old prob = 4.40573353578737e-04
new prob = 0
```

### 3.5 对函数设计的判断

这不是“C++ 把一个正数算成了零”的下溢问题，也不是分母或 `min.cells.sr` 的计数误差。它是：

- **如果 pinned upstream 行为是兼容 oracle：** 新函数存在已确认的兼容性缺陷；
- **如果 `signif()` 是有意选择的数值语义：** 新函数内部设计是自洽的，但这是一次未显式记录的行为变化，不能继续宣称“与 pinned upstream 离散结果逐位一致”。

在维护者明确选择前，不应默默把该差异降级为浮点噪声，也不应直接改动假设检验方法。

---

## 4. 次因：CSC 聚合与 BLAS 聚合的求和顺序不同

### 4.1 upstream 的两段归约

旧路径本质上是：

```r
X <- prob %*% cell.type.mat
num <- Matrix::crossprod(cell.type.mat, X)
```

第一段先按目标组形成每个细胞的部分和，第二段再按发送组归约。第二段归约顺序由 Matrix/BLAS 实现决定。

### 4.2 新 kernel 的单遍归约

`sc_boot_block_avg` 在 `src/SpatialChat_Rcpp.cpp:275-280` 中按 CSC 列顺序直接做：

```cpp
num[f] += x[t];
den[f] += 1.0;
```

两边使用的是同一批 cell-level 权重，但加法顺序不同。共同支撑的全量结果：

- 10,204 条共同边；
- 7,901 条 `prob` 尾数不同；
- 最大绝对差 `3.6049722235143e-18`；
- 最大相对差 `6.05822e-14`；
- 没有共同边超过 `1e-12`。

代表性 p 值变化：

```text
TGFB1_TGFBR1_TGFBR2, 1 -> 3
prob old = 4.42935009099750e-05
prob new = 4.42935009099750e-05
|Δprob|  = 6.77626357803440e-21
pval      0.08 -> 0.09
```

该 LR 层不属于 36 个支撑差异层，说明至少这一类 p 值变化可以在 gate 相同的情况下由浮点求和顺序触发。

### 4.3 对函数设计的判断

这部分符合已有 Plan A 设计判断：

- 精确复现 BLAS 的第二段归约顺序会牺牲 kernel 化收益，并且结果还会随 BLAS 线程、FMA 和构建环境变化；
- CSC 单遍聚合的数值差异处于已测 `1 ulp` 级别；
- 新线程之间处理不同 boot，输出列互斥，已有 `nthreads=1/4` 逐位一致测试。

因此目前没有理由因为这部分尾差而回退 Plan A。需要把 p 值的“严格零差”承诺与浮点归约边界分开描述，而不是把所有差异归结为同一个 tie 机制。

---

## 5. p 值差异的传播路径

当前两侧都使用：

```r
nReject <- rowSums(Pboot - Pnull > 0)
p <- nReject / nboot
```

新代码位置为 `R/modeling.R:1263-1269`；pinned upstream 也使用相同的严格 `> 0` 比较。

### 5.1 门控引起的完全翻转

如果旧侧观察边通过 gate、新侧不通过：

```text
old Prob > 0
new Prob = 0
```

随后新侧执行：

```r
Pval[Prob == 0] <- 1
```

于是旧侧可能是 `p = 0`，新侧必为 `p = 1`。反方向同理，解释了 `80` 个 `0 -> 1` 与 `37` 个 `1 -> 0` 完全翻转。

17 个发生完全翻转的 LR 层全部属于支撑差异层集合，因此门控传播是这些翻转的直接解释。

### 5.2 浮点边界引起的量化位移

在 gate 相同、观察边和 null 边都保留的情况下，`Pboot` 或 `Pnull` 的一 ulp 差异可能改变某一个 boot 的严格比较结果。`nboot = 100` 时，每一次计数改变对应 `0.01` 的 p 值位移。

这解释共同支撑上的诸如：

```text
0.08 -> 0.09
0.15 -> 0.14
0.92 -> 0.93
```

不能把这些位移修成 `>=` 或引入任意 tolerance；那会改变既有统计定义。若未来需要严格 p 值可复现，应单独讨论统计量的数值稳定性，而不是在本次兼容性问题中混入假设检验修改。

### 5.3 存储形式不是根因

旧 `net$pval` 是 dense array，显式存储 p 值 0；新 `net$group$pval` 是 `SparseChatArray`，零值可以不存储。`group_edge_diff.csv` 中的 `NA` 可能只是“该侧显式零未存储”，不是缺失值。全空间比较已经把两侧恢复到同一坐标语义后得到 1,033 条差异。

因此：

```text
SparseChatArray 的零省略 = 存储约定
不是本次 p 值翻转的原因
```

---

## 6. 新函数设计审视：已确认问题、条件风险和未发现问题

### 6.1 `BUG-COMPAT-01`：min.percent 兼容性缺陷（已确认）

**位置：** `R/modeling.R:87-102`、`R/modeling.R:1178-1184`、`src/SpatialChat_Rcpp.cpp:257-263`。
**现象：** 新函数以 `signif()` 数值语义预计算固定阈值，不能复现 pinned upstream 的 `format()` 字符语义。
**影响：** 全量 1,563 层出现 148/37 条 group 支撑差异，并传播到 p 值。
**性质：** 这是兼容性 BUG；是否也是科学语义 BUG，取决于维护者是否把 upstream 的 `format()` 行为视为应保留的 oracle。`format()` 的字符比较本身明显可疑，但不能在本文件中擅自把它改成新的统计决策。

### 6.2 `TEST-GAP-01`：冻结基线 oracle 没有覆盖 pinned upstream 实际门控（已确认）

`tests_dev/reference-computeAvgCommunProb-baseline.R` 选择当前项目的 `computeAvgCommunProb_LR_Avg`；当前 helper 使用 `signif()`，而 pinned upstream package 使用 `format()`。

因此现有端到端测试证明的是：

```text
新 kernel signif 语义 ≈ 当前本地 reference signif 语义
```

不是：

```text
新 kernel ≈ pinned upstream format 语义
```

建议增加一个小型兼容性 fixture，至少覆盖：

- 整列进入科学计数法的低比例输入；
- `0.096~0.099` 与 `0.1` 的边界输入；
- observed 和 permutation 两种标签状态；
- group support、group prob 和 pval 三层结果。

这属于测试边界修正，不是改变假设检验方法。

### 6.3 `NUM-RISK-01`：严格 p 值比较对归约顺序敏感（已观察，暂不判定为 BUG）

`Pboot - Pnull > 0` 是离散判定，数值上极小的尾差可以改变 reject 计数。该现象在 TGFB1 样例中已观察到。

目前不建议：

- 把 `>` 改成 `>=`；
- 给差值加固定 tolerance；
- 更改 p 值分母或 null 构造。

这些都会改变既有统计方法。当前应把它记录为“同一统计方法在不同浮点归约路径下的可复现性边界”，待 gate 语义统一后再单独评估剩余 p 值差异。

### 6.4 `HARDEN-01`：参数输入校验可以加强（低优先级，未证明为本次 BUG）

`computeAvgCommunProb` 的公开参数包括 `min.percent`、`min.cells.sr`、`nboot` 和 `nthreads`。当前实现对 `nboot < 1` 在置换路径中有显式检查，但对以下非法输入没有形成统一的入口校验：

- `min.percent` 是否为有限值且位于 `[0, 1]`；
- `min.cells.sr` 是否为有限的非负整数；
- `nthreads` 是否为正整数；
- `nboot` 是否为正整数。

这是接口硬化候选，不是本次 MBM05 差异的原因。有效参数路径已有测试，不能因为这个低优先级问题推翻 kernel 设计。

### 6.5 非 BUG：`min.cells.sr` 和 denominator 的整数化

旧代码将 `prob@x` 全部替换为 1，再通过行/列求和计算：

- sender/receiver 的 stored-entry 数量；
- avg 的 denominator（包括显式零存储条目）。

新代码使用 `sr_out`、`sr_in` 和 `den` 直接计数。对于当前稀疏矩阵语义，这些量是整数值，累加顺序不会改变结果。已有测试包含显式零条目；没有证据表明该设计造成了本次差异。

### 6.6 非 BUG：permutation gather 和线程隔离

新 C++ permutation kernel 使用：

```cpp
g[c] = group_int[perm(c, b) - 1];
```

即与 R 的 `group[permutation[, nE]]` 一致，而不是把 permutation 值误当作组码。该问题曾在早期实现中存在，已经修复，并有值域校验和测试。

不同 boot 写入不同输出列，scratch buffer 位于每次 boot 内部；已有 `nthreads=1` 与多线程逐位一致测试。当前证据不支持 race condition 归因。

### 6.7 条件成立时可接受：直接从 signaling assay 构造 support

新实现不再读取旧 `net$tmp$Lavg/Ravg`，而是从 signaling assay 和 complex subunit 构造 `suppL/suppR`。这个优化依赖：

```text
cofactor factor = 1 + non-negative expression >= 1
```

因此 cofactor 不会改变 `> 0` 支持模式；复合物的 geometric mean 大于零等价于所有 subunit 大于零。

在当前 SpatialCellChat 的非负 signaling 数据契约下，这个替换有依据，不应把它误判为本次差异来源。如果未来允许带负值或改变 cofactor 语义，则应重新验证这个不变量。

---

## 7. 保守的后续处理建议

### 7.1 先明确兼容性边界，不改假设检验

当前有两个逻辑上不同的选择：

1. **严格 pinned-upstream compatibility：** 新实现需要复现 upstream 的实际门控行为，包括其 `format()` 字符比较语义；这会保留一个可疑但可观测的旧行为，且要考虑 R 版本/locale 的可移植性。
2. **明确采用 numeric `signif()` semantics：** 认为 upstream 的 `format()` 比较是旧代码 bug，保留新实现的数值门控；这需要把行为变化写入契约，并把 reference/oracle 改为明确的数值语义，而不能继续声称离散结果与旧代码逐位一致。

在维护者裁定前，不建议静默修改一侧，也不建议仅修改 p 值比较符来掩盖门控差异。

### 7.2 最小验证顺序

不需要重新跑 7 小时全流程。已有 post-filter object dump 足够支持以下阶段化验证：

1. 用统一的 gate policy 对同一 post-filter checkpoint 重放 observed average；
2. 用同一 gate policy 重放 permutation average；
3. 重新比较 group support、group prob 和 pval；
4. 再把剩余 p 值差异归类为 gate 传播或浮点归约边界。

只有在 gate 语义统一后，才有意义讨论剩余 p 值是否需要新的数值契约。

### 7.3 建议新增的回归断言

不改变假设检验方法，只增加以下行为护栏：

- `format()` 科学计数法与 numeric `signif()` 的差异 fixture；
- `0.09687` 这类 `signif(..., 1) == 0.1` 的边界 fixture；
- observed/permutation gate 的支持集比较；
- `Prob == 0` 后 pval 强制为 1 的存储与全空间读取比较；
- group support 差异必须能定位到具体 gate，而不能只比较总和。

---

## 8. 复现证据索引

已有全量证据：

- `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/old_postfilter_object.rds`
- `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/new_postfilter_object.rds`
- `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/group_summary.csv`
- `tests_dev/benchmarks/out/mbm05_serial_tensor_comparison/group_edge_diff.csv`
- `tests_dev/benchmarks/analyze_group_diff3.R` 至 `analyze_group_diff6.R`
- `tests_dev/benchmarks/compare_serial_tensors.R`

基础对比命令：

```bash
"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" \
  tests_dev/benchmarks/compare_serial_tensors.R \
  tests_dev/benchmarks/out/mbm05_serial_tensor_comparison

"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" \
  tests_dev/benchmarks/analyze_group_diff6.R
```

当前文档对 `SERIAL_COMPARISON_GROUP_AVG.md` 中旧的“`>=`/`>` tie 语义是主要原因”表述进行了机制层面的修正；原始数值证据不变。

---

## 9. 当前结论与开放风险

**已经确定：**

- 差异起点是 `computeAvgCommunProb`；
- 185 条 support 差异全部与旧/新 gate 方向匹配；
- 共同 support 的 `prob` 尾差来自求和顺序；
- p 值完全翻转与 gate 删除/新增以及 `Prob == 0` 规则直接相关；
- permutation 检验原理、null 构造和旧决策不需要在本次审查中修改。

**尚未裁定：**

- pinned upstream 的 `format()` 行为是必须兼容的历史 oracle，还是应被正式归类为旧代码 bug 并以 numeric semantics 取代；
- gate policy 统一后，剩余 1,033 个 p 值差异中有多少属于浮点归约平局；
- 是否需要把 `nthreads`、`min.percent`、`nboot` 等入口参数校验作为单独的低优先级 hardening 任务。

在上述裁定完成前，Plan A Note 应保持 `proposed`/审查状态，不应把当前全量结果直接标记为“严格离散等价”。
