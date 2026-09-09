# SpatialCellChat 计算优化方案讨论纪要

> 2026-08-09 · 针对 `computeCommunProb` / `computeAvgCommunProb` 的并行策略审计与加速方案
> 关联文档：[REFACTORING_PLAN.md](REFACTORING_PLAN.md)、[VISUALIZATION_AUDIT.md](VISUALIZATION_AUDIT.md)
> 本文档整理自连续三轮讨论：实现细节审计 → 加速方案设计 → KD-tree 正确性论证。所有结论均附实测数据，可复现。

---

## 1. 背景与目标

`computeCommunProb`（modeling.R:29）与 `computeAvgCommunProb`（modeling.R:989）是 CellChat 空间版推理管线的两个耗时主体。当前实现依赖 `future.multisession` 并行，但在 Visium 规模（5k spots、~1000 LR）下**主循环为 O(nLR × nC²)**，置换检验段为 O(nLR × nboot × 稀疏矩阵乘)。本讨论的目标：

1. 审计现有并行策略与实现细节；
2. 用基准实测定位瓶颈；
3. 设计分层加速方案（R 层向量化 → Rcpp → 架构级）；
4. 论证 KD-tree 替换全对距离计算的正确性。

---

## 2. 现状实现细节

### 2.1 computeCommunProb 流程（modeling.R:29-303）

```
主进程（串行）                            并行段（future.multisession）
├─ data.signaling → as.matrix 稠密化
├─ computeCellDistance（Rfast::Dist 稠密全对距离）      ← 启动段瓶颈
├─ createPspatialFrom_dspatial（P.spatial, nC×nC 稀疏）
├─ computeExpr_LR ×2 ──────────────► 内部 my_future_sapply（coreceptor 同理）
└─ my_future_lapply(1:nLR)  ★ 主并行段（每 LR 一个 future）
     ├─ dataLR ← crossprod(1×nC, 1×nC)   ★ nC×nC 稠密外积（全非零）
     ├─ Hill(dataLR@x)                     ★ purrr::map_dbl 逐元素
     ├─ dataLR * P.spatial
     └─ AGAN 分支（computeExpr_agonist/antagonist + 稀疏乘法）
└─ my_as_sparse3Darray(Prob.cell_)  ★ 主进程串行 list→3D 合并
```

**关键事实**：`dataLR = crossprod(x_, y_)`（x_, y_ 为 1×nC 稀疏行）是 **nC×nC 稠密外积矩阵**（25M 元素、全非零），Hill 函数逐元素 `purrr::map_dbl`，再与 P.spatial 逐元素相乘。**每个 LR 的复杂度是 O(nC²)，与空间矩阵稀疏度无关**——这是主循环的设计瓶颈。

### 2.2 computeAvgCommunProb 流程（modeling.R:989-1245）

```
串行: prob.sum 计数（map_dbl，快）→ 写回 Prob/Pval 的 for 循环
并行段1: my_future_sapply(LRsig.use.idx)  组平均
         每 LR: computeAvgCommunProb_LR_Avg
          = model.matrix(~group-1) + 2×aggregate() + 2×crossprod()  ≈ 14 ms
并行段2: my_future_lapply(LRsig.use.idx)  置换检验  ★ 最大热点
         每 LR: 内层 sapply(1:nboot) 串行 100 次同样的 14 ms 计算
```

### 2.3 并行机制（utilities.R:658/707）

- `my_future_lapply/sapply` = `future.apply::future_lapply` + progressr 进度条；
- **未设 `future.scheduling`** → 默认每元素一个 future（nLR 次任务、nLR 次 globals 传输）；
- 闭包自由变量（`P.spatial`、`dataLavg/dataRavg` 稠密、`prob.cell_` 整个 list、`permutation` nboot×nC 矩阵）**随每个 future 完整序列化到 worker**；
- 无嵌套并行（computeExpr_agonist/antagonist 为纯函数）；
- **不一致**：`computeCommunProbPathway` do.cell 段（modeling.R:1574）用 `pbapply::pbsapply`——**默认串行**，与 future 体系脱节。

---

## 3. 瓶颈诊断（实测）

### 3.1 基准 1：computeCommunProb 主循环

合成 Visium 规模：5k spots、~200k 空间边、50 LR、4 workers（`tests_dev/bench_parallel.R`）：

| 场景 | 耗时 | 结论 |
|---|---|---|
| 当前实现，顺序（纯计算） | **64.3 s / 50 LR**（1.28 s/LR） | 稠密外积是主瓶颈 |
| 外积稀疏化 O(nnz) | **0.23 s / 50 LR**（4.6 ms/LR） | **279×**，`identical()` 逐位一致 |
| multisession 默认调度 | 25.97 s | 并行有收益（传输被计算掩盖） |
| multisession + 分块 | 24.97 s | 此规模下无差异 |
| 稀疏版 + 并行默认 | 0.28 s | 任务太轻，调度开销开始主导 |
| 稀疏版 + 并行分块 | 1.85 s | 分块反而慢（任务 4.6ms 级） |

### 3.2 基准 2：组平均 / 置换检验

合成：5k spots、k=10 组、200k 边（`tests_dev/bench_avgcommunprob.R`）：

| 实现 | 单次调用 | 置换总量（nLR=500, nboot=100） |
|---|---|---|
| 当前（model.matrix+aggregate×2+crossprod×2） | 14.0 ms | **~700 s** |
| 稀疏矩阵版（one-hot 预计算） | 13.0 ms | ~680 s |
| Rcpp block-sum 预估 | ~0.1 ms 级 | ~5 s（~140×） |

**关键结论**：组平均/置换段里 aggregate/model.matrix 只占小头，**4 次稀疏矩阵乘法才是成本**——R 层重排（14→13ms）几乎无收益，必须 Rcpp 化（单次遍历 O(nnz) 完成全部块聚合）。

### 3.3 实测外推（Visium 5k spots、1000 LR）

- `computeCommunProb` 主循环：当前 ≈ 21 min → 稀疏化后 ≈ 5 s（单进程）；
- `computeAvgCommunProb` 置换段：当前 ≈ 25 min（500 有信号 LR）→ Rcpp 后 ≈ 5 s。

### 3.4 顺带发现的隐藏行为（需确认是否修复）

`computeAvgCommunProb_LR_Avg` 的 `1 * (format(dataLR_percent, digits = 1) >= min.percent)` 是**字符串比较**而非数值比较。实测：

```
format(0.0999999, digits=1) = "0.1"    # 四舍五入到 1 位有效数字
format(0.0001,    digits=1) = "1e-04"
"1e-04" >= "0.1" → TRUE               # 字符串比较：'1' > '0'
```

即表达率 0.0001（万分之一）的细胞群会被误判为"过表达"（数值语义应为 FALSE）。这是原 CellChat 实现的隐藏 bug/特性。**改为数值比较（signif(x,1) >= min.percent）会改变下游结果**，属行为变更，需拍板。

---

## 4. 加速方案

### 方案 A：computeCommunProb 主循环外积稀疏化（核心，279×，零并行依赖）

**数学基础**：`dataLR[i,j] = L[i]·R[j]`（外积），且 `P1_Pspatial = Hill(dataLR) * P.spatial` 只在 P.spatial 非零处非零。**无需构建稠密 dataLR**——直接在 P.spatial 的 CSC 非零位置计算：

```r
## 替换 crossprod + Hill + 乘法三步（O(nC²) → O(nnz)）
ps <- P.spatial
col <- rep.int(seq_len(ncol(ps)), diff(ps@p))    # CSC 列索引（dgCMatrix 无 @j 槽）
lr  <- L[ps@i + 1L] * R[col]                     # 外积值
ps@x <- ps@x * (lr^n / (Kh^n + lr^n))            # Hill 向量化（同时消灭 map_dbl）
```

- 数值逐位一致（已 `identical()` 验证）；
- 顺带消灭 `purrr::map_dbl` 逐元素 Hill；
- AGAN 分支不受影响（data.agonist 为 1×nC 稠密行，O(nnz) 乘法）；
- **并行从"必须"降级为"可选"**：单进程 4.6ms/LR → 1000 LR ≈ 5s。

### 方案 B：置换检验 Rcpp 化（~140×，第二个大头）

```r
## cpp_block_avg(prob, perm, group, dataLR, min.percent, min.cells.sr)
## 对每个 perm e：遍历 prob 的 nnz 一次，
##   累加 block sum / 二值计数 → 桶 (group[perm[i]], group[perm[j]])
##   另算 dataLR 两列的块均值（O(nC)）
## 输出 nboot × k × k 三件套（Prob.avg、二值计数、percent）
```

- 单次 ~0.1ms 级 → 5e4 次 ≈ 5s（vs 700s）；
- 更激进：一次遍历把 nboot 个累加器全算完（O(nnz×nboot)，内存 nboot×k²=1e4 桶）→ **整个置换检验 1 个 Rcpp 调用 ~0.5s**；
- 可配 RcppParallel 按 perm 并行，或保持单线程（已足够快）；
- **数值一致性注意**：R/BLAS 累加顺序与 Rcpp 逐元素不同，浮点差异 ~1e-16，理论上可能翻转个别 `Pboot - Pnull > 0` 边界比较；`min.percent` 的 `format(digits=1)` 字符串语义需在 Rcpp 中显式复现或改为数值语义（见 3.4）。

### 方案 C：并行配置修正（低风险辅助）

- `my_future_lapply/sapply` 暴露 `future.scheduling` 参数，默认 `ceiling(n / nbrOfWorkers())`——globals 传输从 nLR 次降到 workers 次（真实规模 dataL/R 稠密 80MB+、prob.cell_ GB 级时是硬瓶颈）；
- `computeCommunProbPathway` do.cell 的 `pbsapply` → `my_future_lapply`（统一并行体系）；
- 稀疏化后任务变轻（4.6ms），**并行收益需重测**——可能单进程 + OpenMP 更优，或保留并行但只分块到 workers 数。

### 方案 D：架构级（配合现有重构，REFACTORING_PLAN Phase 3）

- 删除 `my_as_sparse3Darray` 的 list→3D 合并（主进程串行 + 长表内存）；新 schema 的 `SparseChatArray` 就是 list，直接保留；
- pathway 聚合用已实现的 `cpp_sum_layers`（跨层求和），替代 `my_as_sparse3Darray + marginSumsSparse`；
- `computeCellDistance` 换 KD-tree（见下节，与方案 A 正交：A 解决每 LR 的 O(nC²)，KD-tree 解决构建 P.spatial 本身的 O(nC²)）；
- `data.use` 保持稀疏（不 `as.matrix`）。

**落地后全链路复杂度**：`computeCommunProb` 从 O(nLR × nC² + nC²) 变为 O(nLR × nnz + nC log nC)。

---

## 5. KD-tree 正确性论证

### 5.1 原理

KD-tree 把空间点组织成二叉树：反复"选一维、按中位数切分"，直到叶子里只剩少数点。范围搜索（radius search）时沿树下行，**用点到子树超平面的距离做剪枝下界**——该下界超过搜索半径的子树，数学上不可能含结果点。

**标准 KD-tree 的范围搜索是精确算法，不是近似**。剪枝条件有严格几何证明，非概率性。因此**理论上零遗漏**。

### 5.2 实测对抗性验证（tests_dev/kdtree_correctness.R）

构造 n=3070 点：随机点云 + 20 对**距离精确等于阈值**的边界点 + 30 个**近重复坐标**点（间距 1e-12），与暴力全对距离（`dist()`，即原 `Rfast::Dist` 语义）对照：

```
n = 3070  threshold = 10
ref pairs: 134838  | kd pairs: 134838
MISSING in kd : 0
EXTRA in kd   : 0
boundary pairs dropped by kd: 0 of 20     # 边界（<= thr）全部保留
dup-group neighbor counts: 一致            # 近重复点邻居数完全相同
```

结论：在 2D 空间坐标 + 默认 `KmknnParam` 下，`BiocNeighbors::findNeighbors` 与暴力全对**完全一致**，含边界与近重复点。

### 5.3 遗漏风险真实存在的 4 个场景（均为用法问题，非 KD-tree 本身）

| 风险 | 场景 | 后果 | 规避 |
|---|---|---|---|
| 近似参数 | `HnswParam()`、RANN/FNN 的 `eps > 0` | 允许"接近最优"邻居，可能漏 | 用精确实现：`KmknnParam()`/`VptreeParam()`（BiocNeighbors 默认即精确） |
| kNN 冒充半径搜索 | `findKNN(k=固定值)` 近似"距离 ≤ thr 的所有邻居" | **半径内但不在前 k 的点被漏**——最常见坑 | 必须用 `findNeighbors(threshold=)`（半径搜索） |
| 完全重复坐标 | 两个 spot 坐标完全相同（d=0），部分实现只返回重复点之一 | 同坐标第二、三个点可能漏 | 2D 空间数据极少见；且原代码 `d > 0` 本身排除 d=0，语义对齐 |
| 高维退化 | 维度 > ~20 时退化为线性扫描 | 性能问题，非正确性问题 | 2D 坐标无此问题 |

### 5.4 落地要点与验证防线

1. **必须用 `findNeighbors`（半径搜索），不要用 `findKNN`**——这是唯一真正的结构性风险，选对 API 即消失；
2. **保留 `d > 0` 过滤**：`findNeighbors` 会把自身（d=0）计入邻居，原代码 `d.spatial > 0` 排除——落地时过滤 0 距离对齐原语义；
3. **边界语义**：原实现 `<= threshold`，实测 `findNeighbors` 同为包含边界（`distance <= threshold`），一致；
4. **对称性**：`findNeighbors` 按"每点自己的邻居"返回，天然去重（每对一次），顺带解决原实现"每对距离存两次"的问题——结果矩阵可直接 `symmetric = TRUE` 构建；
5. **回归防线**：`tests_dev/kdtree_correctness.R` 为对抗性测试骨架；改 `computeCellDistance` 后追加 tutorial 数据对比新旧 `d.spatial` 的 `nnz` 与 `sum` 完全一致。

### 5.5 改造示意（替换 spatial.R:547-565）

```r
# 现在：Rfast::Dist 全对 + 逐列 which()
d.spatial <- Rfast::Dist(coordinates)
idx <- purrr::map_dfr(1:nC, ~ which(d.spatial[, .x] <= thr & d.spatial[, .x] > 0) ...)

# 改造后：
thr <- (interaction.range + tol) / ratio                 # 阈值换算逻辑保留
res <- BiocNeighbors::findNeighbors(
  coordinates, threshold = thr,
  BNPARAM = BiocNeighbors::KmknnParam(),                 # 内部即 KD-tree，精确
  get.distance = TRUE
)
## res$index[[i]] / res$distance[[i]] → 稀疏矩阵 (i, j, x) 三元组（过滤 d > 0）
```

依赖说明：`BiocNeighbors` 已在 Imports，**无新依赖**；`FNN`/`RANN` 已安装但为备选（RANN 对重复点的索引语义有历史不一致，优先 BiocNeighbors）。

---

## 6. 待拍板讨论点

1. **优先级**：方案 A（279×）与 B（140×）覆盖两个最大热点，且都**不需要并行**。方向选择：先稀疏化/Rcpp 化、并行降级为可选配置？还是保留并行框架 + 逐步加速？
2. **数值一致性契约**：A 可逐位一致；B（Rcpp block-sum）与 R/BLAS 累加顺序不同，p 值结果可能有个位数级边界差异。可接受吗？
3. **`format(digits=1)` 字符串比较 bug**（3.4）：重构时改为数值比较会改变下游结果——修还是保持原语义？
4. **真实数据规模**：合成基准为 Visium（5k spots、200k 边、k=10）。若目标含 slide-seq（5 万 spots），`computeCellDistance` 的稠密瓶颈会先爆——方案 D 的 KD-tree 需提前。
5. **落地方式**：方案 A+B 可直接实现并用 tutorial 数据（`visium_human_psoriasis.RData`）回归。

---

## 7. 复现与验证指引

| 脚本 | 内容 |
|---|---|
| `tests_dev/bench_parallel.R` | computeCommunProb 主循环：稠密外积 vs 稀疏 O(nnz)（279×）、multisession 调度对比、数值一致性断言 |
| `tests_dev/bench_avgcommunprob.R` | 组平均/置换检验：当前 vs 稀疏矩阵版单次调用成本、置换总量外推 |
| `tests_dev/kdtree_correctness.R` | KD-tree 对抗性正确性测试：边界点 + 近重复点，与暴力全对断言零差异 |

运行环境：R 4.5.3；依赖 Matrix、future、future.apply、BiocNeighbors（均已安装）。

---

## 8. 遗留与未验证项

- 方案 B 的 Rcpp 实现未编码（预估数字基于单次调用成本外推，落地后需以 tutorial 数据实测）；
- 稀疏化后并行收益的重新标定未做（4.6ms/LR 任务量下并行可能无益甚至有害）；
- `computeCellDistance` KD-tree 改造未落地（正确性已验证，性能收益未实测）；
- 完全重复坐标（d 恰为 0）场景仅在语义层面论证（原代码排除 d=0，KD 实现过滤 d>0 后一致），未做专门构造测试。

---

## 9. 2026-08-30 十万细胞压力审计补充

> 针对第 8 节"遗留与未验证项"的实测闭合：距离缓存 10 万细胞压力、稀疏化后并行收益重新标定、稠密外积内存上界。

### 9.1 十万细胞距离缓存压力实验（`computeCellDistance`）

| 场景 | N | 阈值 | 耗时 | d.spatial nnz | adj.contact off-diag | 结果内存 |
|---|---|---|---|---|---|---|
| 均匀低密度（0–100000 范围） | 100,000 | 10 | 0.24 s | 354 | 92（+10 万对角） | 26.3 MB |
| 400×250 网格（avg deg 11.9） | 100,000 | 2.1 | 0.78 s | 1,193,504 | 398,700（+10 万对角） | 20.1 MB |

- `queryNeighbors`（VptreeParam, `num.threads=1`）在 10 万细胞上保持亚秒级，结果严格 `0 < d <= threshold`、`adj.contact` 全 1 且对角补 1；`benchmark-computeCellDistance.R` 的稠密参照一致性断言通过。
- **结论：距离缓存本身已满足 10 万细胞（O(nnz) 上界，内存 ~20 MB 级），不再是扩展瓶颈。**

### 9.2 稠密外积内存上界（computeCommunProb 主循环）

n=2000 探针：`crossprod(1×n, 1×n)` 产出 4M nnz 的 dgCMatrix，45.8 MB。外推 n=100,000：**111.8 GiB / 每个 LR 对**。`computeCommunProb` 主循环（modeling.R:202-209）在 10 万细胞下**单 LR 即超出常规内存，不可运行**——方案 A（O(nnz) 外积稀疏化）由"279× 提速"升级为"10 万细胞可运行性的前提"。

### 9.3 稀疏化后的并行重新标定（复跑 `bench_parallel.R`）

| 实现 | 顺序 | multisession 默认 | multisession 分块 |
|---|---|---|---|
| 当前稠密外积 | 68.86 s / 50 LR | 26.28 s（2.6×） | 23.83 s（2.9×） |
| 稀疏 O(nnz) | **0.23 s**（299×） | 0.60 s（**慢 2.6×**） | 2.45 s（**慢 10.7×**） |

- 已核实 `future.apply::future_lapply` 默认 `future.scheduling = 1`（逐元素一个 future、globals 随任务重复传输）；分块到 workers 数在稠密路径仅 +10%。
- **结论：方案 A 落地后并行从"可选增益"变为"有害"（任务 4.6 ms 级，调度/序列化占主导）；应默认顺序执行或只分块到 workers 数。**

### 9.4 组平均/置换检验（复跑 `bench_avgcommunprob.R`）

| 实现 | 单次调用 | 置换总量外推（nLR=500, nboot=100） |
|---|---|---|
| 当前（model.matrix+aggregate×2+crossprod×2） | 14.50 ms | ~725 s |
| 稀疏 one-hot 版 | 13.00 ms | ~694 s |
| Rcpp block-sum 预估（文档第 4 节方案 B） | ~0.1 ms 级 | ~5 s |

- **`identity (no perm): FALSE` 实测复现**：稀疏版的数值 `pct >= min.percent` 与当前 `format(digits=1)` 字符串比较在边界值（如 0.0999999 → "0.1"）语义不同——即文档 3.4 节隐藏行为，任何替换都必须显式拍板保持字符串语义还是修正为数值语义。
- **结论：R 层稀疏重写收益仅 ~4%，置换段必须 Rcpp 化（~140×）才有实质改善。**

### 9.5 其余内存与并行风险复核

| 路径 | 位置 | 10 万细胞量级 | 风险 |
|---|---|---|---|
| `data.use <- as.matrix(data)` 稠密化 | modeling.R:49, 1025 | 2k 基因 × 100k 细胞 = 1.49 GiB/份 | 主进程常驻 ×2（LR 表达计算） |
| `dataLavg/dataRavg` 稠密 nLR×nC | modeling.R:176-177 | 1000 LR × 100k = 0.75 GiB/份 | 随每个 future 传输（scheduling=1） |
| `permutation` nboot×nC | modeling.R:1154 | 100×100k int = 38.1 MB（实测 0.75 s） | 可接受 |
| `my_as_sparse3Darray` 串行长表合并 | modeling.R:279, 1582-1619；analysis.R:4071, 4144；visualization.R:5516, 5810, 6086 | 1000 LR × 1.19M nnz ≈ 22+ GiB 长表 | 主进程串行 + 内存峰值 |
| `makeGridSpatialCellChat` / `computeGridSize` 的点-网格关系 | spatial.R:271, 371 | 100k 点与网格的稠密逻辑阵会随 `n × G` 爆炸 | **已稀疏化（2026-09-02）**：`sgbp` + `lengths`/`tabulate`/`split` 派生，O(nnz)，不分配稠密点-网格矩阵 |
| `computeCommunProbPathway` do.cell 用 `pbsapply` | modeling.R:1574 | — | 与 future 体系不一致，串行执行 |

### 9.6 优化优先级

2. **P0（10 万细胞可运行性前提）**：方案 A 外积稀疏化（299×，逐位一致，含消灭 `map_dbl` Hill）；`makeGridSpatialCellChat` 与 `computeGridSize` 已使用 `sparse = TRUE` 成员关系（**2026-09-02**，见 `tests_dev/test-makeGridSpatialCellChat.R` 与 `tests_dev/test-computeGridSize.R`）；新 schema 下删除 `my_as_sparse3Darray` 串行合并（`SparseChatArray` 即 list，直接保留）。
2. **P1**：方案 B 置换检验 Rcpp 化（~140×，需先拍板 `format(digits=1)` 语义）。
3. **P2**：并行配置修正——`my_future_lapply/sapply` 暴露 `future.scheduling` 默认分块到 workers 数；`computeCommunProbPathway` do.cell 统一为 `my_future_lapply`；方案 A 落地后重测并行取舍（当前证据指向顺序执行）。
4. **P3（决策项）**：`min.percent` 字符串比较语义（修复会改变下游结果）。

### 9.7 2026-09-02 稀疏点-网格规模复验

`tests_dev/test-makeGridSpatialCellChat.R` 的确定性 10k 压测通过：10,000
输入点、9,801 个占用网格、39,204 个点-网格命中，耗时 3.980 s，R
`Ncells_delta` 为 2.6 MB。新增
`tests_dev/test-grid-scale-100k.R` 在同一实现上通过 100,172 个输入点、
99,540 个占用网格、398,160 个命中，耗时 168.070 s，R
`Ncells_delta` 为 18.1 MB；输出对象保持 `dgCMatrix` assay 并通过
`validateSpatialCellChat()`。

这里的 `Ncells_delta` 是 R 垃圾回收器报告的对象节点增量，不是 Windows
进程峰值 RSS。两项规模测试均使用 `sf::st_intersects(..., sparse = TRUE)`
及命中列表派生聚合量；未声明不可移植的峰值 RSS 结论。
