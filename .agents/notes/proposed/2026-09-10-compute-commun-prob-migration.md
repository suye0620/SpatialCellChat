# Agent Note: computeCommunProb 迁移至最终 11-slot schema 与方案 A 稀疏化

Status: proposed
Governance: v1
Date: 2026-09-10
Decision type: architecture
Scope: `R/modeling.R` `computeCommunProb`/`computeAvgCommunProb`/`computeAvgCommunProb_LR_{Avg,Sum}`/`filterProbability`、`net$cell` 内部结构契约、`tests_dev/test-computeCommunProb.R`（新增）
Owner: SpatialCellChat maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

1. **schema 断裂**：`computeCommunProb`（modeling.R:29）仍读写旧 14-slot（`@data.signaling`、`@data.project`、`@options$parameter`、`net$tmp$prob.cell` 旧内部布局），在最终 11-slot 对象上必然报错；它是推断主链（阶段 7-16）第一个未迁移函数，是整段迁移的闸门。
2. **可运行性瓶颈**：主循环 `crossprod(1×nC, 1×nC)` 产出稠密 nC×nC 外积（每 LR 全非零）。实测 10 万细胞下单 LR 111.8 GiB，不可运行（COMPUTATION_OPTIMIZATION.md §9.2）；方案 A 稀疏化实测 299× 加速且逐位一致（`bench_parallel.R` `identical()` 断言）。
3. **并行反优化**：稀疏化后单任务 ~4.6ms，multisession 调度开销主导（实测顺序 0.23s vs 并行 0.60s/分块 2.45s，§9.3）。
4. **min.percent 门控真 bug**：`format(x, digits=1)` 返回字符串，`>=` 退化为字典序比较；实测（2026-09-10 format probe）：值 < 1e-4 触发科学计数法后该列门控**整列失效**——0.09（"9e-02"）、0.0001（"1e-04"）、2.5e-05（"2e-05"）、甚至 0（"0e+00"）全部错误通过 min.percent=0.1 门控。bench 脚本 `identity (no perm): FALSE` 已复现口径分歧。
5. **误导性报错**：scale.distance 校验的 stop 文案自相矛盾（"slightly smaller than 1/d.min" 照做仍 <1 失败）。
6. **死代码**：`d.spatial@x <- d.spatial@x / scale.distance`（modeling.R:130）恢复的局部副本在 L137 即被 `rm`，函数内再无引用；权威距离数据经 `res` 未缩放缓存入 `images$.distance`（L292）。该行可证零影响。

背景：基线代码约两年前部分由 AI 生成（维护者确认，2026-09-10），已对照论文第 2 章（式 2.7-2.17）逐项数学核查——核心模型（Hill 核心、复合物几何均值、SCL=1/d、自通讯 diag=max(1/d)、AGAN 叉积结构、方向性、DᵀPD 组级聚合）实现正确；上述 4-6 为确认的缺陷。

## Proposal

以**当前基线代码为正确计算基准**（维护者 2026-09-10 裁定），迁移到 11-slot schema 并实施方案 A；仅一处有意的行为变更（⑤ 门控修复）。五条基准裁定：

| # | 裁定 | 内容 |
|---|---|---|
| ② | 多基因 cofactor 组合 | 保持基线"逐基因调节后连乘"（`computeExpr_coreceptor/agonist/antagonist` 的 `apply(..., 2, prod)` 分支）为正确基准；单基因分支与论文式 2.8/2.9 逐字一致。论文式 2.8 后半与式 2.9 前文的"平均表达水平"表述按实现勘误（影响面实测：agonist 320/342、antagonist 425/459、co-I 320/486 对含 ≥2 基因成员） |
| ③ | 旁分泌边界容差 | 保持 `(interaction.range + tol)` 搜索阈值（tol=细胞半径，覆盖"扩散极限处恰有细胞"的极端情形，维护者裁定）；论文式 2.10 勘误注记。接触阈值 `contact.range + tol` = 15μm 与论文式 2.11 一致，不动 |
| ④ | 组级双汇聚 | `avg.type="sum"`（DᵀPD，论文式 2.14 原义，实际信号总量口径）与 `"avg"`（连接数平均，消除 group size 影响）双分支并存；默认 `"avg"` 不变，文档注明两口径语义 |
| ⑤ | min.percent 门控 | **唯一有意的行为变更**：两处 `1 * (format(dataLR_temp[,-1], digits=1) >= min.percent)` → `1 * (signif(dataLR_temp[,-1], 1) >= min.percent)`（Avg: L1342、Sum: L1271）。保留 1 位有效数字容差意图（维护者确认原意即此，字符串化为 AI 代写笔误），恢复数值语义。迁移测试量化 MBM05 上被新清零的组对数 |

### Schema I/O 契约

```text
读取：
  assay(object, "signaling")   ← 保持稀疏；全局 max 归一化在函数内（对应旧 @data.signaling 路径）
  LR(object, "LRsig")
  images(object, ".distance")  ← 校验 .parameters 与请求阈值一致（复用 spatial.R 缓存校验模式）；
                                  缺失时自动调 computeCellDistance 并以 L1 info 告知（proposed 默认）
  idents(object) / misc$.datatype

写入：
  net$cell$prob   ← SparseChatArray（nC×nC×nLR，命名 dgCMatrix list 直接构造，无 my_as_sparse3Darray 合并）
  net$cell$pval   ← LR 级不写（空间版以 filterProbability 分位截断替代 p 值；proposed 默认，评审可改）
  misc$.param$communication ← raw.use/Kh/n/scale.distance/use.AGAN/interaction.range/nLR/nLR1/
                              all.contact.dependent/all.diffusible/运行耗时
  .log_operation("computeCommunProb") + validObject（经 assay<-/params<- 访问器天然执行）

删除：@options$parameter 写入、net$tmp$prob.cell 临时层、my_as_sparse3Darray 调用
```

### 算法与性能契约（方案 A）

- 主循环每 LR 仅在 `P.spatial` 非零处计算：`col <- rep.int(seq_len(ncol(ps)), diff(ps@p)); lr <- L[ps@i+1] * R[col]; ps@x <- ps@x * (lr^n/(Kh^n + lr^n))`——与基线 `crossprod → map_dbl Hill → 逐元素乘` 数学等价且逐位一致（无求和归约，`bench_parallel.R` identical 断言背书）
- `data.use` 保持稀疏（删 `as.matrix` 稠密化）；L/R 均值向量循环内按 LR 现算（O(nnz)），不做 nLR×nC 稠密预计算；`computeAvgCommunProb` 同样按 LR 现算 L/R 行（删除对 `net$tmp$Lavg/Ravg` 缓存的依赖）
- AGAN 分支（computeExpr_agonist/antagonist + myElementwiseProduct 叉积）语义原样保留
- 默认顺序执行；`my_future_lapply` 保留为显式选项并暴露 `future.scheduling`（默认分块到 workers 数）
- Rcpp 加速（本期落地，维护者 2026-09-13 裁定）：新增 `cpp_prob_layer` kernel（src/SpatialChat_Rcpp.cpp，与 cpp_sum_layers 同文件），单遍 O(nnz) 融合"gather → Hill（n=1 快路径避 pow）→ AGAN 因子 gather → 紧凑化"，CSC 序直出 (i,p,x) 绕开 triplet 转换；OpenMP 按 nnz 并行（nthreads 参数暴露）；R 向量化实现保留为 @noRd 语义规范，测试中对 kernel 做 identical 逐位对拍。`cpp_sum_layers` 本函数不用（属 pathway/组平均跨层求和）；置换检验 Rcpp 化（方案 B，~140×）属 computeAvgCommunProb 后续
- 全链路复杂度：O(nLR × nnz + nC log nC)；MBM05 全量（967 LR × 29,536 cells）预估 ~4.5s 单进程

### 子函数契约与重写清单（Wave 1a/1b，依赖树核查 2026-09-16）

| 函数 | 处置 | 契约要点 |
|---|---|---|
| `computeCellDistance` / `createCellCellContactMatrixFrom_dspatial` | 保留 | 已迁移 KD-tree，10 万细胞实测 0.24-0.78s |
| `createPspatialFrom_dspatial` | 原地保留 | O(nnz) 取倒数 + 对角填 max(1/d)（自通讯=式 2.12）|
| `computeExpr_LR`/`computeExpr_complex` | 重写为 `.sc_expr_row`（@noRd） | 逐 LR 现算：单基因=稀疏行提取；复合物=子矩阵稠密化 + `exp(colMeans(log))`（与基线 `geometricMean` log0=-Inf→0 逐字节同语义；**禁 sparse log 陷阱**）；去重基因行缓存（471 唯一名 vs 3126 次查找）；缺失基因显式 stop（LRsig 可用性保证合法流不可达，`.sc_lr_feature_state` L1303） |
| `computeExpr_coreceptor`/`_agonist`/`_antagonist` | 重写为 `.sc_cofactor_factor` | 逐 LR 因子 + 按唯一 cofactor 名缓存 + 无 cofactor 跳过；多基因 ∏ 分支用 `apply(..., 2, prod)` 统一（单基因分支 prod 单元素=恒等，逐位一致） |
| `HillFunctionFordataLR`/`myElementwiseProduct`/`crossprod` | kernel 取代 | 单遍融合 + 精确基线乘法分组：v=(x0·h)·(AG_i·AG_j)·(AN_i·AN_j)，contact 掩码在 AGAN 前（基线顺序） |
| `my_as_sparse3Darray` | 删除 | SparseChatArray 即命名 list；层 dimnames 置空（基线 `list(NULL,NULL)` 内存注释语义） |
| `my_future_lapply`/`_sapply` | 本路径停用 | 顺序 progressr 进度条；其余调用方不受影响 |
| `scMatrixTruncation` | 下一波 | filterProbability 用 |

子函数验收：Wave 1a 单测以**基线函数在测试场直接对拍**（dense/sparse 同输入，`identical()`）；主体级对拍用**冻结基线参考实现**（tests_dev 快照脚本，plain 输入）于小 fixture 上对照——合法流逐位一致，非法输入从崩溃变明确报错（等效或更优）。

### 同波迁移约束

`filterProbability`（`tmp$prob.cell`/`@options` 的最小消费者）必须同批适配：读 `unclass(net$cell$prob)`（SparseChatArray 即命名 list，tmp 与 3D 数组合一）与 `misc$.param$communication`；`scMatrixTruncation` 截断后写回 `net$cell$prob`。治理规则：不加旧结构兼容分支。

### CLI 规约

按 2026-09-10 信息提示标准：subheader + 输入摘要 + 关键分支（contact/diffusible 分派）+ 结果统计（dims、k 参数、耗时）+ success；进度条仅经 `my_future_lapply` progressor。

### Spill 分层存储接口（本期仅 rds 实现，BPCells 预留）

大张量内存墙（10 万细胞均值口径 ~34 GB）以**可选 spill** 解决，后端可插拔：

```r
computeCommunProb(..., spill.dir = NULL, spill.backend = c("rds"))
# spill.dir 非空：主循环逐层 spill_write_layer 后释放，内存仅驻留当前层（≤166MB）
# spill.dir 为空：全内存（默认，MBM05 规模 7.3GB 无感）；外推超阈值时报错并建议设置 spill.dir
```

内部接缝为三函数（@noRd）：`spill_write_layer` / `spill_open_layer` / `spill_remove_layer`，以
`backend` 参数分发（本期唯一实现 "rds"：逐层 saveRDS(compress=FALSE)，double 逐位无损）。
后端标记记入 `misc$.param$communication$spill.backend`，`filterProbability` 等消费者按标记分发。
**BPCells 后端为预留扩展点**（维护者 2026-09-14 裁定暂不实现）：准入条件为验证清单全过——
float64 dgCMatrix roundtrip 逐位一致（含 0/Inf/亚正规）、磁盘占用与读写速度实测、
r-universe 安装顺畅度；过门后仅需在 spill_* 三函数加分支 + roundtrip 测试，签名零破坏。
BPCells 候选价值：索引位压缩（磁盘 ~0.5-0.7×）、>2^31 nnz 原生支持（百万细胞路径）、
对象路径指针持久化模式（Seurat 同款）；其真正的最优场景是下一波 assay 整数 counts 层
（位压缩 6-8×），已列为 Phase 2 占位。

## Constraints and invariants

- 方案 A 与基线在相同输入（同 scale.distance、同 ⑤ 修复前门控）下 `identical()` 逐位一致——⑤ 门控只作用于组级，细胞级 Prob 必须与基线逐位一致
- SCL ≤ 1 校验语义保留；`net$cell$prob` 每层为 dgCMatrix 且列名 = 对象细胞名
- 相同 seed 与输入重复运行逐位一致（rsvd 不涉及本函数；AGAN 为纯函数）
- 论文勘误清单（②③④①的论文表述偏差）随 note 附带，供论文修订引用

## Alternatives considered

1. 保留稠密外积 + multisession 并行：10 万细胞不可运行（111.8 GiB/LR）；稀疏化后并行实测慢 2.6-10.7×，拒绝
2. 为 `net$tmp$prob.cell` 建兼容别名：违反 schema 治理"迁移后删除旧结构、不叠加兼容分支"，拒绝
3. 本期引入置换检验 Rcpp：属 computeAvgCommunProb 范围，主循环稀疏化后它才是下一瓶颈；接口预留、实现后置，拒绝本期混入
4. ⑤ 保留字符串语义（"代码为基准"例外讨论）：确认为 AI 代写笔误而非有意设计（维护者 2026-09-10 确认原意即 signif 数值容差），保留将使门控在极小表达率场景整列失效，修复

## Consumer impact

- `R/modeling.R`：`computeCommunProb` 重写（方案 A + schema）、`filterProbability` 同波适配、`computeAvgCommunProb_LR_{Avg,Sum}` 门控行替换、`computeAvgCommunProb` 改为自算 L/R 行
- 下游链（下一波）：`filterCommunication`/`computeCommunProbPathway`/`aggregateNet`/`netAnalysis_computeCentrality`/`computeCommunField`/`subsetCommunication`/`identifyEnrichedInteractions` 按 `net$cell`/`net$group` 新结构迁移
- `tests_dev/test-computeCommunProb.R`（新增）：小 fixture 新旧对拍（`identical`，含 scale=1 与 scale=9 两档）、⑤ 门控差异量化、拒绝路径（无距离缓存且无法自动补、无 LRsig）、params/log 检查
- `tests_dev/bench_parallel.R`：复跑确认顺序执行最优（已有骨架）
- 端到端：`tests_dev/test_data/MBM05Chat_pre_inference.rds` 检查点直接续跑；老脚本（Slide-seq(1).R 语义、scale.distance=9、use.AGAN=T）作为回归基准
- `man/`：computeCommunProb/filterProbability Rd 更新
- 无内部其他调用方（grep 确认 computeCommunProb 无 R/ 内部调用）

## Consequences

- 10 万细胞规模主循环可行（预估 ~5s/千 LR 单进程）；Phase 3 后续函数按 `net$cell`/`net$group` 新结构推进
- ⑤ 修复使部分组对组级通讯被新清零（幽灵通讯消除）——与历史结果的组级数值差异需在论文/文档中注明版本界限
- `net$cell$prob` 成为本 LR 级结果的唯一权威布局；旧序列化对象（legacy 14-slot RDS）经重建壳后不再直接可用，需走 updateObject 类迁移（后续 note）
- 方案 A 后主循环 ~4.6ms/LR，置换检验（computeAvgCommunProb，外推 ~25min）成为新瓶颈——下一 note 的输入

## Evidence

- COMPUTATION_OPTIMIZATION.md §3（瓶颈实测）、§4（方案 A 代码与 279×/299×）、§9.2（10 万细胞 111.8 GiB/LR）、§9.3（并行有害实测）、§9.4（identity FALSE 复现）
- `tests_dev/bench_parallel.R`（方案 A identical 断言骨架）、`bench_avgcommunprob.R`、`kdtree_correctness.R`（距离侧已闭合）
- 论文第 2 章式 2.7-2.17（201-2023282010085-苏烨 revised.pdf）与 2026-09-10 五条基准裁定（②③④①⑤，维护者逐条拍板）
- 2026-09-10 format probe 实测（科学计数法整列失效证据）
- 基线代码逐行核对记录（modeling.R:29-319、337-351、369-375、419-436、2316-2352、2491-2512、2527-2548、1265-1370）
- `tests_dev/test_data/MBM05Chat_pre_inference.rds`（阶段 1-6.5 检查点，967 信号基因 × 29,536 细胞，1563 LR 对）
- 2026-09-13 nnz 标定（双源）：SpatialChat_1/2 已算结果每层 nnz/nC² 中位 2.8e-2（post-filter，r̂≈4.6e-2 ↔ 0.735² 表达密度交叉验证交集模型）；MBM05 精确模式预测（150 LR 采样）：填补后 suppL 中位 4.1% / suppR 3.3%，每层 nnz 中位 3.1e4 / 均值 4.2e5 / 最大 1.4e7（r̂ 重尾 1.4e-3/1.9e-2/0.645，均值:中位 = 13.7）；MBM05 输出峰值 6.57e8 nnz = 7.3 GB（pre-filter），filterProbability 后 0.37-0.73 GB；10 万细胞投影（度恒定 735、nLR=2000）均值口径 ~34 GB（±1 数量级：3.4-340 GB）——重尾由广泛表达对主导（top: CDH3-CDH3 supp 79.7%/79.7%，count=1.4e7），运行时外推守卫为必要防线
- 2026-09-19 filterProbability v2 优化（模式行计数 hoist + d 侧模式索引跨层缓存 + P 侧手工 CSR）：MBM05 全量 1563 层 1198.2s -> 54.9s（21.8x），v2 输出与 v1 存档逐位一致（1563/1563，@x/@i/@p 全等）；SpatialChat_2 三方（v2/v1/基线复刻）35/35 逐位一致；fixture 回归 41 PASS。未来全局检测：sequential plan 下 future.globals=FALSE（缓存 env 被引用图三重计数 36.46 GiB 误报上限）

## Acceptance criteria

- [x] 11-slot 对象上 `computeCommunProb` 全流程可运行：assay 层读取、`net$cell$prob` SparseChatArray 写入、`misc$.param$communication` 记录、日志与 validObject 通过
- [ ] 细胞级 Prob 与基线在相同参数下 `identical()` 逐位一致（含 scale=9 与 scale=1 两档、AGAN 两分支、contact-dependent 分派）
- [ ] 方案 A 路径无稠密 nC×nC 外积、无 nLR×nC 稠密预计算、无 my_as_sparse3Darray；MBM05 全量实测耗时记录在案
- [ ] `cpp_prob_layer` kernel 落地（src/SpatialChat_Rcpp.cpp）：与 R 向量化语义规范在全部 fixture 上 `identical()` 逐位一致；nthreads 参数生效（多线程结果与单线程逐位一致）；10 万细胞合成数据实测耗时与峰值内存记录在案
- [ ] ⑤ 两处门控行替换为 `signif` 数值语义（已替换，fixture 验证 ghost 微速率清零、正常速率不变）；MBM05 上量化新清零组对数并记录（待全量运行）
- [x] `filterProbability` 在新结构上运行：读 `net$cell$prob`/`misc$.param$communication`，spill 模式读盘消费并回填（文件清理 + spill params 清空），截断结果与旧结构相同输入行为逐位一致（fixture 4 层对拍 + 零值捷径与基线 quantile==0 分支等价验证）
- [ ] scale.distance 校验逻辑逐位保留；L130 死代码已删除；报错文案给出可执行指引
- [ ] 默认顺序执行；`future.scheduling` 暴露且分块默认 `ceiling(nLR/workers)`
- [ ] spill 接口（rds 唯一实现后端）：spill.dir 模式与全内存模式在同输入下 `identical()` 逐位一致；spill 目录逐层即读即删、流程结束清空；守卫超阈值报错文案包含"设置 spill.dir"指引；`spill.backend` 参数定型为 match.arg 单选 "rds"（未来后端纯 additive 扩展）
- [ ] BPCells 后端为预留扩展点（2026-09-14 裁定暂不实现）：spill_* 三函数 backend 参数化就位；验证清单（float64 roundtrip 逐位一致、磁盘/速度实测、r-universe 安装）作为后端准入条件保留于本 note，落地时不改主流程
- [ ] CLI 消息符合分级标准；`tests_dev/test-computeCommunProb.R` 全绿；`bench_parallel.R` 复跑通过
