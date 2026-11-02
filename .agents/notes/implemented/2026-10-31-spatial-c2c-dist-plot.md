# Agent Note: merge the communication distance plots into a single `spatialC2CDistPlot()`

Status: implemented
Governance: v1
Date: 2026-10-31 (implemented 2026-11-01)
Decision type: api
Scope: R/spatial.R (`spatialCommDistPlot` → `spatialC2CDistPlot`), R/visualization.R (`spatialCCCDistPlot`), NAMESPACE, man/, tests_dev/test-comm-distance.R
Owner: maintainer
Impact: medium (one baseline function becomes a deprecated alias; the working-tree `spatialCommDistPlot` is renamed before its first commit)
Supersedes: §2.5/2.6 of `2026-10-29-b-group-remaining-migration.md` (distance plots as standalone migration)
Superseded by: none

## 1. 背景与动机

- 工作区有一份**未提交**的 `spatialCommDistPlot()` 实现（已合并 `communicationDistPlot` / `communicationDistPlot2`，旧名转 `.Deprecated`，测试 21 项全过）。
- 复核发现其 `enriched.only = FALSE` 分支与既有 `spatialCCCDistPlot()`（基线即有、已导出）**能力完全重叠**：同一批数据（`d.spatial@x` 与 `d.spatial * adj.contact`）、同一横轴区间，仅样式与对角线处理不同。
- 用户拍板：**三合一**，单一函数命名 `spatialC2CDistPlot()`（cell-to-cell）。本笔记是整合的文档计划，实施待确认。

## 2. 单一入口

```r
spatialC2CDistPlot(
  object,
  split.by          = c("annotation", "distance"),
  enriched.only     = TRUE,
  density.alpha     = 0.5,
  interaction.range = NULL,   # 仅在距离缓存缺失、需要现算时使用
  contact.range     = NULL,
  tol               = NULL
)
```

## 3. 核心契约：split.by × enriched.only 的四类行为

| # | `split.by` | `enriched.only` | 上面板 | 下面板 | 面板关系 | 需要 `computeCommunProb()`？ |
|---|---|---|---|---|---|---|
| 1 | `"annotation"` | `TRUE` | 可扩散层（Secreted / ECM / Non-protein） | 接触层（Cell-Cell Contact） | **互斥** | 是 |
| 2 | `"distance"` | `TRUE` | 全部推断记录 | 其中落在 `contact.range` 内 | **嵌套** | 是 |
| 3 | `"annotation"` | `FALSE` | 全部细胞对 | 接触范围内细胞对 | 嵌套 | **否** |
| 4 | `"distance"` | `FALSE` | 同 3（输出完全相同） | 同 3 | 嵌套 | **否** |

约定：

- **3 ≡ 4**：`enriched.only = FALSE` 时没有层信息，`split.by` 无从发挥，两个取值落到同一段几何面板代码。保留两个取值只为参数矩阵完整；文档与实现都要明说"此时 `split.by` 被忽略"。
- 行为 1/2 的记录单位是"(细胞对 × 非零 LR 层)"，一个对带几层就计几次；行为 3/4 的单位是"细胞对"。同一张图内不混两种单位。
- "非零"指概率 > 0，**不是统计显著**（不设阈值，`filterProbability()` 属后续步骤）。
- 行为 1 按 `LR$LRsig$annotation` 逐层分类，不用旧实现按位置切的 `nLR1` 阈值。

## 4. 前置条件矩阵（消除能力重复的关键）

| 数据 | 行为 1/2（enriched） | 行为 3/4（几何） |
|---|---|---|
| `net$cell$prob` | 必须（spill 模式 NULL 即报错） | 不读 |
| `images$.distance` | 必须存在（`computeCommunProb()` 必然留下）；缺失即报错 | **缺失或过期 → 用 `images$spatial.factors` + 三个显式参数经 `computeCellDistance()` 现算**（复用现成的 `.spatial_distance_cache_matches()` 校验，与基线 CCC 同一套逻辑） |
| `misc$.param$communication` | 不再需要（见下） | 不需要 |
| `DB` / `LR$LRsig` | 仅行为 1 分桶需要 | 不读 |

**xmax 来源统一**：横轴上限一律读 `images$.distance$.parameters$interaction.range`（+10）。该字段由 `computeCellDistance()` 写入、`computeCommunProb()` 校验，两条路径都保证存在。现实现里读 `misc$.param$communication` 的分支删除——几何模式由此彻底摆脱"必须先跑推理"的隐藏依赖，行为 3/4 在**未推理**对象上可用，这正是基线 CCC 的定位。

**自配对**：所有面板剔除 `distance <= 0` 的对角线。基线 CCC 会把 900 个 0 距离自配对画进分布，属无意行为；统一实现修正它，并写进文档作为有意差异。

## 5. 函数名的去向

最终决定（2026-11-01，用户明确要求"只保留实现好的，其他全部不要"）：**只留一个名字，其余全部删除，不设弃用别名。**

| 名字 | 去向 |
|---|---|
| `spatialC2CDistPlot` | 唯一实现，唯一导出 |
| `spatialCommDistPlot` | 直接改名，从未入库，无残留 |
| `spatialCCCDistPlot` | **删除**（含 `R/visualization.R` 中的定义、`man/` 帮助页、`export()`） |
| `communicationDistPlot` | **删除** |
| `communicationDistPlot2` | **删除** |

原先设计的 `.Deprecated()` 转发层被取消：三个旧名在包里完全不存在，外部旧脚本调用会直接报 "could not find function"。这是有意的破坏性变更，比保留一年别名更彻底；`docs/` 下的历史审计文档（`VISUALIZATION_AUDIT.md`、`FUNCTION_INVENTORY.md`）保留原样，作为归档记录。

## 6. 命名风险（知情项）

`spatialC2CDistPlot` 与 `spatialCCCDistPlot` 仅两字之差，过渡期两者同时在 NAMESPACE，极易打混。缓解：CCC 立即弃用（调用即警告并指向新名）；两个 Rd 互相 `@seealso`；下一轮删除 CCC 别名时风险自然消失。若嫌易混，备选是区分度更高的名字——默认按你指定的 `spatialC2CDistPlot` 执行。

## 7. 文档计划（man/）

- 新 `man/spatialC2CDistPlot.Rd`：
  - `\details{}` 开头即第 3 节的四类行为表；
  - 前置条件矩阵（第 4 节），重点标注"几何模式无需推理、缓存可自算"；
  - enriched 计数单位、非零≠显著两段（从现 `spatialCommDistPlot.Rd` 原样迁移）；
  - `@seealso` 链接三个弃用名。
- 三个弃用名的 Rd 改为 stub：一句话说明 + 指向新名 + 参数表保留。

## 8. 测试计划（扩展 tests_dev/test-comm-distance.R，现 21 项）

新增断言：

1. 2×2 矩阵逐格：返回类型（`patchwork`）、面板名、面板数；
2. **行为 3 ≡ 行为 4**：两个调用的面板标题与数据行数完全一致（显式断言，防止将来有人"顺手"分化它们）；
3. **几何模式零依赖推理**：构造一个从未跑 `computeCommunProb` 的对象，`enriched.only = FALSE` 直接出图；
4. 缓存缺失自算路径：删掉 `images$.distance` 后几何模式仍可画，且 `.parameters` 被重建；
5. 三个弃用别名各自发出弃用警告、转发正确（`communicationDistPlot2(signaling.type = "Contact")` 仍返回单面板 `ggplot`）；
6. CCC 与 C2C 同时被导出，CCC 调用触发警告。

回归：现有 12 个测试脚本全量重跑，不得有新增失败。

## 9. 待确认的四个决策

| | 问题 | 建议 | 另一面 |
|---|---|---|---|
| **D6** | `enriched.only` 默认值 | `TRUE`（推断结果图是主用例） | `FALSE` 作默认则函数在任何对象上都能跑，但最常用场景要多传一个参 |
| **D7** | 行为 3/4 下 `split.by` 被忽略，要不要警告 | **警告**（与 `asSeurat()` 的参数门控风格一致） | 静默忽略更安静，但用户排查"为什么图没变"时缺线索 |
| **D8** | 行为 1 中未知 annotation 的记录 | **归入 `Other` 兜底面板**，不再静默丢弃 | 丢弃版更"纯"，但自定义 DB 一接入就会出现看不见的数据洞 |
| **D9** | `spatialCCCDistPlot` 何时转弃用 | **立即**（行为 3/4 完全覆盖它） | 保留双导出一个周期 = 把命名混淆期拉长 |

## 10. 实施与提交

实施顺序：

1. `R/spatial.R`：`spatialCommDistPlot` 改名 `spatialC2CDistPlot`；FALSE 分支改为"缓存校验（复用 `.spatial_distance_cache_matches()`）→ 缺失则 `computeCellDistance()` 现算"；xmax 改读缓存 `.parameters`；D8 若通过则加 `Other` 面板；
2. `R/visualization.R`：`spatialCCCDistPlot` 函数体替换为弃用转发；
3. `NAMESPACE`：新增 `export(spatialC2CDistPlot)`（其余四个名字的 export 行不动）；
4. `man/`：新 Rd + 三个 stub；
5. `tests_dev/test-comm-distance.R`：按第 8 节扩展。

提交：由于 `spatialCommDistPlot` 从未入库，整个整合是**一个提交**，`feat: integrate the communication distance plots into spatialC2CDistPlot`。demo 脚本 `.workbuddy/tmp/commdist/demo_slide_seq_dist.R` 同步改用新名（gitignore，不进 git）。

## 11. 实施结果（2026-11-01，与上述计划的差异）

四个决策按卡片结果落地：**D6 = `FALSE`**（其余 D7 警告 / D8 Other 面板 / D9 立即弃用，与建议一致）。

计划之外补的三处：

1. **几何模式恒定输出两个面板**（`all` + `contact.range`），即使窗口内没有接触对也用 `geom_blank()` 占位——这是为了不破坏基线 `spatialCCCDistPlot(signaling.type = "both")` 的两面板外观；推断模式则相反，空面板会被丢弃。
2. **`enriched.only = FALSE` 的 `split.by` 警告只在显式传参时触发**（`missing(split.by)` 判定），否则默认调用会每次都报警。三个弃用别名走内部的 `.sc_c2c_panels()`，不触发这条警告，避免"弃用警告 + 忽略警告"叠加。
3. **`interaction.range` / `contact.range` / `tol` 三个参数只在几何模式需要现算缓存时使用**，推断模式读缓存、不读它们。

实现后的结构：

| 函数 | 职责 |
|---|---|
| `spatialC2CDistPlot` | 唯一公开入口，做参数门控后转调 `.sc_c2c_panels()` |
| `.sc_c2c_panels` | 备好缓存 + 出面板（公开函数与三个别名共用） |
| `.sc_c2c_inference_cache` / `.sc_c2c_geometry_cache` | 两条取缓存路径 |
| `.sc_comm_distance_frame` / `.sc_comm_dist_panels` / `.sc_comm_dist_panel` | 取数 - 分桶 - 单面板 |

验证：`tests_dev/test-comm-distance.R` **26 项**全过（2×2 逐格、3≡4、未推理对象可画、缓存自算、显式 range 生效、"四个旧名全部不存在"与"只剩一个 DistPlot 导出"）。另外把 `tests_dev/test-SpatialCellChat-spatial-image.R`（3 处）和 `tests_dev/audit_visualization_findings.R`（1 处）里的 `spatialCCCDistPlot` 调用改为新名。回归 seurat-bridge 41、projectData 19、computeCommunProb 61、computeAvgCommunProb 162、filterCommunication 22、normalizeData 9、class 22、accessors 22、cpp-prob-layer 20、helpers 16、pre-inference 15、computeCellDistance 8、grid-visualization 2 全绿；四个 Rd `checkRd` 通过；14 个 R 文件 parse 通过。

Slide-seq 实图（900 细胞窗口 × 120 LR，`scale.distance = 9`）用新 API 重出三张：几何图在**推理之前**就画出来了，验证了"几何模式彻底不需要 inference"这一目标。

## 12. 验收标准

- `spatialCommDistPlot` 名字在 `R/` 下消失；
- 四类行为各有测试且全绿；行为 3 ≡ 行为 4 有显式断言；
- 几何模式在**未推理**对象上可用（对基线 CCC 能力的保真）；
- 三个弃用别名签名不变、警告正确；
- 12 个回归脚本无新增失败；全部 R 文件 parse 通过；`R CMD INSTALL` 仍只报那 4 个 undefined exports。
