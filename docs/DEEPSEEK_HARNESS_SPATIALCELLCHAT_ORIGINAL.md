我完整看了 DeepSeek Harness 根目录、`AGENTS.md`、架构、测试策略、CI 门禁、Agent Notes、skills 和 hooks。结论：

> 最值得借鉴的不是 “Everything is a plugin” 本身，而是它把 **契约、决策、验证、生命周期和门禁** 做成了可执行的工程系统。

## 一、最值得借鉴的实践

### 1. 把“权威事实源”明确下来

DeepSeek 有一条很强的原则：

> Model-visible ⟺ logged

凡是进入模型上下文的内容，都必须能从持久化 session log 重建。不能一份存在内存、一份存在 UI、一份存在缓存，最后彼此漂移。

对应文件：

- [architecture.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/architecture.md)
- [AGENTS.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/AGENTS.md)

对 SpatialCellChat 的可迁移版本：

| 类型 | 建议权威来源 |
|---|---|
| 表达矩阵 | `assay()` |
| 空间信息 | `images()` |
| LR 元数据 | `LR()` |
| 通讯结果 | `communication()` / `net` / `netP` |
| 对象参数 | `params()` / `misc$.param` |
| 修改历史 | `misc$.log` |
| 绘图数据 | 每次由对象状态重新生成，不持久化多份副本 |
| 缓存 | 必须绑定输入参数和坐标版本 |

这对你当前重构特别重要：

> 不要让 `spatialFeaturePlot()` 自己拼一套数据访问规则，同时 `plotly_spatialLRpairPlot()` 再拼另一套。所有图都应该从同一个 accessor 和标准化 plotting table 派生。

---

### 2. “能力接缝”设计：定义、提供者、消费者分开

DeepSeek 把一个可替换能力定义为三个角色：

```text
Service Definition
Service Provider
Consumer
```

例如文件系统、Shell、LLM、Session、Subagent 都通过稳定接口连接，而不是让某个消费者反向决定整个接口。

来源：

- [architecture.md — Capability seams](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/architecture.md)
- [packages/AGENTS.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/packages/AGENTS.md)

对应到 SpatialCellChat：

```text
数据访问契约
  ├── assay / images / meta / LR / communication
  └── accessor 层校验

计算提供者
  ├── Gi / Lee
  ├── distance / contact
  ├── communication aggregation
  └── SparseChatArray operations

绘图消费者
  ├── ggplot2
  ├── Plotly 2D
  ├── Plotly 3D
  └── Plotly 3D stack
```

关键原则：

- 绘图函数不应该直接决定对象内部数据结构；
- 统计函数不应该依赖某一个绘图函数的数据整理方式；
- Plotly 不应该重新实现一套与 ggplot2 不同的坐标语义；
- accessor 的接口应服务当前所有消费者，而不是为某一个函数临时定制。

这正好支持你已经开始的 `assay()`、`images()`、`communication()` 访问层重构。

---

### 3. 在真正执行操作的边界上强制规则

DeepSeek 明确反对只在表面 schema、wrapper 或 UI 层做限制：

> 规则必须在真正执行决定的操作中 enforced。

来源：

- [packages/AGENTS.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/packages/AGENTS.md)

对 SpatialCellChat 的典型应用：

错误做法：

```r
spatialGiPlot(..., do.binary = FALSE)
# 只在绘图层解释参数
```

更可靠的做法：

```text
spatialGiPlot()
  -> 计算函数
      -> 明确处理 do.binary
          -> 结果带上实际统计语义
              -> 绘图只负责显示
```

同样适用于：

- `input.assay` 与 `normalize`；
- Gi 的二值化；
- Lee 的负值截断；
- 距离缓存失效；
- 坐标系转换；
- `SparseChatArray` 的稀疏聚合；
- Plotly 3D stack 的全局颜色范围。

验证不能只测“参数被接受”，要测“直接绕过上层调用后，底层仍然拒绝非法状态”。

---

### 4. 语义不变量，而不是“slot 存在性”检查

DeepSeek 的 package invariant 不是检查：

```text
某个 service 存不存在
某个字段是不是非 NULL
某个函数是否被注册
```

而是检查真实的运行关系：

- 事件是否对应真实状态；
- 注册对象是否在 dispose 后被清除；
- durable log 是否能重放出当前状态；
- provider 和 consumer 是否真实连通。

来源：

- [packages/AGENTS.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/packages/AGENTS.md)
- [testing.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/testing.md)

SpatialCellChat 应把 `validateSpatialCellChat()` 的重点放在语义不变量：

```text
assay$norm 的列名 == meta 的行名 == images$coordinates 的行名
idents 覆盖全部细胞且没有非法 NA
SparseChatArray 每一层维度一致
net 与 netP 的 layer names 和元数据一致
images$.distance 的参数匹配当前坐标和阈值
raster 坐标能通过 transform 对齐分析坐标
```

这比单纯检查：

```r
!is.null(object@images)
!is.null(object@net)
```

更有价值。

---

### 5. “验证世界”，不要相信 Agent 自己的报告

DeepSeek 测试规则里最值得借鉴的一句是：

> Verify the world, not the self-report.

也就是不要只检查 Agent 说“完成了”，而要：

- 重新运行命令；
- 重新读取文件；
- 检查外部状态；
- 检查实际输出；
- 检查未修改文件是否保持不变。

来源：

- [testing.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/testing.md)

对 SpatialCellChat：

不够强的测试：

```r
expect_s3_class(plot, "ggplot")
```

更强的测试：

```r
p <- spatialFeaturePlot(chat, features = "CXCL12")
built <- ggplot2::ggplot_build(p)

expect_equal(built$data[[1]]$x, expected_x)
expect_equal(built$data[[1]]$y, expected_y)
expect_true(any(built$data[[1]]$colour == expected_color))
```

Plotly 也应该检查：

```text
trace 数量
trace 类型
colorbar 数量
全局 cmin/cmax
z 轴 ticktext
hover 文本
空 highlight 是否产生空 trace
```

3D stack 原型已经验证了“可以构建 Plotly 对象”，但生产验证还需要检查：

- 单一 colorbar；
- 各层共用色域；
- z 轴层名；
- 空值和全零层；
- 大规模 spot 的渲染路径。

---

### 6. 测试分层，而不是所有改动都跑全套

DeepSeek 的测试层级比较成熟：

```text
unit
coverage
real-api e2e
snapshot
browser snapshot
built-artifact smoke
```

来源：

- [testing.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/testing.md)
- [dsh-pre-push-checks](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/skills/dsh-pre-push-checks/SKILL.md)

对 SpatialCellChat 可以采用较轻的版本：

| 改动类型 | 最小验证 |
|---|---|
| accessor / class | `tests_dev/test-SpatialCellChat-class.R` |
| normalization | `tests_dev/test-normalizeData.R` |
| spatial image | `test-SpatialCellChat-spatial-image.R` |
| distance | `test-computeCellDistance.R` |
| Gi / Lee | 指标边界测试 + plot smoke |
| ggplot2 | `ggplot_build()` + 实际打印 |
| Plotly | `plotly_build()` + trace/scene 检查 |
| 3D stack | 多层、空层、统一色域和 z tick 检查 |
| 包公开 API | 从 package source 或安装包实际调用 |
| 文档和 Notes | 文档链接、格式、路径检查 |

核心不是照搬测试数量，而是建立：

> 一个改动对应一条能失败的验证路径。

---

### 7. 本地 hooks 要快，完整门禁交给 CI

DeepSeek 的 `lefthook.yml` 有明确分工：

本地 hook 只做：

- staged lint；
- whitespace；
- vendor manifest；
- translation pairing；
- 基础 typecheck。

完整测试、coverage、构建、snapshot、兼容性矩阵由 CI 承担。

来源：

- [lefthook.yml](https://github.com/deepseek-ai/deepseek-harness/blob/master/lefthook.yml)
- [development.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/development.md)

对 SpatialCellChat 的建议：

### 本地提交前

只做：

```text
R 语法解析
NAMESPACE / roxygen 基本检查
新增测试脚本可 source
changed file 相关 smoke test
git diff --check
```

### CI 或人工完整检查

再做：

```text
devtools::check()
全部 tests_dev
大型真实对象 fixture
plotly / ggplot2 全部可视化 smoke
性能 benchmark
```

不要把大型 psoriasis fixture、全包 check 或所有可视化都挂到每一次提交上，否则开发反馈会变慢，最终大家会绕过 hook。

---

### 8. 把门禁写成依赖图，而不是一串命令

DeepSeek 的 `scripts/run-gates.ts` 不是简单地：

```text
先跑 A，再跑 B，再跑 C
```

而是把每个 gate 定义为：

```text
id
label
command
needs
allowFailure
```

然后根据依赖关系并行执行，并限制并发数。

来源：

- [scripts/run-gates.ts](https://github.com/deepseek-ai/deepseek-harness/blob/master/scripts/run-gates.ts)

SpatialCellChat 可以借鉴成：

```text
class-parse
  └── class-tests

accessor-parse
  └── accessor-tests

spatial-parse
  └── spatial-image-tests
  └── distance-tests

visualization-parse
  └── ggplot-smoke
  └── plotly-smoke

all-targeted-tests
  └── devtools-check
```

独立分支可以并行，但应限制并发，因为：

- R 进程启动成本高；
- 大矩阵会占用大量内存；
- Plotly / raster 测试可能产生大量临时对象；
- `future` 并行和测试并行叠加会造成过量并发。

这和你现有的 `my_future_lapply()`、`setEnvironment()` 设计可以互补：一个控制包内计算，一个控制测试门禁。

---

### 9. Agent Notes 不只是记录结论，还要记录“放弃了什么”

DeepSeek Agent Notes 的核心不在目录，而在：

```text
Problem
Decision / Proposal
Alternatives considered
Consequences
```

来源：

- [Agent Notes README](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/notes/README.md)

你已经建立了轻量版三目录。下一步最值得借鉴的两个细节是：

1. 增加类别子目录：

```text
.agents/notes/
├── implemented/
│   ├── architecture/
│   ├── process/
│   └── testing/
├── proposed/
│   ├── feature/
│   ├── architecture/
│   └── testing/
└── rejected/
    ├── architecture/
    └── simplification/
```

2. 增加一个非常小的格式检查器，检查：

```text
文件路径状态 == Status
文件名包含日期
implemented 必须有 Decision
proposed 必须有 Proposal
所有 Note 必须有 Alternatives considered
```

不建议现在引入：

- 双语 sidecar；
- hash manifest；
- 完整 archived 冻结系统；
- 强制每个机械改动写 Note。

---

### 10. 只保留一个事实的权威归属

DeepSeek 的文档规则要求“一条事实只有一个家”：

| 内容 | 所属位置 |
|---|---|
| Agent 必须遵守的规则 | 根 `AGENTS.md` |
| 子目录特殊规则 | 子目录 `AGENTS.md` |
| 架构地图 | `docs/architecture.md` |
| 决策理由 | Agent Notes |
| 事故过程 | postmortem |
| 操作教程 | cookbook |
| API 契约 | roxygen / package README |
| 测试策略 | testing 文档或测试 Note |
| 生成内容 | generator source |

来源：

- [docs/AGENTS.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/AGENTS.md)

对 SpatialCellChat，当前存在一些事实在多个地方重复：

- `docs/DATA_STRUCTURE.md`
- `docs/REFACTORING_PLAN.md`
- `R/SpatialCellChat_class.R` roxygen；
- `.agents/notes/implemented/`；
- 测试脚本注释。

建议按这个原则收敛：

```text
最终 schema 的正式定义 → DATA_STRUCTURE.md
重构阶段顺序 → REFACTORING_PLAN.md
为什么这样设计 → implemented Agent Note
函数参数和返回契约 → Roxygen
可执行不变量 → tests_dev
工作方式 → AGENTS.md
```

Notes 不应变成第二份完整数据结构文档，只保留决策理由和代价。

---

## 二、特别适合你的三个“闪光点”

### A. “显式优于隐式”

DeepSeek 明确要求：

- 默认值在 `resolve()` 阶段显式解析；
- 不把默认值藏在执行函数中；
- 配置错误尽早失败；
- 不静默跳过缺失引用；
- 不在多个层重复推断。

这对你的 API 很直接：

```r
createSpatialCellChat(
  input.assay = "raw",
  normalize = TRUE
)
```

比：

```r
createSpatialCellChat(data)
```

然后内部猜测输入类型更可靠。

同样地，绘图 API 应显式表达：

```r
coord.system = "analysis"
value.scale = "global"
z.axis.space = 1
image = TRUE
```

而不是让每个函数偷偷使用不同默认行为。

### B. “简化必须有消费者证据”

DeepSeek 的 simplification skill 不接受：

```text
这段代码看起来复杂，所以删掉
```

它要求先分类所有消费者：

```text
生产代码
测试
示例
脚本
文档
动态加载
配置入口
```

再判断是否可以删除。

来源：

- [dsh-find-simplifications](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/skills/dsh-find-simplifications/SKILL.md)

这特别适合 SpatialCellChat，因为旧 API 和旧 slot 很多。删除 `@data.signaling`、`my_as_sparse3Darray` 或旧 helper 前，必须检查：

- `R/` 生产调用；
- `tests_dev/`；
- tutorial；
- man 文档；
- NAMESPACE；
- 动态字符串访问；
- 示例数据脚本。

你已经在 Notes 中记录“不保留旧 schema 兼容层”，但每次实际删除仍应按消费者证据执行。

### C. “防御逻辑要围绕生命周期和所有权”

DeepSeek 的 defensive patterns 不是泛泛的安全建议，而是针对真实 bug 类：

- 独立结果不能嵌套丢失；
- 异步状态不等于同步状态；
- dispose 必须真正 quiescent；
- callback 异常不能破坏 dispatcher；
- 缓存必须有明确 owner；
- 不可信输出不能继承环境和可预测路径。

来源：

- [defensive-patterns.md](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/defensive-patterns.md)

SpatialCellChat 可迁移的对应问题：

| DeepSeek 问题 | SpatialCellChat 对应物 |
|---|---|
| async lifecycle | `future` workers、parallel computation |
| durable state | `misc$.log`、派生结果 |
| cache ownership | `images$.distance`、raster cache |
| callback cleanup | Shiny Plotly observer |
| bounded output | feature 数量、Plotly trace、hover 文本 |
| independent outcome | Gi score、p-value、highlight 状态分别处理 |
| quiescent disposal | future plan、临时文件、Plotly/Shiny 连接清理 |

尤其是你审计出的 Shiny observer 累积问题，正是“生命周期 owner 不清晰”的典型。

---

## 三、暂时不要照搬的内容

### 1. “Everything is a plugin”

对 DeepSeek Harness，这是核心架构；对 SpatialCellChat 不应直接复制成几十个 R package 或大量 plugin。

可以借鉴：

- 稳定接口；
- provider / consumer 分离；
- 可替换的计算后端；
- 明确的扩展边界。

不必复制：

- Cordis；
- 全局事件总线；
- 每个绘图函数一个 plugin；
- 为了抽象而抽象的运行时容器。

### 2. 每文件 100% coverage

DeepSeek 的 per-file 100% coverage 适合它的 TypeScript harness，但对 SpatialCellChat：

- R 的大量分支依赖外部对象；
- 可视化边界情况多；
- 大型真实 fixture 成本高；
- 纯行覆盖不等于统计结果正确。

建议借鉴它的思想：

> 未覆盖代码可能是死代码。

但不建议机械设置全包 100%。更适合按核心模块设置：

- class/accessor：高覆盖；
- normalization/distance：高覆盖；
- Gi/Lee：关键边界和数值基准；
- ggplot/Plotly：行为 smoke 与结构断言；
- 旧兼容路径：不存在即测试失败。

### 3. 双语文档和完整 archive 机制

DeepSeek 的双语 pairing、sidecar、归档哈希和冻结 manifest 很强，但对当前项目投入过大。

你现在的轻量版：

```text
implemented/
proposed/
rejected/
```

已经足够。优先增加：

- 类别子目录；
- status/路径检查；
- 相对 Markdown 链接；
- 一个事实一个归属。

### 4. 大型 CI failover 和多平台 runner 管理

DeepSeek 的 CI 包含：

- 多 Node 版本；
- 大型 self-hosted runner；
- failover；
- 浏览器 snapshot；
- real API；
- 多套 build consumer。

这些是产品规模和团队规模的结果，不适合当前直接复制。

---

## 四、给 SpatialCellChat 的推荐落地顺序

### P0：现在就可以采用

1. 在项目 `AGENTS.md` 写入少量硬规则：
   - 最终 11-slot schema 是唯一目标；
   - 公开 API 不保留旧 alias；
   - 修改核心模块必须运行对应 targeted test；
   - Notes 与当前代码冲突时，以代码和权威 schema 为准。
2. 给 Notes 增加 `architecture/process/testing/feature` 子目录。
3. 为 `validateSpatialCellChat()` 增加语义不变量，而不是只做 slot 检查。
4. 统一绘图数据预处理和坐标语义。
5. 测试中检查实际绘图数据、Plotly trace 和输出结构。

### P1：下一阶段

1. 写一个轻量 `run-checks.R` 或 Rscript gate runner。
2. 每个 gate 有：
   - id；
   - command；
   - needs；
   - allow failure；
   - 适用文件范围。
3. 默认只运行 changed-surface checks。
4. 大型测试和完整 `devtools::check()` 放到 CI 或显式命令。

### P2：成熟后再考虑

1. Notes 状态迁移检查；
2. archived 目录；
3. 自动检查 Note 与代码路径的一致性；
4. Plotly/ggplot2 输出 snapshot；
5. 真实大型空间数据的 built/package smoke。

## 最终判断

DeepSeek Harness 最有价值的工程经验可以浓缩成五句话：

1. **每个事实只有一个权威来源。**
2. **每条规则必须在真正执行的边界上生效。**
3. **每个抽象都要有当前消费者和明确 owner。**
4. **验证外部世界，不验证 Agent 自己的描述。**
5. **本地快速反馈，CI 承担完整矩阵。**

对 SpatialCellChat，优先采用这五条，比直接复制 Cordis、插件树、复杂 CI 或双语文档基础设施更有收益。