# DeepSeek Harness Engineering：Agent Coding 的成功经验与可借鉴实践

> 研究备份稿｜面向后续微信推文整理
>
> 本文基于 DeepSeek Harness 仓库公开的 `AGENTS.md`、架构文档、测试策略、Agent Notes、skills、CI 门禁和 Git hooks 整理。文中示例均使用电商订单、文件处理器、测试系统和通用数据分析软件，避免引用未公开项目。

## 一句话结论

DeepSeek Harness 最值得借鉴的地方，不是它使用了某个特定模型，也不是“所有东西都是插件”这句架构口号，而是它把 Agent Coding 需要的几件事做成了一个可执行系统：明确的契约、稳定的事实来源、可追溯的决策记录、真实世界验证、分层测试、受控门禁和清晰的生命周期所有权。

这意味着 Agent 不再只是“会改代码的聊天机器人”，而是被放进了一套能够约束、验证和复盘的软件工程流程中。

## 一、先解释几个容易混淆的设计概念

### 1. Consumer：消费者

Consumer 可以理解为“使用某个能力的一方”。例如，一个订单系统需要发送短信，它就是短信能力的消费者；一个网页需要读取用户信息，它就是用户信息服务的消费者；一个绘图函数需要读取分析结果，它也是数据服务的消费者。

Consumer 不应该决定底层能力的全部实现细节。它只应该依赖一份稳定的接口，例如：

```text
sendMessage(to, content)
```

至于消息由短信、邮件还是企业微信发送，应该由 Provider 决定。

### 2. Provider：提供者

Provider 是“真正提供能力的一方”。同一个接口可以有多个 Provider：

```text
MessageProvider
├── SmsProvider
├── EmailProvider
└── FeishuProvider
```

消费者只依赖 `MessageProvider`，不直接绑定某个具体厂商。这样替换服务商时，不需要重写业务逻辑。

### 3. Service Definition：服务定义

Service Definition 是消费者和提供者之间的契约，规定：

- 能调用什么操作；
- 参数是什么；
- 返回什么；
- 失败如何表达；
- 生命周期如何结束。

它不是一个“把所有可能功能都塞进去”的大接口，而是当前消费者真正需要的最小稳定边界。

### 4. Capability Seam：能力接缝

DeepSeek Harness 把 Service Definition、Provider 和 Consumer 组成的可替换边界称为 capability seam，可以翻译成“能力接缝”。

接缝的意义是：系统内部可以替换实现，但外部调用方式保持稳定。

例如文件读取能力：

```text
FileService 定义读取文件的契约
LocalFileProvider 从本地磁盘读取
SandboxFileProvider 从沙箱读取
ToolConsumer 把读取能力暴露给 Agent
```

这比让每个工具自己读取文件、自己处理权限、自己处理错误更容易维护。

### 5. Source of Truth：权威事实源

Source of Truth 指“某类事实到底以哪里为准”。如果订单状态同时存在数据库、缓存、前端状态和日志中，就必须明确哪一个是权威来源，其余都是派生视图。

否则系统很容易出现：数据库显示“已支付”，缓存显示“待支付”，前端又显示“处理中”。

### 6. Invariant：不变量

Invariant 是系统在任何合法状态下都必须保持的关系。

例如购物车系统中的不变量：

```text
购物车中的商品数量不能小于 0
订单中的商品必须来自当前商品目录
已支付订单不能重新变成待支付
订单总额必须等于明细之和
```

不变量不是“某个字段存在”，而是字段之间的真实关系。

### 7. Lifecycle Owner：生命周期所有者

一个异步任务、缓存、临时文件或监听器都应该有明确的 owner：谁创建、谁负责结束、谁负责清理、谁决定失败。

如果一个任务由 A 创建、由 B 取消、由 C 记录结果、由 D 清理资源，就容易出现孤儿进程、重复回调和状态悬挂。

### 8. Gate Dependency DAG：门禁依赖图

质量检查不一定是一串线性命令。很多检查可以并行，但有些检查依赖前置产物。

例如：

```text
源代码解析
├── 单元测试
├── 文档检查
└── 静态分析

构建
├── 安装包 smoke test
└── 发布产物检查
```

这就是一个有依赖关系的有向无环图。它比“所有命令顺序执行”更快，也比“全部同时启动”更可靠。

## 二、DeepSeek Harness 最值得借鉴的实践

### 1. 把系统真正看到的内容记录下来

DeepSeek Harness 的架构文档有一个非常重要的原则：凡是进入模型请求的内容，都必须能够从 session log 重新构建。模型看到的上下文不能只存在于某个暂时的内存变量、UI 状态或未记录的缓存中。

这解决的是 Agent 系统最容易被忽视的问题：如果一个 Agent 得出了错误结论，开发者能不能复现它当时到底看到了什么？

通用例子：文件处理 Agent 接到任务“整理这个目录”。如果它读到了目录列表、文件内容、用户补充指令和工具返回结果，这些内容都应该能从会话记录中恢复。否则下一次复现时，可能只有最终回复，却没有当时的输入环境。

可借鉴原则：

- 模型可见内容必须可重建；
- 工具调用和工具结果必须进入可追溯记录；
- UI 展示、压缩摘要和缓存都应从权威记录派生；
- 不要让多个地方分别维护“差不多一样”的状态。

### 2. 新能力应该挂在扩展点上，而不是随意修改主循环

DeepSeek Harness 使用插件和事件扩展能力。模型适配器、工具注册、会话记录、Agent Loop、沙箱和系统提示词都被放在可组合的扩展结构中。

这并不意味着普通项目必须复制 Cordis，而是说明一个很重要的设计思想：

> 新功能应该接入已经定义好的扩展点，而不是每次都直接改核心执行循环。

通用例子：一个代码 Agent 想新增“生成测试报告”能力。比较稳妥的方式是新增一个报告工具，接入已有工具注册和结果记录机制；不应该为了一个新工具，直接修改 Agent 主循环、消息格式和所有调用方。

这样做可以降低 blast radius，也就是降低一次改动影响整个系统的范围。

### 3. 接口由所有当前消费者共同决定

一个服务接口不应为了某个单独消费者加入大量特殊方法。

例如，通知服务当前只有两个消费者：注册流程和密码重置流程。如果只有密码重置流程需要“发送带验证码的模板消息”，不一定要把模板逻辑硬塞进所有通知 Provider 的公共接口。可以把这个能力作为私有的、局部的适配器或调用闭包提供给密码重置流程。

DeepSeek 的实践强调：

- 先列出当前真实消费者；
- 找出它们共有的最小能力；
- 把消费者专属逻辑留在消费者或专用 Provider 中；
- 不为猜测中的未来需求提前扩大公共 API。

这是一种非常有效的反过度抽象方法。

### 4. 显式优于隐式

DeepSeek 的规则要求，在能力边界上显式解析请求和默认值，而不是把默认行为隐藏在执行函数深处。

例如一个报告导出服务，不应该根据文件名猜测用户想要 PDF 还是 CSV：

```text
exportReport(request)
```

内部偷偷猜格式，会让调用方和测试都难以理解。

更清晰的形式是：

```text
resolveExportSpec(request)
  -> format = "pdf"
  -> includeCharts = true
  -> locale = "zh-CN"

runExport(spec)
```

这样“用户请求是什么”和“最终执行规格是什么”分成两个阶段。默认值、冲突配置和非法配置都有明确位置可以验证。

适用场景包括：

- 配置解析；
- 文件格式选择；
- 权限模式；
- 重试策略；
- 缓存策略；
- 数据导入模式；
- 可视化颜色范围和坐标系。

### 5. 在真正执行的位置强制规则

只在 UI、schema 或 wrapper 层声明一条规则还不够。真正执行操作的函数必须再次保证它不会被其他调用路径绕过。

例如支付系统要求“库存不足时不能创建支付订单”。如果这个检查只写在网页按钮上，API 调用者仍然可以绕过它。规则应该进入创建支付订单的核心操作，并对直接调用和间接调用都测试。

这条原则也适合 Agent 工具：

- 工具 schema 可以限制参数；
- Prompt 可以提醒 Agent；
- 但最终执行器仍然必须拒绝非法路径。

否则只要出现另一个调用入口，之前的安全规则就可能失效。

### 6. 语义不变量比字段存在更重要

一个优秀的 invariant 检查不会满足于“对象里有这个字段”。它会检查对象内部的事实关系。

通用数据分析软件可以有这样的不变量：

```text
数据矩阵的列名必须和元数据的行名一致
分组标签必须覆盖全部观测
结果矩阵的维度必须和输入样本集合一致
缓存必须对应当前数据版本和参数
```

这些关系比下面的检查更有价值：

```text
!is.null(result)
"metadata" %in% names(object)
```

DeepSeek 的 package invariant 也强调，不能通过“服务存在”“方法注册了”“插件加载了”来假装系统关系正确，必须检查真实运行数据或事件关系。

### 7. 验证真实世界，不要相信 Agent 自己的总结

这是 DeepSeek 测试策略中最值得传播的一条经验。

弱验证是：

```text
Agent 说“文件已经生成”
```

强验证是：

```text
重新检查文件是否存在
读取文件内容
验证格式是否正确
检查旧文件是否被意外修改
```

通用例子：Agent 被要求把 CSV 转换为 JSON。测试不应该只检查 Agent 回复中出现了“成功”，而应该：

1. 在临时目录创建真实 CSV；
2. 启动真正的命令或工具入口；
3. 读取生成的 JSON；
4. 检查字段、记录数和编码；
5. 确认源文件没有被修改；
6. 测试输入错误时命令确实失败。

这类测试能捕获“单元测试全绿，但产品入口坏了”的问题。

### 8. 测试按照证据层级分层

DeepSeek Harness 把测试划分为多个层次：

```text
Unit Test
Coverage Gate
Real API E2E
Snapshot Test
Browser Snapshot
Built Artifact Smoke Test
```

每一层回答不同问题：

| 层级 | 主要问题 |
|---|---|
| 单元测试 | 这个函数在局部输入下是否正确？ |
| Coverage | 是否存在未执行的实现分支？ |
| E2E | 真实外部服务是否能完成任务？ |
| Snapshot | 用户或模型可见输出是否发生意外变化？ |
| Browser Snapshot | 实际界面是否变化？ |
| Built Smoke | 构建出来的真正产物能否运行？ |

这比所有改动都运行整个测试套件更有效。正确做法是：根据改动表面选择最小但足够的证据；只有跨越整个系统时才运行完整门禁。

### 9. 本地反馈要快，CI 承担完整检查

DeepSeek 的 Git hooks 被有意控制在较小范围：staged lint、格式检查、配对文件检查、基础 typecheck。完整 coverage、构建、snapshot、兼容性矩阵和发布产物检查交给 CI。

这背后的工程判断很现实：如果提交前检查每次都运行几十分钟，开发者最终会关闭它们。

一个实用的分工是：

```text
本地 hook：秒级或分钟级反馈
CI：完整、昂贵、跨平台、发布前检查
```

本地检查应关注“明显错误能否尽早发现”；CI 关注“最终产品是否完整可靠”。

### 10. 把门禁做成有依赖关系的调度器

DeepSeek 的 `run-gates.ts` 为每个检查定义 id、命令、依赖和失败策略，然后根据依赖图调度并发执行。

例如：

```text
parse
├── unit-test
├── lint
└── docs

build
├── built-smoke
└── package-hygiene
```

独立检查可以并行；依赖构建产物的检查必须等待 build；高内存任务需要限制并发。

这比简单写一个巨大的 shell 脚本更容易解释失败原因，也更容易在 CI 中拆分成不同 lane。

### 11. Agent Notes 保存“为什么”，而不是重复代码

DeepSeek 的 Agent Notes 用来记录代码和普通文档无法完整表达的内容：

- 为什么采用当前方案；
- 放弃了哪些替代方案；
- 当前决策带来什么代价；
- 未来什么条件下可以重新考虑。

它还使用 `proposed`、`implemented`、`rejected` 和 `archived` 管理生命周期，并要求已实施 Note 与当前代码保持一致。

这类记录可以避免团队重复争论同一个问题，也能防止后续维护者因为看不懂某个约束而把它删掉。

但不应把每个 typo 修复都写成设计文档。适合写 Note 的通常是：

- 数据结构；
- 公共 API；
- 算法选择；
- 测试策略；
- 依赖选择；
- 兼容性取舍；
- 明确拒绝的设计。

### 12. 简化必须有消费者证据

DeepSeek 的 simplification 工作流不会因为一段代码“看起来复杂”就建议删除。它会先区分：

```text
生产代码调用
测试调用
示例调用
文档引用
动态加载
配置入口
```

然后才判断一个 API、helper、配置项或兼容分支是否真的没有价值。

这对于大型重构尤其重要。删除一个看似无用的函数，可能会破坏示例、插件加载、动态字符串调用或发布脚本。

可靠的简化应回答：

1. 当前谁在使用它？
2. 它保护了哪个契约？
3. 删除后真正减少了什么？
4. 测试、文档、配置和生成文件是否一起清理？
5. 是否只是把复杂度移动到另一个 wrapper？

### 13. 防御代码围绕生命周期和所有权

DeepSeek 的 defensive patterns 来自实际 bug 经验，特别强调：

- 不同结果事实要独立表达，例如超时、信号和退出码不能互相覆盖；
- 异步状态不能简单当成同步状态；
- dispose 不只是“发出取消请求”，而是要等任务真正停止；
- 一个 listener 抛错不能破坏整个 dispatcher；
- 外部命令不能继承不必要的密钥环境变量；
- 临时文件不能使用可预测且公开可写的路径。

这些原则不只适用于 Agent Harness。

例如一个后台文件压缩任务：

```text
任务创建者拥有任务
任务控制器拥有取消状态
输出文件由任务控制器负责提交
任务结束后才发布完成事件
失败、取消、超时分别记录
```

如果压缩任务还没结束，系统就先把“完成”状态发布出去，后面的 UI、缓存和通知都会看到错误事实。

## 三、哪些做法不应该直接照搬

### 1. 不必因为 Harness 使用插件，就把普通应用拆成几十个插件

可以借鉴能力边界和可替换 Provider，但不必复制 Cordis、全局事件系统或复杂运行时容器。抽象应该由当前消费者和实际替换需求驱动。

### 2. 不必机械追求每个文件 100% coverage

高覆盖率可以帮助发现死代码，但行覆盖率不等于行为正确。数据分析、图形渲染和外部服务更需要数值基准、边界测试、真实入口 smoke 和结果断言。

### 3. 不必一开始复制完整 CI 基础设施

多版本 Node、浏览器 snapshot、self-hosted failover、真实 API 矩阵和大规模 runner 管理，都是大型项目的结果。小项目应先建立少量高价值门禁，再逐步扩展。

### 4. 不必照搬双语和归档元数据系统

DeepSeek 的双语配对、sidecar 文件、哈希冻结和自动归档很成熟，但维护成本也很高。小型项目先实现：

```text
决策记录
状态分类
必需章节
相对链接
当前代码一致性
```

已经能覆盖大部分收益。

## 四、对 Agent Coding 最实用的落地清单

### 规则层

- 明确项目的权威数据结构；
- 明确哪些旧 API 不再兼容；
- 明确修改后必须运行哪些检查；
- 把必须每次执行的规则放在 `AGENTS.md`。

### 设计层

- 为每个能力识别 Definition、Provider 和 Consumer；
- 让公共接口由多个真实消费者共同决定；
- 把默认值和配置冲突显式解析；
- 为缓存、异步任务和派生状态指定 owner；
- 为关键对象定义语义不变量。

### 测试层

- 先写能失败的回归测试；
- 测试真实入口，而不是只测试手工拼装的内部对象；
- 重新读取外部输出验证结果；
- 对用户可见输出使用 snapshot 或结构化断言；
- 按改动表面选择最小充分检查。

### 记忆层

- 重大决策写 Agent Note；
- Note 记录替代方案和代价；
- 已实施 Note 使用当前时态描述真实代码；
- 已过时 proposal 转为 rejected，而不是伪装成 implemented；
- 不让 `AGENTS.md`、README、设计文档和 Note 重复承担同一事实。

### 简化层

- 删除前搜索所有消费者；
- 区分生产调用和测试/文档调用；
- 优先删除没有当前 owner 的抽象；
- 不要用 wrapper 把复杂度从一个地方搬到另一个地方；
- 把小型清理留给 TODO，把结构性取舍写成 Note。

## 五、最后的判断

DeepSeek Harness 的真正经验可以浓缩成五句话：

1. 每类事实只有一个权威来源。
2. 每条规则必须在真正执行的位置生效。
3. 每个抽象都要有当前消费者和明确 owner。
4. 验证外部世界，不验证 Agent 自己的总结。
5. 本地保持快速反馈，CI 承担完整矩阵。

Agent Coding 的质量，最终不取决于 Agent 能不能一次生成大量代码，而取决于系统能不能让错误尽快暴露、让正确结果可复现、让设计理由可以延续。

这也是 Harness Engineering 和普通“给模型加几个工具”之间的区别：前者是在设计一个可观察、可约束、可验证的软件生产系统。

## 参考资料

- [DeepSeek Harness README](https://github.com/deepseek-ai/deepseek-harness/blob/master/README.md)
- [AGENTS.md：项目工作规则与工程约束](https://github.com/deepseek-ai/deepseek-harness/blob/master/AGENTS.md)
- [Architecture：插件、事件、会话日志与能力接缝](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/architecture.md)
- [Testing Policy：分层测试与真实入口验证](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/testing.md)
- [Defensive Patterns：生命周期、并发和安全边界](https://github.com/deepseek-ai/deepseek-harness/blob/master/docs/defensive-patterns.md)
- [Agent Notes README：决策记录与生命周期](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/notes/README.md)
- [Packages AGENTS.md：Consumer、Provider、Invariant 与真实组合测试](https://github.com/deepseek-ai/deepseek-harness/blob/master/packages/AGENTS.md)
- [dsh-pre-push-checks：最小充分验证](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/skills/dsh-pre-push-checks/SKILL.md)
- [dsh-find-simplifications：基于消费者证据的简化](https://github.com/deepseek-ai/deepseek-harness/blob/master/.agents/skills/dsh-find-simplifications/SKILL.md)
- [run-gates.ts：带依赖关系和并发上限的门禁调度](https://github.com/deepseek-ai/deepseek-harness/blob/master/scripts/run-gates.ts)
- [lefthook.yml：快速本地 hooks 与 CI 分工](https://github.com/deepseek-ai/deepseek-harness/blob/master/lefthook.yml)
