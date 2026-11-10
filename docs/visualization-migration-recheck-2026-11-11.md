# 可视化迁移复测：当前生产代码与三对象

本轮从生产代码提交 `449c181` 重绘，运行时 `R/` 和 `src/` 没有未提交修改；测试脚本为本轮工作区版本。独立目录为 `tests_dev/figures/recheck-449c181/`，避免覆盖此前图像和汇总。环境和部分包版本见该目录的 `environment.txt`，每个阶段的控制台输出保存在 `log-*.txt`。

`PASS` 表示所选调用成功执行并写出图，不代表已经完成人工视觉或生物学验收；`VISUAL-REVIEW` 表示已写图、存在已知展示问题或仍需重点检查。三对象采用的具体通路、LR、参数和警告见各 CSV 与日志，不代表所有参数组合均经过测试。

## R6-001：重跑会覆盖历史产物，且筛选为空时没有明确结论

位置：`tests_dev/render_viz_figures.R` 和 `tests_dev/render_p0_fresh_slideseq.R`。

原因与影响：原脚本固定写入 `tests_dev/figures`；重跑某一阶段会覆盖历史 PDF 和 CSV，难以区分不同代码版本。筛选条件不匹配时主脚本也没有明确的“无测试”错误。

处理：两个脚本支持 `--output-dir=...`，默认路径保持原样；主脚本对无匹配项明确报错。核验 `log-no-matching-figure.txt` 记录了预期的非零退出和错误。Slide-seq 距离图还明确记录它使用的是经缓存一致性检查的既有预过滤检查点，而非本次重新推断；全部 29,536 个细胞及 1,563 个 LR 层参与本次过滤和绘图。

## R6-002：待复核的场图曾被汇总成普通通过

位置：`render_viz_figures.R` 的 Slide-seq 场图状态判定。默认 Grid、Flow 已被标为 `VISUAL-REVIEW`，但调整参数后的两张对照图此前被写作 `PASS`，与 P3 记录中“仍需视觉复核”不一致。旧汇总还把 Flow 的积分问题写进了 Grid 的描述。

处理：默认与调整参数的 Grid 和 Flow 都标记为 `VISUAL-REVIEW`，分别描述相关显示问题；不改变绘图、场计算或默认参数。脚本的终端摘要分别列出通过、视觉复核和错误；为保持既有自动化语义，只要还有视觉复核或错误仍返回非零退出，因此非零退出不能直接解释为绘图报错。

本轮第一次 Slide-seq P1 运行时恰好发生测试脚本的并发编辑，进程在写出结果后出现解析错误；该次退出不作为证据。脚本稳定后 P1 已完整重跑，以本目录最终日志和汇总为准。

## R6-003：隔离目录中的交互依赖与浏览器验证

位置：`tests_dev/render_viz_p4.R`、`tests_dev/test-p4-browser-cdp.js`。P4 生成器支持通过 `SC_P4_OUT_ROOT` 指定隔离目录，该目录的每个数据集共享本地依赖目录；浏览器脚本通过 `SC_P4_ROOT` 指向同一批新生成页面。默认路径和原有页面目录布局保留。

`htmlwidgets` 复制 Plotly JavaScript 时未带许可证，本轮先由安装的 Plotly 包补入许可证，并让生成器今后自动复制；`log-P4-license-test.txt` 对独立临时目录进行了生成验证。三对象 15 个新页面经过 Chrome CDP 的 Plotly 初始化、点/细胞覆盖、3D 场景、堆叠层和一次 `relayout` 测试，详见 `P4/browser-summary-2026-11-10.csv`。这些测试不覆盖全部浏览器、全部用户交互或生物学解释。

## 剩余边界

- 本批主脚本 P0–P3 的结果与 PDF 由 `audit-visualization-artifacts.py` 校验数量、PDF 结构与哈希；P4 检查页面本地直接依赖和浏览器汇总。此审计不检查所有间接请求或图像科学正确性。
- P2 Visium 高密度 chord、P3 Slide-seq Grid 与 Flow 的显示效果仍需复核；所有其他 `PASS` 图同样需要按用户要求进行视觉及生物学验收。
- merged 真实独立样本比较仍缺满足共同基因轴和可解释群注释的输入；真实 Visium 自身对照不是独立样本。旧的统计阈值与方向变化语义按既定要求暂不修改。
- 旧对象、依赖恢复、其他参数组合及完整包的独立环境验收不因本轮重绘而自动通过。独立对象和可视化产物的历史记录保留。

## 2026-11-11 安装包门禁

当前源码以 `R CMD INSTALL --no-multiarch --with-keep.source` 安装到临时库后，`test-streamline-installed.R` 通过了公开 `netVisual_CommunFlow()` 的非零对角场、全零场、箭头倍率、序列化和新 R 进程重载检查。新增的 `test-installed-visualization-smoke.R` 使用已安装命名空间读取真实模拟对象，并成功生成 geometry 距离图、heatmap、表达图和 `plotly_spatialStackPlot()` HTML；结果在 `tests_dev/figures/recheck-449c181/installed-smoke/`。

安装测试曾先暴露测试脚本使用未导出的 `assay()` 泛型；测试已改用对象的公开数据结构重新执行通过。该发现不要求修改生产包 API，但说明源码直载测试不能替代安装命名空间测试。安装日志、安装后流线日志和 smoke 日志均保存在复测目录。该门禁仍使用当前机器的 renv 依赖，尚未证明在全新机器上 `renv::restore()` 后无需额外依赖安装。
