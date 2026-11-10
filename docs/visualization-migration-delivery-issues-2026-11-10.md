# 可视化迁移交付问题记录

## DELIVERY-001：测试汇总引用的图像和交互资源尚未全部提交

位置：`tests_dev/figures/` 下的 P0–P4 历史产物。原因：先前提交了部分测试汇总和截图，但若干 PDF、HTML 及其依赖仍为未跟踪文件。影响：本机可查看不代表远端检出后可查看，尤其 HTML 在缺少旁边的 JavaScript/CSS 时无法正常运行。

处理：根据三对象 P0–P3 主汇总、重新推断的距离图汇总、额外 P2 分析汇总、P3 通讯流与 P4 交互图汇总建立产物清单，包含明确待复核的图。补充 HTML 直接引用的本地资源，以及安装包提供的 Plotly 和 typedarray 许可证；`.gitattributes` 将 PDF、生成的 HTML 和 widget 依赖按二进制字节保留，避免换行转换导致清单哈希失效。清单不自动包含工作区所有其他未跟踪文件。

验证：`python tests_dev/audit-visualization-artifacts.py --check-tracked` 检查 PDF 头尾、HTML 直接依赖、本地文件 SHA-256 及 Git 索引字节一致性；2026-11-10 的清单包含 128 PDF、18 HTML、23 预览/浏览器截图、147 直接依赖与 30 许可证。此检查不证明历史图采用当前最新代码、不执行所有交互，也不替代视觉或生物学验收。已记录的 `VISUAL-REVIEW` 不因纳入 Git 而变为通过。
