# 可视化迁移基线问题记录：P0

日期：2026-11-10

## P0-001：旧测试对象缺少距离缓存来源，通信距离图拒绝绘制

位置：`spatialC2CDistPlot(mode = "communication")` 调用 `.sc_c2c_inference_cache()`，检查推断时的坐标、空间参数和当前距离缓存是否来自同一次 `computeCommunProb()`。

现象：P0 回归中，模拟、Visium、完整 Slide-seq 的旧缓存对象均报距离缓存或概率来源过期；geometry 模式可自行重算，已绘图。已直接检查的模拟与 Visium RDS 的 `.distance` 只有 `d.spatial`、`adj.contact`、`.parameters`，缺少 `.coordinates`，通信推断参数也没有 `distance.coordinates` 和 `distance.parameters` 来源记录。

原因和本质：这些缓存生成于距离缓存来源校验引入之前，无法证明当前细胞坐标、距离和已保存概率是一致的。没有证据证明它们必然已错，但也不能证明它们正确；目前的显式报错是防止误用旧概率所必需的安全行为，不应直接补写元数据绕过检查。

修复方向：使用原始输入和原参数重新运行 `computeCommunProb()` 及必要的下游分析，生成具备真实来源的新对象；保留原 RDS，分别验证通信与几何图。完整 Slide-seq 推断成本较高，且现有 fixture 优先读取旧 benchmark 的后处理对象，不能当成已经重新推断。影响三对象旧缓存的通信距离图及其他使用该来源校验的消费者；几何模式不受同样的概率来源限制。

重测：模拟和 Visium 均从原始输入重新执行通信推断后成功绘制通信距离图；Visium 还执行了原 fixture 的 group averaging 与通信过滤，但距离图不依赖其后耗时的 group centrality，因此定向测试不再计算该部分。旧对象仍保持原状及原失败记录。Slide-seq 重新推断与绘图结果单独记录在 `summary-slideseq-P0-spatialC2CDistPlot-fresh.csv`。

## P0-002：Slide-seq fixture 的过滤函数参数过期

全量 Slide-seq 首次新推断成功：29,536 个细胞、1,563 个 LR 层、407,065,015 个非零值，耗时约 335 秒；随后定向脚本调用 `filterProbability(verbose = TRUE)` 报 unused argument，没有进入绘图。原 fixture 的无 benchmark 分支也传入 `verbose = FALSE`，但当前 `filterProbability()` 没有该参数。去除无效参数并重跑；新增仅在不存在时写入的本地新推断检查点，后续使用前必须检查来源，不覆盖原始 RDS。

第二次推断完成并保存新检查点，但过滤时 future 导出的全局对象约 10.31 GiB，超过默认 500 MiB 限额。定向脚本使用与全量 fixture 一致的 20 GiB 上限，从检查点继续；这只是解除大小限制，不解决闭包的内存和传输开销。

最终重测：从新检查点成功完成 1,563 层过滤并绘制完整 Slide-seq 的通信距离 PDF；检查点为未过滤对象，仅保存在本地测试数据目录，不覆盖原对象。PDF 预览中的两个 annotation 面板均有实际密度曲线。该结果不消除旧 RDS 无来源记录的限制。
