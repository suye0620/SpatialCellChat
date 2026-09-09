# SpatialCellChat 数据结构设计文档

> 版本: 2026-08-03
> 基于合作者确认的最终数据结构整理
>
> **当前权威说明**：本文档的“最终数据结构”节以合作者最新口头确认为准：`.mode` 和 `.datatype` 归入 `misc`，不是顶层 S4 slot。因此 SpatialCellChat 的最终目标是 **11 个顶层 slot**。文档中其他历史讨论、13-slot 表述和旧命名仅作设计演进记录，不作为当前实现依据。

---

## 基石：SparseChatArray

SpatialCellChat 的核心数据是 3D 通讯概率矩阵：**细胞 × 细胞 × LR对**（nC × nC × N）。传统 3D 稠密数组在单细胞规模（几万细胞 × 几千 LR 对）下会膨胀到不可接受的内存占用，必须使用稀疏表示。

### `SparseChatArray` 的设计

`SparseChatArray` 是一个轻量 S3 类，本质是**命名列表 + 类属性**，每个元素是一个 `dgCMatrix`（Matrix 包标准稀疏列压缩矩阵），所有矩阵维度相同。

```r
SparseChatArray <- function(x) {
  # x: 等维 dgCMatrix 列表
  structure(x, class = "SparseChatArray")
}
```

三维索引通过 S3 泛型 `[.SparseChatArray` 实现，语义与原生 R 数组一致：

- `arr[i, j, k]` — 按三维下标取值，`k` 选出单个 LR 对层，`i/j` 对该层矩阵做行列索引
- `arr[i, j, ]` 或 `arr[i, , ]` — 跨层索引，返回 `SparseChatArray`
- `dim(arr)` → `c(nrow, ncol, nlayers)`
- `dimnames(arr)` → `list(rownames, colnames, layernames)`

### 为什么选择 `list_of_dgCMatrix` 而非其他方案？

| 方案 | 缺点 |
|------|------|
| `spatstat.sparse::sparse3Darray` | 旧 CellChat v2 使用。构造慢（需要 C 级索引重排），依赖 `spatstat.sparse` 庞大生态链，单层访问需全局解压 |
| `DelayedArray` + `HDF5Array` | 引入 Bioconductor 重型依赖，小/中数据集反向降低性能 |
| 3D `dgTMatrix` / `COO` 平铺 | 按 LR 对聚合统计复杂，行/列切片效率低于列压缩格式 |
| 原生 3D 稠密数组 | nC × nC × N 在稀疏场景下完全不可接受 |

### `list_of_dgCMatrix` 的优势

1. **单层 O(1) 访问** — `x[[k]]` 直接取出第三维任意切片，无需像 `sparse3Darray` 那样从全局三元组索引重建
2. **按层利用不同稀疏度** — 每个 LR 对的连接模式不同（有的广泛分泌，有的只在少数细胞对间接触），各层独立压缩，不做无效填充
3. **复用 Matrix 包优化** — `dgCMatrix` 的列压缩存储、`[` 索引、矩阵乘法和列/行切片全部由 C 代码实现，经过数十年工业级验证
4. **GPU 友好** — 单层 `dgCMatrix` 可直接传输到 GPU（如 cuML / torch），无需额外格式转换
5. **无外部依赖** — 仅依赖 `Matrix`（CRAN 核心包），不引入 `spatstat.sparse` 或 `DelayedArray`

### 与 `future` 并行策略的配合：双重存储模式

实践中的核心模式是**同时在 `net` 中保存两种视图**：

```r
# modeling.R:287-288
Tmp <- list(prob.cell = Prob.cell_, ...)  # list of dgCMatrix（并行视图）
net <- list(prob.cell = Prob.cell,        # sparse3Darray（索引视图）
            tmp = Tmp)
```

| 存储位置 | 格式 | 用途 |
|----------|------|------|
| `net$prob.cell` | 3D array (sparse3Darray / SparseChatArray) | 三维索引、可视化、subset |
| `net$tmp$prob.cell` | **named list of dgCMatrix** | **`my_future_lapply` 并行迭代** |

后续函数的典型用法：

```r
# 从 net 中取出 list 视图
prob.cell_ <- object@net$tmp$prob.cell

# 按 LR 对并行遍历——每个 worker 拿到一个独立的 dgCMatrix
my_future_lapply(
  X = seq_along(prob.cell_),
  FUN = function(i) {
    prob.cell.i <- prob.cell_[[i]]  # 单个 dgCMatrix，直接操作
    # ... 独立计算 ...
  }
)
```

这一模式的关键好处：

- **无需 3D→2D 解压** — `prob.cell_` 本身就是 list，`[[i]]` 直达第 i 层矩阵，没有索引重排开销
- **future worker 天然隔离** — 每个 worker 拿到的 dgCMatrix 是独立 S4 对象，不存在共享内存竞争
- **按需取子集** — 只迭代感兴趣的 LR 对子集时，直接 `prob.cell_[idx]` 过滤 list，不用从 3D array 切片
- **旧版兼容** — `net$prob.cell` 保留 `sparse3Darray` 格式（旧 CellChat 兼容），新代码逐步迁移到 `SparseChatArray`

总结：list_of_dgCMatrix 的底层格式同时满足了**三维语义索引**（通过 S3 方法）和**按层独立并行**（通过 `net$tmp$prob.cell`）两个需求，是 `future_lapply` 并行策略与 3D 稀疏存储的自然结合。

### 与旧 CellChat v2 的 `sparse3Darray` 对比

旧 CellChat v2 使用 **`spatstat.sparse` 的 `sparse3Darray`** 存储 `prob.cell`，构造时需将列表通过 `my_as_sparse3Darray` 合并为一个全局三元组（i, j, k, x）。这一转换：

- 在 LR 对较多时导致显著的构造开销
- 单层访问 `prob.cell[,,k]` 需要从全局索引中筛选出 k 层对应的行，不如 `list[[k]]` 直接
- 跨层操作如 `marginSumsSparse` 依赖 `spatstat.sparse` 内部优化，如果只需操作特定几层仍需解压全局数组

`SparseChatArray` 放弃全局三元组压缩，返回"按层分散"的策略。代价是跨第三维的聚合操作（如对所有 LR 对的概率求和）需要逐层迭代，但这类操作在分析流程中远少于单层访问。



> 问题：counts weights 维度求和是否还原，是否会产生开销

### 稀疏度预期

通讯概率矩阵天然高度稀疏：一个 LR 对通常只在一小部分细胞对之间有非零连接概率。对于 nC = 10000 细胞、N = 2000 LR 对的典型数据集：

- 稠密数组：10000 × 10000 × 2000 × 8 bytes ≈ **1.5 TB**
- `SparseChatArray` 单层稀疏度通常 0.01%~1%，内存占用降到 **MB~GB 级**

`print.SparseChatArray` 会在打印时计算整体稀疏度，便于监控。`SparseChatArray` 也可以无损转换为长格式 DataFrame，便于下游分析和可视化。

### 未来扩展

`SparseChatArray` 的接口是 S3 泛型，可在不改变类定义的情况下添加新方法：
- `as.data.frame.SparseChatArray` — 转为长格式（source, target, layer, value）
- `as.matrix.SparseChatArray` — 按第三维折叠/展开
- `Ops.SparseChatArray` — 逐层算术运算
- `marginSums.SparseChatArray` — 沿任意维度求和

---

## 设计背景

SpatialCellChat（CellChat v3）的核心创新是从 2D 基因表达空间（cell × gene）进入 3D 通讯空间（cell × cell × LR/pathway）。这一维度提升对数据结构提出了新的要求——不再仅仅是"表达矩阵 + 元数据"的扁平结构，而是需要同时管理两个空间的数据及其元数据。

---

## 设计原则（历史讨论版，术语已由最终结构修正）

> 本节保留早期 v3 讨论中的 AnnData 对照。当前实现目标不再使用 `colData`、`rowData`、`reducedDims`、`signalingData`、`netData` 作为顶层 slot；请以本文“最终数据结构”节中的 `meta`、`features`、`dr`、`LR`、`net`、`netP` 为准。

1. **借鉴 AnnData 的哲学，但完全独立实现** — 采用分层表达数据、元数据和降维结构，但 S4 class 独立，不继承任何框架
2. **骨架与内容分离** — 最终 S4 class 定义 11 个顶层 slot，各 slot 内部结构由分析函数按需填充，不预先在 class 定义中硬编码
3. **单样本与多样本使用同一结构** — 避免 single/merged 模式下 slot 结构不一致导致下游函数写大量 `if (mode == "merged")` 分支
4. **只存核心产出，不存中间状态** — `tmp` 字段被清理，计算中间结果用完即弃
5. **高维预留** — 从 2D 空间位置到 3D 空间位置，数据结构天然兼容，无需结构性改动

---

## 与 AnnData 的对应关系（历史映射参考）

> 下表用于解释设计来源，不代表当前顶层 slot 命名。当前实施以“最终数据结构”节为准：`meta`、`features`、`dr`、`LR`、`net`、`netP` 分别承担旧讨论版中部分 `colData`、`rowData`、`reducedDims`、`signalingData`、`netData` 的职责。

| AnnData (Python) | SpatialCellChat (R) | 说明 |
|---|---|---|
| `X` / `layers` | `assay` | 表达矩阵列表 |
| `obs` | `colData` | 细胞级元数据 |
| `var` | `rowData` | 基因级元数据 |
| — | `signalingData` | **通讯条目级元数据** — AnnData 没有的概念，对应 3D 空间的第三维度 |
| `obsm` | `reducedDims` | 降维/嵌入坐标 |
| `uns` | `misc` | 杂项数据 |
| — | `net` | **通讯核心数据** — 3D 稀疏数组 `SparseChatArray` |
| — | `netData` | **通讯衍生数据** — pval/centr/count/weight |
| — | `DB` | **L-R 数据库** — 独立 slot，和 assay 平级 |

---

## 关键讨论与决策（含历史方案演进）

> 本节记录从 v3 讨论版到最终版的设计过程。若本节与“最终数据结构”节存在差异，以最终数据结构节及合作者最新确认的 11-slot schema 为准。

### 1. 是否继承 Seurat / SingleCellExperiment？

**结论：不继承，独立 S4 class。**

- 继承 Seurat 丧失独立性，论文需引用，且 Seurat API 变动会牵制开发
- 继承 SCE 更轻量，但仍引入 Bioconductor 依赖链
- 独立设计最灵活，用户需要时自行写 `as()` 转换

### 2. 合并模式如何设计？

**结论：`@net` 始终并列存放，结果物归原位，不改变 slot 结构。`@idents` 保持为直接的 factor，不包装为 list。**

Phase 0 不变量：`idents` 始终是覆盖全部细胞的 factor，保持与旧版对象的直接访问方式兼容。样本级信息和样本级通讯结果归档到 `misc$.datasets`；`misc$.mode` 只取 `single` 或 `merged`，顶层 `net` / `netP` 表示当前 joint 结果。

- 错误的直觉：合并是把两个样本的细胞合到一起"加总"
- 正确的理解：合并是创建一个**对比容器**，各样本的通讯结果并排放置，用于差异比较
- "表达数据"是加和（rbind），"通讯结果"是并列（list）
- 单样本：`net$cell$prob` → SparseChatArray nC×nC×N
- 多样本：`net$cell$prob` 为空，各样本结果存 `misc$.datasets[[name]]`，其内部镜像保存对应的 `net`、`netP`、`dr` 和 `misc` 数据
- 跑完 joint 分析后，`net$group$prob` 存放统一的 K×K×N 结果
- 对比分析（如"处理组 vs 对照组"）是函数层能力，不在 slot 层预设

### 3. prob / pval / count / weight 的层级关系

```
net                  → LR-level 通讯结果（prob, pval, count, weight, centr, field）
netP                 → pathway-level 通讯结果（prob, pval, count, weight, centr, field）
  ├── cell$pval      → SparseChatArray  nC×nC×N（3D 置信度掩码）
  ├── cell$count     → dgCMatrix        nC×nC（汇总计数）
  ├── cell$weight    → dgCMatrix        nC×nC（汇总强度）
  ├── cell$centr     → list（细胞级出入度）
  ├── cell$field     → list（向量场）
  ├── group$pval     → SparseChatArray  K×K×N
  ├── group$count    → matrix           K×K
  ├── group$weight   → matrix           K×K
  └── group$centr    → list（群级出入度）
```

`LR` 是 LR-level 通讯条目的元数据，`netP` 保存 pathway-level 结果；`pval` 是置信度标注，`count`/`weight` 是降维汇总，`centr` 是网络中心性。三维 `prob` / `pval` 统一使用 `SparseChatArray`。

### 4. `signalingData` 方案（历史讨论，已废止）

早期方案曾把 LR 和 pathway 视为同一通讯空间的不同**粒度**，并计划通过 `signalingData$type` 区分：

- 任何函数通过 `signalingData` 的行号索引，不区分 `net` 与 `netP`
- `computeCommunProbPathway` 在 `signalingData` 尾部追加 pathway 行，同时追加 `prob` 的第三维
- 将来加其他粒度（如单个基因、复合物）无需改结构

该方案已由最终 schema 替代。当前使用独立的 `LR`、`net` 和 `netP` 结构，见后文第 5 节。

### 5. `net` 与 `netP` 的最终关系

早期讨论曾提出通过 `signalingData$type` 合并 LR 和 pathway 结果；该方案已被最终设计替代。当前保留两个独立 slot：

- `net`：LR-level 通讯结果；
- `netP`：pathway-level 通讯结果。

两者均采用 `$cell` / `$group` 结构，三维 `prob` / `pval` 使用 `SparseChatArray`。`signalingData` 不再是最终 schema 的顶层 slot。

### 6. 为什么 `DB` 是独立 slot？

CellChatDB 是核心输入资源，和 `assay` 平级。它不属于任何表达层或通讯层，且在整个分析流程中多次被引用（subsetDB、searchPair、computeCommunProb），独立 slot 最自然。

### 7. `misc` slot 设计

沿袭 Seurat 的 `misc`（miscellaneous）命名。根据合作者最终确认，`.mode` 和 `.datatype` 保留在 `misc` 内部，而不是提升为顶层 slot。内部设计：

- `.param` — 重要、复用、提示性参数（如 `datatype`、`spatial.factors`）
- `.log` — 函数操作日志（记录每次对 object 的修改）
- `.var.features` — 高变基因
- `.distance` — 已废止的旧位置；当前空间距离矩阵和接触邻接矩阵统一缓存于 `images$.distance`
- `.datasets` — 多 rep 存档（内部用，不对外暴露）

加 `.` 前缀标记"内部字段"，和用户可设字段区分。

### 8. 合并后的 `misc$.datasets` 内部结构

> 以下镜像结构按最终 schema 使用 `net`、`netP`、`dr` 和 `misc`；早期草稿中的 `communication`、`netData`、`reducedDims` 名称均已废止。

```r
misc$.datasets$sample1
  └─ net$cell$prob     → SparseChatArray (原样归档)
  └─ net$group$prob    → SparseChatArray (原样归档)
  └─ net$cell$...       → LR-level 各衍生数据
  └─ net$group$...      → LR-level 各衍生数据
  └─ netP$cell$...      → pathway-level 各衍生数据
  └─ netP$group$...     → pathway-level 各衍生数据
  └─ dr                 → 原样本的降维/空间坐标等
  └─ misc$.param       → 原样本的参数

misc$.datasets$sample2 → 同上
```

`misc$.datasets` 的各元素结构和单样本对象的 `net`/`netP`/`dr`/`misc` 完全一致，形成镜像关系。`misc$.mode` 与 `misc$.datatype` 仍按最终 schema 保存在 `misc` 内部。

---

## 最终数据结构

> **最终确认版（2026-08-03）**：`.mode` 和 `.datatype` 属于 `misc` 内部字段，不是顶层 slot。以下 11-slot 结构是当前后续开发唯一采用的目标 schema。

SpatialCellChat S4 class 共 **11 个顶层 slot**：

```text
assay, images, meta, idents, features, LR, dr, net, netP, DB, misc
```

```
SpatialCellChat (S4)
│
├── assay ───────────── 表达层 (list)
│   ├── raw       = dgCMatrix | NULL   原始计数
│   ├── norm      = dgCMatrix          标准化后（默认主矩阵）
│   ├── scale     = dgCMatrix | NULL   缩放后
│   ├── smooth    = dgCMatrix | NULL   平滑后（原 project 更名）
│   └── signaling = dgCMatrix | NULL   信号基因子集
│
├── images ──────────── 空间信息 (list)
│   ├── coordinates      通讯分析坐标矩阵 (nC × 2/3)
│   ├── coordinate.system 坐标单位、校准状态和轴方向
│   ├── spatial.factors  坐标→微米转换参数 (ratio, tol)
│   ├── rasters           命名 raster 记录列表（可选）
│   │   └── <name>
│   │       ├── image             raster array
│   │       ├── scale.factors     spot/fiducial/hires/lowres 等
│   │       ├── spot.coordinates  图像像素坐标 data.frame
│   │       ├── spot.radius       spot 半径（可选）
│   │       ├── coordinate.system raster 坐标单位与方向
│   │       ├── transform         analysis → raster_pixel 3×3 仿射变换
│   │       └── source             图像来源（可选）
│   ├── fov                视野几何、分子和边界（可选）
│   ├── .distance          空间距离/接触矩阵（缓存）
│   └── .grid              网格化空间缓存（可选）
│
├── meta ─────────────── 细胞级元数据 (data.frame)
│   ├── ...              用户自定义（不含 ident，见 idents slot）
│   └── samples          样本来源（多样本时）
│
├── idents ───────────── cell group 标签 (factor)
│                       始终覆盖全部细胞；保持 object@idents 直接访问兼容
│
├── features ─────────── 基因级元数据 (data.frame)
│   └── ...              基因注释
│
├── LR ───────────────── LR 通讯条目信息 (list)
│                       保持与旧 CellChat 一致的格式，
│                       包含 LRsig、pairwiseRank 等
│
├── dr ───────────────── 降维/嵌入 (list)
│   ├── pca              matrix (nC × ...)  PCA
│   └── umap             matrix (nC × ...)  UMAP
│   └── ...              其他降维方法
│
├── net ──────────────── LR-level 通讯结果 (list)
│   ├── cell$prob        SparseChatArray    nC × nC × N
│   ├── cell$pval        SparseChatArray    nC × nC × N
│   ├── cell$count       Matrix             nC × nC
│   ├── cell$weight      Matrix             nC × nC
│   ├── cell$centr       list               细胞级中心性
│   ├── cell$field       list               向量场
│   ├── group$prob       SparseChatArray    K × K × N
│   ├── group$pval       SparseChatArray    K × K × N
│   ├── group$count      matrix             K × K
│   ├── group$weight     matrix             K × K
│   └── group$centr      list               群级中心性
│
├── netP ─────────────── pathway-level 通讯结果 (list)
│                       结构同 net，所有 cell/group 字段一致
│                       维度第三维从 N (LR 对) 变为 P (通路)
│
├── DB ───────────────── L-R 数据库 (list)
│   ├── interaction      data.frame
│   ├── complex          data.frame
│   ├── cofactor         data.frame
│   ├── geneInfo         data.frame
│   │
│   DB 按物种分：DB.human, DB.mouse, DB.zebrafish 等
│
└── misc ─────────────── 杂项 (list)
    ├── .mode            "single" | "merged"
    ├── .datatype        "RNA" | "spatial"
    ├── .param           其他重要/复用参数
    ├── .log             操作日志（仿 Seurat：记录函数名、时间、参数）
    ├── .var.features    高变基因
    └── .datasets        多 rep 存档（内部用）
```

### v3 讨论版 → 最终版的变更说明

> 本表记录历史方案如何收敛到当前最终设计；“最终版”一列是当前实施口径。

| 变化 | v3 讨论版 | 最终版 | 理由 |
|------|-----------|--------|------|
| slot 数量 | 9 个 | **11 个** | 恢复部分语义独立的 slot；`.mode` 和 `.datatype` 保留在 `misc` 内 |
| `assay$project` | project | **smooth** | 更准确描述数据性质 |
| `colData` → `meta` | colData | **meta** | 与旧 CellChat 一致，降低用户迁移成本 |
| `rowData` → `features` | rowData | **features** | 同上 |
| `signalingData` 撤销 | 独立 slot | 归入 **LR** | 保持与旧版兼容，避免破坏现有分析流程 |
| `reducedDims` → `dr` | reducedDims | **dr** | 与旧 CellChat 一致 |
| `images` 恢复 | 拆入 reducedDims + misc | **独立 slot** | `images` 语义强，独立更自然 |
| `idents` 恢复 | colData$ident | **独立 factor slot** | 核心标识符，保持旧版直接访问兼容 |
| `netData` 撤销 | 独立 slot | 归入 `net` / `netP` | prob/pval/centr 共存同一 list，不强行拆分 |
| `netP` 恢复 | 并入 signalingData$type | **独立 slot** | pathway 与 LR 分析流程不同，分离更实用 |
| `misc$.mode` | `.param` 内部 | **`misc$.mode`** | 合作者确认保留在 `misc`，用于区分 single/merged |
| `misc$.datatype` | `.param` 内部 | **`misc$.datatype`** | 合作者确认保留在 `misc`，用于区分 RNA/spatial |
| `misc$.distance` 转移 | `misc` 中 | 移入 **`images$.distance`** | 距离矩阵与空间数据关系更紧密 |

### slot 分类

| 类别 | slot | 数量 |
|------|------|------|
| **表达数据** | `assay` | 1 |
| **空间数据** | `images` | 1 |
| **元数据** | `meta`, `idents`, `features`, `LR` | 4 |
| **降维** | `dr` | 1 |
| **通讯结果** | `net`, `netP` | 2 |
| **知识库** | `DB` | 1 |
| **配置/日志** | `misc` | 1 |

---

## 设计优势

### 1. 维度提升语义清晰
3D 通讯空间中的所有数据（prob、pval）用同样的前两维索引（cell × cell），第三维通过 `LR` slot 提供语义标注。LR 和 pathway 两种粒度分别在 `net` / `netP` 中管理。

### 2. single/merged 结构一致
`net` 和 `netP` 永远是 `$cell` / `$group` 两层，不因单样本或多 rep 而改变。各 rep 的原始结果存 `misc$.datasets`，和主结构镜像，避免了下游函数写大量 `if (.mode == "merged")` 分支。

### 3. 聚合函数的降维结果自然归位
`aggregateNet` 对 `net$cell$prob` 第三维求和，结果 `count` 和 `weight` 存入 `net$cell` 和 `net$group`。降维汇总后 2D 矩阵从 3D 数组中自然导出，存放位置明确。

### 4. 对 3D 空间位置天然兼容
`images$coordinates` 只需增加一列（x, y, z），无需任何结构变更。对于 Visium，原始染色切片图像和图像像素坐标保存在命名的 `images$rasters` 记录中；`images$coordinates` 仍专门用于空间分析和距离计算。

### 5. 骨架稳定，内容灵活
S4 class 只定义 11 个顶层 slot，各 slot 内部的"形状"由分析函数运行时决定。`setClass` 不用为未来的 `centr`、`field`、各种衍生指标预留名字，它们自然出现在 `net$cell` 和 `net$group` 下；`.mode` 和 `.datatype` 作为 `misc` 内部字段保存。

---

## 与旧 CellChat v2 结构对比

| 旧结构 | 新结构 | 改进 |
|--------|--------|------|
| `@data` / `@data.raw` / `@data.signaling` / `@data.scale` / `@data.project` 5 个独立 slot | `@assay` 一个 list slot | 减少冗余 slot，表达数据生命周期清晰 |
| `@options` 散装参数 | `@misc$.param` 统一，`.mode` / `.datatype` 保留在 `misc` | 参数和模式信息集中管理 |
| `@net$tmp` 中间状态堆积 | 清理，中间结果参数间传递 | 减少对象体积 |
| single/merged 模式下 slot 结构不一致 | `@net` / `@netP` 始终 `$cell` / `$group` 两层 | 函数无需写 mode 分支 |
| 无 `@features` | `@features` 新增 | 基因级元数据独立存放 |
| 无 `SparseChatArray` | 3D 通讯数据用 `SparseChatArray`（list_of_dgCMatrix） | 替代 `spatstat.sparse::sparse3Darray` |
| 无操作日志 | `@misc$.log` | 可追溯全流程 |
| `@net$prob.cell` 用 sparse3Darray | 逐步迁移到 `SparseChatArray` | 无外部依赖，GPU 友好 |
