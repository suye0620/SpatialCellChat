# Agent Note: SpatialCellChat 最终 11-slot 数据结构

Status: implemented
Date: 2026-08-15
Decision type: architecture

## Problem

SpatialCellChat 同时承载表达矩阵、空间图像、细胞元数据、LR 通讯结果和 pathway 通讯结果。旧版 slot 命名和多样本分支会让下游函数直接依赖多个历史路径，增加迁移成本和结构不一致风险。

## Decision

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

最终目标对象使用 11 个顶层 S4 slot：

```text
assay, images, meta, idents, features, LR, dr, net, netP, DB, misc
```

约束如下：

- `assay` 保存 `raw`、`norm`、`scale`、`smooth`、`signaling` 表达层；`assay$project` 更名为 `assay$smooth`。
- `images` 独立保存空间坐标、坐标系统、空间因子、raster、FOV 和距离/接触缓存。
- `meta` 与 `features` 取代旧的 `colData` 与 `rowData`。
- `idents` 在单个切片（rep）是覆盖全部细胞的直接 factor，不包装为 `joint` 或 `samples` 子字段。
- `LR` 保存 LR 条目元数据；`net` 保存 LR-level 结果；`netP` 保存 pathway-level 结果。
- `net` 与 `netP` 均使用 `$cell` / `$group` 子结构，分别容纳 `prob`、`pval`、`count`、`weight`、`centr`、`field` 等结果。
- `DB` 作为独立的 LR 数据库输入资源保存。
- `.mode` 与 `.datatype` 放入 `misc`，不是顶层 slot；多样本归档放入 `misc$.datasets`。
- `misc$.param` 保存重要参数，`misc$.log` 保存对象修改日志。

对象写入后必须通过 `validObject()` 或 `validateSpatialCellChat()` 检查；表达、元数据、空间坐标和通讯矩阵的名称必须保持对齐。

## Alternatives considered

- 继续保留旧的 14-slot 结构：会延长旧路径和兼容分支的生命周期。
- 使用 v3 讨论版的 9-slot 结构：与合作者最终确认的命名和职责不一致。
- 把 `.mode`、`.datatype` 提升为顶层 slot：增加顶层结构而没有提供必要的独立对象语义。

## Consequences

这是后续 `spatial.R`、`analysis.R`、`modeling.R` 和 `visualization.R` 迁移的唯一目标 schema。旧 slot 访问应在迁移完成后删除，而不是继续添加新的兼容分支。

## Evidence

- `docs/DATA_STRUCTURE.md`
- `docs/REFACTORING_PLAN.md`
- `R/SpatialCellChat_class.R`
