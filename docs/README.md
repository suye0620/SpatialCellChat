# CellChat: Inference and analysis of cell-cell communication

## 当前开发文档入口

当前数据结构与重构实施以以下文档为准：

- [DATA_STRUCTURE.md](DATA_STRUCTURE.md)：合作者确认的最终 11-slot schema；`.mode`、`.datatype` 位于 `misc` 内部。
- [REFACTORING_PLAN.md](REFACTORING_PLAN.md)：当前实施顺序、SparseChatArray 落点和迁移边界。
- [FUNCTION_INVENTORY.md](FUNCTION_INVENTORY.md)：`R/` 下 300 个函数的逐函数档案（迁移标签、功能、流程位置、上下游、renv 依赖、死代码），用于评估重构进度。
- [ANALYSIS_SEMANTICS.md](ANALYSIS_SEMANTICS.md)：下游分析的科学问题、约简/分解、返回契约、信息损失与待裁定事项；不是迁移完成报告。
- [VISUALIZATION_SEMANTICS.md](VISUALIZATION_SEMANTICS.md)：逐函数可视化语义卡，明确表达、通讯、角色、空间统计的视觉编码与解释边界。
- [PP1_DOWNSTREAM_REFACTORING_PLAN.md](PP1_DOWNSTREAM_REFACTORING_PLAN.md)：以真实 PP1 脚本为调用方的 IO/参数矩阵、上游前置契约、分批迁移与分层验收方案。
- [DATA_STRUCTURE_11SLOT_INTRO.docx](DATA_STRUCTURE_11SLOT_INTRO.docx)：2026-10-21 生成的 11-slot 结构与运行时派生位置 Word 介绍快照；行号证据随源码演进，以源码与 DATA_STRUCTURE.md 为准。
- [DATA_STRUCTURE_11SLOT_SIMPLE.docx](DATA_STRUCTURE_11SLOT_SIMPLE.docx)：同一快照的树状简明版（结构树 + 运行时写入标注），适合快速查阅。
- [UPDATE_PLAN.md](UPDATE_PLAN.md)：早期 v3 讨论稿，仅供历史参考，不作为当前实现依据。
- [PERMUTATION_TEST_AUDIT.md](PERMUTATION_TEST_AUDIT.md)：组级汇聚置换检验的原理、空间自相关混淆问题与改进方案分层（2026-09-22）。
- [COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md](COMPUTE_AVG_NUMERIC_DIVERGENCE_REVIEW.md)：`computeAvgCommunProb` 全量数值差异归因、门控兼容性边界与新函数设计审查。
- `refactoring-plan.html`：更早的可视化讨论稿，仅供历史参考；当前 `docs/` 下未提供该 HTML，不作为可访问文档入口。

## Important update!!

## Capabilities

## Installation

SpatialChat is a powerful R package in the way of inference and analysis of cell-cell communication. In order to fully utilize its capabilities and ensure the reproducibility of your analyzing code, we recommend that you do not use SpatialChat in the base R environment. Instead, you should use `renv` to manage the R environment for working with SpatialChat.

The `renv` package helps you create reproducible environments for your R projects. Please refer to [the vignette](https://rstudio.github.io/renv/articles/renv.html) which introduces you to the basic nouns and verbs of renv, like the user and project libraries, and key functions.

Firstly, install the `renv` package using `install.packages('renv')`.

Then, you can change the working directory to the project directory you want and use `renv::init(bare = T)` to initialize the project.

``` r
#### Preparations ####
renv::init(bare = T)
# Refer to: https://rstudio.github.io/renv/articles/package-sources.html
```

Put the following code into the `.Rprofile` file under the current project directory and source it. This will help you to find the right place to install some local packages (such as some private or recompiled R packages). See details in the [renv cellar documentation](https://rstudio.github.io/renv/articles/package-sources.html?q=cellar#the-package-cellar).

``` r
Sys.getenv("RENV_PATHS_CELLAR")  
Sys.setenv(RENV_PATHS_CELLAR="renv/cellar/")
```

Secondly, install necessary dependencies using scripts below.

``` r
#### Necessary ####
renv::install("devtools") # or install.packages("devtools")
renv::install("BiocManager")

# install some necessary dependencies
renv::install("bioc::NMF")
renv::install("bioc::S4Vectors") # Linux
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
renv::install('bioc::BiocNeighbors')
devtools::install_github("JEFworks-Lab/MERINGUE")
devtools::install_github("KlugerLab/ALRA")
# If you are using  Ubuntu, you might need to install terra and RcppGSL.
# https://github.com/rspatial/terra
# https://installati.one/install-r-cran-rcppgsl-ubuntu-22-04/

# install Seurat-4.4.0
# renv::remove("Matrix")
# renv::install("Matrix") # recommend: 1.7
devtools::install_github("satijalab/seurat-object@v5.0.2")
devtools::install_github("satijalab/seurat@v4.4.0")
renv::install("Rfast2")

# Locally install SpatialChat
renv::install("SpatialChat",rebuild = T)

# Remotely install SpatialChat
......
```

Thirdly, SpatialChat also depends on some Python packages. Please create a new environment for Python and install the following packages.

``` bash
conda create -n Chat python=3.10
conda activate Chat
conda install pandas
conda install seaborn
conda env list
```

Lastly, you can also install some other R packages to do some analysis.

``` r
#### Optional ####
renv::install("readr")
renv::install("readxl")

# install VisualizeChat
devtools::install_github("zdebruine/RcppML")     # install github version
renv::install("spdep") 
devtools::install_github("rcannood/SCORPIUS")
renv::install("VisualizeChat")

# install SCP
renv::install("bioc::S4Arrays")
renv::install("bioc::rhdf5")
renv::install("bioc::XVector")
renv::install("bioc::Rsamtools")
renv::install("bioc::HDF5Array")
devtools::install_github("zhanghao-njmu/SCP")

devtools::install_github("bschilder/scKirby")  
# Refer to: https://github.com/zdebruine/singlet/tree/51646589b4475cdc4bb0b58144884ce3d651c226
# devtools::install_github("SydneyBioX/scMerge")     # install scMerge
 # devtools::install_github("immunogenomics/harmony")

# renv::install("bioc::fgsea")
# renv::install("bioc::limma")
# devtools::install_github("zdebruine/singlet")    # install singlet


renv::install("RcppPlanc", repos = "https://welch-lab.r-universe.dev")
renv::install("rliger")
```

### Installation of other dependencies

-   Install [NMF (\>= 0.23.0)](http://renozao.github.io/NMF/devel/PAGE-INSTALLATION.html) using `install.packages('NMF')`. Please check [here](https://github.com/sqjin/CellChat/issues/16) for other solutions if you encounter any issue. You might can set `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)` if it throws R version error.
-   Install [circlize (\>= 0.4.12)](https://github.com/jokergoo/circlize) using `devtools::install_github("jokergoo/circlize")` if you encounter any issue.
-   Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.
-   Install UMAP python pacakge for dimension reduction: `pip install umap-learn`. Please check [here](https://github.com/lmcinnes/umap) if you encounter any issue.

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

-   **Installation on Mac OX with R \> 3.6**: Please re-install [Xquartz](https://community.rstudio.com/t/imager-package-does-not-work-in-r-3-6-1/38119).
-   **Installation on Windows, Linux and Centos**: Please check the solution for [Windows](https://github.com/sqjin/CellChat/issues/5) and [Linux](https://github.com/sqjin/CellChat/issues/131).

## Tutorials

## Web-based “CellChat Explorer”

## Help, Suggestion and Contribution

If you have any question, comment or suggestion, please use github issue tracker to report coding related [issues](https://github.com/sqjin/CellChat/issues) of CellChat. I will answer you timely, and please remind me again if you have not received response more than three days.

### Before reporting an issue

-   First **check the GitHub [issues](https://github.com/sqjin/CellChat/issues)** to see if the same or a similar issues has been reported and resolved. This relieves the developers from addressing the same issues and helps them focus on adding new features!
-   The best way to figure out the issues is **running the sources codes** of the specific functions by yourself. This will also relieve the developers and helps them focus on the common issues! I am sorry, but I have to say I have no idea on many errors except that I can reproduce the issues.
-   Minimal and **reproducible example** are required when filing a GitHub issue. In certain cases, please share your CellChat object and related codes to reproduce the issues.
-   Users are encouraged to discuss issues and bugs using the github [issues](https://github.com/sqjin/CellChat/issues) instead of email exchanges.

### Contribution

CellChat is an open source software package and any contribution is highly appreciated!

We use GitHub's [Pull Request](https://github.com/sqjin/CellChat/pulls) mechanism for reviewing and accepting submissions of any contribution. Issue a pull request on the GitHub website to request that we merge your branch's changes into CellChat's master branch. Be sure to include a description of your changes in the pull request, as well as any other information that will help the CellChat developers involved in reviewing your code.

## System Requirements

-   Hardware requirements: CellChat package requires only a standard computer with enough RAM to support the in-memory operations.
-   Software requirements: This package is supported for macOS, Windows and Linux. The package has been tested on macOS: Mojave (10.14.5) and Windows 10. Dependencies of CellChat package are indicated in the Description file, and can be automatically installed when installing CellChat pacakge. CellChat can be installed on a normal computer within few mins.

## How to cite?

Suoqin Jin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Raul Ramos, Chen-Hsiang Kuan, Peggy Myung, Maksim V. Plikus, Qing Nie. Inference and analysis of cell-cell communication using CellChat. Nature Communications, 12:1088 (2021). <https://www.nature.com/articles/s41467-021-21246-9>
