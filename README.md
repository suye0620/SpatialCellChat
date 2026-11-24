# SpatialCellChatDev

SpatialCellChatDev is a development fork of SpatialCellChat (CellChat v3) for inferring, analyzing, and visualizing cell–cell communication from spatial transcriptomics. It introduces an 11-slot object and accelerates inference by operating on sparse spatial edges and ligand–receptor (LR) layers. Depending on the assay, a node represents a cell or a spot. [Upstream SpatialCellChat](https://github.com/jinworks/SpatialCellChat)

The R package is named `SpatialCellChatDev` and can be installed alongside `SpatialCellChat`. Both packages define an S4 class named `SpatialCellChat` with different layouts. Run complete upstream and Dev workflows in separate R sessions when comparing results.

[Installation](#installation) · [Visium quick start](#quick-start-with-the-bundled-visium-data) · [Migrate an existing object](#migrate-a-saved-upstream-object)

## What changed from upstream

The upstream layout below refers to the 15-slot baseline used for migration testing. Dev groups expression layers in `assay`, spatial information in `images`, and communication results by both signaling level and node resolution.

| Data | Upstream | SpatialCellChatDev |
| --- | --- | --- |
| Normalized expression | `@data` | `@assay$norm` |
| Raw, scaled, and projected expression | `@data.raw`, `@data.scale`, `@data.project` | `@assay$raw`, `@assay$scale`, `@assay$smooth` |
| Signaling-gene expression | `@data.signaling` | `@assay$signaling` |
| Cell-level LR probabilities | `@net$prob.cell` | `@net$cell$prob` |
| Group-level LR probabilities and p-values | `@net$prob`, `@net$pval` | `@net$group$prob`, `@net$group$pval` |
| Pathway probabilities | `@netP$prob.cell`, `@netP$prob` | `@netP$cell$prob`, `@netP$group$prob` |

`net` stores LR-level results and `netP` stores pathway-level results. A `SparseChatArray` stores each three-dimensional communication array as named, compressed sparse-column `dgCMatrix` layers. Its axes are source × target × LR pair or pathway. Use `dim(chat@net$cell$prob)` to inspect the full shape and `chat@net$cell$prob[[1L]]` to retrieve the first sparse LR layer. Spatial coordinates, calibration, and the distance cache are stored together in `images`. The distance cache records the inputs used to build it so inference and plotting can detect stale distances. [11-slot object design](docs/DATA_STRUCTURE.md#最终数据结构)

## Performance and agreement with upstream

- `computeCommunProb()` evaluates each LR layer with an Rcpp kernel over the spatial edge mask, avoiding an outer product of ligand and receptor expression across all cells. Per-layer work scales with the number of stored spatial edges, so spot density and distance thresholds affect runtime and memory.
- `computeAvgCommunProb()` uses compiled kernels for group aggregation and permutations. It applies `min.percent` numerically after rounding to one significant digit. The pinned upstream implementation used formatted strings for that comparison, which can change which group-level interactions pass the threshold.
- Grid operations use sparse point-to-grid membership. Communication-field aggregation reduces the component vectors directly, avoiding repeated `sf` row subsetting. [Optimization notes](docs/COMPUTATION_OPTIMIZATION.md) · [Grid aggregation benchmark](tests_dev/figures/P3/grid-aggregation-benchmark-2026-11-10.csv)

The October 2026 **MBM05 serial benchmark** compared the pinned upstream implementation (`be88f30…`) with the refactored source on the same input checkpoint. MBM05 had 29,536 cells, 8 groups, and 1,563 LR layers. Both workflows ran serially with 100 permutations for group aggregation. [Benchmark setup and stage timings](docs/SERIAL_COMPARISON_GROUP_AVG.md#串行阶段耗时serial-vs-serial)

| Stage | Upstream | Refactored | Speedup |
| --- | ---: | ---: | ---: |
| Cell-level inference (`computeCommunProb()`) | 18,489 s | 318 s | 58× |
| Cell-level probability filtering (`filterProbability()`) | 4,365 s | 86 s | 50× |
| Group aggregation (`computeAvgCommunProb()`) | 2,519 s | 73 s | 34× |
| Full benchmark, through centrality | 7.3 h | 23 min | 18.9× |

These timings describe the measured workflow on MBM05. Centrality itself took longer in the refactored run. Cell-level probabilities matched exactly on MBM05, while some group-level significance calls changed. In separate 3cts and 4cts simulation comparisons, cell-level probabilities and the centrality outputs checked agreed with upstream. Five 4cts group p-values differed by `0.01`, with all five above `0.05` in both runs. Replaying Dev with the upstream RNG state immediately before drawing permutations reproduced those p-values. [MBM05 numerical comparison](docs/SERIAL_COMPARISON_GROUP_AVG.md) · [Simulation comparison](tests_dev/results/legacy-inference-2026-11-16/comparison-report.md)

## Installation

Install from this repository using the `renv` project library. The checkout supplies `.Rprofile`, `renv/activate.R` and `renv.lock`. The startup profile activates `renv`, the project library holds installed packages, and the lockfile records the dependency versions used for development. Each combination of platform, R minor version and `renv` profile gets its own library. [renv project workflow](https://rstudio.github.io/renv/articles/renv.html)

### 1. Install R and a compiler, then clone the repository

- **Windows.** Use R 4.4.3 and Rtools44 for the R 4.4 installation below. If you select a different R series, install the matching Rtools. RStudio lets you choose the R installation under **Tools → Global Options → General**. [R for Windows](https://cran.r-project.org/bin/windows/base/old/) · [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS or Linux.** A source build requires a compiler toolchain for the selected R release. On macOS, install the Xcode Command Line Tools with `xcode-select --install` and use the matching Fortran toolchain if a dependency requires it. On Linux, install your distribution's R development tools and system libraries required by source dependencies. For example, `sf` uses GDAL, GEOS and PROJ. [R for macOS tools](https://mac.r-project.org/tools/) · [sf system requirements](https://r-spatial.github.io/sf/)

Clone the repository from a terminal.

```sh
git clone https://github.com/suye0620/SpatialCellChat.git C:/src/SpatialCellChat
cd C:/src/SpatialCellChat
```

### 2. Choose and activate an environment

Choose your environment before starting R from the checkout.

| Your starting point | Environment and dependency setup |
| --- | --- |
| New Windows setup with R 4.4.3 | Activate the `compatible` profile below, then install its dependencies in step 3. |
| Existing checkout with a working project library | Open the library you already use and run the dependency check in step 4. |
| Need the dependency versions recorded by the repository | Open the default library and follow the lockfile restore instructions in step 3. |

#### Activate the R 4.4 `compatible` profile

The checked-in lockfile records R 4.4.3 and Bioconductor 3.22. Bioconductor 3.22 targets R 4.5, and its Windows R 4.4 binary index returns HTTP 404. Bioconductor 3.20 provides an R 4.4 repository. The `compatible` profile below installs R 4.4-compatible packages into its own library and uses the GitHub revisions recorded for `ALRA` and `MERINGUE`. Its package versions can differ from those in the development lockfile. [Bioconductor release and R version table](https://bioconductor.org/about/release-announcements/) · [renv profiles](https://rstudio.github.io/renv/articles/profiles.html)

From a terminal in the checkout, launch `R --vanilla` to set the profile before the project loads. On Windows, if `R` is missing from `PATH` or resolves to another version, the default R 4.4.3 executable can be invoked in PowerShell as follows.

```powershell
& 'C:\Program Files\R\R-4.4.3\bin\x64\R.exe' --vanilla
```

At the R prompt, activate the profile. For subsequent sessions, repeat this block to reopen the same library. In RStudio you can set `RENV_PROFILE=compatible` in a project-local `.Renviron` and restart the session with R 4.4.3. Keep the checkout under a short path such as `C:/src/SpatialCellChat` to avoid Windows path-length failures while deeply nested package headers are unpacked.

```r
Sys.setenv(RENV_PROFILE = "compatible")
cran <- c(CRAN = "https://cloud.r-project.org")
options(repos = cran, timeout = 120, renv.bioconductor.version = "3.20")
source(".Rprofile")
project_lib <- renv::paths$library()
renv::settings$bioconductor.version("3.20")
project_lib
```

#### Activate the default project library

Open `SpatialCellChat.Rproj` in RStudio, or start `R` from the checkout. The repository already contains the `renv` startup files. `.Rprofile` loads `renv/activate.R`, which bootstraps the project's `renv` version if needed. If you started R with `--vanilla`, run `source(".Rprofile")` from the checkout to activate the environment manually. If you opened R outside the checkout, restart R from the checkout directory.

#### Verify the active library

Run the following in the R console before installing packages. Confirm that the first entry in `.libPaths()` matches `renv::paths$library()`. On Windows with R 4.4, the default library is under `renv/library/windows/R-4.4/`, and the `compatible` profile uses a path under `renv/library/profiles/compatible/`.

```r
getwd()
R.version.string
stopifnot(file.exists("DESCRIPTION"), file.exists("renv.lock"))
renv::project()
renv::paths$library()
.libPaths()
```

Start a fresh session when changing profiles or R versions, or before reinstalling a loaded package.

### 3. Install dependencies

#### Populate the R 4.4 `compatible` profile

Run these blocks in the `compatible` profile activated in step 2. You can rerun step 4 directly once this library has the required packages.

Install `BiocManager` from CRAN if needed and select repositories for Bioconductor 3.20. Keep `repos` for both dependency installation steps.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager", lib = project_lib, repos = cran, prompt = FALSE)
}
options(repos = cran, renv.bioconductor.version = "3.20")
repos <- BiocManager::repositories(version = "3.20")
```

Install the package's required dependencies from `DESCRIPTION`. The additional packages below cover the build command, the Visium example and the GitHub packages' dependencies. `renv::install()` writes them into the active profile library. On Windows use binary packages for this step to avoid compiling the whole dependency tree.

```r
fields <- c("Depends", "Imports", "LinkingTo")
depends <- read.dcf("DESCRIPTION", fields = fields)[1L, ]
packages <- trimws(sub("\\s*\\(.*$", "", unlist(strsplit(depends, ",", fixed = TRUE))))
packages <- setdiff(unique(packages), c("R", "methods", "ALRA", "MERINGUE"))
extras <- c(
  "devtools", "jsonlite", "SeuratObject", "rsvd", "RANN", "akima",
  "dynamicTreeCut", "geometry", "DT", "clue", "RcppArmadillo", "RcppProgress"
)
package_type <- if (.Platform$OS.type == "windows") "binary" else getOption("pkgType")
renv::install(
  unique(c(packages, extras)), lib = project_lib, repos = repos,
  type = package_type, dependencies = NA, prompt = FALSE
)
```

Install the two GitHub dependencies after their CRAN dependencies are available. `MERINGUE` compiles C++ code in this step.

```r
renv::install(
  c(
    "github::KlugerLab/ALRA@f34d46570b1221179047d9f99c235ed880cc3bae",
    "github::JEFWorks-Lab/MERINGUE@c0d28fec4659e891dd81db620f58318536ce2d58"
  ),
  repos = repos,
  dependencies = FALSE
)
```

#### Restore the versions in `renv.lock`

Use the default `renv` profile for this route. First unset a profile inherited from `.Renviron`, then start a fresh R session from the repository root. `renv::restore()` installs the versions recorded in the lockfile into the active project library. It can download archived sources and compile packages. The lockfile uses Bioconductor 3.22, whose Windows binaries are published for R 4.5. The lockfile therefore needs an R 4.5 installation for a normal Bioconductor binary restore. `renv` will report that the lockfile was generated with R 4.4.3. Full restoration from this lockfile remains an environment-specific path and has not been completed on this R 4.4.3 machine.

```r
Sys.unsetenv("RENV_PROFILE")
renv::paths$lockfile()  # Expect renv.lock at the repository root
renv::restore()
renv::status()
```

A lockfile record with `Source` set to `URL` pins a specific download location. If restore fails on one of those records, check the package name, version and `RemoteUrl` in `renv.lock`. Changing the CRAN mirror only affects repository lookup. A package source archive for the recorded version can be supplied through the project cellar under `renv/cellar/`. [Restore reference](https://rstudio.github.io/renv/reference/restore.html) · [Package sources and the cellar](https://rstudio.github.io/renv/articles/package-sources.html)

### 4. Check the dependencies, then build and install the package

Run the check in the same R session and `renv` profile you intend to use for analysis. It checks the dependencies from `DESCRIPTION` and the build/example packages. An empty `missing` vector means the namespaces load successfully.

```r
fields <- c("Depends", "Imports", "LinkingTo")
depends <- read.dcf("DESCRIPTION", fields = fields)[1L, ]
required <- trimws(sub("\\s*\\(.*$", "", unlist(strsplit(depends, ",", fixed = TRUE))))
required <- setdiff(unique(c(required, "devtools", "jsonlite", "SeuratObject")), "R")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing or unloadable dependencies: ", paste(missing, collapse = ", "))
}
pkgbuild::check_build_tools(debug = TRUE)
```

Build the source archive into `dist/`, then install it explicitly into the active project library. If you received an existing source archive, set `tarball` to its local path, such as `"dist/SpatialCellChatDev_0.1.0.tar.gz"`, then run the `file.exists()` check and the install command.

```r
dir.create("dist", showWarnings = FALSE)
tarball <- devtools::build(pkg = ".", path = "dist", vignettes = FALSE)
stopifnot(file.exists(tarball))
utils::install.packages(tarball, lib = renv::paths$library(), repos = NULL, type = "source")
```

A successful installation prints `* DONE (SpatialCellChatDev)`. The archive in `dist/` is a local build output. A `.tar.gz` source package still needs its dependency library and a compiler at install time. Run these checks in a fresh session with the same profile, then proceed to the Visium example.

```r
library(SpatialCellChatDev)
packageVersion("SpatialCellChatDev")   # 0.1.0 for this checkout
find.package("SpatialCellChatDev")      # Should be inside the active renv library
stopifnot(is.function(SpatialCellChatDev::updateSpatialCellChat))
```

### 5. Diagnose project library problems

| Symptom | Check and action |
| --- | --- |
| `renv` is unavailable or the package resolves to a global library | Check `getwd()` and `.libPaths()`. Activate the checkout from its root and confirm the selected profile before installing anything. |
| `The project is out-of-sync` | Run `renv::status()` to see which versions differ. For a working library, use the namespace and build checks above to establish whether the difference affects installation. To restore specific packages to their recorded versions in the default profile, use `renv::restore(packages = c("package_name"))`. |
| Bioconductor index URL contains `/3.22/` and `/contrib/4.4/` | Use the R 4.4 `compatible` profile with Bioconductor 3.20, or start R 4.5 for the lockfile's Bioconductor release. |
| A dependency is installed but its namespace will not load | Run `loadNamespace("package_name")` directly to see the dependency or DLL error. Check that the current R version and project library are the ones used to build it. |
| `make` or `g++` is missing | Run `pkgbuild::check_build_tools(debug = TRUE)` and restart R after installing the toolchain for your R version. |
| Windows reports an in-use DLL or a `00LOCK-*` directory | Close other R sessions using the project library. If the prior installation has stopped, remove its stale lock directory from that project library before retrying. |

`renv::install()` adds or updates packages in the active profile library. `renv::restore()` uses an existing lockfile to populate libraries. `renv::snapshot()` writes dependency records from a library back to a lockfile. For repository development, review the resulting `renv.lock` diff when intentionally changing dependencies. Routine package builds use the active library as-is. [renv reference](https://rstudio.github.io/renv/reference/index.html)

## Quick start with the bundled Visium data

Run the following from the repository root. The example analyzes CXCL signaling in a human Visium dataset and runs through group-level centrality. [Visium dataset](tutorial/visium_human_psoriasis.RData)

### Create the object

For your own data, supply a gene-by-cell expression matrix, cell metadata with a grouping column, and a two- or three-column coordinate matrix. Match metadata and coordinate row names to the expression column names. Use `input.assay = "raw"` to normalize counts during construction or `input.assay = "norm"` for previously normalized expression. Set `spatial.factors$ratio` to the number of micrometers per coordinate unit and `spatial.factors$tol` to the distance tolerance in micrometers. [Constructor reference](man/createSpatialCellChat.Rd)

```r
library(SpatialCellChatDev)

visium <- new.env()
load("tutorial/visium_human_psoriasis.RData", envir = visium)
seu <- visium$seu
# In this particular fixture, the Spatial data layer holds raw counts.
counts <- SeuratObject::LayerData(seu, assay = "Spatial", layer = "data")
meta <- seu[[]]
pixels <- Seurat::GetTissueCoordinates(
  seu, scale = NULL, cols = c("imagerow", "imagecol")
)
coordinates <- as.matrix(pixels[, c("imagecol", "imagerow")])
colnames(coordinates) <- c("x", "y")
scale_factors <- jsonlite::fromJSON(file.path(
  "tutorial", "spatial_imaging_data_visium_psoriasisSkin", "scalefactors_json.json"
))
spot_diameter_um <- 65

chat <- createSpatialCellChat(
  counts,
  meta = meta,
  group.by = "final.clusters",
  input.assay = "raw",
  datatype = "spatial",
  coordinates = coordinates,
  coordinate.system = list(unit = "pixel", calibrated = TRUE, y.direction = "down"),
  spatial.factors = list(
    ratio = spot_diameter_um / scale_factors$spot_diameter_fullres,
    tol = spot_diameter_um / 2
  )
)
```

### Infer and summarize communication

Select a CellChatDB for the species and choose a feature-selection method for your analysis. The distance settings below are for this example. `interaction.range` and `contact.range` are in micrometers, and `scale.distance` sets the scale for inverse-distance weighting. [Inference prerequisites](docs/COMPUTE_COMMUN_PROB_PREREQUISITES.md) · [Inference parameters](man/computeCommunProb.Rd)

```r
chat@DB <- subsetDB(CellChatDB.human, search = "CXCL", key = "pathway_name")
chat <- subsetData(chat)
chat <- identifyOverExpressedGenes(chat, selection.method = "wilcox")
chat <- identifyOverExpressedInteractions(chat, variable.both = FALSE)
chat <- computeCommunProb(
  chat, interaction.range = 250, contact.range = 10, scale.distance = 0.2
)

chat <- filterProbability(chat)
chat <- filterCommunication(chat, min.cells = NULL)
chat <- computeAvgCommunProb(chat)
chat <- computeCommunProbPathway(chat)
chat <- aggregateNet(chat)
chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = TRUE)

head(subsetCommunication(chat))
print(spatialC2CDistPlot(chat, mode = "geometry"))
print(spatialC2CDistPlot(chat, mode = "communication", split.by = "annotation"))
```

The geometry plot uses geometric pairs within the cached distance range. The communication plot uses the current cell-level probability layers and counts an unordered pair once per panel if either direction has at least one positive LR probability in that panel. Any statistical filtering comes from upstream workflow steps such as `filterProbability()`. [Distance plot reference](man/spatialC2CDistPlot.Rd)

### Use your own data and save results

| Input | What to supply |
| --- | --- |
| Expression | A matrix with gene identifiers in `rownames()` and cell or spot barcodes in `colnames()`. Supply counts with `input.assay = "raw"` or already normalized values with `input.assay = "norm"`. |
| Metadata | One row per barcode. Put the cell type or cluster labels in the column named by `group.by`. |
| Coordinates | One row per barcode, with two or three numeric coordinate columns. Supply the conversion to micrometers as `ratio`. For coordinates already in micrometers, use `ratio = 1`. Choose `tol` from your assay's spatial resolution and keep `contact.range` and `interaction.range` in micrometers. |
| Database | Use `CellChatDB.human`, `CellChatDB.mouse`, or `CellChatDB.zebrafish` for the corresponding species. Use `subsetDB()` to select pathways or annotation categories. The CXCL subset above keeps the example small. |

With a Seurat object, you can also pass it directly to `createSpatialCellChat()`. Set `assay_name` to the expression assay, use `input.assay` to select counts or normalized data, and pass the spatial coordinates and calibration for that assay. The constructor can take metadata from the Seurat object. [Constructor arguments](man/createSpatialCellChat.Rd)

`computeCommunProb()` creates the cell-by-cell LR layers. `filterProbability()` and `filterCommunication()` apply the cell-level filters. `computeAvgCommunProb()` produces the group-by-group LR probabilities and, by default, group p-values. Run `computeCommunProbPathway()` to summarize LR pairs into pathways and `aggregateNet()` for network counts and weights. With `do.group = TRUE`, `netAnalysis_computeCentrality()` computes group-level centrality. Use `do.group = FALSE` for cell-level centrality.

```r
lr_results <- subsetCommunication(chat, slot.name = "net", thresh = 0.05)
pathway_results <- subsetCommunication(chat, slot.name = "netP")
utils::write.csv(lr_results, "communication-lr.csv", row.names = FALSE)
utils::write.csv(pathway_results, "communication-pathway.csv", row.names = FALSE)
saveRDS(chat, "chat-dev.rds")
```

The LR table reports source and target group, LR pair, probability, and group permutation p-value. The pathway table summarizes pathway-level results and reports `NA` for p-values. To analyze multiple conditions, run inference for each dataset separately, then use `mergeSpatialCellChat()` on the named list of completed objects for downstream comparisons. [Result extraction](man/subsetCommunication.Rd) · [Merging datasets](man/mergeSpatialCellChat.Rd)

## Migrate a saved upstream object

In a session using `SpatialCellChatDev`, read the old object, update it, and save the converted object under a new name.

```r
library(SpatialCellChatDev)

old_chat <- readRDS("upstream-chat.rds")
chat <- SpatialCellChatDev::updateSpatialCellChat(old_chat)
methods::validObject(chat)
saveRDS(chat, "chat-dev.rds")
```

`updateSpatialCellChat()` handles the tested 15-slot `SpatialCellChat` and historical 15-slot `SpatialChat` layouts. It also validates an existing 11-slot object. It carries forward expression, spatial inputs, and compatible network results when their axes can be aligned. Legacy distance caches are discarded. Rerun `computeCommunProb()` and the desired downstream steps if you need the communication distance plot or want to re-infer with different spatial settings. A legacy merged object can be migrated directly if `@data` contains normalized expression. Otherwise, migrate the individual objects before merging them in Dev. [Migration reference](man/updateSpatialCellChat.Rd)

## Documentation and license

See the [project documentation](docs/README.md) for the schema, inference, and visualization notes. The [older Visium tutorial](tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.Rmd) provides additional analysis context and uses earlier function signatures. Use the quick start above for the current API. The package is distributed under the [GPL-3 license](LICENSE).
