# Agent Note: 固定 SpatialCellChat 的 R 4.4.3 renv 基线

Status: implemented
Governance: v1
Date: 2026-10-11
Decision type: dependency
Scope: renv.lock and project R library
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The requested local interpreter is R 4.4.3, while the existing `renv.lock` recorded R 4.5.3. A restore under R 4.4.3 initially failed because the lockfile contained an R-4.5-only `Deriv` build, an optional ALRA `fastRPCA` remote that fails on Windows staging, and dependency cascades from the failed roots. The project library consequently contained only bootstrap packages.
## Decision

Use `C:/Program Files/R/R-4.4.3/bin/Rscript.exe` with the activated project library at `renv/library`, install through renv, and record the working result in `renv.lock`. Set the lockfile R version to 4.4.3, keep Bioconductor 3.22, use the Westlake CRAN mirror, pin the installed GitHub packages to MERINGUE commit `c0d28fec4659e891dd81db620f58318536ce2d58` and SeuratObject commit `58bf437fe058dd78913d9ef7b48008a3e24a306a`, and remove ALRA's optional `fastRPCA` Suggests/Remotes trigger. Record `spatstat.utils` as 3.2-5 because the locked `spatstat.geom` 3.8-1 source requires `spatstat.utils >= 3.2-5`.
## Constraints and invariants

- All project packages remain under `renv/library`; the global R library is not the install target.
- `.Rprofile` continues to set `RENV_PATHS_LIBRARY=renv/library` and source `renv/activate.R`.
- CRAN resolution uses `https://mirrors.westlake.edu.cn/CRAN`; Bioconductor remains release 3.22.
- The 237 package records in `renv.lock` remain loadable in the activated R 4.4.3 session.
- MERINGUE and SeuratObject retain their exact recorded GitHub SHAs.
- `fastRPCA` remains optional and is not a lockfile package or restore target.
## Alternatives considered

- Install R 4.5.3 to match the original lockfile: rejected because the user supplied R 4.4.3 as the available interpreter and the project must restore with that executable.
- Keep the R 4.5.3 lockfile and suppress the restore failures: rejected because it leaves the environment unreproducible under the requested interpreter.
- Install optional `fastRPCA` on Windows: rejected because it is only an ALRA Suggests/Remotes path and repeatedly failed during renv graph staging; removing that optional trigger avoids an unrelated package from blocking the required environment.
## Consumer impact

- Producers: `renv::restore()`/`renv::install()` commands, the Westlake CRAN repository, Bioconductor 3.22 sources, and the pinned GitHub package sources.
- Consumers: `.Rprofile`, `renv/activate.R`, `renv/settings.json`, `renv.lock`, the project renv library, and any R script launched from this project.
- Package consumers: SpatialCellChat's declared CRAN, Bioconductor, and GitHub dependencies, including Seurat, MERINGUE, ALRA, ComplexHeatmap, and the spatstat family.
- Tests and smoke checks: `renv::status()`, the full lockfile package load check, and GitHub SHA assertions.
- No source API or runtime package code was changed.
## Consequences

The project now has a reproducible R 4.4.3-specific lockfile rather than an unusable R 4.5.3 declaration. The lockfile is synchronized and the complete 237-package set loads from the project library. Future restoration must use R 4.4.x or deliberately create a new lockfile for another R release. The optional ALRA fastRPCA path is not available, but ALRA itself remains installed and loadable; code that explicitly requests that optional method must install it as a separate, reviewed dependency.
## Evidence

- `.Rprofile:1-2` sets the project library and activates renv; `renv/settings.json:1-20` defines explicit snapshots and dependency fields.
- `renv.lock:1-32` records R 4.4.3, the Westlake CRAN mirror, and Bioconductor 3.22; `renv.lock:35-68` records ALRA without the fastRPCA trigger.
- The final `renv::status()` run reported `No issues found -- the project is in a consistent state`, `SYNCHRONIZED=TRUE`, 237 package records, and lock R 4.4.3.
- The full smoke command loaded all 237 lockfile packages from `D:/AIProjects/SpatialCellChat/renv/library/windows/R-4.4/x86_64-w64-mingw32` and returned `SMOKE_OK`.
- The smoke command asserted MERINGUE SHA `c0d28fec4659e891dd81db620f58318536ce2d58` and SeuratObject SHA `58bf437fe058dd78913d9ef7b48008a3e24a306a`.
- Failed restore logs documented the original roots: `Deriv`, fastRPCA staging, and downstream dependency cascades; exact source installs and the final metadata snapshot resolved them.
## Acceptance criteria

- `renv::status()$synchronized` is `TRUE` under the supplied R 4.4.3 executable.
- `renv.lock` records R 4.4.3, Westlake CRAN, Bioconductor 3.22, 237 packages, and no `fastRPCA` or `rpca-mkl` entry.
- Every package name in `renv.lock` passes `requireNamespace(..., quietly=TRUE)` from the activated project library.
- MERINGUE and SeuratObject package descriptions report the exact SHAs above.
- The project library path printed by the smoke check is the repository's `renv/library`, not a global library.
- The governance validator passes with no new v1 errors.

