# Agent Note: 固定上游 SpatialCellChat 比较基线

Status: implemented
Governance: v1
Date: 2026-08-23
Decision type: dependency
Scope: renv library, renv.lock, current source comparison
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The project renv library cannot load `SpatialCellChat` or `CellChat`, while `tests_dev/simulated_data/SpatialChat_1.rds` contains a legacy object of class `SpatialCellChat`. The object cannot be validated or upgraded until a compatible upstream package and its hard dependencies are available.

## Decision

Choose the `jinworks/SpatialCellChat` main-branch commit `be88f300464cf970b36a6165ea7d2e41a47a7511` as the comparison baseline and install it plus its hard dependencies into the current project renv library with pak; reject installing the separate `jinworks/CellChat` repository because this task specifies the SpatialCellChat repository.

## Constraints and invariants

- R 4.5.3 and the activated project renv library remain the execution environment.
- pak is the only installer for CRAN, Bioconductor, and GitHub packages; renv only snapshots installed results.
- `dependencies = NA` resolves hard dependencies only; upstream Suggests remain outside this phase.
- The upstream package remains version `0.1.0` with `RemoteSha` exactly `be88f300464cf970b36a6165ea7d2e41a47a7511`.
- The local source package is never installed over the same-library upstream baseline.
- The legacy fixture must read as `SpatialCellChat`, upgrade through `updateSpatialCellChat()`, and pass `methods::validObject()`.

## Alternatives considered

- Install the separate `jinworks/CellChat` repository: rejected because it is not the specified comparison repository and does not provide the SpatialCellChat object upgrade boundary.
- Use `renv::install()` or `devtools::install_github()`: rejected because the approved procedure requires pak semantics and keeps renv limited to activation and lockfile recording.
- Install all upstream Suggests: deferred because visualization, SCE, topic/NMF, and Shiny paths are outside the current core compatibility goal.

## Consumer impact

- Producers: pak installation commands and the upstream fixed commit metadata.
- Consumers: the activated renv library, `renv.lock`, `tests_dev/simulated_data/SpatialChat_1.rds`, the upstream `updateSpatialCellChat()` entry point, and current source-based regression tests.
- Configuration and declarations: `DESCRIPTION`, `renv/settings.json`, `.Rprofile`, and `renv.lock`.
- Documentation/evidence: the pinned upstream `DESCRIPTION`, README, commit API response, current `DESCRIPTION`, `renv::status()` output, and fixture compatibility results.

## Consequences

The renv lockfile gains pak, the upstream comparison package, ten missing hard dependencies, and the explicitly requested already-installed dependencies `sf`, `circlize`, `shape`, and `ks`. The legacy fixture gains a tested load-and-upgrade path. Optional visualization and extension dependencies remain unavailable by design and must be installed in a later phase if those paths are exercised. The same-name upstream/local package boundary requires independent R sessions for comparison versus source tests.

## Evidence

- `DESCRIPTION:1-27` and `renv/settings.json:1-20` define the current package dependency fields and explicit snapshot policy.
- `renv.lock` records R 4.5.3, Bioconductor 3.22, and the current locked package baseline.
- `tests_dev/simulated_data/SpatialChat_1.rds` is the legacy fixture; `tests_dev/test-SpatialCellChat-class.R` and `tests_dev/test-SpatialCellChat-spatial-image.R` are the current source regression entries.
- Upstream fixed-commit `DESCRIPTION`, README, commit API response, and `R/SpatialCellChat_class.R` provide the dependency declaration, baseline identity, and `updateSpatialCellChat()` behavior.
- pak package-source and `pkg_install()` documentation define GitHub refs, library targeting, and hard-dependency resolution.
- Preflight command confirmed R 4.5.3, Bioconductor 3.22, pak 0.10.0, the activated renv library path, and `pkgbuild::has_build_tools() == TRUE`.
- pak install output confirmed CRAN hard dependencies, Bioconductor `ComplexHeatmap` 2.26.1, GitHub `ALRA` and `MERINGUE`, followed by `CORE_OK`; the first MERINGUE download was retried successfully and ALRA was retried after the combined transaction aborted.
- The pinned pak install reported `UPSTREAM_OK`; `packageDescription()` returned version `0.1.0`, `RemoteRepo=SpatialCellChat`, `RemoteUsername=jinworks`, and the exact requested SHA.
- Legacy verification returned `OLD_CLASS=SpatialCellChat`, `UPGRADED_CLASS=SpatialCellChat`, and `LEGACY_OBJECT_OK`; upstream emitted the expected warning adding `meta$samples` as `sample1`.
- `tests_dev/test-SpatialCellChat-class.R` and `tests_dev/test-SpatialCellChat-spatial-image.R` both exited 0 with all applicable checks passing; the spatial image test skipped only optional `scatterpie` proportion mode.
- Explicit snapshot with `options(renv.snapshot.ignore.self=FALSE)` recorded all requested targets, including `SpatialCellChat` as GitHub `be88f300464cf970b36a6165ea7d2e41a47a7511`.
- `renv::status()` with self inclusion reported only unrelated legacy `renv` drift plus target-independent `pak` usage state; no requested target had `installed = n` or `recorded = n`.
- Completed evidence satisfies the ten core dependency, pinned upstream metadata, legacy migration, source regression, lockfile, and governance validation gates.


## Acceptance criteria

- All ten named hard dependencies and upstream `SpatialCellChat` load with `requireNamespace(..., quietly = TRUE)` under R 4.5.3.
- `packageDescription("SpatialCellChat")$Version` is `0.1.0` and `$RemoteSha` equals `be88f300464cf970b36a6165ea7d2e41a47a7511`.
- The fixture passes `readRDS()`, `SpatialCellChat::updateSpatialCellChat()`, class checks, and `methods::validObject()`.
- The two named current source regression scripts exit successfully without loading the upstream package.
- `renv.lock` contains the requested package set and `renv::status()` has no `installed = n` or `recorded = n` for those targets; unrelated legacy drift may remain documented.
- The Agent Notes validator reports no new v1 errors before this note transitions from `proposed` to `implemented`.

