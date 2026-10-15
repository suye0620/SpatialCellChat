# Agent Note: MBM05 upstream-refactored benchmark boundary

Status: superseded
Governance: v1
Date: 2026-10-12
Decision type: performance
Scope: tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R and tests_dev/benchmarks/out/mbm05_upstream_vs_refactored
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: .agents/notes/implemented/2026-10-13-mbm05-rscript-baseline-compatibility.md

## Problem

The MBM05 Slide-seq fixture has 29,536 cells and 1,563 LR pairs. Comparing two independent marker-selection/preprocessing runs would mix algorithm performance with different inference inputs, while loading both implementations in one R process is invalid because they register incompatible `SpatialCellChat` S4 schemas.

## Decision

Compare the pinned upstream `SpatialCellChat` commit `be88f300464cf970b36a6165ea7d2e41a47a7511` and the refactored source in separate R 4.4.3 worker processes. Use `MBM05Chat_pre_inference.rds` as the common pre-inference boundary, measure normalization plus object construction separately, and run `computeCommunProb`, `filterProbability`, `computeAvgCommunProb`, `computeCommunProbPathway`, `aggregateNet`, and both cell-level centrality calls with the Slide-seq parameters. Record elapsed time, retained object size, R Vcell high-water mark, boundary RSS, sparse nnz, output signatures, and a structured failure when an implementation cannot finish a stage.

## Constraints and invariants

- The pinned upstream package remains the exact GitHub SHA recorded in the environment.
- Workers run in independent R processes; no same-name upstream/local S4 class is loaded together.
- Both downstream workers consume the same expression matrices, metadata, spatial cache, database, LR table, and selected features extracted from the checkpoint.
- `scale.distance = 9`, `Kh = 0.5`, `n = 1`, `interaction.range = 250`, contact range `10`, AGAN enabled, 100 filter/average permutations, and centrality `degree.only = TRUE` remain fixed.
- The refactored worker uses sequential execution and one OpenMP thread; the upstream worker uses the Slide-seq-style future multisession plan with two workers, making the legacy duplication cost explicit on Windows.
- The checkpoint was generated from the raw and spatial files named by the user; the benchmark script verifies those files exist. Raw preprocessing is intentionally outside this timing run.

## Alternatives considered

- Run full raw preprocessing independently for both implementations: rejected as the default comparison because the current and pinned APIs select features differently and the old dense preprocessing path can dominate or exhaust memory before the requested downstream functions are compared.
- Load upstream and refactored code in one R session: rejected because both define `SpatialCellChat` with incompatible slot schemas.
- Compare only sparse tensor sizes: rejected because it cannot identify time spent in averaging, pathway aggregation, or centrality and cannot establish whether retained memory is lower.

## Consumer impact

- Producer: the MBM05 benchmark script and the existing pre-inference checkpoint.
- Consumers: benchmark logs, CSV summaries, correctness signatures, and the maintainer's memory/performance conclusion.
- Configuration: the activated project renv library, R 4.4.3, and the pinned upstream package.
- No public package API, source schema, or production runtime behavior changes.

## Consequences

The result isolates downstream algorithmic cost without claiming identical raw-input preprocessing. The refactored worker completes the full pipeline. The upstream worker completes `computeCommunProb` but cannot launch `filterProbability` under its two-worker multisession plan because that stage attempts to export 37.92 GiB of globals; later upstream timing and correctness fields therefore remain explicitly unavailable rather than inferred.

This decision is retained as the historical uniform two-worker experiment. The successor uses the stage-specific four-then-two-worker profile from `Slide-seq.R`; it does not rewrite these measurements or artifacts.

## Evidence

- `tests_dev/test_data/Slide-seq.R:118-212` defines the requested downstream order and parameter defaults.
- `tests_dev/test_data/test-cli-pipeline.R:61-95` defines the checkpoint-producing preparation boundary.
- `tests_dev/test_data/MBM05Chat_pre_inference.rds` validates under the refactored 11-slot class and contains 29,536 cells, 967 signaling genes, 1,563 LR pairs, and the cached spatial distance.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/summary.csv` records the complete refactored run and the two completed upstream stages.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/comparison.csv` records that refactored `computeCommunProb` used 273.66 seconds, a 5.343 GB retained object, 6,680 MB boundary RSS, and 10,878 MB Vcell high-water mark; upstream used 12,349.42 seconds, 14.000 GB, 15,286 MB, and 22,652 MB respectively.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/failures.csv` records the upstream `filterProbability` `FutureError`, including the 37.92 GiB globals export and its 18.93 GiB function, 13.97 GiB object, and 4.72 GiB probability-list components.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/signatures.csv` records the refactored 1,563 cell/group layers, 198 pathway layers, 285,801,026 retained cell nnz, and both centrality dimensions; upstream signature fields are unavailable because its pipeline failed before aggregation.
- `python .agents/skills/agent-notes-governance/scripts/validate_agent_notes.py --root .agents/notes` validates this governed Note with zero errors; unrelated legacy Notes produce warnings.
- `.agents/notes/archived/2026-10-13-mbm05-stage-specific-worker-benchmark.md` 记录其执行剖面与隔离输出目录，其记录由 `.agents/notes/implemented/2026-10-13-mbm05-rscript-baseline-compatibility.md` 的成功全量基线进一步落实。
## Acceptance criteria

- Satisfied: the workers ran in independent R 4.4.3 processes under the activated project renv library.
- Satisfied: `summary.csv`, `comparison.csv`, `failures.csv`, `signatures.csv`, and `comparison.rds` contain the measured output and structured upstream failure.
- Satisfied: the refactored worker completed all named stages; the upstream worker produced a failure artifact for `filterProbability` after completing `computeCommunProb`.
- Satisfied where execution permits: the refactored signatures cover cell, group, pathway, and centrality outputs; upstream fields are explicitly `NA` after its fatal stage failure.
- Satisfied: at the last common completed stage, refactored `computeCommunProb` uses 38.2% of the retained object size, 43.7% of boundary RSS, and 48.0% of the Vcell high-water mark of upstream; no post-filter upstream comparison is claimed.
