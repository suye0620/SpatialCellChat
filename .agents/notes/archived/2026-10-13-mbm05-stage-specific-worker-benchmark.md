# Agent Note: MBM05 stage-specific upstream-refactored benchmark

Status: superseded
Governance: v1
Date: 2026-10-13
Decision type: performance
Scope: tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R and tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers
Owner: project maintainer
Impact: medium
Supersedes: .agents/notes/archived/2026-10-12-mbm05-benchmark-boundary.md
Superseded by: .agents/notes/implemented/2026-10-13-mbm05-rscript-baseline-compatibility.md

## Problem

The original MBM05 Slide-seq analysis does not use one future worker count for every inference stage. `tests_dev/test_data/Slide-seq.R` selects four workers before `computeCommunProb`, switches to two before `filterProbability`, and uses two for later permutation and centrality work. The interrupted single-worker comparison does not reproduce that stage-specific upstream execution profile and did not save an upstream result signature. The retained two-worker benchmark records a `filterProbability` failure while exporting 37.92 GiB of globals.

## Decision

Run the upstream baseline with four future workers for `computeCommunProb` and two workers for `filterProbability`, `computeAvgCommunProb`, and centrality; keep the refactored worker sequential with `nthreads = 1` and `parallel = FALSE`. Preserve the interrupted serial attempt and prior two-worker artifacts in their existing directories, and store this stage-specific comparison in `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers`.

## Constraints and invariants

- Each implementation runs in an independent R 4.4.3 process under the project renv library.
- The refactored implementation remains sequential and uses one compute thread.
- Upstream uses four future workers only for `computeCommunProb`; it switches to two before `filterProbability` and retains two through averaging, pathway/aggregate operations, and both centrality calls.
- Both workers use the same pre-inference checkpoint components, DB, LR table, metadata, coordinates, spatial cache, and selected features.
- `scale.distance = 9`, `Kh = 0.5`, `n = 1`, interaction range 250, contact range 10, AGAN, permutation counts, seeds, and centrality settings remain unchanged.
- The completed two-worker artifacts in `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored` and the current serial-attempt artifacts in `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_serial` remain untouched.
- The upstream two-worker `filterProbability` 37.92 GiB export failure is a known risk; report a structured failure, never infer unmeasured later stages.

## Alternatives considered

- Keep one upstream worker for every stage: rejected because it does not follow the stage-specific worker selections in the requested Slide-seq reference and is taking many hours on this fixture.
- Run four workers through every stage: rejected because the reference explicitly switches to two before filtering, and the observed global-export failure makes a higher downstream worker count an unjustified memory risk.
- Reuse the retained all-stage two-worker result: rejected because it already fails at `filterProbability` and cannot measure later stages.

## Consumer impact

- Producer: `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` and `MBM05Chat_pre_inference.rds`.
- Configuration source: `tests_dev/test_data/Slide-seq.R` stage-level `setEnvironment` calls.
- Consumers: stage-specific logs, worker result RDS files, CSV/TXT comparison summaries, output signatures, and the maintainer's performance conclusion.
- Historical artifacts: the original two-worker output directory and current serial-attempt directory remain available and unchanged.
- No package API, production object schema, inference algorithm, or raw fixture changes.

## Consequences

The measurements reflect the worker counts in the Slide-seq analysis rather than a worker-matched comparison. Refactored `computeCommunProb` completed in 286.90 seconds; four-worker upstream completed in 6,993.96 seconds, so the observed ratio is 0.0410 under intentionally asymmetric execution plans. Upstream two-worker `filterProbability` then reproduced the 37.92 GiB globals-export failure. Later upstream timings and the final upstream signature are unavailable and remain `NA`; no post-filter equivalence or performance claim is valid.

## Evidence

- `tests_dev/test_data/Slide-seq.R:127-147` sets four workers before `computeCommunProb`, then two before `filterProbability`; the count remains two for `computeAvgCommunProb`.
- `tests_dev/test_data/Slide-seq.R:201-207` sets two workers before both centrality calls.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/baseline.log` records upstream `computeCommunProb` at 12,349.42 seconds and the subsequent two-worker `filterProbability` failure with 37.92 GiB of globals.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_serial/new.log` records successful refactored stages and timings.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_serial/baseline.log` records all eight timed stages, but no `baseline_result.rds` exists; this attempt cannot establish final numerical equivalence.
- `.agents/notes/rejected/2026-10-12-mbm05-serial-benchmark-boundary.md` preserves the rejected uniform-serial proposal and its incomplete-result boundary.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/historical_checksums.json` records 14 original/serial artifacts; all 14 checksums matched during the stage-specific run.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/original_filter_export_failure/summary.csv` records eight successful refactored stages, upstream `create_normalize`, and four-worker upstream `computeCommunProb` at 6,993.96 seconds.
- At the last common completed inference stage, refactored versus upstream retained-object size is 5.343 versus 14.000 GB, boundary RSS is 6,702 versus 13,580 MB, and Vcell high-water mark is 10,884 versus 22,652 MB; both contain 1,563 layers and 407,065,015 pre-filter nnz.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/original_filter_export_failure/baseline.log` records the four-worker compute plan, the two-worker filter plan, and the structured `FutureError` at `filterProbability`.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/new_result.rds` contains eight stages and a full signature; `original_filter_export_failure/baseline_result.rds`, recovered identically from the preserved `comparison.rds`, contains two measured stages, the complete configured worker plan, and no fabricated signature.
- The current refactored signature matches the retained serial and original refactored signatures at numeric tolerance `1e-10`.

## Acceptance criteria

- Satisfied: the benchmark script selects the configured upstream worker plan immediately before each corresponding stage and records the complete plan in `baseline_result.rds`.
- Satisfied: the stage-specific run wrote only to its separate output directory; all 14 checksums captured for prior two-worker and serial-attempt artifacts still match.
- Satisfied: both result artifacts report their execution plans; `baseline.log` confirms four workers for `computeCommunProb` and two for `filterProbability`.
- Satisfied: `new_result.rds` contains eight measured stages and a full signature; `baseline_result.rds` contains two measured stages plus a structured `filterProbability` failure and no downstream measurements.
- Satisfied: `comparison.txt`, `summary.csv`, and `comparison.csv` report the asymmetric worker counts; no worker-matched speedup is claimed.
- Satisfied: recomputed SHA-256 checksums for both pre-existing output directories match all 14 captured values after the stage-specific run.
- Satisfied: the Agent Notes validator reports zero errors after the implemented/rejected/superseded lifecycle transitions.
