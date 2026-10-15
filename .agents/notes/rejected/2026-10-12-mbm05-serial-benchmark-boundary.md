# Agent Note: MBM05 serial upstream-refactored benchmark boundary

Status: rejected
Governance: v1
Date: 2026-10-12
Decision type: performance
Scope: tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R and tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_serial
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The retained two-worker upstream experiment completes `computeCommunProb` but fails in `filterProbability` while exporting 37.92 GiB of globals. The requested follow-up must remove multisession export and worker-count asymmetry by comparing upstream sequential execution with refactored sequential execution, without overwriting the completed two-worker artifacts.

## Rejection
Reject the uniform serial profile as the active benchmark decision in favor of the requested Slide-seq four-then-two-worker profile. The serial refactored worker completed all eight stages and wrote its result; the upstream log records all eight timed stages, but that process was stopped before writing `baseline_result.rds` and the parent comparison artifacts. Retain `--baseline-workers=1` as an optional execution mode and preserve the incomplete serial directory; do not claim final numerical equivalence from timings alone.
## Constraints and invariants

- The original `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored` files and their checksums remain unchanged.
- Both implementations run in independent R 4.4.3 processes with `future::sequential`; the refactored compute kernels additionally retain `nthreads = 1` and `parallel = FALSE`.
- Both workers consume the same pre-inference checkpoint components, DB, LR table, metadata, coordinates, spatial cache, and selected features.
- `scale.distance = 9`, `Kh = 0.5`, `n = 1`, interaction range `250`, contact range `10`, AGAN, permutation counts, seeds, and centrality configuration remain unchanged.
- Raw preprocessing remains outside the timed boundary, as in the retained experiment.

## Alternatives considered

- Overwrite the two-worker output directory: rejected because it would destroy the requested historical result and obscure the different execution plan.
- Keep upstream at two workers: rejected because it preserves the 37.92 GiB multisession export failure and violates the requested one-worker comparison.
- Reuse the two-worker upstream `computeCommunProb` measurement: rejected because worker-count-dependent elapsed time and memory are not a serial baseline, and the failed run did not retain its intermediate object.

## Consumer impact

- Producer: `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` and `MBM05Chat_pre_inference.rds`.
- Consumers: the serial output directory's logs, result RDS files, CSV summaries, correctness signatures, and the maintainer's comparison conclusion.
- Historical consumer: the original two-worker output directory and its implemented benchmark Note remain unchanged by this rejected proposal; the separate stage-specific successor owns any later supersession.
- No package API, production object schema, inference algorithm, or raw fixture changes.

## Consequences

The serial attempt avoids multisession export and records worker-count-matched stage timings, but lacks an upstream result and final output signature. No complete serial correctness conclusion is available. The active default instead follows the stage-specific Slide-seq worker profile; selecting `--baseline-workers=1` remains explicit and uses the preserved serial directory.

## Evidence

- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/failures.csv` records the prior upstream `filterProbability` `FutureError` and 37.92 GiB globals export.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored/comparison.csv` records the prior two-worker common-stage measurements.
- `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` accepts `--baseline-workers=1`, selects `future::sequential` for that profile, and isolates serial output.
- SHA-256 values for every retained two-worker artifact were captured before starting the serial run.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_serial/new_result.rds` has `status = "success"` and eight measured stages; the directory has no `baseline_result.rds`.
- The serial upstream log ends after `centrality_netP` at 103.95 seconds without `WORKER_OK`; no upstream final signature was saved.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/historical_checksums.json` captures all 14 artifacts across the original and serial output directories.

## Acceptance criteria

- Rejected: the serial run did not produce both successful worker result artifacts or the complete comparison set.
- Satisfied: the source parses worker counts and rejects non-positive/non-integer values.
- Satisfied: the refactored serial log records all eight timed stages; the upstream serial log records all eight timed stages before interruption.
- Satisfied: the original two-worker directory remained unchanged during the serial attempt.
- Satisfied: the serial output directory remains available as an incomplete historical attempt, not as a complete numerical conclusion.
- Satisfied: the Agent Notes validator reports zero errors after moving this rejected proposal to `rejected/`.
