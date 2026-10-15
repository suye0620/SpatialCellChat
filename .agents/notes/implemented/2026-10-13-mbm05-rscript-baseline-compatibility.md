# Agent Note: MBM05 upstream Rscript baseline compatibility

Status: implemented
Governance: v1
Date: 2026-10-13
Decision type: performance
Scope: tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R, tests_dev/test_data/debug_upstream_filterProbability_rscript.R, and the pinned upstream SpatialCellChat filter/pathway execution boundary
Owner: project maintainer
Impact: medium
Supersedes: .agents/notes/archived/2026-10-13-mbm05-stage-specific-worker-benchmark.md
Superseded by: none

## Problem
The pinned upstream package `be88f300464cf970b36a6165ea7d2e41a47a7511` imports `Matrix` wholesale. Its legacy `filterProbability` uses unqualified `which` inside a nested `future.apply` closure. Under the project R 4.4.3 renv runtime, that resolves to Matrix's `standardGeneric`; the package's `my_future_lapply` wrapper also retains the complete filter frame. On Windows Rscript multisession, the resulting export reports 37.92 GiB of globals and fails while launching `filterProbability`. At `computeCommunProbPathway`, the upstream object retains the full cell-probability tensor and duplicate per-layer matrices in `net$tmp$prob.cell`, creating an avoidable peak-memory copy.

## Decision
Run the pinned upstream package unchanged, but use benchmark-local wrappers at two legacy memory/runtime boundaries.

For `filterProbability`, resolve `which` to `base::which` and replace only the future wrapper with a worker closure parented by the SpatialCellChat namespace. Capture only `d.spatial`, `adj.contact`, `permutation`, `prob.cell_`, `nLR1`, and `quantile.prob`.

For pathway aggregation, replace duplicate `net$prob.cell` with an empty `sparse3Darray` shell retaining its dimensions and dimnames. Call the original `computeCommunProbPathway`, then restore the canonical tensor from `net$tmp$prob.cell`; verify dimensions, dimnames, nnz, and sum.

Close the idle multisession pool before synchronous pathway/aggregate stages, report sequential plans for them, and restore two workers before centrality. Keep the upstream algorithm, parameters, and refactored worker unchanged.

## Constraints and invariants
- The pinned upstream package and its recorded GitHub SHA remain unchanged.
- Both compatibility wrappers are benchmark-local; production source and package API are not modified.
- The upstream plan is one worker for object construction, four for `computeCommunProb`, two for `filterProbability` and `computeAvgCommunProb`, sequential for synchronous `computeCommunProbPathway` and `aggregateNet` after closing idle workers, then two for both centrality calls.
- `base::which` and `Matrix::which` return identical indices for sparse logical rows used by the legacy filter; the isolated worker captures only the six listed data globals.
- The pathway shell preserves tensor dimensions and dimnames; restoration reapplies the original dimnames because list-to-sparse conversion retains values/dimensions but does not reconstruct them. Restoration preserves nnz and sum, and wrapped pathway outputs match the original on a deterministic fixture.
- The refactored worker remains sequential with one compute thread.

## Alternatives considered
- Modify the pinned package source or renv installation: rejected because it would invalidate the recorded upstream SHA and contaminate the oracle.
- Force the entire upstream baseline to sequential execution: rejected because it no longer reproduces the Slide-seq worker profile and was already rejected as the active benchmark decision.
- Only raise `future.globals.maxSize` or disable reference checks: rejected because the failure is a closure/global-capture and connection-serialization problem, not a configured size limit; it would retain the 14 GB object capture.
- Keep the two-worker multisession pool alive through synchronous pathway and aggregate stages: rejected because those stages do not dispatch futures, while idle Windows workers retain process memory during the peak-memory pathway operation.

## Consumer impact
- Producer: `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` and the pinned `SpatialCellChat` package.
- Runtime consumers: the upstream baseline worker, its `baseline.log`, `baseline_result.rds`, comparison CSV/TXT/RDS artifacts, and the stage-specific governance Note.
- Regression consumer: `tests_dev/test_data/debug_upstream_filterProbability_rscript.R`, which exercises runtime diagnostics, sparse-row equivalence, original-versus-cloned filter output, and original-versus-wrapped pathway outputs.
- No production package API, serialized object schema, or refactored worker changes.

## Consequences
The benchmark preserves the upstream algorithms as the oracle while adding Windows/Rscript-specific execution and peak-memory boundaries. The filter clone reduces exported globals and removes the non-exportable S4 generic. The pathway wrapper avoids holding both cell-probability representations while `computeCommunProbPathway` builds pathway matrices; it restores the canonical cell tensor before `aggregateNet`, which consumes that tensor. These wrappers require re-audit if the pinned upstream filter body, pathway function, or `net$tmp$prob.cell` representation changes. Full-run failures remain structured and no unmeasured downstream values are inferred.

## Evidence
- `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` lines containing `make_upstream_filter_probability_rscript_safe` and `[BASELINE DEBUG]` implement the benchmark-local shim.
- `tests_dev/test_data/benchmark_mbm05_upstream_vs_refactored.R` implements both benchmark-local wrappers and records actual per-stage worker plans.
- `tests_dev/test_data/debug_upstream_filterProbability_rscript.R` provides verbose Rscript diagnostics, sparse-row equivalence, original-versus-cloned filter output equality, and original-versus-wrapped pathway output equality.
- `renv/library/windows/R-4.4/x86_64-w64-mingw32/SpatialCellChat/DESCRIPTION` records the pinned SHA; the R 4.4.3 runtime is printed by the diagnostic.
- `renv/library/windows/R-4.4/x86_64-w64-mingw32/SpatialCellChat/NAMESPACE` imports `Matrix` wholesale; runtime probes show package-namespace `which` is `standardGeneric`, while `base::which` is a regular function.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/original_filter_export_failure/` preserves the original stage-specific filter-export failure log, comparison artifacts, and upstream result recovered identically from its saved `comparison.rds`.
- A direct Rscript inspection of the pinned `computeCommunProbPathway` confirms that it reads `net$prob.cell` for dimensions/dimnames and builds cell-level pathway probabilities from `net$tmp$prob.cell`; the wrapper retains this metadata and restores the canonical tensor after aggregation.
- Direct Rscript inspection of `aggregateNet` confirms it reads `net$prob.cell` to calculate cell-level counts and weights, so the wrapper restores the canonical tensor before that stage.
- The focused deterministic fixture runs the original and cloned filter under two-worker multisession, compares outputs, and checks the pathway's restored sparse tensor and pathway outputs.
- `Rscript.exe tests_dev/test_data/debug_upstream_filterProbability_rscript.R` and the same command with `--vanilla` both print `Rscript_FILTER_DIAGNOSTIC_OK` and `Rscript_PATHWAY_MEMORY_DIAGNOSTIC_OK`; both filter/pathway output comparisons and the exact restored-cell-tensor comparison return `TRUE`.
- `tests_dev/benchmarks/out/mbm05_upstream_vs_refactored_stage_workers/pathway_dimnames_failure/baseline_result.rds` preserves the first full rerun: four-worker inference completed in 7,775.69 seconds, two-worker filtering in 6,413.62 seconds, and two-worker averaging in 3,346.51 seconds. It then failed at `computeCommunProbPathway` with `Restored upstream cell tensor differs from its pre-pathway representation`. The wrapper now reapplies saved dimnames before checking invariants; the focused fixture includes non-null cell and layer names and passes.
- SHA-256 verification after resuming the final rerun matched all 14 entries in `historical_checksums.json`; the prior uniform two-worker and serial-attempt artifacts remain unchanged.
- The latest full baseline attempt completed all eight measured stages (`create_normalize`, `computeCommunProb`, `filterProbability`, `computeAvgCommunProb`, `computeCommunProbPathway`, `aggregateNet`, `centrality_net`, and `centrality_netP`) with elapsed times 1.52, 7141.87, 3104.44, 1970.34, 43.39, 37.84, 495.11, and 142.80 seconds, but its result artifact was marked failed while assembling the final signature; `signature_tail_failure/baseline_result.rds` preserves that artifact. The malformed message `unexpected ',' in \"(\\\"WORKER_FAILED:\\",\\\"` was emitted by the generic error-handler path and did not identify a measured stage.
- The benchmark now reads the numeric `x` payload from `unclass(sparse3Darray)[["x"]]` for nnz/sum metrics, logs signature assembly explicitly, and labels that boundary as `result_signature`; `BENCHMARK_PARSE_OK` and `SPARSE_STATISTICS_OK` pass, and the focused filter/pathway diagnostic still prints both success markers.
- The first post-fix full rerun (launched 16:06 via an in-session background job) was killed at 19:46 by a Windows reboot (`LastBootUpTime` 2026-10-14 19:46:33); it died mid-`computeCommunProb` and wrote no result artifact, so no benchmark state changed.
- Two later detached reruns were killed by console-level Ctrl+C broadcast from the operator terminal (`^C` recorded at ~22:46 for the 20:00 launch and ~23:01 for the 22:47 launch); the surviving run launched via `Start-Process -WindowStyle Hidden` with its own console (`launch_baseline_run4.ps1`, PID saved in `baseline_run4.pid`, logs in `baseline_run4.log`/`baseline_run4.err`).
- The completed full upstream baseline run (launched 2026-10-14 22:47, artifacts written 2026-10-15 02:40:27) reports `status = "success"` with no failure field, all eight stages with worker plan 1/4/2/2/1/1/2/2, and elapsed times 3.44, 7060.07, 3069.20, 2093.79, 47.03, 38.55, 495.29, and 191.35 seconds; `baseline_result.rds` supersedes the preserved `signature_tail_failure` artifact.
- The driver rebuilt `summary.csv`, `comparison.csv`, `comparison.rds`, `comparison.txt`, `signatures.csv`, and `signature_comparison.csv` from both completed artifacts with `Worker failures: none`; signatures agree bit-identically on cell level (nnz 285,801,026, sum 14031.85), pathway level (nnz 89,174,034, sum 14031.85), centrality dimensions, and the pathway list, and differ only on `group_nnz` (10241 vs 10352) and `group_sum` (2.636712 vs 2.640221). Because the filter-stage and pathway-stage nnz/sum are identical between workers, the compatibility shim has no numeric effect; the group-average divergence is owned by the proposed Note `2026-10-15-mbm05-group-avg-signature-divergence.md`.

## Acceptance criteria

- The focused diagnostic passes under `Rscript.exe` and `Rscript.exe --vanilla`, printing runtime/package/search-path evidence.
- The compatibility filter clone resolves `which` to `base::which`, exports only the six listed data globals, and produces output identical to the original function on a deterministic fixture.
- The pathway shell/restoration preserves dimension, dimnames, nnz, and sum; the full restored cell tensor and wrapped pathway outputs are identical to the original fixture.
- Stage rows and execution metadata report four workers for `computeCommunProb`, two for filtering/averaging, sequential for pathway/aggregate, and two for centrality.
- The full stage-specific baseline worker completes the downstream stages, or writes a structured failure artifact with its exact stage and message; no unmeasured downstream values are inferred.
- The Agent Notes validator reports zero errors after this Note transitions to implemented or rejected.
