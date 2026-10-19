# Agent Note: computeCommunProb 并行输出接线

Status: implemented
Governance: v1
Date: 2026-10-19
Decision type: migration
Scope: R/modeling.R computeCommunProb parallel result and spill paths
Owner: maintainer
Impact: medium
Supersedes: .agents/notes/archived/2026-10-19-inference-chain-cli.md
Superseded by: none

## Problem

`computeCommunProb(parallel=TRUE)` computes its LR layers into `Prob.cell_`, but the common output stage consumes the separately initialized `layers` list. With `spill.dir`, the parallel branch also returns without writing the per-layer RDS files that `filterProbability()` requires. These are observable source-level contract breaks: in-memory parallel output cannot be assembled, and parallel spill output cannot be consumed.

## Decision


Make the parallel branch use the same output contract as the serial branch without changing the LR kernel, public arguments, layer order, RNG policy, or spill consumer. After `my_future_lapply()` returns, materialize its named result as the in-memory `layers` list; when `spill.dir` is set, write every returned layer to its canonical RDS path, record the count, and release the temporary in-memory list before the existing NULL-output spill path. Apply the existing projected-output guard to the completed parallel result when no spill directory is provided.

## Constraints and invariants

- `object@net$cell$prob` is a named `SparseChatArray` with one bitwise-identical `dgCMatrix` per LR pair when `spill.dir` is NULL.
- Parallel and serial results preserve layer names/order, values, sparsity pattern, object validity, and `misc$.param$communication` fields.
- Parallel spill writes exactly one `<layer-name>.rds` file per LR layer; `filterProbability()` consumes those files, materializes the filtered array, clears files and spill metadata, and returns the same result as the serial spill path.
- `future.seed`, `future.scheduling`, progress policy, warnings/errors, and the public signature remain unchanged.
- Existing serial behavior and guard threshold remain unchanged; no new backend or compatibility alias is added.

## Alternatives considered

The alternative of assigning `layers <- Prob.cell_` only would repair in-memory output but leave parallel spill mode broken because no RDS files would be written. Reusing the serial loop inside future workers would write files concurrently and change ownership/error behavior; writing completed future results once in the caller preserves deterministic paths and avoids concurrent file races.

## Consumer impact

`R/modeling.R` produces `Prob.cell_` and consumes `layers`; `filterProbability()` reads `misc$.param$communication$layer.names` and the per-layer RDS files. `tests_dev/test-computeCommunProb.R` covers the repaired `parallel=TRUE` output branch, parallel spill metadata/files/contents, spill consumption, and exact serial equality. Its sourceCpp-backed fixture uses a sequential future plan because a temporary `sourceCpp` DLL pointer cannot be serialized to Windows workers; `tests_dev/test-parallel.R` independently verifies the shared multisession wrapper contract. `man/computeCommunProb.Rd`, `docs/CLI_DESIGN.md`, and this successor Note remove the stale limitation. The prior real CLI smoke remains unchanged and temporary fixtures are not retained.

## Consequences

Parallel calls now return usable objects and parallel spill mode adds a caller-side serialization pass after computation. Spill mode still computes the full future result list before writing, so this repair restores correctness but does not claim streaming during parallel execution. The existing serial path has no semantic change.

## Evidence

Source evidence: `R/modeling.R:325-420` initializes `layers`, previously wrote only `Prob.cell_` in the parallel branch, and wrote RDS files only in the serial branch; `R/modeling.R:2127-2134` requires those files in `filterProbability()`.
- R 4.4.3: test-computeCommunProb.R passed 54 checks; test-computeAvgCommunProb.R passed 162; test-filterCommunication.R passed 22; test-netpathway-aggregate.R passed 44; test-centrality-subset.R passed 39; test-computeCellDistance.R passed 8; test-parallel.R reported 147 passed and 0 failed. All scripts exited 0.
- The updated source roxygen, `man/computeCommunProb.Rd`, and `docs/CLI_DESIGN.md` state that parallel spill writes after worker completion rather than streaming during worker execution. `tools::checkRd()` passed for the affected help pages.

## Acceptance criteria

1. Satisfied: the `parallel=TRUE` branch returns a valid `SparseChatArray` with bitwise-identical layers under the future execution path, and the repaired serial/spill metadata contract is covered by `test-computeCommunProb.R`.
2. Satisfied: parallel spill writes exactly one RDS per LR, records `spill.count`, `filterProbability()` consumes and removes them, and the filtered result equals the serial result.
3. Satisfied: all affected numerical and shared-wrapper suites exit successfully without CLI regressions.
4. Satisfied: source roxygen/man help, CLI documentation, and the successor Note describe the repaired boundary; no temporary test fixture remains. The Note validator passes after lifecycle transition.
