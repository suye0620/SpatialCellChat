# Agent Note: 11-slot inference chain CLI contract

Status: superseded
Governance: v1
Date: 2026-10-19
Decision type: api
Scope: Recent inference chain, shared CLI/progress helpers and their consumers
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: .agents/notes/implemented/2026-10-19-compute-commun-prob-parallel-output.md

## Problem

The migrated inference functions mix `.cli`, `cat(cli.symbol())`, `message`, progressr and pbapply. `computeCommunProb(verbose=FALSE)` still emits entry/input/completion messages; `computeAvgCommunProb(verbose=FALSE)` still emits its title. Shared future wrappers ignore global quiet and suppress sequential progress. `computeCommunProb` also ties its serial memory guard to the verbose branch. These conflict with `docs/CLI_DESIGN.md` and make logging controls affect more than presentation.

## Decision

Keep public signatures and logical `verbose` parameters. Ordinary messages use `.cli`; progress uses progressr in both sequential and future execution. Global `SpatialCellChat.verbose=0` and an explicit local `verbose=FALSE` each suppress ordinary messages and progress; `TRUE` does not override global quiet. Warnings/errors retain their condition semantics. No new integer-valued per-function verbosity API. Batch progress follows progressr's user policy rather than being forced on. Each loop reports completed items, not sticky per-item text; synchronous C++ calls have only LR-level completion granularity.

## Constraints and invariants

- Preserve probabilities, p-values, thresholds, RNG draws/order, schemas, return types, iteration order, scheduling, and spill ownership.
- Logging must not disable safety checks: consolidate duplicated serial loops so the existing output memory guard applies in quiet mode too.
- No new parallel execution for filterCommunication; replace its pbapply presentation with a serial progressr loop.
- Shared `.sc_progress_enabled(verbose=TRUE)` checks local logical permission plus global verbosity. It does not set options or handlers.
- Future wrappers retain arguments, result naming/simplification and `future.seed`; progress must not affect values or RNG state.
- Only setEnvironment configures global handlers; computation scopes do not overwrite user handlers, progressr.enable, or future plans.
- Do not delete exported cli.symbol or pbapply dependency: unrelated consumers still exist.
- subsetCommunication remains a getter without an entry/completion banner or progress bar; real warning conditions stay warnings.
- Distance builders `computeCellDistance` and `createCellCellContactMatrixFrom_dspatial` leave ordinary diagnostics to their computation caller, so nested local quiet cannot leak low-level banners. Their validation and matrix outputs are unchanged.

## Alternatives considered

- Introduce `verbose=NULL` with integer levels everywhere: rejected for this change because it expands the public API migration beyond output consistency.
- Keep pbapply or switch serial loops to future just for progress: rejected because it creates competing presentation controls or changes computation scheduling/RNG.
- Rewrite every legacy output in R/: rejected; this scope is the nine recently migrated functions and necessary shared helpers only.

## Consumer impact

- R/modeling.R: computeCommunProb, filterProbability, computeAvgCommunProb, filterCommunication, computeCommunProbPathway, aggregateNet, relabelSpatialCellChat (standard 11-slot path; legacy Visium implementation not migrated here).
- R/analysis.R: netAnalysis_computeCentrality and subsetCommunication; unchanged downstream visualization and netMappingDEG consumers keep the same returned data.
- R/spatial.R: the two distance builders stop emitting ordinary diagnostics; computeCommunProb reports working-unit ranges and cache information instead. Direct callers in visualization and distance tests still receive identical matrices/errors, without duplicate banners.
- R/utilities.R: .cli, setEnvironment, my_future_lapply/sapply. Shared wrappers also serve expression helpers, legacy Visium, spatial distance consumers, clustering/NMF and identifyOverExpressedInteractions; their argument/result contracts remain unchanged, but quiet/sequential progress behavior follows the unified contract.
- Existing numerical suites: test-computeCommunProb.R, test-computeAvgCommunProb.R, test-filterCommunication.R, test-netpathway-aggregate.R, test-centrality-subset.R. The two core suites retain local-quiet regressions. test-parallel.R covers actual condition/progress behavior instead of vacuous stdout/length assertions. test-computeCellDistance.R protects the silent helper's matrix contract. A throwaway real pipeline smoke checks output and unchanged results.
- docs/CLI_DESIGN.md and affected roxygen/man help: distinguish shipped chain behavior from remaining legacy migration. NAMESPACE, DESCRIPTION, renv.lock and benchmark data stay unchanged.

## Consequences

Scripts can silence migrated computation output consistently without losing errors. Existing stdout log collectors must capture message/stderr output instead. Sequential long-running loops become progress-capable; noninteractive scripts still opt in via progressr. No inference backend changes are required.

## Evidence

- R/utilities.R .cli, setEnvironment, my_future_lapply/sapply; R/modeling.R migrated function bodies; R/analysis.R centrality/subset paths.
- docs/CLI_DESIGN.md now distinguishes this shipped chain from remaining legacy adoption. Seven affected man pages parse cleanly with tools::checkRd and their usage parameter names/order match live function formals; roxygen2 was not installed, so help was synchronized directly from source roxygen.
- Before this change, tests_dev/test-parallel.R captured stdout for cli messages and contained a tautological length >= 0 assertion; the replacement observes message conditions and completed progress events.
- LSP reference request returned no R language server; consumers were mapped with targeted repository search.
- Before-edit real 20-cell/3-LR smoke: local false compute emitted one stdout line and 11 message conditions; local false averaging emitted three message conditions. Global quiet compute/average/relabel each still emitted three progress updates; filters/pathways/centrality and subset emitted ordinary output. Exact net/netP/idents results were saved for after-edit comparison.
- The pre-existing optional computeCommunProb parallel branch populates `Prob.cell_` but constructs output from `layers`; its result/spill wiring is outside this presentation change. Parallel progress is verified through actual filterProbability/centrality consumers and shared wrappers, not claimed as a repair of that branch.
- R 4.4.3 verification: test-parallel.R 147 passed / 0 failed; test-computeCommunProb.R 43, test-computeAvgCommunProb.R 162, test-filterCommunication.R 22, test-netpathway-aggregate.R 44, test-centrality-subset.R 39, and test-computeCellDistance.R 8 checks passed. Every script exited 0.
- After-edit 20-cell/3-LR smoke: all nine functions under global quiet, both core functions under local false, and two-worker filterProbability/centrality returned zero stdout/message/progression updates. Exact pipeline net/netP/idents/subset results matched the saved pre-edit fixture (excluding time/log metadata); actual warning/error conditions remained observable. Normal/detailed CLI output, forced-contact counts, sequential/multisession progress, and spill consumption were exercised. Parallel filter/centrality results and serial spill results were identical to their in-memory sequential counterparts.
- A separate six-LR smoke temporarily lowered GUARD.BYTES to one byte with trace(), without modifying source: both verbose=FALSE/global=0 and verbose=TRUE/global=1 raised the existing projected-output memory error. This checks control flow, not a real 16 GB allocation. trace() was removed afterward.
- Validation emitted pre-existing renv out-of-sync, sourceCpp missing inst/include, Matrix/igraph deprecation notices and intentional future error-cancellation warnings. These were not silenced or changed by this CLI work.
- Before transition, validate_agent_notes.py --root .agents/notes checked 26 Notes with 0 errors and 13 legacy-governance warnings. Both throwaway smoke scripts and the pre-edit RDS fixture were removed after successful verification.

## Acceptance criteria

1. Satisfied by the real 20-cell/3-LR smoke: all nine functions honor global quiet, both core functions honor local false, and normal/detailed output renders stage summaries.
2. Satisfied by test-parallel.R (147 checks) plus real two-worker filter/centrality calls: sequential/multisession completion, explicit/global disable, progressr batch opt-out, unchanged seeded values/RNG, and observable conditions.
3. Satisfied by all five numerical suites and shared-wrapper suite exiting 0, with exact pre/post pipeline net/netP/idents comparison; time/log metadata excluded. Distance-helper regression also passes.
4. Satisfied by the updated CLI design and seven checked Rd pages, removal of all temporary smoke artifacts, and the pre-transition Note validator result above.
