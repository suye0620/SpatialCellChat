# Legacy inference comparison

- Date: 2026-11-16
- Source fixtures: `chat_3cts.rds`, `chat_4cts.rds`
- Workflow: database replacement, `subsetData`, Moran's I over-expression, LR selection, `computeCommunProb`, `filterProbability`, `filterCommunication`, `computeAvgCommunProb`, final cell-level `filterCommunication`, and the centrality calls present in `simulationChat.R`.
- Upstream and development runs used separate R processes because both packages define a same-name S4 class.

## Result

Both historical `SpatialChat` fixtures migrate successfully without installing the historical package. Cell-level probabilities, filtering results, and the requested cell-level centrality are numerically equivalent. Group probabilities differ only at floating-point roundoff. The only observed output differences are five 4cts group-level p-values across four LR layers, each differing by one of 100 permutation counts.

## 3cts

- LR selection and selected features: identical (6 LR pairs; 11 selected features).
- Communication probabilities after `computeCommunProb`: maximum absolute difference `0`.
- Probabilities after `filterProbability`: maximum absolute difference `0`.
- Probabilities after `filterCommunication`: maximum absolute difference `0`.
- Final cell-level probabilities: maximum absolute difference `0`.
- Group probabilities: maximum absolute difference `5.551115e-17`.
- Group p-values: identical in value.
- Cell-level centrality: identical at tolerance `1e-10`; maximum absolute difference `0`.
- Group-level centrality: identical at tolerance `1e-10`; maximum absolute difference `0`.

## 4cts

- LR selection and selected features: identical (6 LR pairs; 12 selected features).
- Communication probabilities after `computeCommunProb`: maximum absolute difference `0`.
- Probabilities after `filterProbability`: maximum absolute difference `0`.
- Probabilities after `filterCommunication`: maximum absolute difference `0`.
- Final cell-level probabilities: maximum absolute difference `0`.
- Group probabilities: maximum absolute difference `1.110223025e-16`.
- Cell-level centrality: identical at tolerance `1e-10`; maximum absolute difference `0`.
- Group-level centrality: not requested by the 4cts source workflow and absent in both runs.

The group p-value differences are:

| LR | source | target | upstream | development |
|---|---|---:|---:|---:|
| `L03-R01` | C | B | 0.53 | 0.52 |
| `L13-R11` | D | A | 0.73 | 0.72 |
| `L43-R41` | D | B | 0.93 | 0.94 |
| `L43-R41` | C | D | 0.64 | 0.63 |
| `L53-R51` | A | A | 0.72 | 0.71 |

The five p-values differ by one permutation hit out of 100. Both values in each comparison are above 0.05.

### First divergence and controlled replay (2026-11-22)

The first divergence is the random number generator (RNG) state before `replicate(nboot, sample.int(nC, size = nC))` in `computeAvgCommunProb()`. Both engines call `set.seed(1L)`. Upstream then computes the observed group averages through `my_future_sapply()`, whose default `future.seed = TRUE` advances the caller's RNG state. Dev computes the observed averages with its Rcpp kernel and reaches the permutation draw with the state from `set.seed(1L)`. Calling the upstream helper with an RNG-free callback reproduced the recorded upstream state. Passing `future.seed = FALSE` in that control reproduced the Dev state.

Both workflows used `min.percent = 0`, `min.cells.sr = 5`, and `nboot = 100`. The first permutation differs at 473 of 500 cell indices. The other 99 permutations match exactly. At the five reported network positions, the first permutation moved `Pboot - Pnull` across zero by margins of at least `4.02e-4`. Four more sign changes had zero observed probability and received p = 1 in both runs.

A third full Dev run restored the recorded upstream RNG state immediately before drawing the permutations. It reproduced all 100 upstream permutations, every strict `Pboot - Pnull > 0` decision on the tested layers, and all final group p-values. The diagnostic instruments a private copy of the installed function in each R process. Production inference keeps its existing RNG behavior.

The `permutation-trace-2026-11-22/` directory contains the permutation indices, pre-draw RNG states, `Pnull` and `Pboot`, and final group p-values from the three runs. The check script also verifies that the five original p-value differences leave the 0.05 classification unchanged. Run the following from the repository root with the two packages available in `renv`.

```sh
Rscript tests_dev/check-legacy-permutation-trace.R tests_dev/results/legacy-inference-2026-11-16/permutation-trace-2026-11-22
```

To regenerate the traces using the supplied `chat_4cts.rds` fixture, run each command in a separate R process. Replace `trace-dir` with the output directory.

```sh
Rscript tests_dev/compare-old-spatialchat-inference.R upstream 4cts trace-dir --trace-permutations
Rscript tests_dev/compare-old-spatialchat-inference.R dev 4cts trace-dir --trace-permutations
Rscript tests_dev/compare-old-spatialchat-inference.R dev 4cts trace-dir/matched-seed --trace-permutations trace-dir/upstream-4cts-permutation-trace.rds
Rscript tests_dev/check-legacy-permutation-trace.R trace-dir
```

## Stored network in `chat_4cts.rds`

The November 22 follow-up verified that migration preserves the stored cell-level centrality array exactly. A Dev recomputation on the migrated historical network matches the stored array after alignment by measure, cell and LR names. The maximum absolute difference over 5 measures, 500 cells and 6 LR layers was 0. This check is included in `tests_dev/test-updateSpatialCellChat.R`.

Starting from the expression matrix produces a different network from the historical one in `L33-R31`. The stored layer has 491 nonzero entries and the fresh upstream layer has 70, with a maximum absolute value difference of `0.01395126`. The other five LR layers match the fresh run. The two fresh workflows agree on all six cell-level layers and their centrality scores.
