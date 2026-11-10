# Agent Note: Preserve rank-one NMF diagonal dimensions

Status: implemented
Governance: v1
Date: 2026-11-11
Decision type: algorithm
Scope: `R/analysis.R` NMF topic and rank-selection helpers
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`identifyCellTopics()` and `suggestK1()` built the NMF singular-value scaling
matrix with `Matrix::diag(sqrt(d))`. In the current Matrix API, a length-one
numeric argument is interpreted as a matrix order, so a rank-one NMF result
with a non-integer singular value produced a wrongly sized matrix and failed
with a non-conformable matrix multiplication error. The rank-one path was not
covered by the previous migration regression.

## Decision

Construct the scaling matrix with `Matrix::Diagonal(x = sqrt(d))` through the
shared `.sc_nmf_sqrt_diagonal()` helper. Validate that the singular values are
finite and non-negative, then use the same scaled factor construction for
rank-one and higher-rank results in both topic storage and cross-validation.
This preserves the existing NMF factorization and reconstruction formulas;
only the ambiguous diagonal constructor is replaced.

## Constraints and invariants

- A rank-one result must produce a 1 x 1 scaling matrix.
- Rank-k results must preserve the k x k diagonal scaling dimensions.
- The reconstructed matrix and stored cell/signaling factor axes retain the
  existing orientation and names.
- No dense communication tensor or new statistical gate is introduced.
- Invalid singular values fail explicitly rather than producing non-finite
  topic loadings.

## Alternatives considered

Keep `Matrix::diag()` and require callers to avoid rank one: rejected because
rank one is a valid NMF choice and the failure depends on the numeric value.
Use base `diag()` with an explicit `nrow`: rejected because it is less direct
for sparse Matrix results and duplicates the constructor logic at each caller.

## Consumer impact

The producer/consumers are `identifyCellTopics()`, `suggestK1()`, and the
topic plotting paths that read their stored factors. The focused regression is
`tests_dev/test-remaining-migration.R`; the existing P1 topic and P3 topic
visualization tests remain consumers. No public function signature or schema
slot changes.

## Consequences

Rank-one topic analysis and rank-selection cross-validation now work under the
current Matrix version. The helper adds one shared internal validation point.
`suggestK1()` still requires at least two tested ranks and at least two
replicates for a meaningful variance curve; those historical statistical
semantics are unchanged.

## Evidence

- Before the change, `suggestK1()` with `k.test = c(1, 2)` failed in
  `outs_NMF$w %*% Matrix::diag(sqrt(outs_NMF$d))`.
- `tests_dev/test-remaining-migration.R` now executes `suggestK1()` on an
  8 x 8 sparse matrix with ranks 1 and 2 and two replicates for both `p = 0`
  and `p = 0.5`.
- `tests_dev/test-spatial-p1-migration.R` executes
  `identifyCellTopics(topic.k = 1)` and verifies the cell and signaling factor
  axes.
- Both focused tests pass; `suggestK1()` returns the expected ggplot and four
  rank/replicate rows for each branch.

## Acceptance criteria

- Rank-one and rank-two `suggestK1()` calls complete without a dimension error.
- `identifyCellTopics(topic.k = 1)` stores one cell factor and one signaling
  factor with the expected row names.
- The returned `stats` table retains one row per rank/replicate and the
  diagnostic figure remains a ggplot.
- Existing P1/P3 topic tests and the remaining migration test pass.
- The source package Rd and `R CMD check` gates remain valid after the change.
