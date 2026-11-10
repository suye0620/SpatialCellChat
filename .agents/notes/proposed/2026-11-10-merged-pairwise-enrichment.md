# Agent Note: Merged pairwise enrichment results

Status: proposed
Governance: v1
Date: 2026-11-10
Decision type: migration
Scope: identifyEnrichedInteractions and rankNetPairwise merged snapshots
Owner: maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

`rankNetPairwise()` writes merged results to snapshots, but the consumer reads only the top-level net and fails.

## Proposal

For merged objects return a dataset-named list of independently computed tables, matching `subsetCommunication()`. Single objects retain their existing table result and filtering behavior.

## Constraints and invariants

- Use only the matching snapshot; do not pool evidence, alter p-value rules, or skip missing rankings.
- Preserve positional arguments and single-object return type.
- Reject invalid group selections without independently discarding endpoints;
  otherwise valid-looking results can describe a pair the caller never requested.

## Alternatives considered

Implicitly choosing the first dataset discards samples. Requiring a new dataset selector makes all-dataset consumption unlike existing merged table consumers; callers can index the named list.

## Consumer impact

Producer: `rankNetPairwise()` in `R/analysis.R`. Consumer: `identifyEnrichedInteractions()` in `R/modeling.R`; search found no other consumers. Update its roxygen, Rd, and `tests_dev/test-remaining-entrypoints.R`. No configuration or examples found.

## Consequences

Previously failing merged calls return a list; callers index datasets. Work scales with snapshot count without densifying cell-level data. Existing threshold inconsistencies are outside this migration.

## Evidence

`R/analysis.R` writes `misc$.datasets[[dataset]]$net$pairwiseRank`; `R/modeling.R` previously read `object@net$pairwiseRank`. Searched both names in R, tests_dev and man.

## Acceptance criteria

`tests_dev/test-remaining-entrypoints.R` checks single behavior, dataset-aligned merged results, differing rankings, missing snapshot errors, and invalid endpoint rejection. Document return types and threshold caveat. Run that test, package installation, and the Agent Notes validator.
