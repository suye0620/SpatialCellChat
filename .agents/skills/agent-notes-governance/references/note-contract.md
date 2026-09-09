# Agent Note Contract v1

Use this contract for new governed Notes. Existing Notes without `Governance: v1` are legacy records: preserve them, validate their basic metadata, and migrate only when they are revised or superseded.

## Purpose and authority

A Note records a durable decision. It never overrides the current user instruction, executable code/tests, or the project’s declared data/API contract.

One fact has one home:

| Fact | Authority |
|---|---|
| Runtime behavior and error conditions | Source code and tests |
| Public parameters, return values, examples | API documentation |
| Canonical schema and invariants | Schema/architecture contract |
| Decision rationale and rejected alternatives | Agent Note |
| Task sequencing and completion | Issue, plan, or task tracker |

## Lifecycle

```text
proposed ──adopt + verify──> implemented ──replaced──> superseded (archived)
    │                              │
    └──────reject──────────────────┴───────────────> rejected
```

- `proposed/`: reviewable, not adopted.
- `implemented/`: current decision; written in present tense; backed by verification.
- `rejected/`: consciously not adopted; preserve the reason and evidence.
- `archived/`: non-current historical records. A superseded Note must link to its successor.

Never move a Note to `implemented/` merely because code work began. Never edit an old decision into its replacement; create and link a successor.

## Path and metadata

Store each Note under:

```text
.agents/notes/<lifecycle>/<YYYY-MM-DD>-<lowercase-kebab-topic>.md
```

Use these exact metadata keys below the title:

```markdown
# Agent Note: <concise decision title>

Status: proposed | implemented | rejected | superseded
Governance: v1
Date: YYYY-MM-DD
Decision type: architecture | api | schema | algorithm | dependency | migration | security | performance | testing | release | process
Scope: <modules, APIs, users, or repositories affected>
Owner: <accountable maintainer or team>
Impact: low | medium | high | critical
Supersedes: none | <repository-relative Note path>
Superseded by: none | <repository-relative Note path>
```

`Status` maps to directories as follows: `proposed → proposed`, `implemented → implemented`, `rejected → rejected`, `superseded → archived`.

## Required sections

Every governed Note must include these sections:

```markdown
## Problem
## Proposal | Decision | Rejection
## Constraints and invariants
## Alternatives considered
## Consumer impact
## Consequences
## Evidence
## Acceptance criteria
```

Use the status-appropriate main heading:

- `proposed`: `## Proposal`
- `implemented`: `## Decision`
- `rejected`: `## Rejection`
- `superseded`: `## Decision` and a non-`none` `Superseded by`

### Section semantics

- **Problem**: observable defect, opportunity, or conflict; cite facts rather than speculation.
- **Proposal / Decision / Rejection**: one precise choice and its scope.
- **Constraints and invariants**: relationships that must hold; not vague goals.
- **Alternatives considered**: at least one viable alternative and why it lost. High/critical decisions require two unless evidence shows only one feasible option.
- **Consumer impact**: every producer, consumer, test, example, document, configuration, generated artifact, and dynamic lookup affected. State `none found` only after a search.
- **Consequences**: migration cost, compatibility boundary, performance/operational cost, and future obligations.
- **Evidence**: source paths, tests, benchmarks, issues, prototypes, or official references.
- **Acceptance criteria**: observable pass/fail conditions. An implemented Note must link them to actual evidence.

## Decision gate

Create a proposed Note before the first irreversible implementation step when a decision changes a shared or high-impact contract. Block the change when any of the following is missing:

- a source-backed problem statement;
- consumer mapping before deletion or incompatible migration;
- explicit invariants for data/lifecycle/security/algorithm changes;
- acceptance criteria that can fail;
- an accountable owner.

## Review gate

A reviewer assesses decision quality, not prose style. For material decisions, challenge scope, alternatives, consumer completeness, invariants, migration plan, and evidence. A proposal with missing evidence stays proposed.

## Implementation gate

A Note may become implemented only after all named contract changes and acceptance checks are complete. If verification is unavailable, state the gap and keep the Note proposed; do not substitute a promise for evidence.

## Legacy migration

Do not mass-rewrite legacy Notes. Preserve their history. Add `Governance: v1` and missing metadata/sections only when a Note is changed, reviewed, or superseded. Run the validator without `--strict` during the transition; enable `--strict` when the repository intentionally adopts v1 for all active Notes.
