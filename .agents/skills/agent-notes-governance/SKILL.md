---
name: agent-notes-governance
description: Govern durable project decisions through evidence-backed Agent Notes. Use when proposing, accepting, rejecting, or superseding architecture, public API, schema, algorithm, dependency, migration, security, performance, testing, release, or process decisions; when a multi-agent project needs durable decision records, ownership, invariants, consumer impact, and validation gates; or when reviewing whether project progress is constrained by its recorded decisions.
---

# Agent Notes Governance

Use Agent Notes to preserve durable decisions that code, tests, and ordinary task trackers cannot explain. Treat a Note as a decision contract—not as a work diary, status update, or substitute for source documentation.

## Non-negotiable rules

- Record **why**, alternatives, consequences, consumers, invariants, and evidence. Keep current behavior and API documentation in source/docs/tests.
- Apply the authority order: current user instruction → executable code and tests → declared schema/API contracts → implemented Notes → proposed/rejected/archived Notes → transient chat history.
- Do not mark a Note `implemented` until the code, callers, tests, and documentation named by its acceptance criteria have actually changed and been verified.
- Do not write a Note for mechanical renames, local formatting, isolated typo fixes, or a one-off experiment without a durable decision.
- Never use a proposed Note as evidence that a design is already shipped.
- Before deleting or replacing an interface, identify every producer, consumer, test, example, generated artifact, configuration entry, and dynamic lookup that uses it.

Read [references/note-contract.md](references/note-contract.md) before creating, reviewing, transitioning, or superseding a governed Note.

## When a Note is mandatory

Create a `proposed` Note before implementation when a change alters any of these:

- public API, data schema, persistence format, migration, compatibility, or deletion boundary;
- algorithm, statistical meaning, security posture, dependency, performance strategy, cache ownership, or lifecycle model;
- cross-package/module contract, shared test strategy, release policy, or work that several Agents will perform independently;
- a high-impact decision that will be costly to reverse.

For medium-impact changes, create the Note before merging the first irreversible implementation step. For low-impact changes, record only the issue/task and evidence unless a decision remains relevant after the task closes.

## Workflow

### 1. Establish the decision boundary

Read the applicable project rules, the relevant source contract, current tests, and existing Notes. Search for an active Note covering the same topic before creating another.

State the decision in one sentence:

```text
Choose <option> for <scope> because <constraint>; reject <alternative> because <cost/risk>.
```

If the statement cannot be written, the decision is not ready. Investigate rather than invent a Note.

### 2. Map evidence and consumers

For every material decision, record:

- **Evidence**: exact source paths, test names, benchmarks, issue links, prototypes, or official references;
- **Consumers**: all codepaths, users, tools, docs, tests, examples, scripts, and dynamic/configuration consumers affected;
- **Invariants**: relationships that must remain true after the change;
- **Owner**: one accountable decision owner, even if several Agents contribute.

A consumer map is mandatory before removing an API, old schema field, compatibility path, cache, provider, or generated artifact.

### 3. Create the correct lifecycle record

Use the generator. It creates a governed v1 Note under the appropriate directory.

```text
uv run --isolated python .agents/skills/agent-notes-governance/scripts/new_agent_note.py \
  --root .agents/notes --status proposed --type architecture \
  --slug result-cache-owner --title "结果缓存所有权" \
  --scope "R/cache and public API" --owner "maintainer" --impact high
```

Fill every required section before treating the Note as reviewable. Do not leave template comments in a proposed Note that is submitted for review.

### 4. Challenge the proposal

For high-impact decisions, use an independent reviewer when available. The review must test:

- whether the problem is real and evidence-backed;
- whether the chosen interface serves all current consumers, not one caller;
- whether alternatives were compared on migration cost, correctness, performance, and reversibility;
- whether invariants and acceptance criteria can fail in a meaningful test;
- whether the decision duplicates an existing authority source.

A reviewer may reject the proposal, request evidence, or approve it. Approval does not change its status to `implemented`.

### 5. Implement and close the contract

Implement only the decision described by the Note. Update every caller, test, API document, generated artifact, and configuration entry named by the consumer map.

Transition to `implemented` only after the acceptance criteria have direct verification evidence. Transition to `rejected` when the proposal is not adopted. When a shipped decision is replaced, create the successor Note first, cross-link both Notes, then transition the old Note to `superseded` under `archived/`.

Do not silently rewrite an old `implemented` Note to describe a new design.

### 6. Validate the Note set

Run the project-local validator after creating or transitioning governed Notes:

```text
uv run --isolated python .agents/skills/agent-notes-governance/scripts/validate_agent_notes.py \
  --root .agents/notes
```

Use `--strict` only after all legacy Notes have been migrated to `Governance: v1`. Default mode validates governed Notes strictly and reports legacy Notes as warnings, allowing incremental adoption.

## Required delivery record

When this Skill is used, report only verifiable facts in this shape:

```text
Decision: <note path and status>
Scope: <affected contracts and consumers>
Evidence: <tests, commands, source references>
Validation: <validator result and changed-contract checks>
Open risk: <only unresolved material risk, or none>
```

## Anti-patterns

- Writing a Note that merely repeats a ticket, commit, or code diff.
- Marking a future plan as `implemented`.
- Keeping incompatible decisions simultaneously authoritative.
- Deleting an old path because it looks unused without consumer evidence.
- Making a validator pass by weakening invariants or omitting failing consumers.
- Recording private reasoning instead of decision-relevant evidence and consequences.
