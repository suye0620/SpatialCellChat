# Agent Note: 兼容历史 SpatialChat 类的 15-slot 迁移

Status: implemented
Governance: v1
Date: 2026-11-16
Decision type: migration
Scope: R/merging.R updateSpatialCellChat, historical fixtures, and inference validation
Owner: project maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The two real historical fixtures supplied for migration are serialized with the
S4 class marker `SpatialChat`, not `SpatialCellChat`. Their stored fields are the
same 15-slot upstream layout, but the defining `SpatialChat` package is not
installed in the current renv. The previous migration entry therefore rejected
the class before it could inspect the persisted slots, and `methods::slot()` also
attempted to resolve the unavailable historical class definition.

## Decision

Allow `updateSpatialCellChat()` to migrate a class-marked `SpatialChat` object
only when its persisted attributes exactly match the recognized 15-slot legacy
layout. Read legacy values directly from the object's stored attributes rather
than resolving them through the unavailable class definition. Keep the current
11-slot path restricted to `SpatialCellChat`, and reject all unknown class/layout
combinations.

## Constraints and invariants

- The migration boundary is the exact legacy slot set: `data.raw`, `data`,
  `data.signaling`, `data.scale`, `data.project`, `images`, `net`, `netP`,
  `meta`, `idents`, `DB`, `LR`, `var.features`, `dr`, and `options`.
- Historical `SpatialChat` objects are not treated as current objects and never
  bypass conversion or validation.
- Stored expression, spatial, network, LR, metadata, and identity axes retain
  their values and alignment checks.
- No `SpatialChat` package installation is required to read a serialized object;
  the migration still constructs and validates a development 11-slot object.
- Unknown slots, missing essential axes, and mismatched probability/p-value
  axes continue to fail explicitly.

## Alternatives considered

- Require users to reinstall the historical `SpatialChat` package: rejected
  because the available RDS files are valid serialized objects and the package
  is not part of the current renv; this would make migration depend on an
  unavailable historical runtime.
- Rename the class marker before migration: rejected because it hides the
  historical provenance and can make S4 slot resolution silently use the wrong
  definition.
- Accept any object with a matching slot-name set: rejected because an unknown
  producer could accidentally satisfy the layout while violating the intended
  compatibility boundary.

## Consumer impact

The direct consumer is `updateSpatialCellChat()`. The historical producers are
`chat_3cts.rds`, `chat_4cts.rds`, and the workflow in
`tests_dev/test_data/old_version_chat_objs/simulation data/simulation data/simulationChat.R`.
The regression consumer is `tests_dev/test-updateSpatialCellChat.R`; the inference
consumer is `tests_dev/compare-old-spatialchat-inference.R`. The public migration
help in `man/updateSpatialCellChat.Rd` documents the accepted historical class.
Downstream consumers are the 11-slot validators, inference functions, centrality,
and spatial visualization functions.

## Consequences

The migration can consume historical RDS files without installing their original
package, while retaining a strict schema boundary. Reading attributes directly
avoids S4 class-cache failures, but it means the migration must continue to keep
its explicit slot-set detector synchronized with supported legacy producers.
The old and development inference paths remain separately executable because the
packages define incompatible same-name S4 classes.

## Evidence

- Both supplied fixtures have class marker `structure("SpatialChat", package =
  "SpatialChat")` and the exact 15 persisted legacy fields.
- The installed upstream package defines the same 15-slot field layout under
  `SpatialCellChat`; the historical package itself is absent.
- `tests_dev/test-updateSpatialCellChat.R dev-first` and `upstream-first` both
  pass, including migration and validation of both real fixtures.
- Independent upstream and development workflows using the supplied fixtures
  completed through the requested centrality stage for both 3cts and 4cts.
- Inference comparison artifacts are in
  `tests_dev/results/legacy-inference-2026-11-16/`.

## Acceptance criteria

- Both real `SpatialChat` RDS fixtures migrate without installing `SpatialChat`.
- The migrated objects pass `methods::validObject()` and the development
  validator, with 4cts sparse communication layers and stored centrality retained.
- The full 3cts and 4cts workflows complete in isolated upstream and development
  R processes through the centrality stage.
- Cell-level probabilities, filtering results, and cell-level centrality match
  elementwise; group probabilities differ only at floating-point roundoff.
- Any group-level p-value discrepancy is reported with exact coordinates and is
  not silently called identical.


Additional fixture evidence: `chat_4cts.rds` contains a precomputed `L33-R31`
layer with 491 nonzero entries and centrality up to 500, while rerunning the
current upstream workflow from the same raw matrix produces 70 nonzero `L33-R31`
entries and centrality up to 136. The isolated upstream and development reruns
match each other; the preserved historical result is therefore treated as a
fixture snapshot, not as the fresh-workflow reference.
