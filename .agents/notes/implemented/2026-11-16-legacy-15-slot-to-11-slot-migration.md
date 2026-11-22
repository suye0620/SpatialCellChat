# Agent Note: 旧版 SpatialCellChat 迁移到 11-slot 结构

Status: implemented
Governance: v1
Date: 2026-11-16
Decision type: migration
Scope: R/merging.R updateSpatialCellChat and legacy compatibility tests
Owner: project maintainer
Impact: high
Supersedes: none
Superseded by: none

## Problem

The development package exposes an 11-slot `SpatialCellChat` schema, while the
pinned upstream `SpatialCellChat` package persists 15-slot objects with
`data.raw`, `data`, `data.signaling`, `data.scale`, `data.project`, `options`,
and the legacy communication layout. Before this change,
`updateSpatialCellChat()` only calls the final-schema validator, so an upstream
object can be accepted because both packages use the same S4 class name and
then be returned without conversion.  This is an observable compatibility
failure: `methods::is()` is ambiguous when both packages are loaded. Moreover,
`methods::slotNames(object)` can resolve a conflicting loaded class definition
instead of the object's stored fields.

## Decision

Use `updateSpatialCellChat()` as an explicit one-way migration from the
upstream 15-slot schema to the development 11-slot schema.  Detect the schema
from the object's actual stored slot attributes, read legacy slots with `methods::slot()`, normalize
expression, images, metadata, options, and communication results through
dedicated helpers, construct a new 11-slot object explicitly, and validate it
before returning.  Already-final 11-slot objects remain a validation-only
path; unknown slot layouts fail explicitly.

## Constraints and invariants

- `assay$norm` is copied from legacy `data`; raw, scale, smooth, and signaling
  are mapped from `data.raw`, `data.scale`, `data.project`, and
  `data.signaling`, respectively, preserving gene and cell dimnames.
- Cell metadata and identities remain aligned to the normalized assay columns;
  legacy metadata are retained without inventing a `meta$samples` column.
- Spatial coordinates are renamed from `x_cent`/`y_cent` to `x`/`y`, reordered
  to normalized-assay cells, and stale distance/grid caches are discarded when
  their identity cannot be proven.
- Legacy cell- and group-level probability/p-value arrays are converted to
  `SparseChatArray` with preserved cell/group/LR dimnames. The upstream sparse
  triplets are split by layer once rather than rescanned for each LR. Existing count,
  weight, pathway, centrality, and auxiliary network records are retained when
  they can be represented without changing their axes. Probability and p-value
  axes must agree; mismatches are rejected instead of mixing LR layers.
- `net` and `netP` always expose legal `cell` and `group` branches, and the
  returned object passes both `methods::validObject()` and
  `validateSpatialCellChat()`.
- The installed package registers the existing `SparseChatArray` S3 methods
  needed for ordinary indexing, dimensions, dimnames and printing outside its
  namespace; a general public-export audit remains separate work.
- Legacy merged objects are migrated only when their normalized `data` matrix
  exists and dataset cell axes can be aligned; otherwise the error directs the
  user to migrate and re-merge the original single-dataset objects.
- The migration does not change the upstream package or its 15-slot class; the
  two packages remain a same-process compatibility edge requiring schema
  detection by stored slots rather than class identity or class-definition lookup.

## Alternatives considered

- Keep rejecting legacy objects and require users to rebuild them: rejected
  because it discards computed communication results and contradicts the
  public `updateSpatialCellChat()` migration boundary.
- Use the upstream `SpatialCellChat::updateSpatialCellChat()` first and then
  adapt its result: rejected because upstream only normalizes old upstream
  fields and still returns a 15-slot object; it cannot create the development
  schema.
- Preserve legacy slots as aliases in the new object: rejected because the
  final schema contract intentionally removes the old top-level slots and
  duplicate representations would allow inconsistent values.

## Consumer impact

Producers are the pinned upstream `SpatialCellChat` constructor, merge function,
and serialized 15-slot objects.  The direct consumer is `updateSpatialCellChat()`; downstream
consumers include all 11-slot validators, accessors, inference, merge/subset,
and visualization functions.  Compatibility coverage belongs in a dedicated
`tests_dev/test-updateSpatialCellChat.R` script and the public migration help
page `man/updateSpatialCellChat.Rd`; `NAMESPACE` registers the existing sparse
array methods for installed-package consumers. No package identity, lockfile, or
upstream installation metadata is changed by this migration.

## Consequences

The migration allocates a new object and may copy large expression matrices;
this is intentional because sharing an S4 object whose class definition is
ambiguous would be unsafe.  Sparse communication layers remain sparse, but
legacy dense/array/list representations are normalized layer by layer.  Original
image records and options remain available in `misc$.legacy`; stale distance
caches are not reused. Unalignable essential axes cause an explicit error. Future schema
changes must update the slot detector and migration tests rather than adding
implicit aliases.

## Evidence

Evidence collected in the workspace:

- Before this change, `R/merging.R` validated only final-schema objects; the
  same source file now contains the migration and its helpers.
- The installed upstream class has slots
  `data.raw,data,data.signaling,data.scale,data.project,images,net,netP,meta,idents,DB,LR,var.features,dr,options`.
- The development class has slots
  `assay,images,meta,idents,features,LR,dr,net,netP,DB,misc`.
- Loading both packages exposes a same-name S4 cache ambiguity. An object with
  15 stored legacy slots and a development-package class marker resolves
  `slotNames()` to the 11-slot definition, while its slot attributes remain 15.
- `R/SpatialCellChat_class.R` defines the final assay, image, network, and
  object validation contracts.
- `tests_dev/test-updateSpatialCellChat.R dev-first` and the same command with
  `upstream-first` both passed with the installed development package. These
  create upstream RNA, spatial and merged objects and check stored values,
  dimnames, current-object handling, error boundaries, plotting, downstream
  `subsetCommunication()`, and a 40-layer / 10,000-triplet sparse conversion.
- `tests_dev/test-SpatialCellChat-class.R` and
  `tests_dev/test-merge-spatialcellchat.R` passed; the latter reported 37
  checks including the current-schema update path.
- Both S4 classes validate in the independent dual-package tests, and an
  upstream-only R process independently validated its 15-slot object.
- The Notes validator reports this record as valid; existing unrelated Notes
  account for its nonzero exit status.

## Acceptance criteria

- A 15-slot upstream RNA object migrates to an 11-slot object with aligned
  assay, metadata, idents, features, DB, LR, dr, and misc fields.
- A 15-slot upstream spatial object migrates with finite `x`/`y` or `x`/`y`/`z`
  coordinates, preserved ratio/tol, and no invalid stale distance cache.
- Legacy `prob.cell`, `pval.cell`, group `prob`/`pval`, `sum`/`weight`, and
  pathway layers retain values and dimnames after conversion where present.
- Character and list forms of `var.features` both migrate successfully.
- A final 11-slot object remains valid and an unknown slot layout produces an
  explicit schema error.
- A merged object with complete `data` and alignable dataset axes retains its
  dataset snapshots; one whose upstream merge omitted `data` fails explicitly.
- In an independent R process that loads both packages, every migration case
  passes `methods::validObject()` and `validateSpatialCellChat()`; the upstream
  package can still construct and validate its own 15-slot object.
- The Agent Notes validator reports no new v1 errors.
