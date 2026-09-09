# Agent Note: Communication-flow widget asset resolution

Status: implemented
Governance: v1
Date: 2026-09-03
Decision type: dependency
Scope: `netVisual_CommunFlowWidget()`, `inst/htmlwidgets`, widget rendering tests
Owner: SpatialCellChat maintainer
Impact: medium
Supersedes: none
Superseded by: none

## Problem

The communication-flow binding is registered through the installed package's `inst/htmlwidgets/spatialcellchat-commun-flow.yaml`. When the R sources are loaded directly during package development and the installed `SpatialCellChat` copy lacks the new JavaScript or CSS, `htmlwidgets` cannot discover the binding through the package path. The widget object is valid but an RStudio Viewer page can omit its local renderer assets.

## Decision

Resolve a complete installed `SpatialCellChat/htmlwidgets` directory first. If it is unavailable or lacks either required asset, resolve `inst/htmlwidgets` under the current package root and attach one explicit `htmltools::htmlDependency` to a widget created under the `htmlwidgets` package. Reject calls when neither complete asset directory exists.

## Constraints and invariants

- Installed packages with both assets use the existing YAML registration path; they do not receive a duplicate dependency.
- Direct-source fallback requires both the JavaScript and CSS from the same local `inst/htmlwidgets` directory.
- Every rendered development widget has exactly one dependency named `spatialcellchat-commun-flow` with version `0.1.0` and both asset files.
- The dependency change does not alter the R payload, communication-field semantics, or public widget parameters.
- A source-loaded widget must render through normal `htmlwidgets::saveWidget()` behavior; validation must not manually inject JavaScript or CSS into the HTML.

## Alternatives considered

1. Rely only on package YAML registration. Rejected: a stale installed copy has no binding assets while source development is active.
2. Always attach a local dependency. Rejected: this duplicates installed package dependencies and requires a package-root working directory in normal installed use.
3. Embed script and stylesheet strings into every generated page. Rejected: it bypasses `htmlwidgets` dependency handling, duplicates assets, and prevents normal cacheable library output.

## Consumer impact

- `R/visualization_widget.R` resolves the asset directory and constructs the widget dependency.
- `inst/htmlwidgets/spatialcellchat-commun-flow.yaml` remains the installed-package registration authority.
- `inst/htmlwidgets/spatialcellchat-commun-flow.js` and `.css` are copied by normal `htmlwidgets` dependency rendering in source development.
- `tests_dev/test-commun-flow-widget.R` verifies the dependency object and saved source-development HTML.
- The real `SpatialChat_2` rendering script verifies a normal saved page without manual asset injection.
- No public API, caller, generated package artifact, or data schema changes.

## Consequences

Developers must run source-loaded widgets from the package root unless a complete installed package copy is available. Installed use remains unchanged. `htmltools` becomes a direct runtime import because the fallback creates an `htmlDependency` explicitly.

## Evidence

- `R/visualization_widget.R:57-91` implements installed-first asset resolution and the source-development `htmlDependency` fallback.
- `inst/htmlwidgets/spatialcellchat-commun-flow.yaml` registers the same binding version and asset files for installed packages.
- `tests_dev/test-commun-flow-widget.R` passes the deterministic payload checks and verifies that normal `saveWidget()` copies and references both local assets.
- `C:/Users/Administrator/AppData/Local/Temp/render-spatialchat2-real-lr.R` generated the real `SpatialChat_2` page with standard `saveWidget()` only; no manual script/style injection remains.
- Headless Chromium loaded the generated local JS/CSS, rendered the real `dL7-dR5` LR page with two populated canvases, z-index cells `1` below trails `2`, and controls for Visual speed, Particles, Point opacity, and Zoom.
- Browser lifecycle API verification showed a stable particle snapshot while paused and changed positions after resume; screenshot: `C:/Users/Administrator/AppData/Local/Temp/omp-sshots-15774d54fdc1a57f.webp`.
- `node --check inst/htmlwidgets/spatialcellchat-commun-flow.js`, R/Rd parsing, and Agent Notes validation passed; validation reported only pre-existing legacy Note warnings.

## Acceptance criteria

- [x] The asset resolver prefers a complete installed asset directory and otherwise returns a complete source asset directory.
- [x] The source-loaded contract widget has an `htmlwidgets` package attribute and one matching `html_dependency`.
- [x] A normal non-self-contained `htmlwidgets::saveWidget()` page references and copies the dependency JavaScript and CSS without manual injection.
- [x] The package YAML path still names the same binding version and asset files.
- [x] Widget contract tests, R/Rd parsing, JavaScript syntax validation, browser smoke, and Agent Notes validation pass.
