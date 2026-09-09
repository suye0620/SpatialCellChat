#!/usr/bin/env python3
"""Create a governed Agent Note without external dependencies."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

STATUS_DIRECTORIES = {
    "proposed": "proposed",
    "implemented": "implemented",
    "rejected": "rejected",
    "superseded": "archived",
}
DECISION_TYPES = (
    "architecture", "api", "schema", "algorithm", "dependency", "migration",
    "security", "performance", "testing", "release", "process",
)
IMPACTS = ("low", "medium", "high", "critical")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".agents/notes", type=Path)
    parser.add_argument("--status", required=True, choices=STATUS_DIRECTORIES)
    parser.add_argument("--type", required=True, choices=DECISION_TYPES)
    parser.add_argument("--slug", required=True,
                        help="lowercase kebab-case filename topic, without date or extension")
    parser.add_argument("--title", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--impact", required=True, choices=IMPACTS)
    parser.add_argument("--date", default=dt.date.today().isoformat())
    parser.add_argument("--supersedes", default="none")
    parser.add_argument("--superseded-by", default="none")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", args.slug):
        print("error: --slug must be lowercase kebab-case", file=sys.stderr)
        return 2
    try:
        dt.date.fromisoformat(args.date)
    except ValueError:
        print("error: --date must be YYYY-MM-DD", file=sys.stderr)
        return 2
    if args.status == "superseded" and args.superseded_by == "none":
        print("error: superseded Notes require --superseded-by", file=sys.stderr)
        return 2

    directory = args.root / STATUS_DIRECTORIES[args.status]
    path = directory / f"{args.date}-{args.slug}.md"
    if path.exists() and not args.force:
        print(f"error: {path} already exists; use --force to overwrite", file=sys.stderr)
        return 1

    main_heading = {
        "proposed": "Proposal",
        "implemented": "Decision",
        "rejected": "Rejection",
        "superseded": "Decision",
    }[args.status]
    content = f"""# Agent Note: {args.title}

Status: {args.status}
Governance: v1
Date: {args.date}
Decision type: {args.type}
Scope: {args.scope}
Owner: {args.owner}
Impact: {args.impact}
Supersedes: {args.supersedes}
Superseded by: {args.superseded_by}

## Problem

<!-- REQUIRED: State the observable problem and cite its source evidence. -->

## {main_heading}

<!-- REQUIRED: State one precise choice, not a list of possibilities. -->

## Constraints and invariants

<!-- REQUIRED: List relationships that must remain true. -->

## Alternatives considered

<!-- REQUIRED: Compare viable alternatives and explain their rejection. -->

## Consumer impact

<!-- REQUIRED: List every producer, consumer, test, doc, script, configuration, generated artifact, and dynamic lookup affected. -->

## Consequences

<!-- REQUIRED: State migration, compatibility, performance, operational, and maintenance effects. -->

## Evidence

<!-- REQUIRED: Cite paths, tests, benchmarks, issues, prototypes, or official references. -->

## Acceptance criteria

<!-- REQUIRED: Write observable pass/fail conditions and the validation command or scenario. -->
"""
    directory.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")
    print(path.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
