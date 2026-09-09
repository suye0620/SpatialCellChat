#!/usr/bin/env python3
"""Validate Agent Note metadata, lifecycle paths, and governed Note contracts."""

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
REQUIRED_V1_METADATA = (
    "status", "governance", "date", "decision type", "scope", "owner", "impact",
    "supersedes", "superseded by",
)
COMMON_SECTIONS = {
    "problem", "constraints and invariants", "alternatives considered", "consumer impact",
    "consequences", "evidence", "acceptance criteria",
}
MAIN_SECTIONS = {
    "proposed": "proposal",
    "implemented": "decision",
    "rejected": "rejection",
    "superseded": "decision",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".agents/notes", type=Path)
    parser.add_argument("--strict", action="store_true",
                        help="fail Notes that have not adopted Governance: v1")
    return parser.parse_args()


def metadata_and_sections(text: str) -> tuple[dict[str, str], set[str]]:
    metadata: dict[str, str] = {}
    for key, value in re.findall(r"(?m)^([A-Za-z][A-Za-z ]+):\s*(.*?)\s*$", text):
        metadata[key.lower()] = value.strip()
    sections = {section.strip().lower() for section in re.findall(r"(?m)^##\s+(.+?)\s*$", text)}
    return metadata, sections


def validate_reference(value: str, repo_root: Path, label: str, problems: list[str]) -> None:
    if value.lower() == "none":
        return
    candidate = repo_root / value
    if not candidate.is_file():
        problems.append(f"{label} does not resolve: {value}")


def validate_note(path: Path, root: Path, strict: bool) -> tuple[list[str], list[str]]:
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    warnings: list[str] = []
    if not re.search(r"(?m)^# Agent Note: .+\S\s*$", text):
        errors.append("missing '# Agent Note: <title>'")

    metadata, sections = metadata_and_sections(text)
    status = metadata.get("status")
    if status not in STATUS_DIRECTORIES:
        errors.append("invalid or missing Status")
        return errors, warnings

    if path.parent.name != STATUS_DIRECTORIES[status]:
        errors.append(f"Status '{status}' must live in {STATUS_DIRECTORIES[status]}/")
    if not re.match(r"^\d{4}-\d{2}-\d{2}-[a-z0-9]+(?:-[a-z0-9]+)*\.md$", path.name):
        errors.append("filename must be YYYY-MM-DD-lowercase-kebab-topic.md")

    try:
        dt.date.fromisoformat(metadata.get("date", ""))
    except ValueError:
        errors.append("missing or invalid Date (expected YYYY-MM-DD)")
    if not metadata.get("decision type"):
        errors.append("missing Decision type")

    governed = metadata.get("governance") == "v1"
    if not governed:
        message = "legacy Note without Governance: v1"
        (errors if strict else warnings).append(message)
        return errors, warnings

    for key in REQUIRED_V1_METADATA:
        if not metadata.get(key):
            errors.append(f"missing {key.title()} metadata")
    if metadata.get("owner", "").lower() == "unassigned":
        errors.append("Owner cannot be unassigned for Governance: v1")
    if metadata.get("impact") not in {"low", "medium", "high", "critical"}:
        errors.append("Impact must be low, medium, high, or critical")
    if "<!--" in text:
        errors.append("contains unfilled template comment")

    missing_sections = COMMON_SECTIONS - sections
    if missing_sections:
        errors.append("missing sections: " + ", ".join(sorted(missing_sections)))
    main_section = MAIN_SECTIONS[status]
    if main_section not in sections:
        errors.append(f"missing status-specific section: {main_section.title()}")

    repo_root = root.parent.parent
    validate_reference(metadata.get("supersedes", "none"), repo_root, "Supersedes", errors)
    validate_reference(metadata.get("superseded by", "none"), repo_root, "Superseded by", errors)
    if status == "superseded" and metadata.get("superseded by", "none").lower() == "none":
        errors.append("superseded Note requires Superseded by")
    return errors, warnings


def main() -> int:
    args = parse_args()
    if not args.root.is_dir():
        print(f"error: Notes root does not exist: {args.root}", file=sys.stderr)
        return 2

    notes = sorted(
        path for directory in {"proposed", "implemented", "rejected", "archived"}
        for path in (args.root / directory).glob("*.md") if (args.root / directory).is_dir()
    )
    error_count = warning_count = 0
    for path in notes:
        errors, warnings = validate_note(path, args.root, args.strict)
        relative = path.as_posix()
        for message in errors:
            print(f"ERROR {relative}: {message}")
        for message in warnings:
            print(f"WARN  {relative}: {message}")
        if not errors and not warnings:
            print(f"OK    {relative}")
        error_count += len(errors)
        warning_count += len(warnings)

    print(f"Checked {len(notes)} Note(s): {error_count} error(s), {warning_count} warning(s)")
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
