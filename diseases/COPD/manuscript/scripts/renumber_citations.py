#!/usr/bin/env python3
"""Renumber numeric citations and reference entries by first manuscript appearance."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "COPD_regulatory_genomics_manuscript.md"
CITATION = re.compile(r"\[((?:\d+(?:-\d+)?)(?:,\d+(?:-\d+)?)*)\]")
REFERENCE = re.compile(r"^(\d+)\.\s+(.+)$")


def expand(group: str) -> list[int]:
    values: list[int] = []
    for item in group.split(","):
        if "-" in item:
            start, end = (int(value) for value in item.split("-", 1))
            values.extend(range(start, end + 1))
        else:
            values.append(int(item))
    return values


def compress(values: list[int]) -> str:
    ordered = sorted(set(values))
    parts: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:] + [None]:
        if value is not None and value == previous + 1:
            previous = value
            continue
        parts.append(str(start) if start == previous else f"{start}-{previous}")
        if value is not None:
            start = previous = value
    return ",".join(parts)


def main() -> None:
    text = MANUSCRIPT.read_text()
    marker = "\n## References\n"
    if text.count(marker) != 1:
        raise SystemExit("Expected exactly one References heading")
    body, reference_block = text.split(marker)

    entries: dict[int, str] = {}
    suffix_lines: list[str] = []
    in_references = True
    for line in reference_block.lstrip("\n").splitlines():
        match = REFERENCE.match(line) if in_references else None
        if match:
            number = int(match.group(1))
            if number in entries:
                raise SystemExit(f"Duplicate reference number: {number}")
            entries[number] = match.group(2)
        else:
            in_references = False
            suffix_lines.append(line)

    first_use: list[int] = []
    seen: set[int] = set()
    for match in CITATION.finditer(body):
        for old_number in expand(match.group(1)):
            if old_number not in entries:
                raise SystemExit(f"Citation [{old_number}] has no reference entry")
            if old_number not in seen:
                seen.add(old_number)
                first_use.append(old_number)

    uncited = sorted(set(entries) - seen)
    if uncited:
        raise SystemExit(f"Uncited reference entries: {uncited}")
    mapping = {old: new for new, old in enumerate(first_use, start=1)}

    def replace(match: re.Match[str]) -> str:
        return "[" + compress([mapping[value] for value in expand(match.group(1))]) + "]"

    renumbered_body = CITATION.sub(replace, body)
    ordered_references = "\n".join(
        f"{mapping[old]}. {entries[old]}" for old in first_use
    )
    suffix = "\n".join(suffix_lines).lstrip("\n")
    output = renumbered_body + marker + "\n" + ordered_references + "\n\n" + suffix.rstrip() + "\n"
    MANUSCRIPT.write_text(output)
    print(f"Renumbered {len(first_use)} references in {MANUSCRIPT}")


if __name__ == "__main__":
    main()
