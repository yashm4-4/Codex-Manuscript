#!/usr/bin/env python3
"""Build the checksum manifest for the static COPD manuscript data bundle.

The live whole-workflow validation table and log are intentionally excluded:
the validator rewrites them after checking this manifest, so including them
would create a checksum cycle.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
COPD = REPO / "diseases" / "COPD"
OUTPUT = COPD / "manuscript" / "supplementary_file_manifest.tsv"

RESULT_DIRS = {
    "Supplementary Data S1": COPD / "01_background" / "results",
    "Supplementary Data S2": COPD / "02_gwas" / "results",
    "Supplementary Data S3-S4": COPD / "03_regulatory_landscape" / "results",
    "Supplementary Data S5": COPD / "04_modeling" / "results",
    "Supplementary Data S6": COPD / "05_computational_validation" / "results",
    "Supplementary Data S7": COPD / "06_experimental_validation" / "results",
    "Figure source data": COPD / "manuscript" / "figure_data",
    "Figures": COPD / "manuscript" / "figures",
}

EXTRA_FILES = {
    "Search strategy": [
        COPD / "01_background" / "data" / "search_log.tsv",
        COPD / "03_regulatory_landscape" / "data" / "COPD-S3-R004_literature_search_log.tsv",
        COPD / "03_regulatory_landscape" / "data" / "encode_peaks" / "manifest.tsv",
    ],
    "Traceability and reproducibility": [
        COPD / "README.md",
        COPD / "results_register.tsv",
        COPD / "activity_log.tsv",
        COPD / "decisions.tsv",
        COPD / "sources.tsv",
        COPD / "manuscript" / "COPD_regulatory_genomics_manuscript.md",
        COPD / "manuscript" / "README.md",
        COPD / "manuscript" / "result_to_manuscript_map.tsv",
        COPD / "manuscript" / "software_environment.tsv",
        COPD / "manuscript" / "software_version_capture.txt",
        COPD / "manuscript" / "scripts" / "build_figures.py",
        COPD / "manuscript" / "scripts" / "build_supplementary_manifest.py",
        COPD / "manuscript" / "scripts" / "renumber_citations.py",
        COPD / "manuscript" / "scripts" / "validate_complete_workflow.py",
    ],
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    entries: list[tuple[str, Path]] = []
    for group, directory in RESULT_DIRS.items():
        entries.extend((group, path) for path in sorted(directory.iterdir()) if path.is_file())
    for group, paths in EXTRA_FILES.items():
        entries.extend((group, path) for path in paths)

    seen: set[Path] = set()
    rows = []
    for group, path in entries:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not path.is_file():
            raise SystemExit(f"Missing supplied file: {path}")
        rows.append(
            {
                "supplement_group": group,
                "path_relative_to_copd_root": path.relative_to(COPD).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} checksummed files to {OUTPUT}")


if __name__ == "__main__":
    main()
