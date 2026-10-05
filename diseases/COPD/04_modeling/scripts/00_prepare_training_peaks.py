#!/usr/bin/env python3
"""Combine donor-matched ENCODE lung-lobe peaks for TREDNet input building."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
S3 = ROOT / "diseases" / "COPD" / "03_regulatory_landscape"
S4 = ROOT / "diseases" / "COPD" / "04_modeling"
MANIFEST = S3 / "data" / "encode_peaks" / "manifest.tsv"
OUT = S4 / "data" / "training_peaks"
OUT.mkdir(parents=True, exist_ok=True)
CANONICAL = {f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]}


def chrom_key(chrom: str) -> tuple[int, int | str]:
    value = chrom.removeprefix("chr")
    return (0, int(value)) if value.isdigit() else (1, value)


def combine(paths: list[Path], output: Path) -> tuple[int, int]:
    rows: set[tuple[str, ...]] = set()
    total = 0
    for path in paths:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = tuple(line.rstrip("\n").split("\t"))
                total += 1
                if fields[0] in CANONICAL:
                    rows.add(fields)
    ordered = sorted(rows, key=lambda x: (chrom_key(x[0]), int(x[1]), int(x[2]), x))
    with gzip.open(output, "wt") as handle:
        for fields in ordered:
            handle.write("\t".join(fields) + "\n")
    return total, len(ordered)


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    if set(manifest["donor_accessions"]) != {"ENCDO520EJG"}:
        raise RuntimeError("Training peaks are not restricted to ENCDO520EJG")
    definitions = {
        "accessibility": "COPD_ENCDO520EJG_ATAC_all_lobes.narrowPeak.gz",
        "H3K27ac": "COPD_ENCDO520EJG_H3K27ac_all_lobes.narrowPeak.gz",
        "H3K27me3": "COPD_ENCDO520EJG_H3K27me3_all_lobes.narrowPeak.gz",
    }
    rows = []
    for mark, filename in definitions.items():
        selected = manifest.loc[manifest["mark"].eq(mark)]
        paths = [Path(value) for value in selected["path"]]
        output = OUT / filename
        raw_rows, unique_rows = combine(paths, output)
        rows.append(
            {
                "donor": "ENCDO520EJG",
                "disease_context": "severe emphysema",
                "mark": mark,
                "n_lobes": selected["lobe"].nunique(),
                "n_experiments": selected["experiment"].nunique(),
                "raw_peak_rows": raw_rows,
                "unique_canonical_peak_rows": unique_rows,
                "path": str(output.resolve()),
                "source_manifest": str(MANIFEST.resolve()),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(S4 / "data" / "training_peak_manifest.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
