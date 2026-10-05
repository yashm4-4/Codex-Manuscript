#!/usr/bin/env python3
"""Normalize shared GRCh38 regulatory references for COPD Section 3.

The shared source files are treated as immutable.  This script writes compact,
six-column, canonical-chromosome BED files beneath the COPD workspace and an
auditable inventory with source/output checksums.  Coordinates are always
0-based, half-open GRCh38 intervals.

Normalized BED columns:
    chrom, start, end, feature_id, feature_type, source_detail
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import date, datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Callable, Iterable, Iterator, Optional

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION = SCRIPT.parents[1]
PROJECT = SECTION.parents[2]
SHARED = PROJECT / "data"
DEFAULT_OUTPUT = SECTION / "data" / "regulatory_references"
DEFAULT_RESULTS = SECTION / "results"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
BED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "feature_id",
    "feature_type",
    "source_detail",
]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--access-date", default="2026-09-21")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: str) -> str:
    """Make a value safe for one BED field without altering its meaning."""
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def canonical_chrom(value: str) -> str:
    return value if value.startswith("chr") else f"chr{value}"


def text_lines(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from handle


ParsedRow = tuple[str, int, int, str, str, str]
Parser = Callable[[Path], Iterable[Optional[ParsedRow]]]


def parse_screen(path: Path) -> Iterator[Optional[ParsedRow]]:
    for line in text_lines(path):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            yield None
            continue
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            yield None
            continue
        yield (
            canonical_chrom(fields[0]),
            start,
            end,
            fields[3],
            fields[5],
            f"linked_accession={fields[4]}",
        )


def parse_gff_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
    return result


def parse_ensembl(path: Path) -> Iterator[Optional[ParsedRow]]:
    record = 0
    for line in text_lines(path):
        if line.startswith("#"):
            continue
        record += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            yield None
            continue
        try:
            start, end = int(fields[3]) - 1, int(fields[4])
        except ValueError:
            yield None
            continue
        attrs = parse_gff_attributes(fields[8])
        feature_id = attrs.get("ID", f"ensembl_record_{record}")
        details = []
        for key in ("gene_id", "gene_name", "gene_biotype"):
            if attrs.get(key):
                details.append(f"{key}={attrs[key]}")
        yield (
            canonical_chrom(fields[0]),
            start,
            end,
            feature_id,
            fields[2],
            ";".join(details) if details else ".",
        )


def parse_fantom5(path: Path) -> Iterator[Optional[ParsedRow]]:
    for line in text_lines(path):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            yield None
            continue
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            yield None
            continue
        score = fields[4] if len(fields) > 4 else "."
        yield (
            canonical_chrom(fields[0]),
            start,
            end,
            fields[3],
            "enhancer",
            f"score={score}",
        )


def parse_repeatmasker(path: Path) -> Iterator[Optional[ParsedRow]]:
    for record, line in enumerate(text_lines(path), start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 17:
            yield None
            continue
        try:
            start, end = int(fields[6]), int(fields[7])
        except ValueError:
            yield None
            continue
        row_id = fields[16] or str(record)
        yield (
            canonical_chrom(fields[5]),
            start,
            end,
            f"rmsk_{row_id}",
            fields[11],
            f"name={fields[10]};family={fields[12]};strand={fields[9]}",
        )


def parse_blacklist(path: Path) -> Iterator[Optional[ParsedRow]]:
    for record, line in enumerate(text_lines(path), start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            yield None
            continue
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            yield None
            continue
        reason = fields[3] if len(fields) > 3 else "blacklist"
        yield (
            canonical_chrom(fields[0]),
            start,
            end,
            f"blacklist_{record}",
            "blacklist",
            f"reason={reason}",
        )


def write_deterministic_gzip(path: Path, rows: Iterable[ParsedRow]) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as out:
                for row in rows:
                    out.write("\t".join(clean(str(value)) for value in row) + "\n")
    temporary.replace(path)


def project_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT))
    except ValueError:
        return str(path)


def normalize(
    reference_id: str,
    source_name: str,
    release: str,
    source_path: Path,
    source_url: str,
    parser: Parser,
    output_path: Path,
    access_date: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if not source_path.is_file():
        raise FileNotFoundError(source_path)

    raw_records = 0
    noncanonical = 0
    invalid = 0
    retained = 0
    feature_counts: Counter[str] = Counter()

    def accepted_rows() -> Iterator[ParsedRow]:
        nonlocal raw_records, noncanonical, invalid, retained
        for parsed in parser(source_path):
            raw_records += 1
            if parsed is None:
                invalid += 1
                continue
            chrom, start, end, feature_id, feature_type, detail = parsed
            if chrom not in CANONICAL:
                noncanonical += 1
                continue
            if start < 0 or end <= start or not feature_id or not feature_type:
                invalid += 1
                continue
            retained += 1
            feature_counts[feature_type] += 1
            yield chrom, start, end, feature_id, feature_type, detail

    write_deterministic_gzip(output_path, accepted_rows())
    verification_rows = sum(1 for _ in text_lines(output_path))
    if verification_rows != retained:
        raise RuntimeError(
            f"{reference_id}: wrote {verification_rows} rows, expected {retained}"
        )
    if retained == 0:
        raise RuntimeError(f"{reference_id}: normalization retained no intervals")

    manifest = {
        "reference_id": reference_id,
        "source": source_name,
        "release": release,
        "genome_build": "GRCh38",
        "coordinate_system": "0-based half-open BED",
        "chromosome_scope": "chr1-22,X,Y",
        "source_path": project_path(source_path),
        "source_url": source_url,
        "source_accessed_on": access_date,
        "source_sha256": sha256(source_path),
        "normalized_path": project_path(output_path),
        "normalized_sha256": sha256(output_path),
        "normalized_columns": ";".join(BED_COLUMNS),
        "raw_data_records": raw_records,
        "retained_canonical_records": retained,
        "noncanonical_records_excluded": noncanonical,
        "invalid_records_excluded": invalid,
        "prepared_on": date.today().isoformat(),
    }
    counts = [
        {
            "reference_id": reference_id,
            "feature_type": feature_type,
            "n_intervals": count,
        }
        for feature_type, count in sorted(feature_counts.items())
    ]
    return manifest, counts


def main() -> None:
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    specifications = [
        (
            "screen_ccre_v3",
            "ENCODE SCREEN candidate cis-regulatory elements",
            "Registry V3",
            SHARED / "screen_ccres" / "GRCh38-cCREs.bed",
            "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed",
            parse_screen,
            args.output_dir / "screen_ccre_v3.grch38.canonical.bed.gz",
        ),
        (
            "ensembl_regulatory_v116",
            "Ensembl Regulatory Build",
            "release 116",
            SHARED
            / "ensembl_regulation"
            / "Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz",
            "https://ftp.ensembl.org/pub/release-116/regulation/homo_sapiens/GRCh38/annotation/",
            parse_ensembl,
            args.output_dir / "ensembl_regulatory_v116.grch38.canonical.bed.gz",
        ),
        (
            "fantom5_cage_enhancers",
            "FANTOM5 CAGE enhancer atlas",
            "hg38_latest accessed 2026-09-21",
            SHARED / "fantom5" / "F5.hg38.enhancers.bed.gz",
            "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/",
            parse_fantom5,
            args.output_dir / "fantom5_cage_enhancers.grch38.canonical.bed.gz",
        ),
        (
            "ucsc_repeatmasker",
            "UCSC RepeatMasker table",
            "hg38 accessed 2026-09-21",
            SHARED / "repeats" / "rmsk.txt.gz",
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
            parse_repeatmasker,
            args.output_dir / "ucsc_repeatmasker.grch38.canonical.bed.gz",
        ),
        (
            "encode_blacklist_v2",
            "ENCODE blacklist",
            "hg38 v2",
            SHARED / "blacklist" / "hg38-blacklist.v2.bed.gz",
            "https://github.com/Boyle-Lab/Blacklist",
            parse_blacklist,
            args.output_dir / "encode_blacklist_v2.grch38.canonical.bed.gz",
        ),
    ]

    manifests: list[dict[str, object]] = []
    counts: list[dict[str, object]] = []
    for specification in specifications:
        manifest, feature_counts = normalize(*specification, args.access_date)
        manifests.append(manifest)
        counts.extend(feature_counts)
        print(
            f"{manifest['reference_id']}\t"
            f"{manifest['retained_canonical_records']} retained\t"
            f"{manifest['noncanonical_records_excluded']} noncanonical excluded\t"
            f"{manifest['invalid_records_excluded']} invalid excluded"
        )

    inventory = pd.DataFrame(manifests)
    count_frame = pd.DataFrame(counts)
    inventory_path = args.results_dir / "COPD-S3-R001_regulatory_reference_inventory.tsv"
    counts_path = args.results_dir / "COPD-S3-R001_regulatory_reference_feature_counts.tsv"
    inventory.to_csv(inventory_path, sep="\t", index=False)
    count_frame.to_csv(counts_path, sep="\t", index=False)
    run_manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exact_replay_command": shlex.join([sys.executable, str(SCRIPT), *sys.argv[1:]]),
        "working_directory": str(Path.cwd()),
        "script": project_path(SCRIPT),
        "script_sha256": sha256(SCRIPT),
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "genome_build": "GRCh38",
        "coordinate_system": "0-based half-open BED",
        "chromosome_scope": "chr1-22,X,Y",
        "normalized_bed_columns": BED_COLUMNS,
        "references_prepared": len(inventory),
        "inventory": project_path(inventory_path),
        "feature_counts": project_path(counts_path),
    }
    manifest_path = args.results_dir / "COPD-S3-R001_preparation_manifest.json"
    manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n")
    print(f"inventory\t{project_path(inventory_path)}")
    print(f"feature_counts\t{project_path(counts_path)}")
    print(f"run_manifest\t{project_path(manifest_path)}")


if __name__ == "__main__":
    main()
