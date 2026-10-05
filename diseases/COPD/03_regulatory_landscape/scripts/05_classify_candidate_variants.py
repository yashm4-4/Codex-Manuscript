#!/usr/bin/env python3
"""Classify the frozen COPD tag plus LD-proxy set against regulatory features.

The primary output retains all candidate records from COPD-S2-R006E, including
the three tags without BED-eligible coordinates.  BED-eligible records are
intersected with GENCODE v50 CDS intervals, the severe-emphysema donor element
unions, SCREEN V3, Ensembl Regulatory Build 116, FANTOM5 enhancers, UCSC
RepeatMasker, and the ENCODE hg38 blacklist.

All coordinates are GRCh38 and all overlap tests require at least one base.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, Iterator, Optional

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION = SCRIPT.parents[1]
PROJECT = SECTION.parents[2]
GWAS = PROJECT / "diseases" / "COPD" / "02_gwas"
SHARED = PROJECT / "data"
RESULTS = SECTION / "results"
DATA_OUT = SECTION / "data" / "candidate_classification"
INTERSECTIONS = SECTION / "data" / "candidate_intersections"
TMP = SECTION / "tmp"
CANONICAL_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CANONICAL = set(CANONICAL_ORDER)
CHROM_RANK = {chrom: rank for rank, chrom in enumerate(CANONICAL_ORDER)}
PANELS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
NOT_ELIGIBLE = "NA_not_bed_eligible"


@dataclass(frozen=True)
class FeatureSpec:
    key: str
    flag_column: str
    path: Path
    bed_columns: int
    parser_kind: str
    source_name: str
    source_id: str
    definition: str


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates",
        type=Path,
        default=GWAS / "results" / "COPD-S2-R006E_candidate_variants_grch38.tsv.gz",
    )
    parser.add_argument(
        "--candidate-bed",
        type=Path,
        default=GWAS / "results" / "COPD-S2-R006E_candidate_variants_grch38.bed",
    )
    parser.add_argument(
        "--gencode-gtf",
        type=Path,
        default=SHARED / "gencode" / "gencode.v50.annotation.gtf.gz",
    )
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--data-dir", type=Path, default=DATA_OUT)
    parser.add_argument("--intersection-dir", type=Path, default=INTERSECTIONS)
    parser.add_argument(
        "--bedtools",
        type=Path,
        default=Path("/usr/local/apps/bedtools/2.31.1/bin/bedtools"),
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_path(path: Path) -> str:
    absolute = path if path.is_absolute() else Path.cwd() / path
    try:
        return str(absolute.relative_to(PROJECT))
    except ValueError:
        return str(absolute)


def write_deterministic_gzip_rows(
    path: Path, rows: Iterable[Iterable[object]]
) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as out:
                for row in rows:
                    out.write("\t".join(str(value) for value in row) + "\n")
    temporary.replace(path)


def write_deterministic_dataframe(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as out:
                frame.to_csv(out, sep="\t", index=False, lineterminator="\n")
    temporary.replace(path)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def prepare_cds(gtf: Path, output: Path) -> int:
    """Write unique gene-aware GENCODE CDS intervals as a six-column BED."""
    records: set[tuple[str, int, int, str, str, str]] = set()
    with gzip.open(gtf, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS" or fields[0] not in CANONICAL:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"{gtf}:{line_number}: invalid CDS coordinates")
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            gene_type = attrs.get("gene_type", "")
            if not gene_id or not gene_name or not gene_type:
                raise ValueError(f"{gtf}:{line_number}: incomplete CDS gene annotation")
            records.add((fields[0], start, end, gene_id, gene_name, gene_type))
    ordered = sorted(
        records,
        key=lambda row: (CHROM_RANK[row[0]], row[1], row[2], row[3], row[4], row[5]),
    )
    write_deterministic_gzip_rows(output, ordered)
    return len(ordered)


def read_candidates(candidate_path: Path, bed_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(
        candidate_path, sep="\t", compression="gzip", dtype=str, keep_default_na=False
    )
    expected_columns = {
        "candidate_record_id",
        "chromosome_grch38",
        "position_grch38",
        "ref",
        "alt",
        "variant_class",
        "tag_supported_panels",
        "tag_polymorphic_panels",
        "source_focal_tags",
        "ld_panels",
        "is_gws_tag",
        "is_ld_proxy",
        "candidate_origin",
        "bed_eligible",
    }
    missing = expected_columns - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table lacks columns: {sorted(missing)}")
    if candidates["candidate_record_id"].duplicated().any():
        raise ValueError("Candidate record identifiers are not unique")
    if len(candidates) != 15_389:
        raise ValueError(f"Expected 15,389 candidate records, found {len(candidates)}")

    bed = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["bed_chrom", "bed_start", "bed_end", "candidate_record_id", "score", "strand"],
        dtype={"bed_chrom": str, "bed_start": int, "bed_end": int, "candidate_record_id": str},
    )
    if len(bed) != 15_386 or bed["candidate_record_id"].duplicated().any():
        raise ValueError("Candidate BED must contain 15,386 unique records")
    eligible_ids = set(candidates.loc[candidates["bed_eligible"].eq("True"), "candidate_record_id"])
    if eligible_ids != set(bed["candidate_record_id"]):
        raise ValueError("BED identifiers do not equal the BED-eligible candidate set")
    if int(candidates["bed_eligible"].eq("False").sum()) != 3:
        raise ValueError("Expected exactly three BED-ineligible candidate records")

    candidate_by_id = candidates.set_index("candidate_record_id")
    for row in bed.itertuples(index=False):
        source = candidate_by_id.loc[row.candidate_record_id]
        expected_chrom = f"chr{source['chromosome_grch38']}"
        expected_start = int(source["position_grch38"]) - 1
        expected_span = max(1, len(source["ref"]))
        if (
            row.bed_chrom != expected_chrom
            or row.bed_start != expected_start
            or row.bed_end != expected_start + expected_span
        ):
            raise ValueError(f"BED coordinate mismatch for {row.candidate_record_id}")
    return candidates, bed


def resolve_bedtools(requested: Path) -> Path:
    if requested.is_file() and requested.stat().st_mode & 0o111:
        return requested
    discovered = shutil.which("bedtools")
    if discovered:
        return Path(discovered)
    raise FileNotFoundError(f"bedtools not executable: {requested}")


def feature_specs(cds: Path) -> list[FeatureSpec]:
    donor = SECTION / "data" / "copd_donor_regulatory_elements"
    refs = SECTION / "data" / "regulatory_references"
    return [
        FeatureSpec(
            "coding_CDS",
            "coding_CDS",
            cds,
            6,
            "cds",
            "GENCODE v50 comprehensive gene annotation",
            "COPD-SRC-043",
            "At least 1 bp overlap of a GENCODE v50 CDS interval",
        ),
        FeatureSpec(
            "donor_preliminary_enhancer",
            "donor_preliminary_enhancer",
            donor / "donor_union.enhancer.preliminary_histone_peak.grch38.bed.gz",
            8,
            "donor",
            "ENCODE severe-emphysema donor ENCDO520EJG",
            "COPD-SRC-018",
            "At least 1 bp overlap of the three-lobe union of H3K27ac peaks",
        ),
        FeatureSpec(
            "donor_preliminary_silencer",
            "donor_preliminary_silencer",
            donor / "donor_union.silencer.preliminary_histone_peak.grch38.bed.gz",
            8,
            "donor",
            "ENCODE severe-emphysema donor ENCDO520EJG",
            "COPD-SRC-018",
            "At least 1 bp overlap of the three-lobe union of H3K27me3 peaks",
        ),
        FeatureSpec(
            "donor_refined_enhancer",
            "donor_refined_enhancer",
            donor
            / "donor_union.enhancer.refined_accessibility_histone_overlap.grch38.bed.gz",
            8,
            "donor",
            "ENCODE severe-emphysema donor ENCDO520EJG",
            "COPD-SRC-018",
            "At least 1 bp overlap of the three-lobe union of same-lobe ATAC plus H3K27ac elements",
        ),
        FeatureSpec(
            "donor_refined_silencer",
            "donor_refined_silencer",
            donor
            / "donor_union.silencer.refined_accessibility_histone_overlap.grch38.bed.gz",
            8,
            "donor",
            "ENCODE severe-emphysema donor ENCDO520EJG",
            "COPD-SRC-018",
            "At least 1 bp overlap of the three-lobe union of same-lobe ATAC plus H3K27me3 elements",
        ),
        FeatureSpec(
            "screen_ccre",
            "screen_ccre",
            refs / "screen_ccre_v3.grch38.canonical.bed.gz",
            6,
            "generic",
            "ENCODE SCREEN cCRE Registry V3",
            "COPD-SRC-038",
            "At least 1 bp overlap of a Registry V3 cCRE",
        ),
        FeatureSpec(
            "ensembl_regulatory",
            "ensembl_regulatory",
            refs / "ensembl_regulatory_v116.grch38.canonical.bed.gz",
            6,
            "generic",
            "Ensembl Regulatory Build 116",
            "COPD-SRC-039",
            "At least 1 bp overlap of an Ensembl 116 regulatory feature",
        ),
        FeatureSpec(
            "fantom5_enhancer",
            "fantom5_enhancer",
            refs / "fantom5_cage_enhancers.grch38.canonical.bed.gz",
            6,
            "generic",
            "FANTOM5 CAGE enhancer atlas",
            "COPD-SRC-040",
            "At least 1 bp overlap of a FANTOM5 hg38 CAGE enhancer",
        ),
        FeatureSpec(
            "repeatmasker",
            "repeatmasker",
            refs / "ucsc_repeatmasker.grch38.canonical.bed.gz",
            6,
            "generic",
            "UCSC RepeatMasker hg38 table",
            "COPD-SRC-041",
            "At least 1 bp overlap of a UCSC hg38 RepeatMasker interval",
        ),
        FeatureSpec(
            "encode_blacklist",
            "encode_blacklist",
            refs / "encode_blacklist_v2.grch38.canonical.bed.gz",
            6,
            "generic",
            "ENCODE hg38 blacklist v2",
            "COPD-SRC-042",
            "At least 1 bp overlap of an ENCODE hg38 blacklist v2 interval",
        ),
    ]


def parse_feature(fields: list[str], spec: FeatureSpec) -> tuple[str, str, str]:
    if spec.parser_kind == "cds":
        gene_id, gene_name, gene_type = fields[3], fields[4], fields[5]
        return gene_id, "CDS", f"gene_name={gene_name};gene_type={gene_type}"
    if spec.parser_kind == "donor":
        feature_id, element_type = fields[3], fields[4]
        detail = (
            f"definition={fields[5]};context={fields[6]};"
            f"source_accessions={fields[7]}"
        )
        return feature_id, element_type, detail
    return fields[3], fields[4], fields[5]


def run_intersection(
    spec: FeatureSpec,
    candidate_bed: Path,
    bedtools: Path,
    output_dir: Path,
    temp_dir: Path,
) -> tuple[dict[str, list[tuple[str, str, str, str, int, int, int]]], int]:
    """Run bedtools and retain a normalized raw-hit table for one feature set."""
    hits: dict[str, list[tuple[str, str, str, str, int, int, int]]] = defaultdict(list)
    normalized_rows: list[tuple[object, ...]] = []
    with tempfile.NamedTemporaryFile(mode="w+", dir=temp_dir, prefix=f"{spec.key}.") as raw:
        completed = subprocess.run(
            [
                str(bedtools),
                "intersect",
                "-wo",
                "-a",
                str(candidate_bed),
                "-b",
                str(spec.path),
            ],
            text=True,
            stdout=raw,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"bedtools failed for {spec.key}: {completed.stderr}")
        raw.flush()
        raw.seek(0)
        for line_number, line in enumerate(raw, start=1):
            fields = line.rstrip("\n").split("\t")
            expected = 6 + spec.bed_columns + 1
            if len(fields) != expected:
                raise ValueError(
                    f"{spec.key} intersection line {line_number}: {len(fields)} fields, expected {expected}"
                )
            candidate_fields = fields[:6]
            feature_fields = fields[6 : 6 + spec.bed_columns]
            overlap_bp = int(fields[-1])
            feature_id, feature_type, detail = parse_feature(feature_fields, spec)
            candidate_id = candidate_fields[3]
            record = (
                feature_id,
                feature_type,
                detail,
                feature_fields[0],
                int(feature_fields[1]),
                int(feature_fields[2]),
                overlap_bp,
            )
            hits[candidate_id].append(record)
            normalized_rows.append(
                (
                    candidate_id,
                    candidate_fields[0],
                    candidate_fields[1],
                    candidate_fields[2],
                    feature_id,
                    feature_type,
                    feature_fields[0],
                    feature_fields[1],
                    feature_fields[2],
                    detail,
                    overlap_bp,
                    spec.key,
                )
            )

    header = (
        "candidate_record_id",
        "candidate_chrom",
        "candidate_start",
        "candidate_end",
        "feature_id",
        "feature_type",
        "feature_chrom",
        "feature_start",
        "feature_end",
        "source_detail",
        "overlap_bp",
        "feature_set",
    )
    output = output_dir / f"COPD-S3-R003_{spec.key}_intersections.tsv.gz"
    write_deterministic_gzip_rows(output, [header, *normalized_rows])
    return hits, len(normalized_rows)


def semicolon(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def annotate_candidates(
    candidates: pd.DataFrame,
    bed: pd.DataFrame,
    specs: list[FeatureSpec],
    feature_hits: dict[str, dict[str, list[tuple[str, str, str, str, int, int, int]]]],
) -> pd.DataFrame:
    annotated = candidates.copy()
    bed_lookup = bed.set_index("candidate_record_id")
    annotated["bed_chrom"] = annotated["candidate_record_id"].map(bed_lookup["bed_chrom"])
    annotated["bed_start"] = (
        annotated["candidate_record_id"]
        .map(bed_lookup["bed_start"])
        .map(lambda value: "" if pd.isna(value) else str(int(value)))
    )
    annotated["bed_end"] = (
        annotated["candidate_record_id"]
        .map(bed_lookup["bed_end"])
        .map(lambda value: "" if pd.isna(value) else str(int(value)))
    )
    annotated["bed_chrom"] = annotated["bed_chrom"].fillna("")

    eligible = annotated["bed_eligible"].eq("True")
    for spec in specs:
        by_candidate = feature_hits[spec.key]
        flags: list[str] = []
        counts: list[object] = []
        identifiers: list[str] = []
        types: list[str] = []
        details: list[str] = []
        for candidate_id, is_eligible in zip(annotated["candidate_record_id"], eligible):
            if not is_eligible:
                flags.append(NOT_ELIGIBLE)
                counts.append("")
                identifiers.append("")
                types.append("")
                details.append("")
                continue
            records = by_candidate.get(candidate_id, [])
            unique_records = sorted(set(records))
            flags.append("True" if unique_records else "False")
            counts.append(len(unique_records))
            identifiers.append(semicolon(record[0] for record in unique_records))
            types.append(semicolon(record[1] for record in unique_records))
            details.append(semicolon(record[2] for record in unique_records))
        annotated[spec.flag_column] = flags
        annotated[f"{spec.key}_n_intervals"] = counts
        annotated[f"{spec.key}_feature_ids"] = identifiers
        annotated[f"{spec.key}_feature_types"] = types
        annotated[f"{spec.key}_source_details"] = details

    annotated["coding_CDS_gene_ids"] = annotated["coding_CDS_feature_ids"]
    annotated["coding_CDS_gene_names"] = annotated["coding_CDS_source_details"].map(
        lambda text: semicolon(
            item.split("=", 1)[1]
            for detail in text.split(";")
            for item in [detail]
            if item.startswith("gene_name=")
        )
    )
    annotated["coding_CDS_gene_types"] = annotated["coding_CDS_source_details"].map(
        lambda text: semicolon(
            item.split("=", 1)[1]
            for detail in text.split(";")
            for item in [detail]
            if item.startswith("gene_type=")
        )
    )

    def derived_flag(row: pd.Series, columns: list[str], mode: str = "any") -> str:
        if row["bed_eligible"] != "True":
            return NOT_ELIGIBLE
        values = [row[column] == "True" for column in columns]
        value = any(values) if mode == "any" else all(values)
        return "True" if value else "False"

    annotated["donor_preliminary_any"] = annotated.apply(
        derived_flag,
        axis=1,
        columns=["donor_preliminary_enhancer", "donor_preliminary_silencer"],
    )
    annotated["donor_refined_any"] = annotated.apply(
        derived_flag,
        axis=1,
        columns=["donor_refined_enhancer", "donor_refined_silencer"],
    )
    annotated["known_regulatory_any"] = annotated.apply(
        derived_flag,
        axis=1,
        columns=["screen_ccre", "ensembl_regulatory", "fantom5_enhancer"],
    )
    biological = [
        "coding_CDS",
        "donor_preliminary_enhancer",
        "donor_preliminary_silencer",
        "donor_refined_enhancer",
        "donor_refined_silencer",
        "screen_ccre",
        "ensembl_regulatory",
        "fantom5_enhancer",
        "repeatmasker",
    ]
    annotated["any_tested_biological_annotation"] = annotated.apply(
        derived_flag, axis=1, columns=biological
    )
    annotated["unannotated_across_tested_biological_features"] = annotated.apply(
        lambda row: (
            NOT_ELIGIBLE
            if row["bed_eligible"] != "True"
            else ("True" if row["any_tested_biological_annotation"] == "False" else "False")
        ),
        axis=1,
    )

    def hierarchy(row: pd.Series, version: str) -> str:
        if row["bed_eligible"] != "True":
            return "bed_ineligible"
        if row["coding_CDS"] == "True":
            return "coding_CDS"
        if version == "framework_preliminary":
            if row["donor_preliminary_enhancer"] == "True":
                return "preliminary_enhancer"
            if row["donor_preliminary_silencer"] == "True":
                return "preliminary_silencer"
            if row["repeatmasker"] == "True":
                return "repetitive_only_in_hierarchy"
            return "other"
        if version == "framework_refined":
            if row["donor_refined_enhancer"] == "True":
                return "refined_enhancer"
            if row["donor_refined_silencer"] == "True":
                return "refined_silencer"
            if row["repeatmasker"] == "True":
                return "repetitive_only_in_hierarchy"
            return "other"
        if row["donor_refined_enhancer"] == "True":
            return "refined_enhancer"
        if row["donor_refined_silencer"] == "True":
            return "refined_silencer"
        if row["donor_preliminary_enhancer"] == "True":
            return "preliminary_enhancer"
        if row["donor_preliminary_silencer"] == "True":
            return "preliminary_silencer"
        if row["known_regulatory_any"] == "True":
            return "other_known_regulatory"
        if row["repeatmasker"] == "True":
            return "repetitive_only_in_hierarchy"
        return "other"

    for version in ("framework_preliminary", "framework_refined", "comprehensive"):
        annotated[f"exclusive_class_{version}"] = annotated.apply(
            hierarchy, axis=1, version=version
        )
    return annotated


def panel_member(series: pd.Series, panel: str) -> pd.Series:
    return series.map(lambda value: panel in value.split(";") if value else False)


def strata(frame: pd.DataFrame) -> list[tuple[str, str, pd.Series]]:
    groups: list[tuple[str, str, pd.Series]] = [
        ("overall", "all_candidate_records", pd.Series(True, index=frame.index)),
    ]
    for origin in ("gws_tag_only", "gws_tag_and_ld_proxy", "ld_proxy_only"):
        groups.append(("candidate_origin", origin, frame["candidate_origin"].eq(origin)))
    groups.extend(
        [
            ("role", "gws_tag_record", frame["is_gws_tag"].eq("True")),
            ("role", "ld_proxy_record", frame["is_ld_proxy"].eq("True")),
        ]
    )
    for variant_class in ("SNV", "indel_or_complex", "unresolved_tag"):
        groups.append(
            ("variant_class", variant_class, frame["variant_class"].eq(variant_class))
        )
    groups.extend(
        [
            (
                "reference_status",
                "reference_matched_record",
                frame["reference_matched"].eq("True"),
            ),
            (
                "reference_status",
                "reference_unmatched_tag_record",
                frame["reference_matched"].eq("False"),
            ),
        ]
    )
    for panel in PANELS:
        groups.append(("ld_panel", panel, panel_member(frame["ld_panels"], panel)))
    for panel in PANELS:
        groups.append(
            ("tag_supported_panel", panel, panel_member(frame["tag_supported_panels"], panel))
        )
    for panel in PANELS:
        groups.append(
            (
                "tag_polymorphic_panel",
                panel,
                panel_member(frame["tag_polymorphic_panels"], panel),
            )
        )
    return groups


def summarize_nonexclusive(frame: pd.DataFrame) -> pd.DataFrame:
    feature_definitions = [
        ("coding_CDS", "GENCODE v50 CDS"),
        ("donor_preliminary_enhancer", "donor-union H3K27ac"),
        ("donor_preliminary_silencer", "donor-union H3K27me3"),
        ("donor_refined_enhancer", "donor-union ATAC plus H3K27ac"),
        ("donor_refined_silencer", "donor-union ATAC plus H3K27me3"),
        ("donor_preliminary_any", "either preliminary donor element"),
        ("donor_refined_any", "either refined donor element"),
        ("screen_ccre", "SCREEN V3 cCRE"),
        ("ensembl_regulatory", "Ensembl Regulatory Build 116"),
        ("fantom5_enhancer", "FANTOM5 CAGE enhancer"),
        ("known_regulatory_any", "SCREEN, Ensembl, or FANTOM5"),
        ("repeatmasker", "UCSC RepeatMasker"),
        ("encode_blacklist", "ENCODE hg38 blacklist v2"),
        ("any_tested_biological_annotation", "coding, donor, known regulatory, or repeat"),
        (
            "unannotated_across_tested_biological_features",
            "no coding, donor, known-regulatory, or repeat overlap",
        ),
    ]
    rows: list[dict[str, object]] = []
    for stratum_type, stratum, mask in strata(frame):
        subset = frame.loc[mask]
        n_total = len(subset)
        n_eligible = int(subset["bed_eligible"].eq("True").sum())
        for feature, definition in feature_definitions:
            n_overlap = int(subset[feature].eq("True").sum())
            rows.append(
                {
                    "stratum_type": stratum_type,
                    "stratum": stratum,
                    "feature": feature,
                    "definition": definition,
                    "n_total_records": n_total,
                    "n_bed_eligible": n_eligible,
                    "n_bed_ineligible": n_total - n_eligible,
                    "n_overlapping": n_overlap,
                    "percent_of_all_records": round(n_overlap / n_total * 100, 6)
                    if n_total
                    else 0.0,
                    "percent_of_bed_eligible": round(n_overlap / n_eligible * 100, 6)
                    if n_eligible
                    else 0.0,
                }
            )
    return pd.DataFrame(rows)


HIERARCHY_CATEGORIES = {
    "framework_preliminary": [
        "bed_ineligible",
        "coding_CDS",
        "preliminary_enhancer",
        "preliminary_silencer",
        "repetitive_only_in_hierarchy",
        "other",
    ],
    "framework_refined": [
        "bed_ineligible",
        "coding_CDS",
        "refined_enhancer",
        "refined_silencer",
        "repetitive_only_in_hierarchy",
        "other",
    ],
    "comprehensive": [
        "bed_ineligible",
        "coding_CDS",
        "refined_enhancer",
        "refined_silencer",
        "preliminary_enhancer",
        "preliminary_silencer",
        "other_known_regulatory",
        "repetitive_only_in_hierarchy",
        "other",
    ],
}


def summarize_exclusive(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for stratum_type, stratum, mask in strata(frame):
        subset = frame.loc[mask]
        n_total = len(subset)
        n_eligible = int(subset["bed_eligible"].eq("True").sum())
        for hierarchy, categories in HIERARCHY_CATEGORIES.items():
            column = f"exclusive_class_{hierarchy}"
            observed_total = 0
            for rank, category in enumerate(categories):
                count = int(subset[column].eq(category).sum())
                observed_total += count
                rows.append(
                    {
                        "stratum_type": stratum_type,
                        "stratum": stratum,
                        "hierarchy": hierarchy,
                        "precedence_rank": rank,
                        "exclusive_class": category,
                        "n_total_records": n_total,
                        "n_bed_eligible": n_eligible,
                        "n_records": count,
                        "percent_of_all_records": round(count / n_total * 100, 6)
                        if n_total
                        else 0.0,
                        "percent_of_bed_eligible": (
                            ""
                            if category == "bed_ineligible"
                            else round(count / n_eligible * 100, 6)
                            if n_eligible
                            else 0.0
                        ),
                    }
                )
            if observed_total != n_total:
                raise RuntimeError(
                    f"Exclusive hierarchy {hierarchy} does not sum for {stratum_type}:{stratum}"
                )
    return pd.DataFrame(rows)


def make_tag_identifier_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Expand shared reference records back to the 660 Catalog tag labels."""
    tags = frame.loc[frame["is_gws_tag"].eq("True")].copy()
    tags["gws_tag_id"] = tags["gws_tag_ids"].str.split(";")
    tags = tags.explode("gws_tag_id", ignore_index=True)
    if len(tags) != 660 or tags["gws_tag_id"].duplicated().any():
        raise RuntimeError("GWS tag expansion did not produce 660 unique tag identifiers")
    columns = ["gws_tag_id"] + [column for column in tags.columns if column != "gws_tag_id"]
    return tags[columns]


def identifier_level_summaries(
    tag_identifiers: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    nonexclusive = summarize_nonexclusive(tag_identifiers)
    nonexclusive = nonexclusive.loc[nonexclusive["stratum_type"].eq("overall")].copy()
    nonexclusive["stratum_type"] = "identifier_level"
    nonexclusive["stratum"] = "gws_tag_identifiers"
    exclusive = summarize_exclusive(tag_identifiers)
    exclusive = exclusive.loc[exclusive["stratum_type"].eq("overall")].copy()
    exclusive["stratum_type"] = "identifier_level"
    exclusive["stratum"] = "gws_tag_identifiers"
    return nonexclusive, exclusive


def hierarchy_definition_rows() -> list[dict[str, object]]:
    descriptions = {
        "framework_preliminary": (
            "bed ineligible is separated; eligible precedence is coding CDS, preliminary "
            "enhancer, preliminary silencer, RepeatMasker, other"
        ),
        "framework_refined": (
            "bed ineligible is separated; eligible precedence is coding CDS, refined "
            "enhancer, refined silencer, RepeatMasker, other"
        ),
        "comprehensive": (
            "bed ineligible is separated; eligible precedence is coding CDS, refined "
            "enhancer, refined silencer, preliminary enhancer, preliminary silencer, "
            "other known regulatory, RepeatMasker, other"
        ),
    }
    rows: list[dict[str, object]] = []
    for hierarchy, categories in HIERARCHY_CATEGORIES.items():
        for rank, category in enumerate(categories):
            rows.append(
                {
                    "hierarchy": hierarchy,
                    "precedence_rank": rank,
                    "exclusive_class": category,
                    "hierarchy_definition": descriptions[hierarchy],
                }
            )
    return rows


def count_lines(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return sum(1 for _ in handle)


def main() -> None:
    args = arguments()
    for directory in (args.results_dir, args.data_dir, args.intersection_dir, TMP):
        directory.mkdir(parents=True, exist_ok=True)
    bedtools = resolve_bedtools(args.bedtools)
    bedtools_version = subprocess.run(
        [str(bedtools), "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()

    candidates, candidate_bed = read_candidates(args.candidates, args.candidate_bed)
    cds_bed = args.data_dir / "gencode_v50_cds.grch38.canonical.bed.gz"
    cds_interval_count = prepare_cds(args.gencode_gtf, cds_bed)
    specs = feature_specs(cds_bed)
    for spec in specs:
        if not spec.path.is_file():
            raise FileNotFoundError(
                f"Missing {spec.path}; run the Section 3 reference and donor-element preparation"
            )

    all_hits: dict[
        str, dict[str, list[tuple[str, str, str, str, int, int, int]]]
    ] = {}
    definition_rows: list[dict[str, object]] = []
    for spec in specs:
        hits, raw_rows = run_intersection(
            spec, args.candidate_bed, bedtools, args.intersection_dir, TMP
        )
        all_hits[spec.key] = hits
        raw_output = args.intersection_dir / f"COPD-S3-R003_{spec.key}_intersections.tsv.gz"
        definition_rows.append(
            {
                "feature": spec.flag_column,
                "source_name": spec.source_name,
                "source_id": spec.source_id,
                "definition": spec.definition,
                "genome_build": "GRCh38",
                "overlap_minimum_bp": 1,
                "feature_bed_path": project_path(spec.path),
                "feature_bed_sha256": sha256(spec.path),
                "n_feature_intervals": count_lines(spec.path),
                "raw_intersection_path": project_path(raw_output),
                "raw_intersection_sha256": sha256(raw_output),
                "n_raw_intersection_rows": raw_rows,
                "n_candidate_records_overlapping": len(hits),
            }
        )
        print(f"{spec.key}\t{len(hits)} candidate records\t{raw_rows} raw overlaps")

    annotated = annotate_candidates(candidates, candidate_bed, specs, all_hits)
    nonexclusive = summarize_nonexclusive(annotated)
    exclusive = summarize_exclusive(annotated)
    tag_identifiers = make_tag_identifier_table(annotated)
    tag_nonexclusive, tag_exclusive = identifier_level_summaries(tag_identifiers)
    nonexclusive = pd.concat([nonexclusive, tag_nonexclusive], ignore_index=True)
    exclusive = pd.concat([exclusive, tag_exclusive], ignore_index=True)
    definitions = pd.DataFrame(definition_rows)
    hierarchy_definitions = pd.DataFrame(hierarchy_definition_rows())

    if len(annotated) != 15_389 or annotated["candidate_record_id"].duplicated().any():
        raise RuntimeError("Annotated table does not preserve all unique input records")
    ineligible = annotated.loc[annotated["bed_eligible"].eq("False")]
    expected_ineligible = {
        "unmatched_tag:15q25.1",
        "unmatched_tag:rs139284640",
        "unmatched_tag:rs751872749",
    }
    if set(ineligible["candidate_record_id"]) != expected_ineligible:
        raise RuntimeError("Unexpected BED-ineligible record set")
    for spec in specs:
        if not ineligible[spec.flag_column].eq(NOT_ELIGIBLE).all():
            raise RuntimeError(f"Ineligible rows are not explicit for {spec.flag_column}")

    output_frames = {
        "COPD-S3-R003_nonexclusive_summary.tsv": nonexclusive,
        "COPD-S3-R003_exclusive_summary.tsv": exclusive,
        "COPD-S3-R003_feature_definitions.tsv": definitions,
        "COPD-S3-R003_hierarchy_definitions.tsv": hierarchy_definitions,
    }
    for filename, frame in output_frames.items():
        path = args.results_dir / filename
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    annotation_output = args.results_dir / "COPD-S3-R003_candidate_classification.tsv.gz"
    write_deterministic_dataframe(annotation_output, annotated)
    tag_output = args.results_dir / "COPD-S3-R003_gws_tag_identifier_classification.tsv.gz"
    write_deterministic_dataframe(tag_output, tag_identifiers)

    output_metadata: list[dict[str, object]] = []
    for path, rows in [
        (annotation_output, len(annotated)),
        (tag_output, len(tag_identifiers)),
        *[(args.results_dir / name, len(frame)) for name, frame in output_frames.items()],
    ]:
        output_metadata.append(
            {
                "path": project_path(path),
                "rows": rows,
                "sha256": sha256(path),
            }
        )

    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exact_replay_command": shlex.join([sys.executable, str(SCRIPT), *sys.argv[1:]]),
        "working_directory": str(Path.cwd()),
        "script": project_path(SCRIPT),
        "script_sha256": sha256(SCRIPT),
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "bedtools": str(bedtools),
        "bedtools_version": bedtools_version,
        "genome_build": "GRCh38",
        "overlap_rule": "at least 1 bp",
        "candidate_input": project_path(args.candidates),
        "candidate_input_sha256": sha256(args.candidates),
        "candidate_input_rows": len(candidates),
        "candidate_bed": project_path(args.candidate_bed),
        "candidate_bed_sha256": sha256(args.candidate_bed),
        "candidate_bed_rows": len(candidate_bed),
        "bed_ineligible_record_ids": sorted(expected_ineligible),
        "gencode_gtf": project_path(args.gencode_gtf),
        "gencode_gtf_sha256": sha256(args.gencode_gtf),
        "gencode_release": "v50",
        "gencode_cds_unique_intervals": cds_interval_count,
        "feature_definitions": definition_rows,
        "hierarchies": {
            hierarchy: categories for hierarchy, categories in HIERARCHY_CATEGORIES.items()
        },
        "stratum_note": (
            "Candidate summaries use reference-record units except the explicit 660-row "
            "gws_tag_identifiers stratum. LD-panel, tag-supported-panel, and "
            "tag-polymorphic-panel strata are multi-membership and must not be summed "
            "across panels"
        ),
        "outputs": output_metadata,
    }
    manifest_path = args.results_dir / "COPD-S3-R003_analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    overall = nonexclusive.loc[
        (nonexclusive["stratum_type"] == "overall")
        & (nonexclusive["stratum"] == "all_candidate_records")
    ]
    print(f"candidate_records\t{len(annotated)}")
    print(f"bed_eligible\t{int(annotated['bed_eligible'].eq('True').sum())}")
    print(f"bed_ineligible\t{len(ineligible)}")
    for row in overall.itertuples(index=False):
        print(
            f"summary\t{row.feature}\t{row.n_overlapping}\t"
            f"{row.percent_of_bed_eligible:.6f}% eligible"
        )
    print(f"manifest\t{project_path(manifest_path)}")


if __name__ == "__main__":
    main()
