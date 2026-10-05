#!/usr/bin/env python3
"""Validate and summarize COPD ancestry-aware 1000 Genomes LD expansion.

This step retains all 660 genome-wide-significant tags in the focal audit and
candidate table, including tags that cannot be matched to the local reference.
Only strictly matched, panel-polymorphic focal records can contribute LD
proxies.  The final GRCh38 candidate table and BED include tag and proxy
provenance for downstream regulatory overlap analyses.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import math
import platform
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
WORK = SECTION_ROOT / "data" / "ld_work"
RESULTS = SECTION_ROOT / "results"
LOGS = SECTION_ROOT / "logs"
REFERENCE_RELEASE = "1000 Genomes Phase 3, 20190312 GRCh38 remap"
PANEL_ORDER = ("AFR", "AMR", "EAS", "EUR", "SAS")
VARIANT_ID_RE = re.compile(r"^([^:]+):(\d+):([^:]+):([^:]+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--focal", type=Path, default=WORK / "focal_variants.tsv"
    )
    parser.add_argument(
        "--assignments", type=Path, default=WORK / "focal_panel_assignments.tsv"
    )
    parser.add_argument(
        "--match-audit",
        type=Path,
        default=RESULTS / "COPD-S2-R006A_focal_match_audit.tsv",
    )
    parser.add_argument(
        "--job-status", type=Path, default=WORK / "job_status.tsv"
    )
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--window-bp", type=int, default=500_000)
    parser.add_argument("--r2", type=float, default=0.8)
    parser.add_argument(
        "--log", type=Path, default=LOGS / "COPD-S2-ld_summary.log"
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2-LD-summary")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ\t%(levelname)s\t%(message)s")
    formatter.converter = lambda *args: datetime.now(timezone.utc).timetuple()
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_variant_id(value: object) -> tuple[str, int, str, str]:
    match = VARIANT_ID_RE.fullmatch(str(value))
    if match is None:
        raise ValueError(f"invalid reference variant ID: {value!r}")
    return match.group(1), int(match.group(2)), match.group(3), match.group(4)


def variant_class(ref: str, alt: str) -> str:
    return "SNV" if len(ref) == 1 and len(alt) == 1 else "indel_or_complex"


def current_job_paths(job_status: pd.DataFrame, filename: str) -> list[Path]:
    paths = [Path(value) / filename for value in job_status["job_dir"]]
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} current-job files are missing; first: {missing[0]}"
        )
    return paths


def read_focal_frequencies(job_status: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for path in current_job_paths(job_status, "focal_frequency.frq"):
        panel = path.parent.name.rsplit("_", 1)[-1]
        frame = pd.read_csv(path, sep=r"\s+")
        require_columns(frame, ["SNP", "MAF", "NCHROBS"], str(path))
        frame["panel"] = panel
        rows.append(frame[["panel", "SNP", "MAF", "NCHROBS"]])
    if not rows:
        raise FileNotFoundError("no focal_frequency.frq files found")
    frequencies = pd.concat(rows, ignore_index=True)
    if frequencies.duplicated(["panel", "SNP"]).any():
        duplicate = frequencies.loc[
            frequencies.duplicated(["panel", "SNP"], keep=False)
        ].head()
        raise AssertionError(f"duplicate panel-frequency rows:\n{duplicate}")
    return frequencies


def build_panel_audit(
    assignments: pd.DataFrame,
    match_audit: pd.DataFrame,
    frequencies: pd.DataFrame,
) -> pd.DataFrame:
    selected_columns = [
        "normalized_variant_id",
        "matched",
        "match_status",
        "match_basis",
        "position_delta_bp",
        "panel_variant_id",
        "panel_chromosome",
        "panel_position_grch38",
        "panel_ref",
        "panel_alt",
        "panel_variant_class",
    ]
    panel_audit = assignments.merge(
        match_audit[selected_columns],
        on="normalized_variant_id",
        how="left",
        validate="many_to_one",
    )
    panel_audit = panel_audit.merge(
        frequencies,
        left_on=["panel", "panel_variant_id"],
        right_on=["panel", "SNP"],
        how="left",
        validate="many_to_one",
    ).drop(columns=["SNP"])

    def status(row: pd.Series) -> str:
        if not bool(row["matched"]):
            return "not_run_global_focal_match_failed"
        if pd.isna(row["MAF"]):
            return "matched_record_missing_from_panel_subset"
        if float(row["MAF"]) <= 0 or float(row["NCHROBS"]) <= 0:
            return "matched_monomorphic_in_panel"
        return "matched_polymorphic_in_panel"

    panel_audit["panel_evaluation_status"] = panel_audit.apply(status, axis=1)
    panel_audit["reference_release"] = REFERENCE_RELEASE
    return panel_audit.sort_values(
        ["panel", "normalized_variant_id"], kind="stable"
    ).reset_index(drop=True)


def read_ld_pairs(job_status: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    rows: list[pd.DataFrame] = []
    raw_count = 0
    for path in current_job_paths(job_status, "ld.ld"):
        panel = path.parent.name.rsplit("_", 1)[-1]
        frame = pd.read_csv(path, sep=r"\s+")
        required = ["CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"]
        require_columns(frame, required, str(path))
        raw_count += len(frame)
        if not frame.empty:
            frame["panel"] = panel
            rows.append(frame)
    if not rows:
        return pd.DataFrame(
            columns=[
                "CHR_A",
                "BP_A",
                "SNP_A",
                "CHR_B",
                "BP_B",
                "SNP_B",
                "R2",
                "panel",
            ]
        ), raw_count
    return pd.concat(rows, ignore_index=True), raw_count


def validate_and_expand_pairs(
    raw: pd.DataFrame,
    panel_audit: pd.DataFrame,
    r2_threshold: float,
    window_bp: int,
) -> pd.DataFrame:
    if raw.empty:
        return raw
    raw = raw.copy()
    raw["R2"] = pd.to_numeric(raw["R2"], errors="raise")
    raw["BP_A"] = pd.to_numeric(raw["BP_A"], errors="raise").astype(int)
    raw["BP_B"] = pd.to_numeric(raw["BP_B"], errors="raise").astype(int)
    if (raw["R2"] < r2_threshold - 1e-12).any() or (raw["R2"] > 1 + 1e-12).any():
        raise AssertionError("PLINK LD output violates requested r2 threshold/range")
    if (raw["CHR_A"].astype(str) != raw["CHR_B"].astype(str)).any():
        raise AssertionError("interchromosomal row found in local-window LD output")
    raw["distance_bp"] = (raw["BP_B"] - raw["BP_A"]).abs()
    if (raw["distance_bp"] > window_bp).any():
        raise AssertionError("PLINK LD output exceeds requested window")
    # PLINK 1.9 emits focal self-pairs (r2=1) in this mode.  They are useful
    # for confirming focal presence but are not LD proxies and are excluded.
    raw = raw.loc[~raw["SNP_A"].eq(raw["SNP_B"])].copy()
    if raw.duplicated(["panel", "SNP_A", "SNP_B"]).any():
        raise AssertionError("duplicate raw PLINK LD pairs")

    focal_map = panel_audit.loc[
        panel_audit["panel_evaluation_status"].eq("matched_polymorphic_in_panel"),
        ["panel", "panel_variant_id", "normalized_variant_id", "MAF", "NCHROBS"],
    ].rename(
        columns={
            "panel_variant_id": "SNP_A",
            "normalized_variant_id": "focal_tag_variant",
            "MAF": "focal_panel_maf",
            "NCHROBS": "focal_panel_n_observed_alleles",
        }
    )
    expanded = raw.merge(
        focal_map,
        on=["panel", "SNP_A"],
        how="left",
        validate="many_to_many",
    )
    if expanded["focal_tag_variant"].isna().any():
        missing = expanded.loc[expanded["focal_tag_variant"].isna(), "SNP_A"].unique()
        raise AssertionError(
            "PLINK emitted focal IDs absent from polymorphic audit: "
            + ", ".join(map(str, missing[:10]))
        )

    parsed_proxy = expanded["SNP_B"].map(parse_variant_id)
    expanded["proxy_chromosome"] = parsed_proxy.map(lambda item: item[0])
    expanded["proxy_position_grch38"] = parsed_proxy.map(lambda item: item[1])
    expanded["proxy_ref"] = parsed_proxy.map(lambda item: item[2])
    expanded["proxy_alt"] = parsed_proxy.map(lambda item: item[3])
    expanded["proxy_variant_class"] = parsed_proxy.map(
        lambda item: variant_class(item[2], item[3])
    )
    parsed_focal = expanded["SNP_A"].map(parse_variant_id)
    expanded["focal_chromosome"] = parsed_focal.map(lambda item: item[0])
    expanded["focal_position_grch38"] = parsed_focal.map(lambda item: item[1])
    expanded["focal_ref"] = parsed_focal.map(lambda item: item[2])
    expanded["focal_alt"] = parsed_focal.map(lambda item: item[3])

    columns = [
        "panel",
        "focal_tag_variant",
        "SNP_A",
        "focal_chromosome",
        "focal_position_grch38",
        "focal_ref",
        "focal_alt",
        "focal_panel_maf",
        "focal_panel_n_observed_alleles",
        "SNP_B",
        "proxy_chromosome",
        "proxy_position_grch38",
        "proxy_ref",
        "proxy_alt",
        "proxy_variant_class",
        "R2",
        "distance_bp",
    ]
    expanded = expanded[columns].rename(
        columns={
            "SNP_A": "focal_panel_variant_id",
            "SNP_B": "proxy_panel_variant_id",
            "R2": "r2",
        }
    )
    expanded.insert(0, "phenotype_scope", "core_copd")
    expanded["reference_release"] = REFERENCE_RELEASE
    return expanded.sort_values(
        [
            "panel",
            "focal_chromosome",
            "focal_position_grch38",
            "focal_tag_variant",
            "proxy_position_grch38",
            "proxy_panel_variant_id",
        ],
        kind="stable",
    ).reset_index(drop=True)


def expansion_by_focal(
    panel_audit: pd.DataFrame, pairs: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    grouped = {
        key: group
        for key, group in pairs.groupby(["panel", "focal_tag_variant"], sort=False)
    }
    for audit_row in panel_audit.itertuples(index=False):
        key = (audit_row.panel, audit_row.normalized_variant_id)
        group = grouped.get(key)
        if group is None or group.empty:
            n_proxy = n_snv = n_indel = 0
            block_start = block_end = block_length = pd.NA
            max_r2 = pd.NA
        else:
            unique = group.drop_duplicates("proxy_panel_variant_id")
            n_proxy = len(unique)
            n_snv = int(unique["proxy_variant_class"].eq("SNV").sum())
            n_indel = int(unique["proxy_variant_class"].eq("indel_or_complex").sum())
            focal_pos = int(unique["focal_position_grch38"].iloc[0])
            block_start = min(focal_pos, int(unique["proxy_position_grch38"].min()))
            block_end = max(focal_pos, int(unique["proxy_position_grch38"].max()))
            block_length = block_end - block_start + 1
            max_r2 = float(unique["r2"].max())
        rows.append(
            {
                "phenotype_scope": "core_copd",
                "panel": audit_row.panel,
                "focal_tag_variant": audit_row.normalized_variant_id,
                "panel_evaluation_status": audit_row.panel_evaluation_status,
                "focal_match_status": audit_row.match_status,
                "focal_panel_variant_id": audit_row.panel_variant_id,
                "focal_panel_maf": audit_row.MAF,
                "n_ld_proxies": n_proxy,
                "n_snv_proxies": n_snv,
                "n_indel_or_complex_proxies": n_indel,
                "max_r2": max_r2,
                "block_start_grch38": block_start,
                "block_end_grch38": block_end,
                "block_length_bp_inclusive": block_length,
                "reference_release": REFERENCE_RELEASE,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["panel", "focal_tag_variant"], kind="stable"
    ).reset_index(drop=True)


def build_candidates(
    match_audit: pd.DataFrame,
    panel_audit: pd.DataFrame,
    pairs: pd.DataFrame,
) -> pd.DataFrame:
    candidate: dict[str, dict[str, object]] = {}

    def new_reference_record(variant_id: str) -> dict[str, object]:
        chrom, pos, ref, alt = parse_variant_id(variant_id)
        return {
            "candidate_record_id": variant_id,
            "panel_variant_id": variant_id,
            "chromosome_grch38": chrom,
            "position_grch38": pos,
            "ref": ref,
            "alt": alt,
            "variant_class": variant_class(ref, alt),
            "gws_tag_ids": set(),
            "tag_match_statuses": set(),
            "tag_supported_panels": set(),
            "tag_polymorphic_panels": set(),
            "source_focal_tags": set(),
            "ld_panels": set(),
            "r2_values": [],
            "is_gws_tag": False,
            "is_ld_proxy": False,
            "reference_matched": True,
        }

    panels_by_tag: dict[str, set[str]] = defaultdict(set)
    polymorphic_by_tag: dict[str, set[str]] = defaultdict(set)
    for row in panel_audit.itertuples(index=False):
        panels_by_tag[row.normalized_variant_id].add(row.panel)
        if row.panel_evaluation_status == "matched_polymorphic_in_panel":
            polymorphic_by_tag[row.normalized_variant_id].add(row.panel)

    for row in match_audit.itertuples(index=False):
        tag = str(row.normalized_variant_id)
        if bool(row.matched):
            key = str(row.panel_variant_id)
            item = candidate.setdefault(key, new_reference_record(key))
        else:
            key = f"unmatched_tag:{tag}"
            if pd.isna(row.chromosome_catalog_grch38):
                chrom = ""
            else:
                chrom_value = str(row.chromosome_catalog_grch38)
                chrom = (
                    str(int(float(chrom_value)))
                    if re.fullmatch(r"\d+(?:\.0)?", chrom_value)
                    else chrom_value
                )
            pos = (
                pd.NA
                if pd.isna(row.position_catalog_grch38)
                else int(float(row.position_catalog_grch38))
            )
            item = {
                "candidate_record_id": key,
                "panel_variant_id": "",
                "chromosome_grch38": chrom,
                "position_grch38": pos,
                "ref": "",
                "alt": "",
                "variant_class": "unresolved_tag",
                "gws_tag_ids": set(),
                "tag_match_statuses": set(),
                "tag_supported_panels": set(),
                "tag_polymorphic_panels": set(),
                "source_focal_tags": set(),
                "ld_panels": set(),
                "r2_values": [],
                "is_gws_tag": False,
                "is_ld_proxy": False,
                "reference_matched": False,
            }
            candidate[key] = item
        item["is_gws_tag"] = True
        item["gws_tag_ids"].add(tag)
        item["tag_match_statuses"].add(str(row.match_status))
        item["tag_supported_panels"].update(panels_by_tag.get(tag, set()))
        item["tag_polymorphic_panels"].update(polymorphic_by_tag.get(tag, set()))

    for row in pairs.itertuples(index=False):
        key = str(row.proxy_panel_variant_id)
        item = candidate.setdefault(key, new_reference_record(key))
        item["is_ld_proxy"] = True
        item["source_focal_tags"].add(str(row.focal_tag_variant))
        item["ld_panels"].add(str(row.panel))
        item["r2_values"].append(float(row.r2))

    rows: list[dict[str, object]] = []
    for item in candidate.values():
        is_tag = bool(item["is_gws_tag"])
        is_proxy = bool(item["is_ld_proxy"])
        if is_tag and is_proxy:
            origin = "gws_tag_and_ld_proxy"
        elif is_tag:
            origin = "gws_tag_only"
        else:
            origin = "ld_proxy_only"
        position = item["position_grch38"]
        bed_eligible = bool(item["chromosome_grch38"]) and not pd.isna(position)
        r2_values = item.pop("r2_values")
        row = {
            **item,
            "gws_tag_ids": ";".join(sorted(item["gws_tag_ids"])),
            "n_gws_tag_ids": len(item["gws_tag_ids"]),
            "tag_match_statuses": ";".join(sorted(item["tag_match_statuses"])),
            "tag_supported_panels": ";".join(
                panel for panel in PANEL_ORDER if panel in item["tag_supported_panels"]
            ),
            "tag_polymorphic_panels": ";".join(
                panel for panel in PANEL_ORDER if panel in item["tag_polymorphic_panels"]
            ),
            "source_focal_tags": ";".join(sorted(item["source_focal_tags"])),
            "n_source_focal_tags": len(item["source_focal_tags"]),
            "ld_panels": ";".join(
                panel for panel in PANEL_ORDER if panel in item["ld_panels"]
            ),
            "n_ld_panels": len(item["ld_panels"]),
            "max_r2_across_links": max(r2_values) if r2_values else pd.NA,
            "candidate_origin": origin,
            "bed_eligible": bed_eligible,
            "reference_release": REFERENCE_RELEASE,
        }
        # Remove the set-valued internal fields replaced above.
        for field in (
            "gws_tag_ids",
            "tag_match_statuses",
            "tag_supported_panels",
            "tag_polymorphic_panels",
            "source_focal_tags",
            "ld_panels",
        ):
            if isinstance(row[field], set):
                raise AssertionError(f"internal set escaped candidate serialization: {field}")
        rows.append(row)

    result = pd.DataFrame(rows)
    result["_chrom_sort"] = result["chromosome_grch38"].map(
        lambda value: int(value) if str(value).isdigit() else 99
    )
    result["_pos_sort"] = pd.to_numeric(
        result["position_grch38"], errors="coerce"
    ).fillna(10**18)
    result = result.sort_values(
        ["_chrom_sort", "_pos_sort", "candidate_record_id"], kind="stable"
    ).drop(columns=["_chrom_sort", "_pos_sort"])
    result.insert(0, "phenotype_scope", "core_copd")
    return result.reset_index(drop=True)


def write_bed(candidates: pd.DataFrame, path: Path) -> int:
    eligible = candidates.loc[candidates["bed_eligible"]].copy()
    rows: list[tuple[str, int, int, str, int, str]] = []
    for row in eligible.itertuples(index=False):
        pos = int(row.position_grch38)
        ref_length = max(1, len(str(row.ref))) if row.ref else 1
        rows.append(
            (
                f"chr{row.chromosome_grch38}",
                pos - 1,
                pos - 1 + ref_length,
                str(row.candidate_record_id),
                0,
                ".",
            )
        )
    with path.open("w", encoding="utf-8") as handle:
        for fields in rows:
            handle.write("\t".join(map(str, fields)) + "\n")
    return len(rows)


def build_summaries(
    match_audit: pd.DataFrame,
    panel_audit: pd.DataFrame,
    raw_pairs: pd.DataFrame,
    raw_pair_count: int,
    pairs: pd.DataFrame,
    expansion: pd.DataFrame,
    candidates: pd.DataFrame,
    bed_rows: int,
    r2_threshold: float,
    window_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    proxy_reference_ids = set(pairs["proxy_panel_variant_id"]) if not pairs.empty else set()
    proxy_classes = (
        pairs.drop_duplicates("proxy_panel_variant_id")["proxy_variant_class"].value_counts()
        if not pairs.empty
        else pd.Series(dtype=int)
    )
    metrics: list[tuple[str, object, str]] = [
        ("gws_tag_variants", len(match_audit), "All unique P <= 5e-8 COPD tags"),
        (
            "focal_panel_assignments",
            len(panel_audit),
            "Assignments supported by initial-stage GWAS Catalog ancestry",
        ),
        (
            "strictly_matched_tag_variants",
            int(match_audit["matched"].sum()),
            "Matched to local GRCh38 reference under documented rules",
        ),
        (
            "unresolved_tag_variants",
            int((~match_audit["matched"].astype(bool)).sum()),
            "Retained in audit/candidate table but not LD-expanded",
        ),
        (
            "exact_coordinate_and_allele_matches",
            int(
                match_audit["match_status"]
                .eq("matched_exact_position_and_reported_allele")
                .sum()
            ),
            "Strict exact matches",
        ),
        (
            "exact_coordinate_no_reported_allele_matches",
            int(
                match_audit["match_status"]
                .eq("matched_exact_position_no_reported_sequence_allele")
                .sum()
            ),
            "Unique record at coordinate; allele unavailable",
        ),
        (
            "adjacent_indel_anchor_matches",
            int(
                match_audit["match_status"]
                .eq("matched_adjacent_indel_anchor_and_reported_allele")
                .sum()
            ),
            "Unique +/-1-bp VCF indel anchor with reported allele",
        ),
        (
            "unique_matched_reference_records",
            match_audit.loc[match_audit["matched"], "panel_variant_id"].nunique(),
            "Multiple Catalog tag labels at one VCF record are collapsed",
        ),
        (
            "polymorphic_focal_panel_assignments",
            int(
                panel_audit["panel_evaluation_status"]
                .eq("matched_polymorphic_in_panel")
                .sum()
            ),
            "Eligible for LD calculation",
        ),
        (
            "monomorphic_focal_panel_assignments",
            int(
                panel_audit["panel_evaluation_status"]
                .eq("matched_monomorphic_in_panel")
                .sum()
            ),
            "Matched but no variation in assigned reference panel",
        ),
        (
            "raw_plink_ld_pairs",
            raw_pair_count,
            "Panel-record rows including focal self-pairs",
        ),
        (
            "raw_plink_focal_self_pairs_excluded",
            int(raw_pairs["SNP_A"].eq(raw_pairs["SNP_B"]).sum()),
            "PLINK r2=1 focal self-pairs are not counted as proxies",
        ),
        (
            "tag_specific_ld_pairs",
            len(pairs),
            "Rows after assigning each focal record to Catalog tag label(s)",
        ),
        (
            "unique_ld_proxy_reference_records",
            len(proxy_reference_ids),
            "Union across focal tags and ancestry panels",
        ),
        (
            "unique_ld_proxy_snvs",
            int(proxy_classes.get("SNV", 0)),
            "Unique proxy records",
        ),
        (
            "unique_ld_proxy_indel_or_complex",
            int(proxy_classes.get("indel_or_complex", 0)),
            "Unique proxy records",
        ),
        (
            "expanded_candidate_records",
            len(candidates),
            "Union of GWS tag records and LD proxies; unresolved tags retained",
        ),
        (
            "candidate_bed_records",
            bed_rows,
            "Candidate records with a valid GRCh38 coordinate",
        ),
        ("ld_r2_threshold", r2_threshold, "Inclusive PLINK threshold"),
        ("ld_window_bp", window_bp, "Absolute focal-proxy distance limit"),
    ]
    summary = pd.DataFrame(metrics, columns=["metric", "value", "definition"])
    summary.insert(0, "phenotype_scope", "core_copd")
    summary["reference_release"] = REFERENCE_RELEASE

    panel_rows: list[dict[str, object]] = []
    for panel in PANEL_ORDER:
        audit = panel_audit.loc[panel_audit["panel"].eq(panel)]
        pair_subset = pairs.loc[pairs["panel"].eq(panel)]
        unique_proxy = pair_subset.drop_duplicates("proxy_panel_variant_id")
        expansion_subset = expansion.loc[
            expansion["panel"].eq(panel)
            & expansion["panel_evaluation_status"].eq(
                "matched_polymorphic_in_panel"
            )
        ]
        panel_rows.append(
            {
                "phenotype_scope": "core_copd",
                "panel": panel,
                "n_assigned_tags": len(audit),
                "n_globally_matched_tags": int(audit["matched"].sum()),
                "n_polymorphic_tag_assignments": int(
                    audit["panel_evaluation_status"]
                    .eq("matched_polymorphic_in_panel")
                    .sum()
                ),
                "n_monomorphic_tag_assignments": int(
                    audit["panel_evaluation_status"]
                    .eq("matched_monomorphic_in_panel")
                    .sum()
                ),
                "n_tag_specific_ld_pairs": len(pair_subset),
                "n_unique_proxy_records": len(unique_proxy),
                "n_unique_proxy_snvs": int(
                    unique_proxy["proxy_variant_class"].eq("SNV").sum()
                ),
                "n_unique_proxy_indel_or_complex": int(
                    unique_proxy["proxy_variant_class"].eq("indel_or_complex").sum()
                ),
                "median_proxies_per_polymorphic_tag": (
                    float(expansion_subset["n_ld_proxies"].median())
                    if not expansion_subset.empty
                    else pd.NA
                ),
                "reference_release": REFERENCE_RELEASE,
            }
        )
    return summary, pd.DataFrame(panel_rows)


def main() -> int:
    args = parse_args()
    logger = setup_logging(args.log)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    focal = pd.read_csv(args.focal, sep="\t", low_memory=False)
    assignments = pd.read_csv(args.assignments, sep="\t", low_memory=False)
    match_audit = pd.read_csv(args.match_audit, sep="\t", low_memory=False)
    job_status = pd.read_csv(args.job_status, sep="\t", low_memory=False)
    require_columns(focal, ["normalized_variant_id"], "focal variants")
    require_columns(
        assignments, ["normalized_variant_id", "panel"], "panel assignments"
    )
    require_columns(
        match_audit,
        [
            "normalized_variant_id",
            "matched",
            "match_status",
            "panel_variant_id",
        ],
        "match audit",
    )
    require_columns(job_status, ["job", "job_dir", "status"], "job status")
    if job_status.empty or not job_status["status"].eq("complete").all():
        raise AssertionError("not every current chromosome-panel LD job is complete")
    if job_status["job"].duplicated().any():
        raise AssertionError("duplicate current LD job names")
    if len(focal) != len(match_audit):
        raise AssertionError("focal and match-audit row counts differ")
    if set(focal["normalized_variant_id"]) != set(match_audit["normalized_variant_id"]):
        raise AssertionError("focal and match-audit variant IDs differ")

    frequencies = read_focal_frequencies(job_status)
    panel_audit = build_panel_audit(assignments, match_audit, frequencies)
    missing_subsets = panel_audit["panel_evaluation_status"].eq(
        "matched_record_missing_from_panel_subset"
    ).sum()
    if missing_subsets:
        raise AssertionError(
            f"{missing_subsets} globally matched focal assignments missing from panel subset"
        )

    raw_pairs, raw_pair_count = read_ld_pairs(job_status)
    pairs = validate_and_expand_pairs(
        raw_pairs, panel_audit, args.r2, args.window_bp
    )
    expansion = expansion_by_focal(panel_audit, pairs)
    candidates = build_candidates(match_audit, panel_audit, pairs)

    panel_audit_path = args.results_dir / "COPD-S2-R006B_focal_panel_audit.tsv"
    pairs_path = args.results_dir / "COPD-S2-R006C_ld_pairs.tsv.gz"
    expansion_path = args.results_dir / "COPD-S2-R006D_ld_expansion_by_focal.tsv"
    candidate_path = args.results_dir / "COPD-S2-R006E_candidate_variants_grch38.tsv.gz"
    bed_path = args.results_dir / "COPD-S2-R006E_candidate_variants_grch38.bed"
    summary_path = args.results_dir / "COPD-S2-R006F_summary.tsv"
    panel_summary_path = args.results_dir / "COPD-S2-R006F_panel_summary.tsv"

    panel_audit.to_csv(panel_audit_path, sep="\t", index=False)
    pairs.to_csv(
        pairs_path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    expansion.to_csv(expansion_path, sep="\t", index=False)
    candidates.to_csv(
        candidate_path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    bed_rows = write_bed(candidates, bed_path)
    summary, panel_summary = build_summaries(
        match_audit,
        panel_audit,
        raw_pairs,
        raw_pair_count,
        pairs,
        expansion,
        candidates,
        bed_rows,
        args.r2,
        args.window_bp,
    )
    summary.to_csv(summary_path, sep="\t", index=False)
    panel_summary.to_csv(panel_summary_path, sep="\t", index=False)

    # Cross-output invariants.
    represented_tags = {
        tag
        for value in candidates.loc[candidates["is_gws_tag"], "gws_tag_ids"]
        for tag in str(value).split(";")
        if tag
    }
    if represented_tags != set(match_audit["normalized_variant_id"].astype(str)):
        raise AssertionError("candidate table does not retain every GWS tag")
    if len(panel_audit) != len(assignments):
        raise AssertionError("focal-panel audit lost assignments")
    if len(expansion) != len(assignments):
        raise AssertionError("per-focal expansion table lost assignments")

    output_rows = {
        panel_audit_path: len(panel_audit),
        pairs_path: len(pairs),
        expansion_path: len(expansion),
        candidate_path: len(candidates),
        bed_path: bed_rows,
        summary_path: len(summary),
        panel_summary_path: len(panel_summary),
    }
    manifest_rows: list[tuple[str, object]] = [
        ("run_utc", datetime.now(timezone.utc).isoformat()),
        ("python", platform.python_version()),
        ("pandas", pd.__version__),
        ("script", str(SCRIPT)),
        ("script_sha256", sha256_file(SCRIPT)),
        ("focal_input_sha256", sha256_file(args.focal)),
        ("assignments_input_sha256", sha256_file(args.assignments)),
        ("match_audit_input_sha256", sha256_file(args.match_audit)),
        ("job_status_input_sha256", sha256_file(args.job_status)),
        ("reference_release", REFERENCE_RELEASE),
        ("r2_threshold", args.r2),
        ("window_bp", args.window_bp),
    ]
    for path, row_count in output_rows.items():
        manifest_rows.extend(
            [
                (f"output_{path.name}_rows", row_count),
                (f"output_{path.name}_sha256", sha256_file(path)),
            ]
        )
    pd.DataFrame(manifest_rows, columns=["field", "value"]).to_csv(
        LOGS / "COPD-S2-ld_summary_manifest.tsv", sep="\t", index=False
    )

    metric_map = dict(zip(summary["metric"], summary["value"]))
    logger.info(
        "Strict focal matches: %s/%s; unresolved retained: %s",
        metric_map["strictly_matched_tag_variants"],
        metric_map["gws_tag_variants"],
        metric_map["unresolved_tag_variants"],
    )
    logger.info(
        "LD pairs: %s raw PLINK, %s tag-specific; unique proxies: %s",
        metric_map["raw_plink_ld_pairs"],
        metric_map["tag_specific_ld_pairs"],
        metric_map["unique_ld_proxy_reference_records"],
    )
    logger.info(
        "Expanded candidates: %s records; BED-eligible: %s",
        metric_map["expanded_candidate_records"],
        metric_map["candidate_bed_records"],
    )
    logger.info("Wrote COPD-S2-R006A through R006F under %s", args.results_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
