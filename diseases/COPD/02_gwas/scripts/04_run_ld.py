#!/usr/bin/env python3
"""Run ancestry-aware COPD LD expansion against local 1000 Genomes GRCh38.

The reference is the 1000 Genomes Phase 3 20190312 biallelic phased SNV/indel
release remapped to GRCh38.  Focal records are matched before LD calculation by
Catalog GRCh38 coordinate and reported allele.  A unique adjacent (plus or
minus one base) indel whose REF or ALT is the reported allele is accepted as a
VCF anchor-normalization match; neighboring SNVs and ambiguous records are not.

For each chromosome and ancestry panel, bcftools subsets the union of 500-kb
focal windows and PLINK 1.9 calculates pairwise r2.  Intermediate VCF/BED,
frequency, and LD files remain under data/ld_work for a complete audit trail.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
REPO_ROOT = SECTION_ROOT.parents[2]
DEFAULT_SHARED_DATA = REPO_ROOT / "data"
WORK = SECTION_ROOT / "data" / "ld_work"
RESULTS = SECTION_ROOT / "results"
LOGS = SECTION_ROOT / "logs"

VCF_STEM = (
    "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_"
    "27022019.GRCh38.phased.vcf.gz"
)
REFERENCE_RELEASE = "1000 Genomes Phase 3, 20190312 GRCh38 remap"
VALID_CHROMS = {str(value) for value in range(1, 23)} | {"X"}
BASES_RE = re.compile(r"^[ACGTN]+$", flags=re.IGNORECASE)


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--focal",
        type=Path,
        default=WORK / "focal_variants.tsv",
    )
    parser.add_argument(
        "--assignments",
        type=Path,
        default=WORK / "focal_panel_assignments.tsv",
    )
    parser.add_argument(
        "--thousand-genomes-dir",
        type=Path,
        default=DEFAULT_SHARED_DATA / "1000genomes",
    )
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--jobs", type=int, default=5)
    parser.add_argument("--window-bp", type=int, default=500_000)
    parser.add_argument("--r2", type=float, default=0.8)
    parser.add_argument(
        "--log", type=Path, default=LOGS / "COPD-S2-ld_expansion.log"
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2-LD")
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


def executable(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(
            f"required executable {command!r} is not on PATH; load bcftools/plink"
        )
    return resolved


def command_version(command: str, arguments: Sequence[str]) -> str:
    result = subprocess.run(
        [command, *arguments], capture_output=True, text=True, check=True
    )
    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0] if text else "unknown"


def chrom_string(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).replace("chr", "").strip()
    if re.fullmatch(r"\d+\.0", text):
        text = str(int(float(text)))
    if text == "23":
        text = "X"
    return text if text in VALID_CHROMS else None


def position_int(value: object) -> int | None:
    if pd.isna(value):
        return None
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0 or not numeric.is_integer():
        return None
    return int(numeric)


def reported_alleles(value: object) -> tuple[set[str], str]:
    """Return literal/underscore-normalized sequence alleles and an audit string."""
    if pd.isna(value):
        return set(), ""
    accepted: set[str] = set()
    raw_items: list[str] = []
    for item in str(value).split(";"):
        allele = item.strip().upper()
        if not allele or allele in {"?", "NR", "N/A", "NAN", "NONE"}:
            continue
        raw_items.append(allele)
        if BASES_RE.fullmatch(allele):
            accepted.add(allele)
        collapsed = allele.replace("_", "")
        if collapsed != allele and BASES_RE.fullmatch(collapsed):
            accepted.add(collapsed)
    return accepted, ";".join(raw_items)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def write_query_regions(
    frame: pd.DataFrame, work_dir: Path
) -> dict[str, Path]:
    query_dir = work_dir / "match_query"
    query_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for chrom, group in frame.groupby("query_chrom", dropna=True):
        intervals = []
        for pos in group["query_pos"].dropna().astype(int).unique():
            # BED is zero-based, half-open; include one base on either side of
            # the one-based Catalog position.
            intervals.append((max(0, pos - 2), pos + 1))
        path = query_dir / f"chr{chrom}.bed"
        with path.open("w", encoding="utf-8") as handle:
            for start, end in merge_intervals(intervals):
                handle.write(f"{chrom}\t{start}\t{end}\n")
        paths[str(chrom)] = path
    return paths


def query_records(
    focal: pd.DataFrame,
    regions: dict[str, Path],
    kg_dir: Path,
    bcftools: str,
) -> dict[str, dict[int, list[VcfRecord]]]:
    records: dict[str, dict[int, list[VcfRecord]]] = {}
    for chrom, region_path in sorted(
        regions.items(), key=lambda item: (item[0] == "X", int(item[0]) if item[0].isdigit() else 99)
    ):
        vcf = kg_dir / VCF_STEM.format(chrom=chrom)
        if not vcf.exists() or not Path(str(vcf) + ".tbi").exists():
            raise FileNotFoundError(f"missing reference VCF or index: {vcf}")
        result = subprocess.run(
            [
                bcftools,
                "query",
                "-R",
                str(region_path),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\n",
                str(vcf),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        by_position: dict[int, list[VcfRecord]] = {}
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 4:
                raise ValueError(f"unexpected bcftools query row: {line!r}")
            record_chrom, pos_text, ref, alt_text = fields
            alts = alt_text.split(",")
            if len(alts) != 1:
                raise AssertionError(
                    f"reference claimed biallelic but has ALT={alt_text}: {line}"
                )
            record = VcfRecord(record_chrom, int(pos_text), ref, alts[0])
            by_position.setdefault(record.pos, []).append(record)
        records[chrom] = by_position
    return records


def validate_reference_samples(
    panels: Iterable[str],
    chromosomes: Iterable[str],
    work_dir: Path,
    kg_dir: Path,
    bcftools: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Intersect panel metadata with VCF headers and audit every exclusion."""
    chromosome_list = sorted(
        set(chromosomes),
        key=lambda value: (value == "X", int(value) if value.isdigit() else 99),
    )
    if not chromosome_list:
        raise ValueError("no valid reference chromosomes were requested")

    reference_samples: list[str] | None = None
    for chrom in chromosome_list:
        vcf = kg_dir / VCF_STEM.format(chrom=chrom)
        result = subprocess.run(
            [bcftools, "query", "-l", str(vcf)],
            capture_output=True,
            text=True,
            check=True,
        )
        samples = [item for item in result.stdout.splitlines() if item]
        if len(samples) != len(set(samples)):
            raise AssertionError(f"duplicate samples in chr{chrom} VCF header")
        if reference_samples is None:
            reference_samples = samples
        elif samples != reference_samples:
            raise AssertionError(
                f"sample order/content differs between reference chromosomes and chr{chrom}"
            )
    assert reference_samples is not None
    reference_set = set(reference_samples)

    audit_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    for panel in sorted(set(panels)):
        requested_path = work_dir / f"samples_{panel}.txt"
        requested = [
            item.strip()
            for item in requested_path.read_text(encoding="utf-8").splitlines()
            if item.strip()
        ]
        if len(requested) != len(set(requested)):
            raise AssertionError(f"duplicate requested sample IDs for {panel}")
        effective = [sample for sample in requested if sample in reference_set]
        missing = [sample for sample in requested if sample not in reference_set]
        if not effective:
            raise ValueError(f"no requested {panel} samples occur in the reference VCF")
        effective_path = work_dir / f"samples_{panel}_reference.txt"
        effective_path.write_text("\n".join(effective) + "\n", encoding="utf-8")
        audit_rows.append(
            {
                "panel": panel,
                "n_panel_metadata_samples": len(requested),
                "n_reference_vcf_samples": len(reference_samples),
                "n_effective_panel_samples": len(effective),
                "n_panel_samples_absent_from_vcf": len(missing),
                "effective_sample_file": str(effective_path.resolve()),
            }
        )
        for sample in missing:
            missing_rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "reason": "listed_in_phase3_panel_metadata_but_absent_from_GRCh38_VCF_header",
                }
            )
    return pd.DataFrame(audit_rows), pd.DataFrame(
        missing_rows, columns=["panel", "sample", "reason"]
    )


def choose_match(
    row: pd.Series, by_position: dict[int, list[VcfRecord]]
) -> dict[str, object]:
    chrom = row["query_chrom"]
    pos = row["query_pos"]
    base = {
        "match_status": "",
        "matched": False,
        "match_basis": "",
        "position_delta_bp": pd.NA,
        "panel_variant_id": "",
        "panel_chromosome": "",
        "panel_position_grch38": pd.NA,
        "panel_ref": "",
        "panel_alt": "",
        "panel_variant_class": "",
        "n_exact_records": 0,
        "n_neighbor_indel_allele_matches": 0,
    }
    if chrom is None or pos is None:
        base["match_status"] = "invalid_or_missing_catalog_coordinate"
        return base

    alleles, reported_text = reported_alleles(row["reported_risk_alleles"])
    exact = list(by_position.get(int(pos), []))
    base["n_exact_records"] = len(exact)

    def allele_match(record: VcfRecord) -> bool:
        return bool(alleles & {record.ref.upper(), record.alt.upper()})

    exact_allele = [record for record in exact if allele_match(record)]
    match: VcfRecord | None = None
    status = ""
    basis = ""

    if alleles and len(exact_allele) == 1:
        match = exact_allele[0]
        status = "matched_exact_position_and_reported_allele"
        basis = "Catalog GRCh38 coordinate plus reported sequence allele"
    elif alleles and len(exact_allele) > 1:
        status = "ambiguous_exact_position_allele_match"
    elif not alleles and len(exact) == 1:
        match = exact[0]
        status = "matched_exact_position_no_reported_sequence_allele"
        basis = "unique record at Catalog GRCh38 coordinate"
    elif not alleles and len(exact) > 1:
        status = "ambiguous_exact_position_without_allele"
    else:
        neighboring = [
            record
            for neighbor_pos in (int(pos) - 1, int(pos) + 1)
            for record in by_position.get(neighbor_pos, [])
            if record.is_indel and allele_match(record)
        ]
        base["n_neighbor_indel_allele_matches"] = len(neighboring)
        if len(neighboring) == 1:
            match = neighboring[0]
            status = "matched_adjacent_indel_anchor_and_reported_allele"
            basis = (
                "unique +/-1-bp VCF indel anchor with reported sequence allele"
            )
        elif len(neighboring) > 1:
            status = "ambiguous_adjacent_indel_anchor_match"
        elif exact:
            status = "exact_position_reported_allele_discordant"
        elif alleles:
            status = "no_position_and_reported_allele_match"
        elif reported_text:
            status = "reported_allele_not_a_resolvable_sequence"
        else:
            status = "no_record_at_catalog_position"

    base["match_status"] = status
    if match is not None:
        base.update(
            {
                "matched": True,
                "match_basis": basis,
                "position_delta_bp": match.pos - int(pos),
                "panel_variant_id": match.variant_id,
                "panel_chromosome": match.chrom,
                "panel_position_grch38": match.pos,
                "panel_ref": match.ref,
                "panel_alt": match.alt,
                "panel_variant_class": (
                    "SNV"
                    if len(match.ref) == 1 and len(match.alt) == 1
                    else "indel_or_complex"
                ),
            }
        )
    return base


def build_match_audit(
    focal: pd.DataFrame,
    records: dict[str, dict[int, list[VcfRecord]]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in focal.iterrows():
        chrom = row["query_chrom"]
        result = choose_match(row, records.get(chrom, {}))
        rows.append(
            {
                "phenotype_scope": "core_copd",
                "normalized_variant_id": row["normalized_variant_id"],
                "identifier_type": row["identifier_type"],
                "chromosome_catalog_grch38": chrom or "",
                "position_catalog_grch38": row["query_pos"],
                "reported_risk_alleles": row["reported_risk_alleles"],
                "gws_study_accessions": row["gws_study_accessions"],
                "supported_panels": row["supported_panels"],
                **result,
                "reference_release": REFERENCE_RELEASE,
            }
        )
    return pd.DataFrame(rows)


def write_ld_regions(
    positions: Iterable[int], chrom: str, window_bp: int, path: Path
) -> None:
    intervals = [
        (max(0, int(pos) - window_bp - 1), int(pos) + window_bp)
        for pos in set(positions)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for start, end in merge_intervals(intervals):
            handle.write(f"{chrom}\t{start}\t{end}\n")


def run_checked(
    command: list[str], stdout_path: Path, stderr_path: Path
) -> None:
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        result = subprocess.run(
            command, stdout=stdout_handle, stderr=stderr_handle, text=True
        )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)


def run_job(
    panel: str,
    chrom: str,
    group: pd.DataFrame,
    args: argparse.Namespace,
    bcftools: str,
    plink: str,
) -> dict[str, object]:
    job_name = f"chr{chrom}_{panel}"
    job_dir = args.work_dir / "jobs" / job_name
    job_dir.mkdir(parents=True, exist_ok=True)
    regions = job_dir / "regions.bed"
    focal_ids = job_dir / "focal_variant_ids.txt"
    subset_vcf = job_dir / "reference_subset.vcf.gz"
    bed_prefix = job_dir / "reference_subset"
    freq_prefix = job_dir / "focal_frequency"
    ld_prefix = job_dir / "ld"
    source_vcf = args.thousand_genomes_dir / VCF_STEM.format(chrom=chrom)
    sample_file = args.work_dir / f"samples_{panel}_reference.txt"

    write_ld_regions(
        group["panel_position_grch38"].astype(int),
        chrom,
        args.window_bp,
        regions,
    )
    unique_ids = sorted(set(group["panel_variant_id"].astype(str)))
    focal_ids.write_text("\n".join(unique_ids) + "\n", encoding="utf-8")

    # Stream uncompressed BCF between bcftools processes to avoid a duplicate
    # intermediate VCF while retaining separate stderr logs.
    view_command = [
        bcftools,
        "view",
        "-R",
        str(regions),
        "-S",
        str(sample_file),
        "-Ou",
        str(source_vcf),
    ]
    annotate_command = [
        bcftools,
        "annotate",
        "--set-id",
        "%CHROM:%POS:%REF:%ALT",
        "-Oz",
        "-o",
        str(subset_vcf),
    ]
    with (job_dir / "bcftools_view.stderr.log").open(
        "w", encoding="utf-8"
    ) as view_err, (job_dir / "bcftools_annotate.stderr.log").open(
        "w", encoding="utf-8"
    ) as annotate_err:
        view_process = subprocess.Popen(
            view_command, stdout=subprocess.PIPE, stderr=view_err
        )
        assert view_process.stdout is not None
        annotate_process = subprocess.Popen(
            annotate_command,
            stdin=view_process.stdout,
            stdout=subprocess.DEVNULL,
            stderr=annotate_err,
        )
        view_process.stdout.close()
        annotate_return = annotate_process.wait()
        view_return = view_process.wait()
    if view_return != 0:
        raise subprocess.CalledProcessError(view_return, view_command)
    if annotate_return != 0:
        raise subprocess.CalledProcessError(annotate_return, annotate_command)

    run_checked(
        [bcftools, "index", "-t", "-f", str(subset_vcf)],
        job_dir / "bcftools_index.stdout.log",
        job_dir / "bcftools_index.stderr.log",
    )
    run_checked(
        [
            plink,
            "--vcf",
            str(subset_vcf),
            "--vcf-half-call",
            "missing",
            "--keep-allele-order",
            "--make-bed",
            "--out",
            str(bed_prefix),
        ],
        job_dir / "plink_make_bed.stdout.log",
        job_dir / "plink_make_bed.stderr.log",
    )
    run_checked(
        [
            plink,
            "--bfile",
            str(bed_prefix),
            "--keep-allele-order",
            "--extract",
            str(focal_ids),
            "--freq",
            "--out",
            str(freq_prefix),
        ],
        job_dir / "plink_frequency.stdout.log",
        job_dir / "plink_frequency.stderr.log",
    )
    run_checked(
        [
            plink,
            "--bfile",
            str(bed_prefix),
            "--keep-allele-order",
            "--r2",
            "--ld-snp-list",
            str(focal_ids),
            "--ld-window-kb",
            str(args.window_bp // 1000),
            "--ld-window",
            "999999",
            "--ld-window-r2",
            str(args.r2),
            "--out",
            str(ld_prefix),
        ],
        job_dir / "plink_ld.stdout.log",
        job_dir / "plink_ld.stderr.log",
    )

    bim_rows = sum(1 for _ in (Path(str(bed_prefix) + ".bim")).open())
    ld_path = Path(str(ld_prefix) + ".ld")
    ld_rows = max(0, sum(1 for _ in ld_path.open()) - 1) if ld_path.exists() else 0
    return {
        "job": job_name,
        "panel": panel,
        "chromosome": chrom,
        "status": "complete",
        "n_focal_tags": group["normalized_variant_id"].nunique(),
        "n_unique_focal_panel_records": len(unique_ids),
        "n_reference_variants": bim_rows,
        "n_ld_pairs": ld_rows,
        "job_dir": str(job_dir.resolve()),
        "error": "",
    }


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.window_bp < 1:
        raise ValueError("--window-bp must be positive")
    if not (0 < args.r2 <= 1):
        raise ValueError("--r2 must be in (0, 1]")
    logger = setup_logging(args.log)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    bcftools = executable(args.bcftools)
    plink = executable(args.plink)
    bcftools_version = command_version(bcftools, ["--version"])
    plink_version = command_version(plink, ["--version"])
    logger.info("Reference: %s", REFERENCE_RELEASE)
    logger.info("bcftools: %s", bcftools_version)
    logger.info("PLINK: %s", plink_version)

    focal = pd.read_csv(args.focal, sep="\t", low_memory=False)
    assignments = pd.read_csv(args.assignments, sep="\t", low_memory=False)
    require_columns(
        focal,
        [
            "normalized_variant_id",
            "identifier_type",
            "chromosome",
            "position_grch38_catalog",
            "reported_risk_alleles",
            "gws_study_accessions",
            "supported_panels",
        ],
        "focal variants",
    )
    require_columns(
        assignments, ["normalized_variant_id", "panel"], "panel assignments"
    )
    if focal["normalized_variant_id"].duplicated().any():
        raise AssertionError("focal variants are not unique")
    if assignments.duplicated(["normalized_variant_id", "panel"]).any():
        raise AssertionError("focal-panel assignments are not unique")

    focal["query_chrom"] = focal["chromosome"].map(chrom_string)
    focal["query_pos"] = focal["position_grch38_catalog"].map(position_int)
    query_regions = write_query_regions(focal, args.work_dir)
    sample_audit, missing_samples = validate_reference_samples(
        assignments["panel"].astype(str),
        query_regions,
        args.work_dir,
        args.thousand_genomes_dir,
        bcftools,
    )
    sample_audit.to_csv(
        args.work_dir / "reference_sample_audit.tsv", sep="\t", index=False
    )
    missing_samples.to_csv(
        args.work_dir / "reference_sample_exclusions.tsv", sep="\t", index=False
    )
    for row in sample_audit.itertuples(index=False):
        logger.info(
            "%s reference samples: %d/%d retained; %d panel-metadata IDs absent",
            row.panel,
            row.n_effective_panel_samples,
            row.n_panel_metadata_samples,
            row.n_panel_samples_absent_from_vcf,
        )
    records = query_records(
        focal, query_regions, args.thousand_genomes_dir, bcftools
    )
    audit = build_match_audit(focal, records)
    audit_path = args.results_dir / "COPD-S2-R006A_focal_match_audit.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)

    status_counts = audit["match_status"].value_counts()
    logger.info("Focal matching: %d/%d matched", int(audit["matched"].sum()), len(audit))
    for status, count in status_counts.items():
        logger.info("Match status %s: %d", status, count)

    matched_assignments = assignments.merge(
        audit.loc[audit["matched"]],
        on="normalized_variant_id",
        how="inner",
        validate="many_to_one",
    )
    jobs: list[tuple[str, str, pd.DataFrame]] = []
    for (panel, chrom), group in matched_assignments.groupby(
        ["panel", "panel_chromosome"], sort=True
    ):
        jobs.append((str(panel), str(chrom), group.copy()))
    logger.info(
        "Running %d chromosome-panel jobs with concurrency %d", len(jobs), args.jobs
    )

    job_results: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {
            executor.submit(
                run_job, panel, chrom, group, args, bcftools, plink
            ): (panel, chrom, group)
            for panel, chrom, group in jobs
        }
        for future in as_completed(future_to_job):
            panel, chrom, group = future_to_job[future]
            job_name = f"chr{chrom}_{panel}"
            try:
                result = future.result()
                job_results.append(result)
                logger.info(
                    "%s complete: %d focal tags, %d reference variants, %d LD pairs",
                    job_name,
                    result["n_focal_tags"],
                    result["n_reference_variants"],
                    result["n_ld_pairs"],
                )
            except Exception as error:  # preserve every completed job for audit
                job_results.append(
                    {
                        "job": job_name,
                        "panel": panel,
                        "chromosome": chrom,
                        "status": "failed",
                        "n_focal_tags": group["normalized_variant_id"].nunique(),
                        "n_unique_focal_panel_records": group[
                            "panel_variant_id"
                        ].nunique(),
                        "n_reference_variants": pd.NA,
                        "n_ld_pairs": pd.NA,
                        "job_dir": str(
                            (args.work_dir / "jobs" / job_name).resolve()
                        ),
                        "error": repr(error),
                    }
                )
                logger.error("%s failed: %r", job_name, error)

    job_frame = pd.DataFrame(job_results).sort_values(
        ["panel", "chromosome"], kind="stable"
    )
    job_frame.to_csv(args.work_dir / "job_status.tsv", sep="\t", index=False)

    manifest_rows: list[tuple[str, object]] = [
        ("run_utc", datetime.now(timezone.utc).isoformat()),
        ("python", platform.python_version()),
        ("pandas", pd.__version__),
        ("script", str(SCRIPT)),
        ("script_sha256", sha256_file(SCRIPT)),
        ("focal_input", str(args.focal.resolve())),
        ("focal_input_sha256", sha256_file(args.focal)),
        ("assignments_input", str(args.assignments.resolve())),
        ("assignments_input_sha256", sha256_file(args.assignments)),
        ("reference_directory", str(args.thousand_genomes_dir.resolve())),
        ("reference_release", REFERENCE_RELEASE),
        ("reference_build", "GRCh38"),
        ("reference_vcf_pattern", VCF_STEM),
        ("bcftools", bcftools_version),
        ("plink", plink_version),
        ("window_bp", args.window_bp),
        ("r2_threshold", args.r2),
        ("parallel_jobs", args.jobs),
        ("n_focal_tags", len(focal)),
        ("n_strictly_matched_focal_tags", int(audit["matched"].sum())),
        ("n_focal_panel_assignments", len(assignments)),
        ("n_matched_focal_panel_assignments", len(matched_assignments)),
        ("n_chromosome_panel_jobs", len(jobs)),
        ("n_failed_jobs", int(job_frame["status"].eq("failed").sum())),
        (
            "n_panel_metadata_samples_absent_from_reference_vcf",
            len(missing_samples),
        ),
    ]
    for row in sample_audit.itertuples(index=False):
        manifest_rows.append(
            (f"{row.panel}_effective_reference_samples", row.n_effective_panel_samples)
        )
    for chrom in sorted(
        set(matched_assignments["panel_chromosome"]),
        key=lambda value: (value == "X", int(value) if str(value).isdigit() else 99),
    ):
        vcf = args.thousand_genomes_dir / VCF_STEM.format(chrom=chrom)
        manifest_rows.extend(
            [
                (f"reference_chr{chrom}_path", str(vcf.resolve())),
                (f"reference_chr{chrom}_size_bytes", vcf.stat().st_size),
                (
                    f"reference_chr{chrom}_index_size_bytes",
                    Path(str(vcf) + ".tbi").stat().st_size,
                ),
            ]
        )
    pd.DataFrame(manifest_rows, columns=["field", "value"]).to_csv(
        LOGS / "COPD-S2-ld_expansion_manifest.tsv", sep="\t", index=False
    )

    failed = job_frame.loc[job_frame["status"].eq("failed")]
    if not failed.empty:
        raise RuntimeError(
            f"{len(failed)} LD jobs failed; inspect data/ld_work/job_status.tsv"
        )
    logger.info("All %d chromosome-panel jobs completed", len(jobs))
    logger.info("Focal match audit: %s", audit_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
