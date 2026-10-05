#!/usr/bin/env python3
"""Build the portable, offline COPD Section 2 core-GWAS result tables.

Only rows explicitly labelled ``core_copd`` by step 01 are analysed.  The
script uses the local GWAS Catalog extracts, GENCODE v50, and (when present)
the local 1000 Genomes GRCh38 VCFs.  It performs no network requests and does
not perform LD expansion.

The stable result identifiers produced here are COPD-S2-R001 through R005.
All paths are derived from this script's location unless explicitly
overridden on the command line.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import logging
import math
import platform
import re
import shlex
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
REPO_ROOT = SECTION_ROOT.parents[2]
DEFAULT_SHARED_DATA = REPO_ROOT / "data"
RESULTS = SECTION_ROOT / "results"
LOGS = SECTION_ROOT / "logs"

GWS_THRESHOLD = 5e-8
CORE_SCOPE = "core_copd"

CODING_CONSEQUENCES = {
    "coding_sequence_variant",
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "protein_altering_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "synonymous_variant",
}
SPLICE_CONSEQUENCES = {
    "splice_acceptor_variant",
    "splice_donor_region_variant",
    "splice_donor_variant",
    "splice_region_variant",
}

COUNT_RE = re.compile(
    r"(?<![\d,])(?P<number>\d{1,3}(?:,\d{3})+|\d+)\s+"
    r"(?P<label>[^,;\d]{0,120}?)\b(?P<kind>cases?|controls?)\b",
    flags=re.IGNORECASE,
)
RSID_RE = re.compile(r"rs0*(\d+)", flags=re.IGNORECASE)
COORD_RE = re.compile(
    r"chr:?([0-9]+|X|Y|M|MT):(\d+)", flags=re.IGNORECASE
)
PRIMARY_CHROM_RE = re.compile(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--associations",
        type=Path,
        default=SECTION_ROOT / "data" / "copd_associations.tsv",
    )
    parser.add_argument(
        "--studies",
        type=Path,
        default=SECTION_ROOT / "data" / "copd_studies.tsv",
    )
    parser.add_argument(
        "--ancestries",
        type=Path,
        default=SECTION_ROOT / "data" / "copd_ancestries.tsv",
    )
    parser.add_argument(
        "--gencode",
        type=Path,
        default=DEFAULT_SHARED_DATA / "gencode" / "gencode.v50.annotation.gtf.gz",
    )
    parser.add_argument(
        "--thousand-genomes-dir",
        type=Path,
        default=DEFAULT_SHARED_DATA / "1000genomes",
        help="Local 1000 Genomes GRCh38 VCF directory; no files are downloaded.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Number of variants retained per inferred effect-measure group.",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=LOGS / "COPD-S2-core_analysis.log",
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2")
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


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def read_core_table(path: Path, required: Iterable[str], label: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    require_columns(frame, [*required, "copd_category"], label)
    core = frame.loc[frame["copd_category"].eq(CORE_SCOPE)].copy()
    if core.empty:
        raise ValueError(f"{label} contains no {CORE_SCOPE!r} rows")
    if not core["copd_category"].eq(CORE_SCOPE).all():
        raise AssertionError(f"non-core phenotype leaked into {label}")
    return core


def clean_number(value: Any) -> float:
    if pd.isna(value):
        return math.nan
    return pd.to_numeric(str(value).strip(), errors="coerce")


def parse_case_control(text: Any) -> tuple[int, int, int]:
    """Return cases, controls, and number of recognized count phrases."""
    if pd.isna(text):
        return 0, 0, 0
    cases = controls = matches = 0
    for match in COUNT_RE.finditer(str(text)):
        number = int(match.group("number").replace(",", ""))
        if match.group("kind").lower().startswith("case"):
            cases += number
        else:
            controls += number
        matches += 1
    return cases, controls, matches


def text_has(text: Any, pattern: str) -> bool:
    return bool(re.search(pattern, "" if pd.isna(text) else str(text), re.I))


def build_study_outputs(
    studies: pd.DataFrame, ancestries: pd.DataFrame, associations: pd.DataFrame
) -> dict[str, pd.DataFrame]:
    studies = studies.copy()
    ancestries = ancestries.copy()

    for column in ("NUMBER OF INDIVIDUALS", "NUMBER OF CASES", "NUMBER OF CONTROLS"):
        ancestries[column] = pd.to_numeric(ancestries[column], errors="coerce")

    ancestry_summary = (
        ancestries.groupby(["STAGE", "BROAD ANCESTRAL CATEGORY"], dropna=False)
        .agg(
            n_cohort_entries=("STUDY ACCESSION", "size"),
            n_studies=("STUDY ACCESSION", "nunique"),
            n_entries_with_individual_count=("NUMBER OF INDIVIDUALS", "count"),
            n_individuals_summed=("NUMBER OF INDIVIDUALS", "sum"),
            median_individuals_per_entry=("NUMBER OF INDIVIDUALS", "median"),
        )
        .reset_index()
        .sort_values(["STAGE", "n_individuals_summed", "BROAD ANCESTRAL CATEGORY"],
                     ascending=[True, False, True])
    )
    ancestry_summary.insert(0, "phenotype_scope", CORE_SCOPE)
    ancestry_summary["interpretation_note"] = (
        "Catalog cohort-entry sum; cohort overlap across studies was not removed"
    )

    for stage, field in (
        ("initial", "INITIAL SAMPLE SIZE"),
        ("replication", "REPLICATION SAMPLE SIZE"),
    ):
        parsed = studies[field].map(parse_case_control)
        studies[f"{stage}_cases_parsed"] = parsed.map(lambda item: item[0])
        studies[f"{stage}_controls_parsed"] = parsed.map(lambda item: item[1])
        studies[f"{stage}_count_phrases_parsed"] = parsed.map(lambda item: item[2])
        studies[f"{stage}_sex_mentioned"] = studies[field].map(
            lambda value: text_has(value, r"\b(?:male|female|men|women)s?\b")
        )
        studies[f"{stage}_age_mentioned"] = studies[field].map(
            lambda value: text_has(value, r"\bage(?:d)?\b")
        )

    by_study_stage = (
        ancestries.groupby(["STUDY ACCESSION", "STAGE"], dropna=False)
        ["NUMBER OF INDIVIDUALS"]
        .sum(min_count=1)
        .unstack("STAGE")
    )
    for stage in ("initial", "replication"):
        if stage in by_study_stage:
            mapping = by_study_stage[stage]
            studies[f"{stage}_individuals_catalog"] = studies["STUDY ACCESSION"].map(mapping)
        else:
            studies[f"{stage}_individuals_catalog"] = math.nan

    studies = studies.sort_values(["DATE", "STUDY ACCESSION"]).reset_index(drop=True)

    sample_rows: list[dict[str, Any]] = []
    for stage, field in (
        ("initial", "INITIAL SAMPLE SIZE"),
        ("replication", "REPLICATION SAMPLE SIZE"),
    ):
        anc_stage = ancestries.loc[ancestries["STAGE"].eq(stage)]
        sample_rows.append(
            {
                "phenotype_scope": CORE_SCOPE,
                "stage": stage,
                "n_studies_with_sample_text": int(studies[field].notna().sum()),
                "n_ancestry_cohort_entries": len(anc_stage),
                "n_studies_in_ancestry_table": anc_stage["STUDY ACCESSION"].nunique(),
                "n_entries_with_individual_count": int(
                    anc_stage["NUMBER OF INDIVIDUALS"].notna().sum()
                ),
                "n_individuals_summed": anc_stage["NUMBER OF INDIVIDUALS"].sum(min_count=1),
                "parsed_cases_summed": int(studies[f"{stage}_cases_parsed"].sum()),
                "parsed_controls_summed": int(studies[f"{stage}_controls_parsed"].sum()),
                "n_studies_with_parsed_case_or_control": int(
                    studies[f"{stage}_count_phrases_parsed"].gt(0).sum()
                ),
                "n_studies_with_sex_mention": int(
                    studies[f"{stage}_sex_mentioned"].sum()
                ),
                "n_studies_with_age_mention": int(
                    studies[f"{stage}_age_mentioned"].sum()
                ),
                "interpretation_note": (
                    "Sums are descriptive Catalog totals; cohorts overlap between studies "
                    "and are not counts of unique people"
                ),
            }
        )
    sample_summary = pd.DataFrame(sample_rows)

    full_stats_yes = studies["FULL SUMMARY STATISTICS"].fillna("").str.lower().eq("yes")
    study_metrics = [
        ("core_study_accessions", studies["STUDY ACCESSION"].nunique(),
         "All study records mapped directly to COPD"),
        ("core_publications", studies["PUBMED ID"].nunique(),
         "Unique PubMed IDs among core studies"),
        ("studies_with_curated_association_rows",
         associations["STUDY ACCESSION"].nunique(),
         "Core studies represented in the extracted association table"),
        ("core_association_rows", len(associations),
         "Rows in copd_associations.tsv after enforcing core_copd"),
        ("studies_with_full_summary_statistics", int(full_stats_yes.sum()),
         "GWAS Catalog FULL SUMMARY STATISTICS equals yes"),
        ("ancestry_cohort_entries", len(ancestries),
         "Initial and replication ancestry rows"),
        ("study_year_min", int(studies["DATE"].astype(str).str[:4].min()),
         "Earliest study publication year"),
        ("study_year_max", int(studies["DATE"].astype(str).str[:4].max()),
         "Latest study publication year"),
    ]
    study_summary = pd.DataFrame(
        study_metrics, columns=["metric", "value", "definition"]
    )
    study_summary.insert(0, "phenotype_scope", CORE_SCOPE)

    return {
        "studies": studies,
        "study_summary": study_summary,
        "ancestry_summary": ancestry_summary,
        "sample_summary": sample_summary,
    }


def normalize_variant_id(value: Any) -> str:
    raw = "" if pd.isna(value) else str(value).strip()
    rs_match = RSID_RE.fullmatch(raw)
    if rs_match:
        return f"rs{int(rs_match.group(1))}"
    coord_match = COORD_RE.fullmatch(raw)
    if coord_match:
        chrom = coord_match.group(1).upper()
        chrom = "MT" if chrom == "M" else chrom
        return f"chr{chrom}:{int(coord_match.group(2))}"
    return re.sub(r"\s+", "", raw)


def identifier_type(variant_id: str) -> str:
    if RSID_RE.fullmatch(variant_id):
        return "rsID"
    if COORD_RE.fullmatch(variant_id):
        return "genomic_coordinate"
    if re.fullmatch(r"(?:[0-9]+|X|Y)[pq]\d+(?:\.\d+)?", variant_id, re.I):
        return "cytogenetic_region"
    if re.search(r";|\bx\b|,", variant_id, re.I):
        return "compound_marker"
    return "other_identifier"


def normalize_chromosome(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().removeprefix("chr")
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    text = text.upper()
    if text == "23":
        return "X"
    if text == "24":
        return "Y"
    if text == "M":
        return "MT"
    return text or None


def risk_allele(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if "-" not in text:
        return None
    allele = text.rsplit("-", 1)[-1].strip().upper()
    return None if allele in {"", "?", "NR", "NA"} else allele


def genes_of(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NR", "NAN"}:
        return []
    genes = [
        part.strip()
        for part in re.split(r"\s+-\s+|\s*[,;]\s*", text)
        if part.strip() and part.strip().upper() not in {"NA", "NR", "NAN"}
    ]
    return list(dict.fromkeys(genes))


def consequence_class(value: Any) -> str:
    if pd.isna(value) or not str(value).strip():
        return "unknown"
    terms = {
        term.strip()
        for term in re.split(r"[;,]", str(value))
        if term.strip()
    }
    if terms & CODING_CONSEQUENCES:
        return "coding"
    if terms & SPLICE_CONSEQUENCES:
        return "splice"
    return "noncoding"


def chromosome_sort_key(chrom: Any) -> tuple[int, str]:
    if pd.isna(chrom):
        return (10_000, "")
    value = str(chrom).removeprefix("chr")
    if value.isdigit():
        return (int(value), value)
    return ({"X": 23, "Y": 24, "MT": 25, "M": 25}.get(value, 999), value)


def build_unique_variants(associations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    assoc = associations.copy()
    assoc["normalized_variant_id"] = assoc["SNPS"].map(normalize_variant_id)
    assoc["identifier_type"] = assoc["normalized_variant_id"].map(identifier_type)
    assoc["P"] = pd.to_numeric(assoc["P-VALUE"], errors="coerce")
    if assoc["P"].isna().any():
        raise ValueError("one or more core association rows have a non-numeric P-VALUE")
    assoc["genome_wide_significant"] = assoc["P"].le(GWS_THRESHOLD)
    assoc["consequence_class"] = assoc["CONTEXT"].map(consequence_class)

    rows: list[dict[str, Any]] = []
    for variant_id, group in assoc.groupby("normalized_variant_id", sort=True):
        ranked = group.sort_values(["P", "STUDY ACCESSION"], kind="stable")
        representative = ranked.iloc[0]
        chrom_values = [normalize_chromosome(value) for value in group["CHR_ID"]]
        chrom = next((value for value in chrom_values if value), None)
        positions = [int(value) for value in group["CHR_POS"] if pd.notna(value)]
        position = positions[0] if positions else None
        coordinate = COORD_RE.fullmatch(variant_id)
        if coordinate:
            chrom = normalize_chromosome(coordinate.group(1))
            position = int(coordinate.group(2))
        contexts = sorted(
            {str(value).strip() for value in group["CONTEXT"] if pd.notna(value)}
        )
        mapped_genes = sorted(
            {gene for value in group["MAPPED_GENE"] for gene in genes_of(value)}
        )
        risk_alleles = sorted(
            {
                allele
                for allele in group["STRONGEST SNP-RISK ALLELE"].map(risk_allele)
                if allele
            }
        )
        rows.append(
            {
                "phenotype_scope": CORE_SCOPE,
                "normalized_variant_id": variant_id,
                "identifier_type": identifier_type(variant_id),
                "source_snps_values": ";".join(sorted(set(group["SNPS"].astype(str)))),
                "chromosome": chrom,
                "position_grch38_catalog": position,
                "reported_risk_alleles": ";".join(risk_alleles),
                "n_association_rows": len(group),
                "n_studies": group["STUDY ACCESSION"].nunique(),
                "n_publications": group["PUBMEDID"].nunique(),
                "study_accessions": ";".join(sorted(set(group["STUDY ACCESSION"].astype(str)))),
                "best_p": group["P"].min(),
                "genome_wide_significant": bool(group["genome_wide_significant"].any()),
                "catalog_consequences": ";".join(contexts),
                "consequence_class": consequence_class(contexts[0] if contexts else None),
                "mapped_genes": ";".join(mapped_genes),
                "representative_strongest_risk_allele": representative[
                    "STRONGEST SNP-RISK ALLELE"
                ],
            }
        )
    variants = pd.DataFrame(rows)
    variants["_chrom_sort"] = variants["chromosome"].map(chromosome_sort_key)
    variants = (
        variants.sort_values(
            ["_chrom_sort", "position_grch38_catalog", "normalized_variant_id"],
            kind="stable",
        )
        .drop(columns="_chrom_sort")
        .reset_index(drop=True)
    )
    return assoc, variants


def sequence_variant_type(ref: str, alts: tuple[str, ...]) -> str:
    if not alts:
        return "unresolved"
    if any(
        alt == "*" or alt.startswith("<") or "[" in alt or "]" in alt
        for alt in alts
    ):
        return "structural_or_symbolic"
    if len(ref) == 1 and all(len(alt) == 1 for alt in alts):
        return "SNV"
    if len(ref) > 1 and all(len(alt) == len(ref) for alt in alts):
        return "MNV"
    if any(len(alt) != len(ref) for alt in alts):
        return "indel"
    return "other_sequence_variant"


def annotate_with_local_1000g(
    variants: pd.DataFrame, directory: Path, logger: logging.Logger
) -> pd.DataFrame:
    variants = variants.copy()
    try:
        import pysam
    except ImportError:
        logger.warning("pysam is unavailable; local 1000 Genomes annotation skipped")
        variants["variant_type_1000g"] = "unresolved"
        variants["local_reference_status"] = "pysam_unavailable"
        variants["local_reference_ref"] = ""
        variants["local_reference_alt"] = ""
        variants["local_reference_n_records"] = 0
        return variants

    handles: dict[str, Any] = {}
    annotations: list[dict[str, Any]] = []
    try:
        for row in variants.itertuples(index=False):
            variant_id = row.normalized_variant_id
            chrom = row.chromosome
            position = row.position_grch38_catalog
            if row.identifier_type == "cytogenetic_region":
                annotations.append(
                    {
                        "variant_type_1000g": "region",
                        "local_reference_status": "not_a_single_sequence_variant",
                        "local_reference_ref": "",
                        "local_reference_alt": "",
                        "local_reference_n_records": 0,
                    }
                )
                continue
            if pd.isna(chrom) or pd.isna(position):
                annotations.append(
                    {
                        "variant_type_1000g": "unresolved",
                        "local_reference_status": "catalog_coordinate_missing",
                        "local_reference_ref": "",
                        "local_reference_alt": "",
                        "local_reference_n_records": 0,
                    }
                )
                continue
            vcf_chrom = "X" if str(chrom) == "23" else str(chrom)
            path = directory / (
                f"ALL.chr{vcf_chrom}.shapeit2_integrated_snvindels_v2a_"
                "27022019.GRCh38.phased.vcf.gz"
            )
            if not path.exists():
                annotations.append(
                    {
                        "variant_type_1000g": "unresolved",
                        "local_reference_status": "local_vcf_missing",
                        "local_reference_ref": "",
                        "local_reference_alt": "",
                        "local_reference_n_records": 0,
                    }
                )
                continue
            if vcf_chrom not in handles:
                handles[vcf_chrom] = pysam.VariantFile(str(path))
            pos = int(position)
            records = [
                record
                for record in handles[vcf_chrom].fetch(vcf_chrom, pos - 1, pos)
                if record.pos == pos
            ]
            if not records:
                annotations.append(
                    {
                        "variant_type_1000g": "unresolved",
                        "local_reference_status": "not_found_at_catalog_position",
                        "local_reference_ref": "",
                        "local_reference_alt": "",
                        "local_reference_n_records": 0,
                    }
                )
                continue

            reported = {
                allele.upper()
                for allele in str(row.reported_risk_alleles).split(";")
                if allele and allele.lower() != "nan"
            }
            allele_matches = [
                record
                for record in records
                if reported
                and reported
                & {record.ref.upper(), *(alt.upper() for alt in (record.alts or ())) }
            ]
            record = (allele_matches or records)[0]
            alts = tuple(record.alts or ())
            status = (
                "matched_position_and_reported_allele"
                if allele_matches
                else "matched_position"
            )
            if len(records) > 1:
                status += "_multiple_records"
            annotations.append(
                {
                    "variant_type_1000g": sequence_variant_type(record.ref, alts),
                    "local_reference_status": status,
                    "local_reference_ref": record.ref,
                    "local_reference_alt": ",".join(alts),
                    "local_reference_n_records": len(records),
                }
            )
    finally:
        for handle in handles.values():
            handle.close()

    annotation_frame = pd.DataFrame(annotations, index=variants.index)
    return pd.concat([variants, annotation_frame], axis=1)


def make_count_table(
    variants: pd.DataFrame, column: str, label: str
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for subset_name, subset in (
        ("all_unique_tag_variants", variants),
        ("gws_unique_tag_variants", variants.loc[variants["genome_wide_significant"]]),
    ):
        counts = subset[column].fillna("unknown").replace("", "unknown").value_counts()
        for value, count in counts.items():
            rows.append(
                {
                    "phenotype_scope": CORE_SCOPE,
                    "subset": subset_name,
                    "classification": label,
                    "category": value,
                    "n_variants": int(count),
                    "fraction_of_subset": count / len(subset),
                }
            )
    return pd.DataFrame(rows)


def consequence_counts(variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for subset_name, subset in (
        ("all_unique_tag_variants", variants),
        ("gws_unique_tag_variants", variants.loc[variants["genome_wide_significant"]]),
    ):
        terms = subset["catalog_consequences"].replace("", "unannotated")
        counts = terms.value_counts()
        for term, count in counts.items():
            rows.append(
                {
                    "phenotype_scope": CORE_SCOPE,
                    "subset": subset_name,
                    "catalog_consequence": term,
                    "consequence_class": consequence_class(
                        None if term == "unannotated" else term
                    ),
                    "n_variants": int(count),
                    "fraction_of_subset": count / len(subset),
                }
            )
    return pd.DataFrame(rows)


def coding_summary(variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for subset_name, subset in (
        ("all_unique_tag_variants", variants),
        ("gws_unique_tag_variants", variants.loc[variants["genome_wide_significant"]]),
    ):
        counts = subset["consequence_class"].value_counts()
        annotated = len(subset) - int(counts.get("unknown", 0))
        rows.append(
            {
                "phenotype_scope": CORE_SCOPE,
                "subset": subset_name,
                "n_unique_variants": len(subset),
                "n_with_catalog_consequence": annotated,
                "n_coding": int(counts.get("coding", 0)),
                "n_splice": int(counts.get("splice", 0)),
                "n_noncoding": int(counts.get("noncoding", 0)),
                "n_unknown": int(counts.get("unknown", 0)),
                "coding_fraction_of_annotated": (
                    counts.get("coding", 0) / annotated if annotated else math.nan
                ),
                "splice_fraction_of_annotated": (
                    counts.get("splice", 0) / annotated if annotated else math.nan
                ),
                "noncoding_fraction_of_annotated": (
                    counts.get("noncoding", 0) / annotated if annotated else math.nan
                ),
                "annotation_coverage": annotated / len(subset) if len(subset) else math.nan,
            }
        )
    return pd.DataFrame(rows)


def build_gene_outputs(gws_associations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expanded = gws_associations.copy().reset_index(drop=True)
    expanded["association_row_id"] = expanded.index
    expanded["gene"] = expanded["MAPPED_GENE"].map(genes_of)
    expanded = expanded.explode("gene").dropna(subset=["gene"])
    expanded = expanded.drop_duplicates(["association_row_id", "gene"])

    rows: list[dict[str, Any]] = []
    for gene, group in expanded.groupby("gene", sort=True):
        rows.append(
            {
                "phenotype_scope": CORE_SCOPE,
                "gene": gene,
                "n_studies": group["STUDY ACCESSION"].nunique(),
                "n_publications": group["PUBMEDID"].nunique(),
                "n_association_rows": group["association_row_id"].nunique(),
                "n_unique_tag_variants": group["normalized_variant_id"].nunique(),
                "best_p": group["P"].min(),
                "variants": ";".join(sorted(set(group["normalized_variant_id"]))),
                "study_accessions": ";".join(sorted(set(group["STUDY ACCESSION"].astype(str)))),
                "replicated_across_studies": group["STUDY ACCESSION"].nunique() >= 2,
                "replicated_across_publications": group["PUBMEDID"].nunique() >= 2,
            }
        )
    all_genes = pd.DataFrame(rows).sort_values(
        ["n_studies", "n_publications", "best_p", "gene"],
        ascending=[False, False, True, True],
        kind="stable",
    )
    replicated = all_genes.loc[all_genes["replicated_across_studies"]].reset_index(drop=True)
    return all_genes.reset_index(drop=True), replicated


def effect_measure(ci_text: Any) -> str:
    text = "" if pd.isna(ci_text) else str(ci_text)
    if re.search(
        r"unit|increase|decrease|\bbeta\b|standard deviation|z[ -]?score|year|month",
        text,
        re.I,
    ):
        return "beta_or_other_unit_effect"
    return "odds_ratio_inferred"


def build_effect_outputs(
    gws_associations: pd.DataFrame, top_n: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = gws_associations.copy()
    effects["effect_value"] = pd.to_numeric(effects["OR or BETA"], errors="coerce")
    effects = effects.loc[effects["effect_value"].notna()].copy()
    effects["effect_measure_inference"] = effects["95% CI (TEXT)"].map(effect_measure)
    effects["ranking_value"] = effects["effect_value"].abs()
    odds_mask = effects["effect_measure_inference"].eq("odds_ratio_inferred")
    valid_odds = odds_mask & effects["effect_value"].gt(0)
    effects.loc[valid_odds, "ranking_value"] = effects.loc[
        valid_odds, "effect_value"
    ].map(lambda value: max(value, 1.0 / value))
    effects["ranking_definition"] = effects["effect_measure_inference"].map(
        {
            "odds_ratio_inferred": "max(OR, 1/OR)",
            "beta_or_other_unit_effect": "absolute reported effect value",
        }
    )
    effects = effects.sort_values(
        ["effect_measure_inference", "ranking_value", "P", "normalized_variant_id"],
        ascending=[True, False, True, True],
        kind="stable",
    )

    columns = [
        "phenotype_scope",
        "normalized_variant_id",
        "STRONGEST SNP-RISK ALLELE",
        "MAPPED_GENE",
        "effect_value",
        "effect_measure_inference",
        "ranking_value",
        "ranking_definition",
        "95% CI (TEXT)",
        "P",
        "FIRST AUTHOR",
        "DATE",
        "STUDY ACCESSION",
        "PUBMEDID",
    ]
    if "phenotype_scope" not in effects.columns:
        effects.insert(0, "phenotype_scope", CORE_SCOPE)
    all_effects = effects[columns].reset_index(drop=True)

    top_frames: list[pd.DataFrame] = []
    for measure, group in effects.groupby("effect_measure_inference", sort=True):
        top = group.drop_duplicates("normalized_variant_id").head(top_n).copy()
        top.insert(0, "rank_within_measure", range(1, len(top) + 1))
        top_frames.append(top)
    top_effects = pd.concat(top_frames, ignore_index=True)
    top_columns = ["rank_within_measure", *columns]
    return all_effects, top_effects[top_columns]


def gene_profile(gws_variants: pd.DataFrame) -> pd.DataFrame:
    expanded = gws_variants.copy()
    expanded["gene"] = expanded["mapped_genes"].map(
        lambda value: [item for item in str(value).split(";") if item]
    )
    expanded = expanded.explode("gene").dropna(subset=["gene"])
    rows: list[dict[str, Any]] = []
    for gene, group in expanded.groupby("gene", sort=True):
        classes = set(group["consequence_class"])
        if "coding" in classes:
            profile = "has_coding"
        elif "splice" in classes:
            profile = "has_splice_no_coding"
        elif classes == {"noncoding"}:
            profile = "only_noncoding"
        elif classes == {"unknown"}:
            profile = "unknown_only"
        else:
            profile = "noncoding_with_unknown"
        rows.append(
            {
                "phenotype_scope": CORE_SCOPE,
                "gene": gene,
                "profile": profile,
                "n_gws_tag_variants": group["normalized_variant_id"].nunique(),
                "consequence_classes": ";".join(sorted(classes)),
                "gws_tag_variants": ";".join(sorted(set(group["normalized_variant_id"]))),
            }
        )
    return pd.DataFrame(rows).sort_values(["profile", "gene"]).reset_index(drop=True)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]*)"', text))


def read_gencode_genes(path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or not PRIMARY_CHROM_RE.fullmatch(fields[0]):
                continue
            attrs = parse_gtf_attributes(fields[8])
            rows.append(
                {
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gencode_gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene": attrs.get("gene_name", ""),
                    "gene_type": attrs.get("gene_type", ""),
                }
            )
    genes = pd.DataFrame(rows)
    genes["locus_length_bp"] = genes["end"] - genes["start"] + 1
    genes["_chrom_sort"] = genes["chromosome"].map(chromosome_sort_key)
    return genes.sort_values(["_chrom_sort", "start", "end", "gene"]).drop(
        columns="_chrom_sort"
    )


def build_locus_table(
    selected: dict[str, set[str]], gencode_path: Path
) -> pd.DataFrame:
    genes = read_gencode_genes(gencode_path)
    by_name = {name: frame for name, frame in genes.groupby("gene", sort=False)}
    flank_pool = genes.loc[genes["gene_type"].isin(["protein_coding", "lncRNA"])].copy()
    flank_by_chrom = {
        chrom: frame.sort_values(["start", "end", "gene"])
        for chrom, frame in flank_pool.groupby("chromosome", sort=False)
    }

    rows: list[dict[str, Any]] = []
    for query_gene in sorted(selected):
        matches = by_name.get(query_gene)
        base: dict[str, Any] = {
            "phenotype_scope": CORE_SCOPE,
            "gene": query_gene,
            "selection_reason": ";".join(sorted(selected[query_gene])),
            "gencode_version": "v50",
            "genome_build": "GRCh38",
        }
        if matches is None or matches.empty:
            rows.append(
                {
                    **base,
                    "gencode_match_status": "not_found_by_exact_gene_name",
                    "gencode_record_count": 0,
                }
            )
            continue
        target = matches.iloc[0]
        chrom_genes = flank_by_chrom.get(target["chromosome"], pd.DataFrame())
        other = chrom_genes.loc[chrom_genes["gene"].ne(query_gene)]
        left = other.loc[other["end"].lt(target["start"])].sort_values(
            ["end", "start", "gene"], ascending=[False, False, True]
        ).head(1)
        right = other.loc[other["start"].gt(target["end"])].sort_values(
            ["start", "end", "gene"], ascending=[True, True, True]
        ).head(1)
        left_row = left.iloc[0] if len(left) else None
        right_row = right.iloc[0] if len(right) else None
        rows.append(
            {
                **base,
                "gencode_match_status": "exact_gene_name",
                "gencode_record_count": len(matches),
                "gencode_gene_id": target["gencode_gene_id"],
                "chromosome": target["chromosome"],
                "start_1based": target["start"],
                "end_1based": target["end"],
                "strand": target["strand"],
                "gene_type": target["gene_type"],
                "locus_length_bp": target["locus_length_bp"],
                "left_flank_gene": "" if left_row is None else left_row["gene"],
                "left_flank_gene_id": "" if left_row is None else left_row["gencode_gene_id"],
                "left_flank_gene_type": "" if left_row is None else left_row["gene_type"],
                "left_flank_distance_bp": (
                    math.nan
                    if left_row is None
                    else int(target["start"] - left_row["end"] - 1)
                ),
                "right_flank_gene": "" if right_row is None else right_row["gene"],
                "right_flank_gene_id": "" if right_row is None else right_row["gencode_gene_id"],
                "right_flank_gene_type": "" if right_row is None else right_row["gene_type"],
                "right_flank_distance_bp": (
                    math.nan
                    if right_row is None
                    else int(right_row["start"] - target["end"] - 1)
                ),
            }
        )
    loci = pd.DataFrame(rows)
    loci["_chrom_sort"] = loci.get("chromosome", pd.Series(index=loci.index)).map(
        chromosome_sort_key
    )
    return (
        loci.sort_values(["_chrom_sort", "start_1based", "gene"], kind="stable")
        .drop(columns="_chrom_sort")
        .reset_index(drop=True)
    )


def manifest_rows(args: argparse.Namespace) -> list[tuple[str, str]]:
    command = shlex.join([sys.executable, str(SCRIPT), *sys.argv[1:]])
    rows = [
        ("run_timestamp_utc", datetime.now(timezone.utc).isoformat()),
        ("exact_replay_command", command),
        ("argv", shlex.join(sys.argv)),
        ("working_directory", str(Path.cwd())),
        ("script", str(SCRIPT)),
        ("python", platform.python_version()),
        ("python_executable", sys.executable),
        ("pandas", package_version("pandas")),
        ("pysam", package_version("pysam")),
        ("gwas_catalog_release", "2026-09-15 (local repository release)"),
        ("gencode_release", "v50 GRCh38"),
        ("1000_genomes_release", "20190312 biallelic SNV/indel GRCh38"),
        ("phenotype_scope", CORE_SCOPE),
        ("gws_threshold", str(GWS_THRESHOLD)),
        ("ld_expansion", "not performed"),
        ("network_calls", "none"),
    ]
    for label, path in (
        ("input_associations", args.associations),
        ("input_studies", args.studies),
        ("input_ancestries", args.ancestries),
        ("reference_gencode", args.gencode),
    ):
        rows.extend(
            [
                (f"{label}_path", str(path.resolve())),
                (f"{label}_sha256", sha256_file(path)),
            ]
        )
    rows.append(("script_sha256_at_run", sha256_file(SCRIPT)))
    return rows


def run(args: argparse.Namespace, logger: logging.Logger) -> dict[str, pd.DataFrame]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    logger.info("phenotype scope: %s only", CORE_SCOPE)
    logger.info("GWS threshold: P <= %.1e", GWS_THRESHOLD)
    logger.info("LD expansion: not performed")
    logger.info("network/API calls: none")

    studies = read_core_table(
        args.studies,
        ["STUDY ACCESSION", "PUBMED ID", "DATE", "INITIAL SAMPLE SIZE",
         "REPLICATION SAMPLE SIZE", "FULL SUMMARY STATISTICS"],
        "studies",
    )
    ancestries = read_core_table(
        args.ancestries,
        ["STUDY ACCESSION", "STAGE", "BROAD ANCESTRAL CATEGORY",
         "NUMBER OF INDIVIDUALS", "NUMBER OF CASES", "NUMBER OF CONTROLS"],
        "ancestries",
    )
    associations = read_core_table(
        args.associations,
        ["STUDY ACCESSION", "PUBMEDID", "SNPS", "P-VALUE", "CONTEXT",
         "MAPPED_GENE", "CHR_ID", "CHR_POS", "STRONGEST SNP-RISK ALLELE",
         "OR or BETA", "95% CI (TEXT)"],
        "associations",
    )
    logger.info(
        "loaded %d studies, %d ancestry rows, %d association rows",
        len(studies), len(ancestries), len(associations),
    )

    output_frames: dict[str, pd.DataFrame] = {}
    study_outputs = build_study_outputs(studies, ancestries, associations)
    output_frames.update(
        {
            "COPD-S2-R001_studies_core.tsv": study_outputs["studies"],
            "COPD-S2-R001_study_summary.tsv": study_outputs["study_summary"],
            "COPD-S2-R001_ancestry_summary.tsv": study_outputs["ancestry_summary"],
            "COPD-S2-R001_sample_summary.tsv": study_outputs["sample_summary"],
        }
    )

    normalized_assoc, variants = build_unique_variants(associations)
    variants = annotate_with_local_1000g(variants, args.thousand_genomes_dir, logger)
    type_counts = pd.concat(
        [
            make_count_table(variants, "identifier_type", "identifier_type"),
            make_count_table(variants, "variant_type_1000g", "sequence_variant_type"),
            make_count_table(variants, "local_reference_status", "local_reference_status"),
        ],
        ignore_index=True,
    )
    cons_counts = consequence_counts(variants)
    coding = coding_summary(variants)

    variant_annotation = variants[
        [
            "normalized_variant_id", "variant_type_1000g", "local_reference_status",
            "local_reference_ref", "local_reference_alt",
        ]
    ]
    normalized_assoc = normalized_assoc.merge(
        variant_annotation, on="normalized_variant_id", how="left", validate="many_to_one"
    )
    normalized_assoc.insert(0, "phenotype_scope", CORE_SCOPE)
    gws_assoc = normalized_assoc.loc[normalized_assoc["genome_wide_significant"]].copy()
    gws_assoc = gws_assoc.sort_values(
        ["P", "normalized_variant_id", "STUDY ACCESSION"], kind="stable"
    ).reset_index(drop=True)
    gws_variants = variants.loc[variants["genome_wide_significant"]].copy()

    output_frames.update(
        {
            "COPD-S2-R002_unique_tag_variants.tsv": variants,
            "COPD-S2-R002_variant_type_counts.tsv": type_counts,
            "COPD-S2-R002_consequence_counts.tsv": cons_counts,
            "COPD-S2-R002_gws_associations.tsv": gws_assoc,
            "COPD-S2-R002_gws_unique_tag_variants.tsv": gws_variants,
        }
    )

    all_genes, replicated_genes = build_gene_outputs(gws_assoc)
    all_effects, top_effects = build_effect_outputs(gws_assoc, args.top_n)
    output_frames.update(
        {
            "COPD-S2-R003A_gene_summary_all_gws.tsv": all_genes,
            "COPD-S2-R003A_genes_multistudy.tsv": replicated_genes,
            "COPD-S2-R003B_gws_effect_sizes_all.tsv": all_effects,
            "COPD-S2-R003B_top_effect_sizes.tsv": top_effects,
        }
    )

    profiles = gene_profile(gws_variants)
    consequence_table = variants[
        [
            "phenotype_scope", "normalized_variant_id", "catalog_consequences",
            "consequence_class", "genome_wide_significant", "mapped_genes",
            "best_p",
        ]
    ].copy()
    output_frames.update(
        {
            "COPD-S2-R004_coding_noncoding_summary.tsv": coding,
            "COPD-S2-R004_variant_consequences.tsv": consequence_table,
            "COPD-S2-R004_gene_coding_profile.tsv": profiles,
        }
    )

    selected: dict[str, set[str]] = defaultdict(set)
    for gene in replicated_genes["gene"]:
        selected[gene].add("gws_in_at_least_2_study_accessions")
    for mapped in top_effects["MAPPED_GENE"]:
        for gene in genes_of(mapped):
            selected[gene].add("top_reported_effect_size")
    loci = build_locus_table(selected, args.gencode)
    output_frames["COPD-S2-R005_locus_table_gencode_v50.tsv"] = loci

    for filename, frame in output_frames.items():
        path = RESULTS / filename
        write_tsv(frame, path)
        logger.info("wrote %s (%d rows)", path.relative_to(SECTION_ROOT), len(frame))

    logger.info(
        "core counts: studies=%d publications=%d tag_variants=%d GWS_rows=%d "
        "GWS_tag_variants=%d replicated_genes=%d",
        studies["STUDY ACCESSION"].nunique(),
        studies["PUBMED ID"].nunique(),
        len(variants),
        len(gws_assoc),
        len(gws_variants),
        len(replicated_genes),
    )
    coding_gws = coding.loc[coding["subset"].eq("gws_unique_tag_variants")].iloc[0]
    logger.info(
        "GWS consequence annotation: %d/%d annotated; coding=%d splice=%d "
        "noncoding=%d unknown=%d; noncoding fraction among annotated=%.6f",
        coding_gws["n_with_catalog_consequence"],
        coding_gws["n_unique_variants"],
        coding_gws["n_coding"],
        coding_gws["n_splice"],
        coding_gws["n_noncoding"],
        coding_gws["n_unknown"],
        coding_gws["noncoding_fraction_of_annotated"],
    )
    logger.info(
        "GENCODE loci: selected=%d exact_matches=%d unresolved=%d",
        len(loci),
        loci["gencode_match_status"].eq("exact_gene_name").sum(),
        loci["gencode_match_status"].ne("exact_gene_name").sum(),
    )
    return output_frames


def main() -> None:
    args = parse_args()
    logger = setup_logging(args.log)
    command = shlex.join([sys.executable, str(SCRIPT), *sys.argv[1:]])
    logger.info("exact replay command: %s", command)
    logger.info("argv: %s", shlex.join(sys.argv))
    logger.info("working directory: %s", Path.cwd())
    logger.info(
        "versions: Python=%s pandas=%s pysam=%s",
        platform.python_version(), package_version("pandas"), package_version("pysam")
    )
    try:
        frames = run(args, logger)
        manifest = manifest_rows(args)
        for filename, frame in frames.items():
            manifest.append((f"output_rows:{filename}", str(len(frame))))
        manifest_path = LOGS / "COPD-S2-core_analysis_manifest.tsv"
        write_tsv(pd.DataFrame(manifest, columns=["key", "value"]), manifest_path)
        logger.info("wrote %s", manifest_path.relative_to(SECTION_ROOT))
        logger.info("completed successfully")
    except Exception:
        logger.exception("analysis failed")
        raise


if __name__ == "__main__":
    main()
