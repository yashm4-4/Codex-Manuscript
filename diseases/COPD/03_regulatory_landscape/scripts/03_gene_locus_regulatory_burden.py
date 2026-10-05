#!/usr/bin/env python3
"""Quantify severe-emphysema lung regulatory elements in COPD GWAS loci.

Definitions follow the disease framework:

* preliminary enhancer: an H3K27ac narrowPeak;
* preliminary silencer: an H3K27me3 narrowPeak;
* refined enhancer: an ATAC narrowPeak overlapping H3K27ac in the same lobe;
* refined silencer: an ATAC narrowPeak overlapping H3K27me3 in the same lobe.

Three lung lobes from donor ENCDO520EJG are analyzed independently.  A donor
union is then formed by merging overlapping or book-ended intervals across the
three lobes.  COPD loci are GENCODE v50 gene bodies plus/minus 100 kb by
default.  Ambiguous gene symbols with multiple GENCODE gene records are
excluded rather than assigning an arbitrary genomic location.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Iterable, Iterator

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION = SCRIPT.parents[1]
PROJECT = SECTION.parents[2]
GWAS = PROJECT / "diseases" / "COPD" / "02_gwas"
CANONICAL_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CANONICAL = set(CANONICAL_ORDER)
CHROM_RANK = {chrom: rank for rank, chrom in enumerate(CANONICAL_ORDER)}
LOBE_ORDER = ["upper_right", "lower_right", "lower_left"]


@dataclass(frozen=True, order=True)
class Interval:
    chrom: str
    start: int
    end: int
    name: str


@dataclass(frozen=True)
class ElementSet:
    context: str
    definition: str
    element_type: str
    intervals: tuple[Interval, ...]
    source_accessions: str
    output_path: Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--locus-table",
        type=Path,
        default=GWAS / "results" / "COPD-S2-R005_locus_table_gencode_v50.tsv",
    )
    parser.add_argument(
        "--replicated-table",
        type=Path,
        default=GWAS / "results" / "COPD-S2-R003A_genes_multistudy.tsv",
    )
    parser.add_argument(
        "--peak-manifest",
        type=Path,
        default=SECTION / "data" / "encode_peaks" / "manifest.tsv",
    )
    parser.add_argument("--flank-bp", type=int, default=100_000)
    parser.add_argument(
        "--element-dir",
        type=Path,
        default=SECTION / "data" / "copd_donor_regulatory_elements",
    )
    parser.add_argument("--results-dir", type=Path, default=SECTION / "results")
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


def natural_key(interval: Interval) -> tuple[int, int, int, str]:
    return CHROM_RANK[interval.chrom], interval.start, interval.end, interval.name


def write_gzip_rows(path: Path, rows: Iterable[tuple[object, ...]]) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as out:
                for row in rows:
                    out.write("\t".join(str(value) for value in row) + "\n")
    temporary.replace(path)


def read_peak_file(path: Path, accession: str) -> tuple[Interval, ...]:
    intervals: list[Interval] = []
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three BED columns")
            chrom = fields[0]
            if chrom not in CANONICAL:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            original_name = fields[3] if len(fields) > 3 and fields[3] else str(line_number)
            intervals.append(Interval(chrom, start, end, f"{accession}:{original_name}"))
    if not intervals:
        raise RuntimeError(f"No canonical intervals in {path}")
    return tuple(sorted(intervals, key=natural_key))


def group_by_chrom(intervals: Iterable[Interval]) -> dict[str, list[Interval]]:
    result: dict[str, list[Interval]] = defaultdict(list)
    for interval in intervals:
        result[interval.chrom].append(interval)
    for chrom in result:
        result[chrom].sort(key=lambda item: (item.start, item.end, item.name))
    return result


def select_overlapping(query: Iterable[Interval], reference: Iterable[Interval]) -> tuple[Interval, ...]:
    """Return query intervals with at least one half-open overlap in reference."""
    query_by_chrom = group_by_chrom(query)
    reference_by_chrom = group_by_chrom(reference)
    selected: list[Interval] = []
    for chrom in CANONICAL_ORDER:
        queries = query_by_chrom.get(chrom, [])
        refs = reference_by_chrom.get(chrom, [])
        left = 0
        for query_interval in queries:
            while left < len(refs) and refs[left].end <= query_interval.start:
                left += 1
            index = left
            overlaps = False
            while index < len(refs) and refs[index].start < query_interval.end:
                if refs[index].end > query_interval.start:
                    overlaps = True
                    break
                index += 1
            if overlaps:
                selected.append(query_interval)
    return tuple(sorted(selected, key=natural_key))


def merge_intervals(intervals: Iterable[Interval], name_prefix: str) -> tuple[Interval, ...]:
    merged: list[Interval] = []
    by_chrom = group_by_chrom(intervals)
    counter = 0
    for chrom in CANONICAL_ORDER:
        chrom_intervals = by_chrom.get(chrom, [])
        if not chrom_intervals:
            continue
        start, end = chrom_intervals[0].start, chrom_intervals[0].end
        for interval in chrom_intervals[1:]:
            if interval.start <= end:
                end = max(end, interval.end)
            else:
                counter += 1
                merged.append(Interval(chrom, start, end, f"{name_prefix}_{counter:07d}"))
                start, end = interval.start, interval.end
        counter += 1
        merged.append(Interval(chrom, start, end, f"{name_prefix}_{counter:07d}"))
    return tuple(merged)


def resolve_peak_path(row: pd.Series) -> Path:
    local = SECTION / "data" / "encode_peaks" / f"{row['file_accession']}.bed.gz"
    if local.is_file():
        return local
    recorded = Path(str(row["path"]))
    if recorded.is_file():
        return recorded
    raise FileNotFoundError(local)


def load_peaks(manifest_path: Path) -> tuple[pd.DataFrame, dict[tuple[str, str], tuple[Interval, ...]]]:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "lobe",
        "mark",
        "file_accession",
        "donor_accessions",
        "donor_health_status",
        "donor_age",
        "donor_age_units",
        "donor_sex",
        "donor_ethnicity",
        "donor_race",
        "assembly",
        "sha256",
    }
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Peak manifest lacks columns: {sorted(missing)}")
    if set(manifest["lobe"]) != set(LOBE_ORDER):
        raise ValueError(f"Unexpected lobe set: {sorted(set(manifest['lobe']))}")
    if set(manifest["mark"]) != {"accessibility", "H3K27ac", "H3K27me3"}:
        raise ValueError(f"Unexpected mark set: {sorted(set(manifest['mark']))}")
    if set(manifest["donor_accessions"]) != {"ENCDO520EJG"}:
        raise ValueError("Peak files are not donor matched to ENCDO520EJG")
    if not manifest["donor_health_status"].str.lower().str.contains("severe emphysema").all():
        raise ValueError("Peak manifest does not consistently record severe emphysema")
    expected_donor_metadata = {
        "donor_age": {"60"},
        "donor_age_units": {"year"},
        "donor_sex": {"male"},
        "donor_ethnicity": {"European"},
    }
    for column, expected in expected_donor_metadata.items():
        observed = set(manifest[column])
        if observed != expected:
            raise ValueError(
                f"Unexpected ENCDO520EJG {column}: {sorted(observed)}; expected {sorted(expected)}"
            )
    if set(manifest["assembly"]) != {"GRCh38"}:
        raise ValueError("Peak manifest is not uniformly GRCh38")
    if manifest.groupby(["lobe", "mark"]).size().ne(1).any():
        raise ValueError("Expected exactly one file for each lobe and mark")

    peaks: dict[tuple[str, str], tuple[Interval, ...]] = {}
    resolved_paths: list[str] = []
    canonical_counts: list[int] = []
    for _, row in manifest.iterrows():
        path = resolve_peak_path(row)
        observed_sha = sha256(path)
        if observed_sha != row["sha256"]:
            raise RuntimeError(
                f"Checksum mismatch for {row['file_accession']}: {observed_sha} != {row['sha256']}"
            )
        intervals = read_peak_file(path, row["file_accession"])
        peaks[(row["lobe"], row["mark"])] = intervals
        resolved_paths.append(project_path(path))
        canonical_counts.append(len(intervals))
    manifest = manifest.copy()
    manifest["resolved_project_path"] = resolved_paths
    manifest["canonical_interval_count"] = canonical_counts
    return manifest, peaks


def prepare_loci(
    locus_table_path: Path, replicated_table_path: Path, flank_bp: int
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    if flank_bp < 0:
        raise ValueError("--flank-bp cannot be negative")
    loci = pd.read_csv(locus_table_path, sep="\t", keep_default_na=False)
    replicated = pd.read_csv(replicated_table_path, sep="\t", keep_default_na=False)
    required = {
        "gene",
        "selection_reason",
        "genome_build",
        "gencode_match_status",
        "gencode_record_count",
        "gencode_gene_id",
        "chromosome",
        "start_1based",
        "end_1based",
        "gene_type",
    }
    missing = required - set(loci.columns)
    if missing:
        raise ValueError(f"Locus table lacks columns: {sorted(missing)}")
    if loci["gene"].duplicated().any():
        raise ValueError("Locus table contains duplicate selected gene symbols")
    if not loci["genome_build"].eq("GRCh38").all():
        raise ValueError("Locus table is not uniformly GRCh38")

    replicated_genes = set(
        replicated.loc[
            replicated["replicated_across_studies"].astype(str).str.lower().eq("true"),
            "gene",
        ]
    )
    selected_genes = set(loci["gene"])
    if not replicated_genes <= selected_genes:
        missing_replicated = sorted(replicated_genes - selected_genes)
        raise ValueError(f"Replicated genes missing from selected table: {missing_replicated}")

    exclusion_rows: list[dict[str, object]] = []
    keep_indices: list[int] = []
    for index, row in loci.iterrows():
        reasons: list[str] = []
        if row["gencode_match_status"] != "exact_gene_name":
            reasons.append(f"gencode_match_status={row['gencode_match_status']}")
        try:
            record_count = int(row["gencode_record_count"])
        except (TypeError, ValueError):
            record_count = 0
        if record_count != 1:
            reasons.append(f"ambiguous_gencode_record_count={record_count}")
        if row["chromosome"] not in CANONICAL:
            reasons.append(f"noncanonical_chromosome={row['chromosome']}")
        try:
            start, end = int(row["start_1based"]), int(row["end_1based"])
            if start < 1 or end < start:
                reasons.append("invalid_gene_coordinates")
        except (TypeError, ValueError):
            reasons.append("missing_gene_coordinates")
        if reasons:
            exclusion_rows.append(
                {
                    "gene": row["gene"],
                    "is_replicated": row["gene"] in replicated_genes,
                    "selection_reason": row["selection_reason"],
                    "exclusion_reason": ";".join(reasons),
                }
            )
        else:
            keep_indices.append(index)

    kept = loci.loc[keep_indices].copy()
    kept["is_replicated"] = kept["gene"].isin(replicated_genes)
    kept["gene_start_0based"] = kept["start_1based"].astype(int) - 1
    kept["gene_end_0based_exclusive"] = kept["end_1based"].astype(int)
    kept["locus_start_0based"] = (kept["gene_start_0based"] - flank_bp).clip(lower=0)
    kept["locus_end_0based_exclusive"] = kept["gene_end_0based_exclusive"] + flank_bp
    kept["locus_span_bp"] = kept["locus_end_0based_exclusive"] - kept["locus_start_0based"]
    kept["flank_bp_each_side"] = flank_bp
    kept["genome_build"] = "GRCh38"
    kept = kept.sort_values(
        ["chromosome", "locus_start_0based", "gene"],
        key=lambda series: series.map(CHROM_RANK) if series.name == "chromosome" else series,
        kind="stable",
    )
    columns = [
        "gene",
        "gencode_gene_id",
        "gene_type",
        "selection_reason",
        "is_replicated",
        "chromosome",
        "gene_start_0based",
        "gene_end_0based_exclusive",
        "locus_start_0based",
        "locus_end_0based_exclusive",
        "locus_span_bp",
        "flank_bp_each_side",
        "genome_build",
    ]
    exclusions = pd.DataFrame(
        exclusion_rows,
        columns=["gene", "is_replicated", "selection_reason", "exclusion_reason"],
    )
    return kept[columns].reset_index(drop=True), exclusions, replicated_genes


def make_element_sets(
    manifest: pd.DataFrame,
    peaks: dict[tuple[str, str], tuple[Interval, ...]],
    output_dir: Path,
) -> tuple[list[ElementSet], pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    sets: list[ElementSet] = []
    metadata_rows: list[dict[str, object]] = []
    definitions = [
        ("enhancer", "H3K27ac"),
        ("silencer", "H3K27me3"),
    ]
    by_key: dict[tuple[str, str, str], ElementSet] = {}

    for lobe in LOBE_ORDER:
        accessibility = peaks[(lobe, "accessibility")]
        atac_accession = manifest.loc[
            (manifest["lobe"] == lobe) & (manifest["mark"] == "accessibility"),
            "file_accession",
        ].iat[0]
        for element_type, mark in definitions:
            histone = peaks[(lobe, mark)]
            histone_accession = manifest.loc[
                (manifest["lobe"] == lobe) & (manifest["mark"] == mark),
                "file_accession",
            ].iat[0]
            options = [
                (
                    "preliminary_histone_peak",
                    histone,
                    histone_accession,
                ),
                (
                    "refined_accessibility_histone_overlap",
                    select_overlapping(accessibility, histone),
                    f"{atac_accession};{histone_accession}",
                ),
            ]
            for definition, intervals, source_accessions in options:
                if not intervals:
                    raise RuntimeError(f"No {definition} {element_type} intervals for {lobe}")
                output = output_dir / f"{lobe}.{element_type}.{definition}.grch38.bed.gz"
                write_gzip_rows(
                    output,
                    (
                        (
                            interval.chrom,
                            interval.start,
                            interval.end,
                            interval.name,
                            element_type,
                            definition,
                            lobe,
                            source_accessions,
                        )
                        for interval in intervals
                    ),
                )
                element_set = ElementSet(
                    lobe,
                    definition,
                    element_type,
                    tuple(intervals),
                    source_accessions,
                    output,
                )
                sets.append(element_set)
                by_key[(lobe, definition, element_type)] = element_set
                metadata_rows.append(
                    {
                        "context": lobe,
                        "definition": definition,
                        "element_type": element_type,
                        "n_intervals": len(intervals),
                        "source_accessions": source_accessions,
                        "genome_build": "GRCh38",
                        "coordinate_system": "0-based half-open BED",
                        "output_path": project_path(output),
                        "sha256": sha256(output),
                    }
                )

    for element_type, _ in definitions:
        for definition in (
            "preliminary_histone_peak",
            "refined_accessibility_histone_overlap",
        ):
            members = [by_key[(lobe, definition, element_type)] for lobe in LOBE_ORDER]
            merged = merge_intervals(
                (interval for member in members for interval in member.intervals),
                f"COPD_{element_type}_{definition}_union",
            )
            sources = ";".join(
                sorted({item for member in members for item in member.source_accessions.split(";")})
            )
            output = output_dir / f"donor_union.{element_type}.{definition}.grch38.bed.gz"
            write_gzip_rows(
                output,
                (
                    (
                        interval.chrom,
                        interval.start,
                        interval.end,
                        interval.name,
                        element_type,
                        definition,
                        "donor_union",
                        sources,
                    )
                    for interval in merged
                ),
            )
            element_set = ElementSet(
                "donor_union", definition, element_type, merged, sources, output
            )
            sets.append(element_set)
            metadata_rows.append(
                {
                    "context": "donor_union",
                    "definition": definition,
                    "element_type": element_type,
                    "n_intervals": len(merged),
                    "source_accessions": sources,
                    "genome_build": "GRCh38",
                    "coordinate_system": "0-based half-open BED",
                    "output_path": project_path(output),
                    "sha256": sha256(output),
                }
            )

    metadata = pd.DataFrame(metadata_rows)
    context_rank = {name: rank for rank, name in enumerate(LOBE_ORDER + ["donor_union"])}
    metadata["_context_rank"] = metadata["context"].map(context_rank)
    metadata = metadata.sort_values(
        ["_context_rank", "definition", "element_type"], kind="stable"
    ).drop(columns="_context_rank")
    return sets, metadata.reset_index(drop=True)


def intervals_from_loci(frame: pd.DataFrame) -> tuple[Interval, ...]:
    return tuple(
        Interval(
            row.chromosome,
            int(row.locus_start_0based),
            int(row.locus_end_0based_exclusive),
            row.gene,
        )
        for row in frame.itertuples(index=False)
    )


def overlap_count_by_locus(
    loci: Iterable[Interval], features: Iterable[Interval]
) -> dict[str, int]:
    feature_by_chrom = group_by_chrom(features)
    result: dict[str, int] = {}
    for locus in loci:
        chrom_features = feature_by_chrom.get(locus.chrom, [])
        starts = [feature.start for feature in chrom_features]
        stop = bisect_left(starts, locus.end)
        result[locus.name] = sum(
            feature.end > locus.start for feature in chrom_features[:stop]
        )
    return result


def count_features_in_union(
    features: Iterable[Interval], merged_loci: Iterable[Interval]
) -> int:
    return len(select_overlapping(features, merged_loci))


def locus_union_bp(loci: Iterable[Interval]) -> int:
    merged = merge_intervals(loci, "locus_union")
    return sum(interval.end - interval.start for interval in merged)


def quantify(
    loci: pd.DataFrame, element_sets: list[ElementSet]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    locus_intervals = intervals_from_loci(loci)
    replicated_names = set(loci.loc[loci["is_replicated"], "gene"])
    per_gene_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for element_set in element_sets:
        counts = overlap_count_by_locus(locus_intervals, element_set.intervals)
        for row in loci.itertuples(index=False):
            count = counts[row.gene]
            per_gene_rows.append(
                {
                    "gene": row.gene,
                    "gencode_gene_id": row.gencode_gene_id,
                    "gene_type": row.gene_type,
                    "selection_reason": row.selection_reason,
                    "is_replicated": row.is_replicated,
                    "chromosome": row.chromosome,
                    "locus_start_0based": row.locus_start_0based,
                    "locus_end_0based_exclusive": row.locus_end_0based_exclusive,
                    "locus_span_bp": row.locus_span_bp,
                    "context": element_set.context,
                    "definition": element_set.definition,
                    "element_type": element_set.element_type,
                    "element_count": count,
                    "elements_per_mb_locus": count / row.locus_span_bp * 1_000_000,
                    "genome_build": "GRCh38",
                }
            )

        for locus_set, names in (
            ("selected_unambiguous", set(loci["gene"])),
            ("replicated_unambiguous", replicated_names),
        ):
            subset_loci = tuple(interval for interval in locus_intervals if interval.name in names)
            merged_loci = merge_intervals(subset_loci, f"{locus_set}_union")
            union_bp = sum(interval.end - interval.start for interval in merged_loci)
            overlapping = count_features_in_union(element_set.intervals, merged_loci)
            values = pd.Series([counts[name] for name in sorted(names)], dtype="int64")
            summary_rows.append(
                {
                    "locus_set": locus_set,
                    "n_gene_loci": len(names),
                    "n_merged_locus_intervals": len(merged_loci),
                    "locus_union_bp": union_bp,
                    "context": element_set.context,
                    "definition": element_set.definition,
                    "element_type": element_set.element_type,
                    "elements_total": len(element_set.intervals),
                    "unique_elements_overlapping_locus_union": overlapping,
                    "fraction_elements_overlapping_locus_union": overlapping
                    / len(element_set.intervals),
                    "elements_per_mb_locus_union": overlapping / union_bp * 1_000_000,
                    "gene_loci_with_at_least_one_element": int((values > 0).sum()),
                    "fraction_gene_loci_with_at_least_one_element": float((values > 0).mean()),
                    "per_gene_min": int(values.min()),
                    "per_gene_q1": float(values.quantile(0.25)),
                    "per_gene_median": float(values.median()),
                    "per_gene_mean": float(values.mean()),
                    "per_gene_q3": float(values.quantile(0.75)),
                    "per_gene_max": int(values.max()),
                    "genome_build": "GRCh38",
                    "flank_bp_each_side": int(loci["flank_bp_each_side"].iat[0]),
                }
            )

    per_gene = pd.DataFrame(per_gene_rows)
    summary = pd.DataFrame(summary_rows)
    context_rank = {name: rank for rank, name in enumerate(LOBE_ORDER + ["donor_union"])}
    for frame in (per_gene, summary):
        frame["_context_rank"] = frame["context"].map(context_rank)
    per_gene = per_gene.sort_values(
        ["_context_rank", "definition", "element_type", "chromosome", "locus_start_0based"],
        key=lambda series: series.map(CHROM_RANK) if series.name == "chromosome" else series,
        kind="stable",
    ).drop(columns="_context_rank")
    summary = summary.sort_values(
        ["locus_set", "_context_rank", "definition", "element_type"], kind="stable"
    ).drop(columns="_context_rank")
    return per_gene.reset_index(drop=True), summary.reset_index(drop=True)


def write_analysis_manifest(
    args: argparse.Namespace,
    peak_manifest: pd.DataFrame,
    element_metadata: pd.DataFrame,
    loci: pd.DataFrame,
    exclusions: pd.DataFrame,
    output: Path,
) -> None:
    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exact_replay_command": shlex.join([sys.executable, str(SCRIPT), *sys.argv[1:]]),
        "working_directory": str(Path.cwd()),
        "script": project_path(SCRIPT),
        "script_sha256": sha256(SCRIPT),
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "genome_build": "GRCh38",
        "coordinate_system": "0-based half-open BED",
        "donor": "ENCDO520EJG",
        "donor_context": "lungs with severe emphysema",
        "donor_metadata": peak_manifest[
            [
                "donor_accessions",
                "donor_health_status",
                "donor_age",
                "donor_age_units",
                "donor_sex",
                "donor_ethnicity",
                "donor_race",
            ]
        ]
        .drop_duplicates()
        .to_dict(orient="records"),
        "lung_lobes": LOBE_ORDER,
        "locus_definition": f"GENCODE v50 gene body plus/minus {args.flank_bp} bp",
        "selected_genes_input": int(len(loci) + len(exclusions)),
        "selected_genes_analyzed": int(len(loci)),
        "selected_genes_excluded": int(len(exclusions)),
        "replicated_genes_analyzed": int(loci["is_replicated"].sum()),
        "preliminary_enhancer_definition": "H3K27ac narrowPeak",
        "preliminary_silencer_definition": "H3K27me3 narrowPeak",
        "refined_enhancer_definition": "ATAC narrowPeak with >=1 bp overlap of same-lobe H3K27ac narrowPeak",
        "refined_silencer_definition": "ATAC narrowPeak with >=1 bp overlap of same-lobe H3K27me3 narrowPeak",
        "donor_union_definition": "merge overlapping or book-ended intervals across all three lobes",
        "blacklist_filter": "not applied to descriptive ENCODE peak counts",
        "locus_table": project_path(args.locus_table),
        "locus_table_sha256": sha256(args.locus_table),
        "replicated_table": project_path(args.replicated_table),
        "replicated_table_sha256": sha256(args.replicated_table),
        "peak_manifest": project_path(args.peak_manifest),
        "peak_manifest_sha256": sha256(args.peak_manifest),
        "peak_files": peak_manifest[
            ["file_accession", "resolved_project_path", "sha256", "canonical_interval_count"]
        ].to_dict(orient="records"),
        "derived_element_sets": element_metadata.to_dict(orient="records"),
    }
    output.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    args = arguments()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.element_dir.mkdir(parents=True, exist_ok=True)

    peak_manifest, peaks = load_peaks(args.peak_manifest)
    loci, exclusions, replicated_genes = prepare_loci(
        args.locus_table, args.replicated_table, args.flank_bp
    )
    element_sets, element_metadata = make_element_sets(
        peak_manifest, peaks, args.element_dir
    )
    per_gene, summary = quantify(loci, element_sets)

    if len(replicated_genes) != 143:
        raise RuntimeError(f"Expected 143 replicated input genes, found {len(replicated_genes)}")
    if len(loci) + len(exclusions) != 155:
        raise RuntimeError("Selected-gene accounting does not sum to 155")
    if len(element_sets) != 16:
        raise RuntimeError(f"Expected 16 lobe/union element sets, found {len(element_sets)}")
    if len(per_gene) != len(loci) * len(element_sets):
        raise RuntimeError("Per-gene burden table has an unexpected number of rows")
    if summary["unique_elements_overlapping_locus_union"].gt(summary["elements_total"]).any():
        raise RuntimeError("A locus-overlap count exceeds its element-set total")

    outputs = {
        "COPD-S3-R002_gene_loci_grch38.tsv": loci,
        "COPD-S3-R002_locus_exclusions.tsv": exclusions,
        "COPD-S3-R002_donor_element_sets.tsv": element_metadata,
        "COPD-S3-R002_lobe_union_summary.tsv": summary,
        "COPD-S3-R002_per_gene_burden.tsv": per_gene,
    }
    for filename, frame in outputs.items():
        path = args.results_dir / filename
        frame.to_csv(path, sep="\t", index=False)
        print(f"{filename}\t{len(frame)} rows")

    write_analysis_manifest(
        args,
        peak_manifest,
        element_metadata,
        loci,
        exclusions,
        args.results_dir / "COPD-S3-R002_analysis_manifest.json",
    )
    print(f"selected genes analyzed\t{len(loci)}")
    print(f"replicated genes analyzed\t{int(loci['is_replicated'].sum())}")
    print(f"ambiguous selected genes excluded\t{len(exclusions)}")
    print("donor\tENCDO520EJG (lungs with severe emphysema)")


if __name__ == "__main__":
    main()
