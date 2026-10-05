#!/usr/bin/env python3
"""Analyze length, chromosome distribution, and conservation of COPD loci.

The analysis set is the 140 replicated COPD genes that resolve to exactly one
GENCODE v50 gene record.  The comparison universe contains all uniquely named
GENCODE v50 genes on canonical chromosomes whose gene biotypes occur in the
COPD set, excluding the COPD genes themselves.  Randomization is stratified by
gene biotype, preserving the observed COPD biotype composition.

Operational loci are gene bodies plus 100 kb on each side, clipped to GRCh38
chromosome bounds.  Conservation is measured across gene bodies with UCSC
100-way phyloP and phastCons bigWigs; mean0 treats bases without a bigWig value
as zero and coverage is reported separately.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import math
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
REPO_ROOT = SECTION_ROOT.parents[2]
SHARED_DATA = REPO_ROOT / "data"
RESULTS = SECTION_ROOT / "results"
LOGS = SECTION_ROOT / "logs"
WORK = SECTION_ROOT / "data" / "locus_properties"

CANONICAL = [f"chr{value}" for value in range(1, 23)] + ["chrX", "chrY"]
CHROM_RANK = {chrom: index for index, chrom in enumerate(CANONICAL)}
ATTRIBUTE_RE = re.compile(r'^\s*([^ ]+)\s+"([^"]*)"\s*$')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--locus-table",
        type=Path,
        default=RESULTS / "COPD-S2-R005_locus_table_gencode_v50.tsv",
    )
    parser.add_argument(
        "--gencode",
        type=Path,
        default=SHARED_DATA / "gencode" / "gencode.v50.annotation.gtf.gz",
    )
    parser.add_argument(
        "--phylop",
        type=Path,
        default=SHARED_DATA / "conservation" / "hg38.phyloP100way.bw",
    )
    parser.add_argument(
        "--phastcons",
        type=Path,
        default=SHARED_DATA / "conservation" / "hg38.phastCons100way.bw",
    )
    parser.add_argument("--bigwig-info", default="bigWigInfo")
    parser.add_argument(
        "--bigwig-average-over-bed", default="bigWigAverageOverBed"
    )
    parser.add_argument("--flank-bp", type=int, default=100_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument(
        "--log", type=Path, default=LOGS / "COPD-S2-locus_properties.log"
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2-locus-properties")
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


def resolve_executable(command: str) -> str:
    path = shutil.which(command)
    if path is None:
        raise FileNotFoundError(
            f"required executable {command!r} is not on PATH; load a UCSC module"
        )
    return path


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        match = ATTRIBUTE_RE.match(item)
        if match:
            result[match.group(1)] = match.group(2)
    return result


def load_gencode(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in CHROM_RANK:
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id", "").split(".", 1)[0]
            gene_name = attributes.get("gene_name", "")
            gene_type = attributes.get("gene_type", "")
            if not gene_id or not gene_name or not gene_type:
                raise ValueError(f"incomplete GENCODE gene attributes: {line[:200]!r}")
            start = int(fields[3])
            end = int(fields[4])
            rows.append(
                {
                    "gencode_gene_id": gene_id,
                    "gene": gene_name,
                    "gene_type": gene_type,
                    "chromosome": fields[0],
                    "start_1based": start,
                    "end_1based": end,
                    "strand": fields[6],
                    "gene_body_length_bp": end - start + 1,
                }
            )
    frame = pd.DataFrame(rows)
    if frame.empty or frame["gencode_gene_id"].duplicated().any():
        raise AssertionError("GENCODE gene records are empty or non-unique by gene ID")
    return frame


def chromosome_lengths(bigwig: Path, bigwig_info: str) -> dict[str, int]:
    result = subprocess.run(
        [bigwig_info, "-chroms", str(bigwig)],
        capture_output=True,
        text=True,
        check=True,
    )
    lengths: dict[str, int] = {}
    line_re = re.compile(r"^\s+(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y))\s+\d+\s+(\d+)$")
    for line in result.stdout.splitlines():
        match = line_re.match(line)
        if match:
            lengths[match.group(1)] = int(match.group(2))
    if set(lengths) != set(CANONICAL):
        raise AssertionError("conservation bigWig lacks canonical GRCh38 chromosome sizes")
    return lengths


def add_operational_loci(
    frame: pd.DataFrame, lengths: dict[str, int], flank_bp: int
) -> pd.DataFrame:
    result = frame.copy()
    result["gene_start_0based"] = result["start_1based"] - 1
    result["gene_end_0based_exclusive"] = result["end_1based"]
    result["locus_start_0based"] = (
        result["gene_start_0based"] - flank_bp
    ).clip(lower=0)
    result["chromosome_length_bp"] = result["chromosome"].map(lengths)
    result["locus_end_0based_exclusive"] = np.minimum(
        result["gene_end_0based_exclusive"] + flank_bp,
        result["chromosome_length_bp"],
    )
    result["operational_locus_length_bp"] = (
        result["locus_end_0based_exclusive"] - result["locus_start_0based"]
    )
    return result


def write_gene_bed(frame: pd.DataFrame, path: Path) -> None:
    ordered = frame.sort_values(
        ["chromosome", "gene_start_0based", "gene_end_0based_exclusive"],
        key=lambda series: (
            series.map(CHROM_RANK) if series.name == "chromosome" else series
        ),
        kind="stable",
    )
    with path.open("w", encoding="utf-8") as handle:
        for row in ordered.itertuples(index=False):
            handle.write(
                f"{row.chromosome}\t{row.gene_start_0based}\t"
                f"{row.gene_end_0based_exclusive}\t{row.gencode_gene_id}\n"
            )


def score_bigwig(
    bigwig: Path,
    bed: Path,
    output: Path,
    tool: str,
    prefix: str,
) -> pd.DataFrame:
    subprocess.run([tool, str(bigwig), str(bed), str(output)], check=True)
    columns = ["gencode_gene_id", "size", "covered", "sum", "mean0", "covered_mean"]
    frame = pd.read_csv(output, sep="\t", names=columns)
    if frame["gencode_gene_id"].duplicated().any():
        raise AssertionError(f"duplicate gene scores from {bigwig}")
    frame[f"{prefix}_coverage_fraction"] = frame["covered"] / frame["size"]
    return frame[
        [
            "gencode_gene_id",
            f"{prefix}_coverage_fraction",
            "mean0",
            "covered_mean",
        ]
    ].rename(
        columns={
            "mean0": f"{prefix}_mean0",
            "covered_mean": f"{prefix}_covered_mean",
        }
    )


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def empirical_test(observed: float, null: np.ndarray) -> dict[str, float]:
    center = float(null.mean())
    deviation = abs(observed - center)
    two_sided = (1 + int((np.abs(null - center) >= deviation).sum())) / (len(null) + 1)
    greater = (1 + int((null >= observed).sum())) / (len(null) + 1)
    less = (1 + int((null <= observed).sum())) / (len(null) + 1)
    sd = float(null.std(ddof=1))
    return {
        "observed": observed,
        "null_expected": center,
        "difference": observed - center,
        "ratio": observed / center if center != 0 else math.nan,
        "percent_difference": 100 * (observed - center) / center if center != 0 else math.nan,
        "standardized_difference": (observed - center) / sd if sd > 0 else math.nan,
        "null_ci025": float(np.quantile(null, 0.025)),
        "null_ci975": float(np.quantile(null, 0.975)),
        "p_two_sided": two_sided,
        "p_greater": greater,
        "p_less": less,
    }


def randomization(
    cases: pd.DataFrame,
    background: pd.DataFrame,
    metric_names: list[str],
    permutations: int,
    seed: int,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    rng = np.random.default_rng(seed)
    case_type_counts = cases["gene_type"].value_counts().sort_index()
    pools: list[tuple[int, np.ndarray, np.ndarray]] = []
    for gene_type, count in case_type_counts.items():
        pool = background.loc[background["gene_type"].eq(gene_type)]
        if len(pool) < count:
            raise ValueError(
                f"background for {gene_type} has {len(pool)} genes, fewer than {count} cases"
            )
        pools.append(
            (
                int(count),
                pool[metric_names].to_numpy(dtype=float),
                pool["chromosome"].map(CHROM_RANK).to_numpy(dtype=int),
            )
        )

    null_mean = np.empty((permutations, len(metric_names)), dtype=float)
    null_median = np.empty((permutations, len(metric_names)), dtype=float)
    null_chrom = np.zeros((permutations, len(CANONICAL)), dtype=np.int16)
    sample = np.empty((len(cases), len(metric_names)), dtype=float)
    sample_chrom = np.empty(len(cases), dtype=int)

    for index in range(permutations):
        offset = 0
        for count, matrix, chroms in pools:
            chosen = rng.choice(len(matrix), size=count, replace=False)
            sample[offset : offset + count] = matrix[chosen]
            sample_chrom[offset : offset + count] = chroms[chosen]
            offset += count
        null_mean[index] = sample.mean(axis=0)
        null_median[index] = np.median(sample, axis=0)
        null_chrom[index] = np.bincount(sample_chrom, minlength=len(CANONICAL))

    nulls: dict[str, np.ndarray] = {}
    for metric_index, metric in enumerate(metric_names):
        nulls[f"{metric}:mean"] = null_mean[:, metric_index]
        nulls[f"{metric}:median"] = null_median[:, metric_index]
    return nulls, null_chrom


def metric_summary(
    cases: pd.DataFrame,
    nulls: dict[str, np.ndarray],
    metrics: list[tuple[str, str, str]],
    family: str,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    positive_direction = "longer" if family == "locus_length" else "higher_conservation"
    negative_direction = "shorter" if family == "locus_length" else "lower_conservation"
    for metric, label, units in metrics:
        values = cases[metric].to_numpy(dtype=float)
        for statistic, observed in (("mean", float(values.mean())), ("median", float(np.median(values)))):
            test = empirical_test(observed, nulls[f"{metric}:{statistic}"])
            rows.append(
                {
                    "phenotype_scope": "core_copd",
                    "analysis_family": family,
                    "metric": metric,
                    "metric_label": label,
                    "units": units,
                    "statistic": statistic,
                    "n_copd_genes": len(cases),
                    "direction": positive_direction if test["difference"] > 0 else negative_direction,
                    **test,
                    "permutations": permutations,
                    "seed": seed,
                    "null_model": "sample without replacement within exact GENCODE gene_type",
                }
            )
    result = pd.DataFrame(rows)
    result["p_two_sided_bh_within_family"] = bh_adjust(result["p_two_sided"])
    result["significant_fdr05"] = result["p_two_sided_bh_within_family"] <= 0.05
    return result


def chromosome_summary(
    cases: pd.DataFrame,
    null_chrom: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_counts = np.array(
        [int(cases["chromosome"].eq(chrom).sum()) for chrom in CANONICAL]
    )
    expected = null_chrom.mean(axis=0)
    rows: list[dict[str, object]] = []
    for index, chrom in enumerate(CANONICAL):
        test = empirical_test(float(observed_counts[index]), null_chrom[:, index].astype(float))
        rows.append(
            {
                "phenotype_scope": "core_copd",
                "chromosome": chrom,
                "observed_genes": observed_counts[index],
                "expected_genes": expected[index],
                "observed_minus_expected": observed_counts[index] - expected[index],
                "enrichment_ratio": (
                    observed_counts[index] / expected[index] if expected[index] > 0 else math.nan
                ),
                "null_count_ci025": test["null_ci025"],
                "null_count_ci975": test["null_ci975"],
                "p_two_sided": test["p_two_sided"],
                "p_enrichment": test["p_greater"],
                "p_depletion": test["p_less"],
                "permutations": permutations,
                "seed": seed,
            }
        )
    per_chrom = pd.DataFrame(rows)
    per_chrom["p_two_sided_bh_24_chromosomes"] = bh_adjust(per_chrom["p_two_sided"])
    per_chrom["significant_fdr05"] = per_chrom["p_two_sided_bh_24_chromosomes"] <= 0.05

    observed_stat = float(((observed_counts - expected) ** 2 / expected).sum())
    null_stats = ((null_chrom - expected) ** 2 / expected).sum(axis=1)
    global_p = (1 + int((null_stats >= observed_stat).sum())) / (permutations + 1)
    global_result = pd.DataFrame(
        [
            {
                "phenotype_scope": "core_copd",
                "test": "stratified Monte Carlo Pearson goodness-of-fit",
                "n_copd_genes": len(cases),
                "chromosomes_tested": len(CANONICAL),
                "statistic": observed_stat,
                "cramers_v_like_sqrt_chisq_over_n": math.sqrt(observed_stat / len(cases)),
                "empirical_p": global_p,
                "permutations": permutations,
                "seed": seed,
                "expected_distribution": (
                    "GENCODE v50 canonical-chromosome frequencies conditional on the "
                    "observed exact gene_type counts"
                ),
                "interpretation_scope": (
                    "chromosome-level over/underrepresentation; not within-chromosome distance clustering"
                ),
            }
        ]
    )
    return per_chrom, global_result


def percentile_within_type(
    cases: pd.DataFrame, background: pd.DataFrame, metric: str
) -> pd.Series:
    values: list[float] = []
    for row in cases.itertuples(index=False):
        pool = background.loc[
            background["gene_type"].eq(row.gene_type), metric
        ].to_numpy(dtype=float)
        value = float(getattr(row, metric))
        values.append((float((pool < value).sum()) + 0.5 * float((pool == value).sum())) / len(pool))
    return pd.Series(values, index=cases.index)


def main() -> int:
    args = parse_args()
    if args.flank_bp < 0 or args.permutations < 999:
        raise ValueError("flank must be nonnegative and permutations at least 999")
    logger = setup_logging(args.log)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    bigwig_info = resolve_executable(args.bigwig_info)
    bigwig_average = resolve_executable(args.bigwig_average_over_bed)

    locus = pd.read_csv(args.locus_table, sep="\t", low_memory=False)
    require_columns(
        locus,
        [
            "gene",
            "selection_reason",
            "gencode_record_count",
            "gencode_gene_id",
            "chromosome",
            "start_1based",
            "end_1based",
            "gene_type",
        ],
        "locus table",
    )
    replicated = locus["selection_reason"].str.contains(
        "gws_in_at_least_2_study_accessions", regex=False, na=False
    )
    exclusions = locus.loc[replicated & locus["gencode_record_count"].ne(1)].copy()
    exclusions = exclusions[
        ["gene", "gencode_record_count", "gencode_gene_id", "chromosome", "gene_type"]
    ]
    exclusions.insert(0, "phenotype_scope", "core_copd")
    exclusions["exclusion_reason"] = "replicated gene name maps to multiple GENCODE v50 records"
    case_locus = locus.loc[replicated & locus["gencode_record_count"].eq(1)].copy()
    if len(case_locus) != 140 or len(exclusions) != 3:
        raise AssertionError(
            f"expected 140 unambiguous and 3 ambiguous replicated genes, observed "
            f"{len(case_locus)} and {len(exclusions)}"
        )

    gencode = load_gencode(args.gencode)
    name_counts = gencode.groupby("gene")["gencode_gene_id"].transform("nunique")
    unique_gencode = gencode.loc[name_counts.eq(1)].copy()
    case_types = sorted(set(case_locus["gene_type"]))
    universe = unique_gencode.loc[unique_gencode["gene_type"].isin(case_types)].copy()
    case_ids = set(case_locus["gencode_gene_id"])
    cases = universe.loc[universe["gencode_gene_id"].isin(case_ids)].copy()
    if set(cases["gencode_gene_id"]) != case_ids or len(cases) != 140:
        raise AssertionError("unambiguous COPD genes do not resolve uniquely in parsed GENCODE")

    # Verify that the precomputed locus coordinates agree exactly with GENCODE.
    verification = case_locus.merge(
        cases,
        on="gencode_gene_id",
        suffixes=("_locus", "_gtf"),
        validate="one_to_one",
    )
    for field in ("gene", "gene_type", "chromosome", "start_1based", "end_1based"):
        if not verification[f"{field}_locus"].astype(str).equals(
            verification[f"{field}_gtf"].astype(str)
        ):
            raise AssertionError(f"locus table disagrees with GENCODE for {field}")

    lengths = chromosome_lengths(args.phylop, bigwig_info)
    cases = add_operational_loci(cases, lengths, args.flank_bp)
    universe = add_operational_loci(universe, lengths, args.flank_bp)
    background = universe.loc[~universe["gencode_gene_id"].isin(case_ids)].copy()

    bed_path = args.work_dir / "gencode_v50_comparable_unique_gene_bodies.bed"
    write_gene_bed(universe, bed_path)
    phy_raw = args.work_dir / "gencode_v50_phyloP100way.tsv"
    phast_raw = args.work_dir / "gencode_v50_phastCons100way.tsv"
    logger.info("Scoring %d comparable unique GENCODE genes with phyloP", len(universe))
    phy = score_bigwig(args.phylop, bed_path, phy_raw, bigwig_average, "phylop100way")
    logger.info("Scoring %d comparable unique GENCODE genes with phastCons", len(universe))
    phast = score_bigwig(
        args.phastcons, bed_path, phast_raw, bigwig_average, "phastcons100way"
    )
    universe = universe.merge(phy, on="gencode_gene_id", validate="one_to_one")
    universe = universe.merge(phast, on="gencode_gene_id", validate="one_to_one")
    cases = universe.loc[universe["gencode_gene_id"].isin(case_ids)].copy()
    background = universe.loc[~universe["gencode_gene_id"].isin(case_ids)].copy()

    metrics = [
        "gene_body_length_bp",
        "operational_locus_length_bp",
        "phylop100way_mean0",
        "phastcons100way_mean0",
    ]
    logger.info(
        "Running %d gene-type-stratified randomizations with seed %d",
        args.permutations,
        args.seed,
    )
    nulls, null_chrom = randomization(
        cases, background, metrics, args.permutations, args.seed
    )
    length_summary = metric_summary(
        cases,
        nulls,
        [
            ("gene_body_length_bp", "GENCODE gene body length", "bp"),
            (
                "operational_locus_length_bp",
                f"gene body plus {args.flank_bp} bp per side, chromosome-clipped",
                "bp",
            ),
        ],
        "locus_length",
        args.permutations,
        args.seed,
    )
    conservation_summary = metric_summary(
        cases,
        nulls,
        [
            ("phylop100way_mean0", "UCSC hg38 phyloP 100-way gene-body mean0", "phyloP"),
            (
                "phastcons100way_mean0",
                "UCSC hg38 phastCons 100-way gene-body mean0",
                "probability",
            ),
        ],
        "evolutionary_conservation",
        args.permutations,
        args.seed,
    )
    chrom_summary, chrom_global = chromosome_summary(
        cases, null_chrom, args.permutations, args.seed
    )

    for metric in metrics:
        cases[f"{metric}_percentile_within_gene_type_background"] = percentile_within_type(
            cases, background, metric
        )
    cases.insert(0, "phenotype_scope", "core_copd")
    cases["replicated_gene_definition"] = "GWS mapped gene in at least 2 GWAS Catalog study accessions"
    cases["genome_build"] = "GRCh38"
    cases["gencode_version"] = "v50"

    background_summary_rows: list[dict[str, object]] = []
    for gene_type, case_count in cases["gene_type"].value_counts().sort_index().items():
        pool = background.loc[background["gene_type"].eq(gene_type)]
        background_summary_rows.append(
            {
                "phenotype_scope": "core_copd",
                "gene_type": gene_type,
                "n_copd_genes": case_count,
                "n_background_genes": len(pool),
                "background_definition": (
                    "unique gene_name, canonical chromosome, same GENCODE v50 gene_type, "
                    "excluding 140 COPD genes"
                ),
            }
        )
    background_summary = pd.DataFrame(background_summary_rows)

    exclusions_path = args.results_dir / "COPD-S2-R007_excluded_ambiguous_replicated_genes.tsv"
    per_gene_path = args.results_dir / "COPD-S2-R007A_locus_conservation_per_gene.tsv"
    background_path = args.results_dir / "COPD-S2-R007A_background_by_gene_type.tsv"
    length_path = args.results_dir / "COPD-S2-R007A_length_tests.tsv"
    chrom_path = args.results_dir / "COPD-S2-R007B_chromosome_tests.tsv"
    chrom_global_path = args.results_dir / "COPD-S2-R007B_chromosome_global_test.tsv"
    conservation_path = args.results_dir / "COPD-S2-R007C_conservation_tests.tsv"
    exclusions.to_csv(exclusions_path, sep="\t", index=False)
    cases.sort_values(
        ["chromosome", "start_1based"],
        key=lambda series: series.map(CHROM_RANK) if series.name == "chromosome" else series,
        kind="stable",
    ).to_csv(per_gene_path, sep="\t", index=False)
    background_summary.to_csv(background_path, sep="\t", index=False)
    length_summary.to_csv(length_path, sep="\t", index=False)
    chrom_summary.to_csv(chrom_path, sep="\t", index=False)
    chrom_global.to_csv(chrom_global_path, sep="\t", index=False)
    conservation_summary.to_csv(conservation_path, sep="\t", index=False)

    output_rows = {
        exclusions_path: len(exclusions),
        per_gene_path: len(cases),
        background_path: len(background_summary),
        length_path: len(length_summary),
        chrom_path: len(chrom_summary),
        chrom_global_path: len(chrom_global),
        conservation_path: len(conservation_summary),
    }
    manifest_rows: list[tuple[str, object]] = [
        ("run_utc", datetime.now(timezone.utc).isoformat()),
        ("python", platform.python_version()),
        ("numpy", np.__version__),
        ("pandas", pd.__version__),
        ("script", str(SCRIPT)),
        ("script_sha256", sha256_file(SCRIPT)),
        ("locus_table", str(args.locus_table.resolve())),
        ("locus_table_sha256", sha256_file(args.locus_table)),
        ("gencode", str(args.gencode.resolve())),
        ("gencode_sha256", sha256_file(args.gencode)),
        ("gencode_version", "v50 / Ensembl 116 / GRCh38"),
        ("phylop", str(args.phylop.resolve())),
        ("phylop_size_bytes", args.phylop.stat().st_size),
        ("phastcons", str(args.phastcons.resolve())),
        ("phastcons_size_bytes", args.phastcons.stat().st_size),
        ("bigWigInfo", bigwig_info),
        ("bigWigAverageOverBed", bigwig_average),
        ("n_replicated_input_genes", int(replicated.sum())),
        ("n_unambiguous_replicated_genes", len(cases)),
        ("n_ambiguous_replicated_genes_excluded", len(exclusions)),
        ("n_comparable_background_genes", len(background)),
        ("flank_bp_each_side", args.flank_bp),
        ("permutations", args.permutations),
        ("seed", args.seed),
    ]
    for path, row_count in output_rows.items():
        manifest_rows.extend(
            [
                (f"output_{path.name}_rows", row_count),
                (f"output_{path.name}_sha256", sha256_file(path)),
            ]
        )
    pd.DataFrame(manifest_rows, columns=["field", "value"]).to_csv(
        LOGS / "COPD-S2-locus_properties_manifest.tsv", sep="\t", index=False
    )

    primary_lengths = length_summary.loc[length_summary["statistic"].eq("mean")]
    for row in primary_lengths.itertuples(index=False):
        logger.info(
            "%s: observed %.3f, expected %.3f, ratio %.4f, p=%.6g",
            row.metric,
            row.observed,
            row.null_expected,
            row.ratio,
            row.p_two_sided,
        )
    logger.info(
        "Chromosome global test: statistic %.4f, empirical p %.6g",
        chrom_global.iloc[0]["statistic"],
        chrom_global.iloc[0]["empirical_p"],
    )
    for row in conservation_summary.loc[
        conservation_summary["statistic"].eq("mean")
    ].itertuples(index=False):
        logger.info(
            "%s: observed %.6f, expected %.6f, ratio %.4f, p=%.6g",
            row.metric,
            row.observed,
            row.null_expected,
            row.ratio,
            row.p_two_sided,
        )
    logger.info("Wrote COPD-S2-R007 outputs to %s", args.results_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
