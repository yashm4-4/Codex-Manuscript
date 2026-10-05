#!/usr/bin/env python3
"""Integrate COPD TREDNet scores, annotations, LD links, and gene loci.

A model-specific predicted causal regulatory call requires all of the following:
the candidate was successfully scored, it is outside the ENCODE blacklist, its
maximum REF/ALT score passes the model's held-out 5% FPR threshold, and its
absolute allele delta is at or above the 95th percentile among eligible COPD
candidates of the same broad variant class (SNV versus indel/complex).
"""

from __future__ import annotations

import hashlib
import json
import platform
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
COPD = ROOT / "diseases" / "COPD"
SECTION = COPD / "04_modeling"
RESULTS = SECTION / "results"

CLASSIFICATION = (
    COPD
    / "03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz"
)
SEQUENCE_AUDIT = RESULTS / "COPD-S4-R002_candidate_sequence_audit.tsv.gz"
SCORES = RESULTS / "COPD-S4-R003_candidate_allele_scores.tsv.gz"
MODEL_THRESHOLDS = RESULTS / "COPD-S4-R001_test_threshold_performance.tsv"
TAGS = COPD / "02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv"
LOCI = COPD / "03_regulatory_landscape/results/COPD-S3-R002_gene_loci_grch38.tsv"

OUT_CANDIDATES = RESULTS / "COPD-S4-R004_prioritized_candidates.tsv.gz"
OUT_CAUSAL = RESULTS / "COPD-S4-R004_predicted_causal_regulatory_variants.tsv"
OUT_SUMMARY = RESULTS / "COPD-S4-R004_variant_summary.tsv"
OUT_DELTA = RESULTS / "COPD-S4-R004_delta_thresholds.tsv"
OUT_LOCI = RESULTS / "COPD-S4-R004_locus_classes.tsv"
OUT_LOCUS_SUMMARY = RESULTS / "COPD-S4-R004_locus_summary.tsv"
OUT_MANIFEST = RESULTS / "COPD-S4-R004_analysis_manifest.json"

MODEL_TYPES = ("enhancer", "silencer")
DELTA_PERCENTILE = 95.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return sorted({part.strip() for part in str(value).split(";") if part.strip()})


def combine_tokens(*values: object) -> str:
    combined: set[str] = set()
    for value in values:
        combined.update(tokens(value))
    return ";".join(sorted(combined))


def load_tag_maps() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    tags = pd.read_csv(TAGS, sep="\t", dtype=str)
    gene_map: dict[str, set[str]] = defaultdict(set)
    study_map: dict[str, set[str]] = defaultdict(set)
    for row in tags.itertuples(index=False):
        tag = str(row.normalized_variant_id)
        gene_map[tag].update(tokens(row.mapped_genes))
        study_map[tag].update(tokens(row.study_accessions))
    return dict(gene_map), dict(study_map)


def mapped_values(tag_string: str, mapping: dict[str, set[str]]) -> str:
    values: set[str] = set()
    for tag in tokens(tag_string):
        values.update(mapping.get(tag, set()))
    return ";".join(sorted(values))


def model_thresholds() -> dict[str, float]:
    table = pd.read_csv(MODEL_THRESHOLDS, sep="\t")
    selected = table[table["target_fpr_percent"] == 5]
    values = dict(zip(selected["model_type"], selected["score_threshold"]))
    if set(values) != set(MODEL_TYPES):
        raise ValueError(f"missing 5% FPR threshold: {values}")
    return {key: float(value) for key, value in values.items()}


def classify_variants(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    thresholds = model_thresholds()
    candidates["variant_class_group"] = np.where(
        candidates["variant_class"].eq("SNV"), "SNV", "indel_or_complex"
    )
    candidates["scored"] = candidates["enhancer_ref_score"].notna() & candidates[
        "silencer_ref_score"
    ].notna()
    candidates["blacklisted"] = as_bool(candidates["encode_blacklist"])
    candidates["causal_call_eligible"] = candidates["scored"] & ~candidates["blacklisted"]

    cutoff_rows = []
    for model_type in MODEL_TYPES:
        ref_col = f"{model_type}_ref_score"
        alt_col = f"{model_type}_alt_score"
        delta_col = f"{model_type}_delta_alt_minus_ref"
        region_col = f"{model_type}_region_score"
        absolute_col = f"{model_type}_abs_delta"
        percentile_col = f"{model_type}_abs_delta_percentile_within_class"
        predicted_col = f"predicted_{model_type}_fpr5"
        causal_col = f"predicted_causal_{model_type}"

        candidates[region_col] = candidates[[ref_col, alt_col]].max(axis=1)
        # Recompute rather than trusting a serialized difference.
        candidates[delta_col] = candidates[alt_col] - candidates[ref_col]
        candidates[absolute_col] = candidates[delta_col].abs()
        candidates[predicted_col] = candidates[region_col].ge(thresholds[model_type])
        candidates[percentile_col] = np.nan
        candidates[causal_col] = False

        for variant_group in ("SNV", "indel_or_complex"):
            eligible = candidates["causal_call_eligible"] & candidates[
                "variant_class_group"
            ].eq(variant_group)
            deltas = candidates.loc[eligible, absolute_col].dropna()
            if deltas.empty:
                cutoff = np.nan
            else:
                cutoff = float(np.percentile(deltas.to_numpy(), DELTA_PERCENTILE))
                candidates.loc[eligible, percentile_col] = (
                    candidates.loc[eligible, absolute_col].rank(method="average", pct=True)
                    * 100.0
                )
                candidates.loc[eligible, causal_col] = (
                    candidates.loc[eligible, predicted_col]
                    & candidates.loc[eligible, absolute_col].ge(cutoff)
                )
            cutoff_rows.append(
                {
                    "model_type": model_type,
                    "variant_class_group": variant_group,
                    "n_eligible_background": int(eligible.sum()),
                    "delta_percentile": DELTA_PERCENTILE,
                    "abs_delta_cutoff": cutoff,
                    "held_out_target_fpr_percent": 5,
                    "region_score_cutoff": thresholds[model_type],
                    "n_predicted_regions": int((eligible & candidates[predicted_col]).sum()),
                    "n_predicted_causal": int(candidates.loc[eligible, causal_col].sum()),
                }
            )

        candidates[f"{model_type}_allelic_direction"] = np.select(
            [candidates[delta_col] > 0, candidates[delta_col] < 0],
            ["gain", "loss"],
            default="no_change_or_unscored",
        )
        candidates[f"{model_type}_threshold_ratio"] = (
            candidates[region_col] / thresholds[model_type]
        )

    candidates["predicted_causal_regulatory"] = candidates[
        [f"predicted_causal_{model_type}" for model_type in MODEL_TYPES]
    ].any(axis=1)
    candidates["causal_model_count"] = candidates[
        [f"predicted_causal_{model_type}" for model_type in MODEL_TYPES]
    ].sum(axis=1)
    candidates["predicted_regulatory_fpr5"] = candidates[
        [f"predicted_{model_type}_fpr5" for model_type in MODEL_TYPES]
    ].any(axis=1)
    candidates["observed_donor_regulatory"] = as_bool(candidates["donor_refined_any"])
    candidates["observed_known_regulatory"] = as_bool(candidates["known_regulatory_any"])
    candidates["observed_regulatory_any"] = (
        candidates["observed_donor_regulatory"]
        | candidates["observed_known_regulatory"]
    )
    candidates["coding"] = as_bool(candidates["coding_CDS"])

    candidates["max_abs_delta_percentile_within_class"] = candidates[
        [f"{model_type}_abs_delta_percentile_within_class" for model_type in MODEL_TYPES]
    ].max(axis=1)
    candidates["max_region_threshold_ratio"] = candidates[
        [f"{model_type}_threshold_ratio" for model_type in MODEL_TYPES]
    ].max(axis=1)
    candidates["priority_tier"] = np.select(
        [
            candidates["causal_model_count"].eq(2),
            candidates["predicted_causal_regulatory"]
            & candidates["observed_donor_regulatory"],
            candidates["predicted_causal_regulatory"]
            & candidates["observed_known_regulatory"],
            candidates["predicted_causal_regulatory"],
            candidates["predicted_regulatory_fpr5"],
            candidates["observed_regulatory_any"],
            candidates["scored"],
        ],
        [1, 2, 3, 4, 5, 6, 7],
        default=8,
    ).astype(int)
    candidates["priority_tier_definition"] = candidates["priority_tier"].map(
        {
            1: "predicted causal in enhancer and silencer models",
            2: "predicted causal in one model with donor refined-element overlap",
            3: "predicted causal in one model with other known regulatory overlap",
            4: "predicted causal in one model without observed regulatory overlap",
            5: "5% FPR predicted region without top-5% matched-class allele delta",
            6: "observed regulatory overlap without a 5% FPR model call",
            7: "scored without regulatory prioritization evidence",
            8: "unscored",
        }
    )
    return candidates, pd.DataFrame(cutoff_rows)


def assign_loci(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    loci = pd.read_csv(LOCI, sep="\t")
    loci["is_replicated"] = as_bool(loci["is_replicated"])
    selected_genes = set(loci["gene"])

    tag_gene_sets = candidates["linked_gwas_genes"].map(lambda value: set(tokens(value)))
    assignments: list[set[str]] = []
    assignment_sources: list[dict[str, set[str]]] = []
    loci_by_chrom = {
        chrom: frame for chrom, frame in loci.groupby("chromosome", sort=False)
    }
    for index, row in candidates.iterrows():
        by_tag = tag_gene_sets.loc[index] & selected_genes
        by_coordinate: set[str] = set()
        if pd.isna(row.chromosome_grch38):
            chrom = ""
        else:
            chrom_value = str(row.chromosome_grch38).removeprefix("chr")
            if chrom_value.endswith(".0") and chrom_value[:-2].isdigit():
                chrom_value = chrom_value[:-2]
            chrom = f"chr{chrom_value}"
        if chrom in loci_by_chrom and pd.notna(row.position_grch38):
            position_zero_based = int(row.position_grch38) - 1
            frame = loci_by_chrom[chrom]
            overlaps = frame[
                frame["locus_start_0based"].le(position_zero_based)
                & frame["locus_end_0based_exclusive"].gt(position_zero_based)
            ]
            by_coordinate.update(overlaps["gene"].astype(str))
        assignments.append(by_tag | by_coordinate)
        assignment_sources.append({"tag": by_tag, "coordinate": by_coordinate})

    candidates["assigned_selected_loci"] = [
        ";".join(sorted(values)) for values in assignments
    ]
    candidates["n_assigned_selected_loci"] = [len(values) for values in assignments]

    locus_rows = []
    for locus in loci.itertuples(index=False):
        by_tag = np.fromiter(
            (locus.gene in value["tag"] for value in assignment_sources), dtype=bool
        )
        by_coordinate = np.fromiter(
            (locus.gene in value["coordinate"] for value in assignment_sources), dtype=bool
        )
        assigned = by_tag | by_coordinate
        subset = candidates.loc[assigned]
        has_coding = bool(subset["coding"].any())
        has_causal = bool(subset["predicted_causal_regulatory"].any())
        has_observed = bool(subset["observed_regulatory_any"].any())
        has_predicted = bool(subset["predicted_regulatory_fpr5"].any())
        if has_coding and has_causal:
            category = "coding + predicted causal regulatory"
        elif has_coding:
            category = "coding only"
        elif has_causal:
            category = "predicted causal regulatory only"
        else:
            category = "other"
        causal_ids = subset.loc[
            subset["predicted_causal_regulatory"], "candidate_record_id"
        ].astype(str)
        locus_rows.append(
            {
                "gene": locus.gene,
                "gencode_gene_id": locus.gencode_gene_id,
                "gene_type": locus.gene_type,
                "selection_reason": locus.selection_reason,
                "is_replicated": bool(locus.is_replicated),
                "chromosome": locus.chromosome,
                "locus_start_0based": locus.locus_start_0based,
                "locus_end_0based_exclusive": locus.locus_end_0based_exclusive,
                "n_assigned_candidates": int(assigned.sum()),
                "n_tag_mapped_candidates": int(by_tag.sum()),
                "n_coordinate_mapped_candidates": int(by_coordinate.sum()),
                "n_coding_candidates": int(subset["coding"].sum()),
                "n_observed_regulatory_candidates": int(
                    subset["observed_regulatory_any"].sum()
                ),
                "n_predicted_regulatory_fpr5": int(
                    subset["predicted_regulatory_fpr5"].sum()
                ),
                "n_predicted_causal_regulatory": int(
                    subset["predicted_causal_regulatory"].sum()
                ),
                "predicted_causal_variant_ids": ";".join(causal_ids),
                "has_coding_variant": has_coding,
                "has_observed_regulatory_variant": has_observed,
                "has_predicted_regulatory_fpr5": has_predicted,
                "has_any_regulatory_variant": has_observed or has_predicted,
                "locus_category": category,
                "alternative_context_assessment": (
                    "additional lung cell types, donors, and disease stages should be evaluated"
                    if category == "other"
                    else "not triggered by the prespecified locus rule"
                ),
            }
        )
    return candidates, pd.DataFrame(locus_rows)


def variant_summary(candidates: pd.DataFrame) -> pd.DataFrame:
    n_all = len(candidates)
    rows = []

    def add(metric: str, mask: pd.Series, denominator: int = n_all, note: str = "") -> None:
        count = int(mask.sum())
        rows.append(
            {
                "metric": metric,
                "n": count,
                "denominator": denominator,
                "fraction": count / denominator if denominator else np.nan,
                "note": note,
            }
        )

    add("all_candidate_records", pd.Series(True, index=candidates.index))
    add("scored_candidate_records", candidates["scored"])
    add("unscored_candidate_records", ~candidates["scored"])
    add("encode_blacklist_overlap", candidates["blacklisted"])
    add("coding_CDS", candidates["coding"])
    add("observed_donor_refined_regulatory", candidates["observed_donor_regulatory"])
    add("observed_known_regulatory", candidates["observed_known_regulatory"])
    for model_type in MODEL_TYPES:
        add(f"predicted_{model_type}_fpr5", candidates[f"predicted_{model_type}_fpr5"])
        add(
            f"predicted_causal_{model_type}",
            candidates[f"predicted_causal_{model_type}"],
            note="5% held-out FPR region threshold plus matched-class top-5% absolute allele delta; blacklist excluded",
        )
    add("predicted_regulatory_fpr5_union", candidates["predicted_regulatory_fpr5"])
    add(
        "predicted_causal_regulatory_union",
        candidates["predicted_causal_regulatory"],
        note="union of enhancer and silencer predicted causal calls",
    )
    add("predicted_causal_both_models", candidates["causal_model_count"].eq(2))
    for group in ("SNV", "indel_or_complex"):
        add(
            f"predicted_causal_{group}",
            candidates["predicted_causal_regulatory"]
            & candidates["variant_class_group"].eq(group),
        )
    add(
        "predicted_causal_gws_tag",
        candidates["predicted_causal_regulatory"] & as_bool(candidates["is_gws_tag"]),
    )
    add(
        "predicted_causal_ld_proxy",
        candidates["predicted_causal_regulatory"] & as_bool(candidates["is_ld_proxy"]),
    )
    return pd.DataFrame(rows)


def locus_summary(locus_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for set_name, frame in (
        ("replicated_gwas_genes", locus_table[locus_table["is_replicated"]]),
        ("all_selected_gwas_genes", locus_table),
    ):
        denominator = len(frame)
        definitions = {
            "genes_with_predicted_causal_regulatory": frame[
                "n_predicted_causal_regulatory"
            ].gt(0),
            "genes_with_no_regulatory_variant": ~frame["has_any_regulatory_variant"],
            "coding_plus_predicted_causal_regulatory": frame["locus_category"].eq(
                "coding + predicted causal regulatory"
            ),
            "coding_only": frame["locus_category"].eq("coding only"),
            "predicted_causal_regulatory_only": frame["locus_category"].eq(
                "predicted causal regulatory only"
            ),
            "other": frame["locus_category"].eq("other"),
        }
        for metric, mask in definitions.items():
            count = int(mask.sum())
            rows.append(
                {
                    "analysis_set": set_name,
                    "metric": metric,
                    "n": count,
                    "denominator": denominator,
                    "fraction": count / denominator if denominator else np.nan,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    classification = pd.read_csv(CLASSIFICATION, sep="\t", low_memory=False)
    audit = pd.read_csv(
        SEQUENCE_AUDIT,
        sep="\t",
        usecols=["candidate_record_id", "sequence_status"],
    )
    scores = pd.read_csv(SCORES, sep="\t")
    if classification["candidate_record_id"].duplicated().any():
        raise ValueError("candidate classification contains duplicate record IDs")
    if scores["candidate_record_id"].duplicated().any():
        raise ValueError("candidate score output contains duplicate record IDs")

    candidates = classification.merge(audit, on="candidate_record_id", how="left").merge(
        scores, on="candidate_record_id", how="left", validate="one_to_one"
    )
    if len(candidates) != len(classification):
        raise ValueError("candidate count changed during score integration")

    gene_map, study_map = load_tag_maps()
    candidates["linked_gwas_tag_ids"] = candidates.apply(
        lambda row: combine_tokens(row["gws_tag_ids"], row["source_focal_tags"]), axis=1
    )
    candidates["n_linked_gwas_tags"] = candidates["linked_gwas_tag_ids"].map(
        lambda value: len(tokens(value))
    )
    candidates["linked_gwas_genes"] = candidates["linked_gwas_tag_ids"].map(
        lambda value: mapped_values(value, gene_map)
    )
    candidates["linked_study_accessions"] = candidates["linked_gwas_tag_ids"].map(
        lambda value: mapped_values(value, study_map)
    )

    candidates, delta_table = classify_variants(candidates)
    candidates, locus_table = assign_loci(candidates)
    summary = variant_summary(candidates)
    locus_summary_table = locus_summary(locus_table)

    candidates = candidates.sort_values(
        [
            "priority_tier",
            "causal_model_count",
            "max_abs_delta_percentile_within_class",
            "max_region_threshold_ratio",
            "max_r2_across_links",
            "candidate_record_id",
        ],
        ascending=[True, False, False, False, False, True],
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
    candidates.insert(0, "priority_rank_all_candidates", np.arange(1, len(candidates) + 1))
    causal = candidates[candidates["predicted_causal_regulatory"]].copy()
    causal.insert(0, "predicted_causal_priority_rank", np.arange(1, len(causal) + 1))

    candidates.to_csv(
        OUT_CANDIDATES,
        sep="\t",
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    causal.to_csv(OUT_CAUSAL, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    delta_table.to_csv(OUT_DELTA, sep="\t", index=False)
    locus_table.to_csv(OUT_LOCI, sep="\t", index=False)
    locus_summary_table.to_csv(OUT_LOCUS_SUMMARY, sep="\t", index=False)

    inputs = {
        "candidate_classification": CLASSIFICATION,
        "sequence_audit": SEQUENCE_AUDIT,
        "allele_scores": SCORES,
        "model_thresholds": MODEL_THRESHOLDS,
        "gwas_tags": TAGS,
        "gene_loci": LOCI,
    }
    outputs = {
        "prioritized_candidates": OUT_CANDIDATES,
        "predicted_causal_variants": OUT_CAUSAL,
        "variant_summary": OUT_SUMMARY,
        "delta_thresholds": OUT_DELTA,
        "locus_classes": OUT_LOCI,
        "locus_summary": OUT_LOCUS_SUMMARY,
    }
    manifest = {
        "result_id": "COPD-S4-R004",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "causal_call_definition": {
            "region_threshold": "maximum of REF and ALT scores at held-out 5% FPR",
            "allelic_threshold": "absolute ALT-minus-REF score at or above the 95th percentile within SNV or indel/complex COPD candidates",
            "quality_filter": "successfully scored and no ENCODE blacklist v2 overlap",
            "delta_percentile": DELTA_PERCENTILE,
        },
        "locus_assignment": "union of exact GWAS Catalog mapped-gene links from associated tags and physical overlap with gene body plus 100 kb on each side",
        "candidate_records": len(candidates),
        "scored_candidate_records": int(candidates["scored"].sum()),
        "predicted_causal_records": int(candidates["predicted_causal_regulatory"].sum()),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in inputs.items()
        },
        "outputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in outputs.items()
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"Prioritized {len(candidates):,} candidates; "
        f"{len(causal):,} met the predicted causal regulatory definition."
    )


if __name__ == "__main__":
    main()
