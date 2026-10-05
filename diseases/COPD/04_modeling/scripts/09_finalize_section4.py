#!/usr/bin/env python3
"""Create the auditable COPD experimental-recovery audit and THE LIST.

R009 compares the eight literature evidence rows from Section 3 with exact
variant identities in the frozen candidate set.  LD-source labels are shown
separately and never treated as identity matches.  R010 joins the prespecified
R004 rank with population, provisional target, and allele-specific motif
annotations without changing that rank.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[4]
DISEASE = ROOT / "diseases/COPD"
SECTION = DISEASE / "04_modeling"
RESULTS = SECTION / "results"

POP_ALL = RESULTS / "COPD-S4-R005_population_annotated_candidates.tsv.gz"
POP_CAUSAL = RESULTS / "COPD-S4-R005_predicted_causal_population_genetics.tsv"
TARGETS = RESULTS / "COPD-S4-R006_candidate_target_summary.tsv"
TFBS = RESULTS / "COPD-S4-R008_candidate_tfbs_summary.tsv.gz"
FUNCTIONAL = (
    DISEASE
    / "03_regulatory_landscape/results/COPD-S3-R004_functional_variant_evidence.tsv"
)

OUT_AUDIT = RESULTS / "COPD-S4-R009_literature_functional_variant_recovery.tsv"
OUT_AUDIT_SUMMARY = RESULTS / "COPD-S4-R009_summary.tsv"
OUT_LIST = RESULTS / "COPD-S4-R010_THE_LIST.tsv"
OUT_MANIFEST = RESULTS / "COPD-S4-R010_analysis_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tokens(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {part.strip() for part in str(value).split(";") if part.strip()}


def bool_value(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def exact_identity_tokens(row: pd.Series) -> set[str]:
    result: set[str] = set()
    for column in (
        "gws_tag_ids",
        "ensembl_variation_ids",
        "selected_ensembl_variation_id",
    ):
        result.update(tokens(row.get(column, "")))
    return result


def linked_identity_tokens(row: pd.Series) -> set[str]:
    result: set[str] = set()
    for column in ("source_focal_tags", "linked_gwas_tag_ids"):
        result.update(tokens(row.get(column, "")))
    return result


def miss_reason(frame: pd.DataFrame) -> str:
    if frame.empty:
        return (
            "no exact identity in the frozen candidate set after GWAS-tag and "
            "Ensembl-variation alias resolution; model recovery is not estimable"
        )
    if frame["predicted_causal_regulatory"].map(bool_value).any():
        return "recovered by the prespecified predicted-causal regulatory definition"
    if not frame["scored"].map(bool_value).all():
        return "an exact candidate identity was present but at least one sequence was not scorable"
    enhancer_region = frame["predicted_enhancer_fpr5"].map(bool_value).any()
    silencer_region = frame["predicted_silencer_fpr5"].map(bool_value).any()
    if enhancer_region or silencer_region:
        return (
            "an exact candidate identity passed at least one 5% FPR region threshold "
            "but did not pass that model's matched-class top-5% absolute-delta cutoff"
        )
    return (
        "an exact candidate identity was scorable but both enhancer and silencer "
        "region scores were below their held-out 5% FPR thresholds"
    )


def build_experimental_audit(all_candidates: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    exact_sets = all_candidates.apply(exact_identity_tokens, axis=1)
    linked_sets = all_candidates.apply(linked_identity_tokens, axis=1)
    rows: list[dict[str, object]] = []
    for item in evidence.itertuples(index=False):
        rsids = sorted(set(re.findall(r"rs\d+", str(item.variant_or_haplotype))))
        exact_mask = exact_sets.map(lambda values: bool(values.intersection(rsids)))
        linked_mask = linked_sets.map(lambda values: bool(values.intersection(rsids)))
        exact = all_candidates.loc[exact_mask].sort_values("priority_rank_all_candidates")
        linked = all_candidates.loc[linked_mask].sort_values("priority_rank_all_candidates")
        causal = exact.loc[exact["predicted_causal_regulatory"].map(bool_value)]
        if not causal.empty:
            recovery = "predicted_causal_exact_identity"
        elif not exact.empty:
            recovery = "exact_candidate_not_predicted_causal"
        else:
            recovery = "no_exact_candidate_identity"
        rows.append(
            {
                "variant_evidence_id": item.variant_evidence_id,
                "source_id": item.source_id,
                "variant_or_haplotype": item.variant_or_haplotype,
                "tested_rsids": ";".join(rsids),
                "n_tested_rsids": len(rsids),
                "conclusion_tier": item.conclusion_tier,
                "endogenous_allele_tested": item.endogenous_allele_tested,
                "target_gene": item.target_gene,
                "n_exact_candidate_records": len(exact),
                "exact_candidate_record_ids": ";".join(exact["candidate_record_id"].astype(str)),
                "exact_candidate_priority_ranks": ";".join(
                    exact["priority_rank_all_candidates"].astype(int).astype(str)
                ),
                "n_exact_predicted_causal_records": len(causal),
                "exact_predicted_causal_record_ids": ";".join(
                    causal["candidate_record_id"].astype(str)
                ),
                "n_ld_linked_candidate_records": len(linked),
                "n_ld_linked_predicted_causal_records": int(
                    linked["predicted_causal_regulatory"].map(bool_value).sum()
                ),
                "ld_linked_candidate_record_ids": ";".join(
                    linked["candidate_record_id"].astype(str)
                ),
                "model_recovery_class": recovery,
                "model_recovery_interpretation": miss_reason(exact),
                "identity_rule": (
                    "exact rsID token among GWAS tag IDs or Ensembl variation aliases; "
                    "source-focal and linked-tag labels are reported only as LD context"
                ),
                "literature_caveat": item.caveat,
            }
        )
    return pd.DataFrame(rows)


def build_the_list(causal: pd.DataFrame, targets: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
    if causal["candidate_record_id"].duplicated().any():
        raise ValueError("R005 causal table contains duplicate candidate_record_id values")
    if set(causal["priority_rank_all_candidates"].astype(int)) != set(range(1, 338)):
        raise ValueError("R005 causal ranks are not exactly 1..337")

    target_keep = [
        "candidate_record_id",
        "model_context",
        "nearest_gene",
        "nearest_gene_tss_distance_bp",
        "nearest_protein_coding_gene",
        "nearest_protein_coding_tss_distance_bp",
        "genes_with_tss_within_100kb",
        "target_mapping_limit",
    ]
    result = causal.merge(
        targets[target_keep], on="candidate_record_id", how="left", validate="one_to_one"
    )

    motif_fields = [
        "candidate_record_id",
        "any_motif_compatible_site",
        "disrupts_any_motif_compatible_site",
        "creates_any_motif_compatible_site",
        "retains_any_motif_compatible_site",
        "largest_effect_motif_id",
        "largest_effect_tf_name",
        "largest_signed_delta_relative_score",
        "largest_absolute_delta_relative_score",
    ]
    for model in ("enhancer", "silencer"):
        subset = tfbs.loc[tfbs["model_type"].eq(model), motif_fields].copy()
        subset = subset.rename(
            columns={column: f"{model}_tfbs_{column}" for column in motif_fields if column != "candidate_record_id"}
        )
        result = result.merge(
            subset, on="candidate_record_id", how="left", validate="one_to_one"
        )

    keep = [
        "priority_rank_all_candidates",
        "candidate_record_id",
        "chromosome_grch38",
        "position_grch38",
        "ref",
        "alt",
        "variant_class_group",
        "gws_tag_ids",
        "source_focal_tags",
        "ld_panels",
        "is_gws_tag",
        "is_ld_proxy",
        "max_r2_across_links",
        "linked_gwas_tag_ids",
        "linked_gwas_genes",
        "linked_study_accessions",
        "assigned_selected_loci",
        "exclusive_class_comprehensive",
        "coding_CDS_gene_names",
        "donor_refined_enhancer",
        "donor_refined_silencer",
        "enhancer_ref_score",
        "enhancer_alt_score",
        "enhancer_delta_alt_minus_ref",
        "enhancer_region_score",
        "enhancer_abs_delta_percentile_within_class",
        "predicted_causal_enhancer",
        "enhancer_allelic_direction",
        "silencer_ref_score",
        "silencer_alt_score",
        "silencer_delta_alt_minus_ref",
        "silencer_region_score",
        "silencer_abs_delta_percentile_within_class",
        "predicted_causal_silencer",
        "silencer_allelic_direction",
        "priority_tier_definition",
        "selected_ensembl_variation_id",
        "global_alt_af",
        "AFR_alt_af",
        "AMR_alt_af",
        "EAS_alt_af",
        "EUR_alt_af",
        "SAS_alt_af",
        "frequency_class",
        "strict_population_specific",
        "ancestral_allele",
        "derived_candidate_allele",
        "global_derived_allele_frequency",
        "global_DAF_gt_0_5",
        "model_context",
        "nearest_gene",
        "nearest_gene_tss_distance_bp",
        "nearest_protein_coding_gene",
        "nearest_protein_coding_tss_distance_bp",
        "genes_with_tss_within_100kb",
        "target_mapping_limit",
    ]
    keep.extend([column for column in result.columns if "_tfbs_" in column])
    result = result[keep].rename(
        columns={"priority_rank_all_candidates": "predicted_causal_priority_rank"}
    )
    result["rank_basis"] = (
        "prespecified R004 tier, model count, matched-class absolute-delta percentile, "
        "region-threshold ratio, maximum LD r2, then stable candidate ID"
    )
    result["interpretation_boundary"] = (
        "computational prioritization, motif, population, and proximity evidence; not causal proof"
    )
    result["chromosome_grch38"] = result["chromosome_grch38"].map(
        lambda value: str(int(float(value)))
        if re.fullmatch(r"\d+(?:\.0+)?", str(value).strip())
        else str(value).strip().removeprefix("chr")
    )
    positions = pd.to_numeric(result["position_grch38"], errors="raise")
    if positions.isna().any() or not np.equal(positions, np.floor(positions)).all():
        raise ValueError("THE LIST contains missing or non-integer GRCh38 positions")
    result["position_grch38"] = positions.astype("int64")
    return result.sort_values("predicted_causal_priority_rank").reset_index(drop=True)


def main() -> None:
    all_candidates = pd.read_csv(POP_ALL, sep="\t", low_memory=False)
    causal = pd.read_csv(POP_CAUSAL, sep="\t", low_memory=False)
    targets = pd.read_csv(TARGETS, sep="\t", low_memory=False)
    tfbs = pd.read_csv(TFBS, sep="\t", low_memory=False)
    functional = pd.read_csv(FUNCTIONAL, sep="\t", low_memory=False)

    audit = build_experimental_audit(all_candidates, functional)
    the_list = build_the_list(causal, targets, tfbs)
    if len(audit) != 8 or len(the_list) != 337:
        raise ValueError(f"unexpected output sizes: audit={len(audit)}, list={len(the_list)}")
    if the_list["candidate_record_id"].duplicated().any():
        raise ValueError("THE LIST contains duplicate candidates")

    distinct_rsids = sorted(
        {token for value in audit["tested_rsids"] for token in tokens(value)}
    )
    exact_rsids: set[str] = set()
    causal_rsids: set[str] = set()
    for row in audit.itertuples(index=False):
        rsids = tokens(row.tested_rsids)
        if row.n_exact_candidate_records:
            exact_rsids.update(rsids)
        if row.n_exact_predicted_causal_records:
            causal_rsids.update(rsids)
    summary = pd.DataFrame(
        [
            ("literature_functional_evidence_rows", len(audit), len(audit), 1.0, "Section 3 R004 evidence units"),
            ("distinct_literature_tested_rsids", len(distinct_rsids), len(distinct_rsids), 1.0, ";".join(distinct_rsids)),
            ("distinct_tested_rsids_with_exact_candidate_identity", len(exact_rsids), len(distinct_rsids), len(exact_rsids) / len(distinct_rsids), ";".join(sorted(exact_rsids))),
            ("distinct_tested_rsids_recovered_as_predicted_causal", len(causal_rsids), len(distinct_rsids), len(causal_rsids) / len(distinct_rsids), ";".join(sorted(causal_rsids))),
            ("predicted_causal_candidates_with_exact_literature_tested_identity", audit["exact_predicted_causal_record_ids"].map(tokens).map(len).sum(), len(the_list), audit["exact_predicted_causal_record_ids"].map(tokens).map(len).sum() / len(the_list), "exact identity; excludes LD-source labels"),
        ],
        columns=["metric", "n", "denominator", "fraction", "note"],
    )

    audit.to_csv(OUT_AUDIT, sep="\t", index=False)
    summary.to_csv(OUT_AUDIT_SUMMARY, sep="\t", index=False)
    the_list.to_csv(OUT_LIST, sep="\t", index=False)

    inputs = {
        "population_all": POP_ALL,
        "population_causal": POP_CAUSAL,
        "target_summary": TARGETS,
        "tfbs_summary": TFBS,
        "functional_literature": FUNCTIONAL,
    }
    outputs = {
        "literature_recovery_audit": OUT_AUDIT,
        "literature_recovery_summary": OUT_AUDIT_SUMMARY,
        "the_list": OUT_LIST,
    }
    manifest = {
        "result_ids": ["COPD-S4-R009", "COPD-S4-R010"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "candidate_count": len(the_list),
        "rank_policy": the_list.loc[0, "rank_basis"],
        "literature_identity_policy": audit.loc[0, "identity_rule"],
        "inputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in inputs.items()
        },
        "outputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in outputs.items()
        },
        "script": {
            "path": str(SCRIPT),
            "bytes": SCRIPT.stat().st_size,
            "sha256": sha256(SCRIPT),
        },
        "software": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {len(the_list)} ranked candidates; exact literature identity recovered for "
        f"{len(exact_rsids)}/{len(distinct_rsids)} tested rsIDs and predicted causal for "
        f"{len(causal_rsids)}/{len(distinct_rsids)}."
    )


if __name__ == "__main__":
    main()
