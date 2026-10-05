#!/usr/bin/env python3
"""Build an auditable COPD experimental-validation design and handoff tables.

This script consumes completed Sections 4 and 5.  It proposes experiments; it
does not represent any wet-lab experiment as performed and does not convert a
model prediction, database overlap, or provisional target map into causal proof.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases/COPD/06_experimental_validation"
DATA = SECTION / "data"
RESULTS = SECTION / "results"
LOGS = SECTION / "logs"

S4 = ROOT / "diseases/COPD/04_modeling/results"
S5 = ROOT / "diseases/COPD/05_computational_validation/results"

R004 = S4 / "COPD-S4-R004_predicted_causal_regulatory_variants.tsv"
R005 = S4 / "COPD-S4-R005_predicted_causal_population_genetics.tsv"
R006_SUMMARY = S4 / "COPD-S4-R006_candidate_target_summary.tsv"
R006_EVIDENCE = S4 / "COPD-S4-R006_candidate_target_evidence.tsv"
R008_SUMMARY = S4 / "COPD-S4-R008_candidate_tfbs_summary.tsv.gz"
R008_MOTIFS = S4 / "COPD-S4-R008_allele_specific_motif_scores.tsv.gz"
R010 = S4 / "COPD-S4-R010_THE_LIST.tsv"
S5_INTEGRATED = S5 / "COPD-S5-R005_integrated_candidate_validation.tsv"

REF_FASTA = ROOT / "diseases/COPD/04_modeling/data/COPD_candidate_variants_ref_2001bp.fa"
ALT_FASTA = ROOT / "diseases/COPD/04_modeling/data/COPD_candidate_variants_alt_2001bp.fa"

OUT_SHORTLIST = RESULTS / "COPD-S6-R001_candidate_shortlist.tsv"
OUT_AUDIT = RESULTS / "COPD-S6-R002_selection_audit.tsv.gz"
OUT_CONSTRUCTS = RESULTS / "COPD-S6-R003_MPRA_constructs.tsv"
OUT_CANDIDATE_PLAN = RESULTS / "COPD-S6-R004_candidate_validation_plan.tsv"
OUT_TF_PLAN = RESULTS / "COPD-S6-R004_TF_first_plan.tsv"
OUT_CONTEXTS = RESULTS / "COPD-S6-R004_cell_context_controls.tsv"
OUT_UPSTREAM = DATA / "COPD-S6-R005_upstream_inputs.tsv"
OUT_COLLAB_README = RESULTS / "COPD-S6-R005_collaborator_README.md"
OUT_REPORT = RESULTS / "section_6_experimental_validation.md"
OUT_MANIFEST = RESULTS / "COPD-S6-R006_analysis_manifest.json"

SHORTLIST_SIZE = 12
MPRA_FLANK_BP = 100
SELECTION_POLICY = "COPD-S6-selection-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def clean(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "<na>"} else text


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return clean(value).lower() in {"true", "1", "yes"}


def tokens(value: object) -> List[str]:
    return sorted({part.strip() for part in clean(value).split(";") if part.strip()})


def joined(values: Iterable[object]) -> str:
    return ";".join(sorted({clean(value) for value in values if clean(value)}))


def write_tsv(frame: pd.DataFrame, path: Path, gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if gzip_output else None
    frame.to_csv(path, sep="\t", index=False, na_rep="", compression=compression)


def require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def load_inputs(require_targets: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    required = [
        R010,
        R004,
        R005,
        R008_SUMMARY,
        R008_MOTIFS,
        S5_INTEGRATED,
        REF_FASTA,
        ALT_FASTA,
    ]
    if require_targets:
        required.extend([R006_SUMMARY, R006_EVIDENCE])
    absent = [str(path) for path in required if not path.exists()]
    if absent:
        raise FileNotFoundError("Missing required Section 4/5 inputs:\n" + "\n".join(absent))

    canonical = pd.read_csv(R010, sep="\t", low_memory=False)
    core = pd.read_csv(R004, sep="\t", low_memory=False)
    population = pd.read_csv(R005, sep="\t", low_memory=False)
    validation = pd.read_csv(S5_INTEGRATED, sep="\t", low_memory=False)
    tf_summary = pd.read_csv(R008_SUMMARY, sep="\t", low_memory=False)

    require_columns(
        canonical,
        [
            "predicted_causal_priority_rank",
            "candidate_record_id",
            "variant_class_group",
            "source_focal_tags",
            "linked_gwas_tag_ids",
            "predicted_causal_enhancer",
            "predicted_causal_silencer",
            "rank_basis",
            "interpretation_boundary",
        ],
        "Section 4 canonical R010",
    )
    require_columns(
        core,
        [
            "predicted_causal_priority_rank",
            "candidate_record_id",
            "source_focal_tags",
            "linked_gwas_tag_ids",
            "sequence_status",
            "blacklisted",
            "predicted_causal_enhancer",
            "predicted_causal_silencer",
        ],
        "Section 4 R004",
    )
    require_columns(population, ["candidate_record_id", "frequency_match_status"], "Section 4 R005")
    require_columns(
        validation,
        [
            "candidate_record_id",
            "gtex_lung_exact_significant_eqtl",
            "mprabase_any_element_evidence",
            "integrated_validation_status",
        ],
        "Section 5 R005",
    )

    canonical["candidate_record_id"] = canonical["candidate_record_id"].astype(str)
    core["candidate_record_id"] = core["candidate_record_id"].astype(str)
    population["candidate_record_id"] = population["candidate_record_id"].astype(str)
    validation["candidate_record_id"] = validation["candidate_record_id"].astype(str)
    tf_summary["candidate_record_id"] = tf_summary["candidate_record_id"].astype(str)

    if canonical["candidate_record_id"].duplicated().any():
        raise ValueError("Section 4 canonical R010 candidate IDs are not unique")
    if core["candidate_record_id"].duplicated().any():
        raise ValueError("Section 4 R004 candidate IDs are not unique")
    canonical_ids = set(canonical["candidate_record_id"])
    if canonical_ids != set(core["candidate_record_id"]):
        raise ValueError("Section 4 canonical R010 and R004 candidate sets differ")
    if canonical_ids != set(population["candidate_record_id"]):
        raise ValueError("Section 4 canonical R010 and R005 candidate sets differ")
    if canonical_ids != set(validation["candidate_record_id"]):
        raise ValueError("Section 4 canonical R010 and Section 5 candidate sets differ")

    canonical_ranks = canonical.set_index("candidate_record_id")[
        "predicted_causal_priority_rank"
    ].astype(int)
    r004_ranks = core.set_index("candidate_record_id")[
        "predicted_causal_priority_rank"
    ].astype(int)
    if not canonical_ranks.sort_index().equals(r004_ranks.sort_index()):
        raise ValueError("R004 ranks do not reconcile with canonical Section 4 R010")

    # R010 is the canonical ranked list. R004-R008 only supplement fields not
    # carried in that presentation table, such as sequence eligibility and
    # full interval annotations.
    canonical = canonical.rename(columns={"variant_class_group": "variant_class"})
    core_only = [
        column
        for column in core.columns
        if column not in canonical.columns and column != "candidate_record_id"
    ]
    combined = canonical.merge(
        core[["candidate_record_id"] + core_only],
        on="candidate_record_id",
        how="left",
        validate="one_to_one",
    )
    combined["canonical_section4_list"] = relative(R010)
    combined["canonical_rank_preserved"] = True
    combined["predicted_causal_regulatory"] = True
    if "donor_refined_any" not in combined.columns:
        combined["donor_refined_any"] = combined["donor_refined_enhancer"].map(as_bool) | combined[
            "donor_refined_silencer"
        ].map(as_bool)

    population_only = [
        column
        for column in population.columns
        if column not in combined.columns and column != "candidate_record_id"
    ]
    combined = combined.merge(
        population[["candidate_record_id"] + population_only],
        on="candidate_record_id",
        how="left",
        validate="one_to_one",
    )

    validation_evidence = [
        "candidate_rsids",
        "gtex_lung_exact_significant_eqtl",
        "gtex_lung_significant_pair_count",
        "gtex_lung_significant_gene_count",
        "gtex_lung_significant_genes",
        "gtex_lung_minimum_nominal_p",
        "gtex_lung_gene_matches_candidate_mapping",
        "gtex_lung_matching_candidate_gene_count",
        "gtex_lung_matching_candidate_genes",
        "gtex_lung_gene_matches_provisional_target",
        "gtex_lung_matching_provisional_targets",
        "mprabase_any_element_evidence",
        "mprabase_coordinate_element_evidence",
        "mprabase_exact_rsid_metadata_evidence",
        "mprabase_matched_element_count",
        "mprabase_elements_with_scores",
        "mprabase_matched_dataset_count",
        "mprabase_matched_datasets",
        "open_targets_COPD_associated_gene_count",
        "open_targets_COPD_associated_genes",
        "provisional_section4_target_count",
        "provisional_section4_targets",
        "variant_or_regional_evidence_channels",
        "gene_context_channels",
        "integrated_validation_status",
        "interpretation_scope",
    ]
    validation_evidence = [column for column in validation_evidence if column in validation.columns]
    combined = combined.merge(
        validation[["candidate_record_id"] + validation_evidence],
        on="candidate_record_id",
        how="left",
        validate="one_to_one",
    )

    target_evidence = pd.DataFrame()
    if R006_SUMMARY.exists() and R006_EVIDENCE.exists():
        target_summary = pd.read_csv(R006_SUMMARY, sep="\t", low_memory=False)
        target_evidence = pd.read_csv(R006_EVIDENCE, sep="\t", low_memory=False)
        target_summary["candidate_record_id"] = target_summary["candidate_record_id"].astype(str)
        target_evidence["candidate_record_id"] = target_evidence["candidate_record_id"].astype(str)
        if set(target_summary["candidate_record_id"]) != canonical_ids:
            raise ValueError("Section 4 R006 target summary candidate set differs from canonical R010")
        grouped_targets: List[Dict[str, object]] = []
        for candidate_id, frame in target_evidence.groupby("candidate_record_id", sort=False):
            method_pairs = joined(
                f"{row.mapping_method}:{row.gene_name}"
                for row in frame[["mapping_method", "gene_name"]].itertuples(index=False)
            )
            grouped_targets.append(
                {
                    "candidate_record_id": candidate_id,
                    "R006_provisional_target_count": int(frame["gene_name"].nunique()),
                    "R006_provisional_target_genes": joined(frame["gene_name"]),
                    "R006_target_mapping_methods": joined(frame["mapping_method"]),
                    "R006_target_method_gene_pairs": method_pairs,
                }
            )
        target_grouped = pd.DataFrame(grouped_targets)
        target_columns = [
            "candidate_record_id",
            "nearest_gene",
            "nearest_gene_tss_distance_bp",
            "nearest_protein_coding_gene",
            "nearest_protein_coding_tss_distance_bp",
            "genes_with_tss_within_100kb",
            "target_mapping_limit",
        ]
        target_additions = [
            column
            for column in target_columns
            if column == "candidate_record_id" or column not in combined.columns
        ]
        if len(target_additions) > 1:
            combined = combined.merge(
                target_summary[target_additions],
                on="candidate_record_id",
                how="left",
                validate="one_to_one",
            )
        combined = combined.merge(
            target_grouped, on="candidate_record_id", how="left", validate="one_to_one"
        )
    elif require_targets:
        raise FileNotFoundError("Section 4 R006 target outputs are required")
    else:
        for column in [
            "nearest_gene",
            "nearest_gene_tss_distance_bp",
            "nearest_protein_coding_gene",
            "nearest_protein_coding_tss_distance_bp",
            "genes_with_tss_within_100kb",
            "target_mapping_limit",
            "R006_provisional_target_count",
            "R006_provisional_target_genes",
            "R006_target_mapping_methods",
            "R006_target_method_gene_pairs",
        ]:
            combined[column] = ""

    for model in ("enhancer", "silencer"):
        model_rows = tf_summary[tf_summary["model_type"] == model].copy()
        if model_rows["candidate_record_id"].duplicated().any():
            raise ValueError(f"R008 has duplicate {model} candidate summaries")
        rename = {
            column: f"{model}_tf_{column}"
            for column in model_rows.columns
            if column not in {"candidate_record_id", "model_type"}
        }
        model_rows = model_rows.drop(columns=["model_type"]).rename(columns=rename)
        combined = combined.merge(model_rows, on="candidate_record_id", how="left", validate="one_to_one")

    def exact_target_match(row: pd.Series) -> str:
        gtex = set(tokens(row.get("gtex_lung_significant_genes", "")))
        r006 = set(tokens(row.get("R006_provisional_target_genes", "")))
        return ";".join(sorted(gtex & r006))

    combined["section6_exact_GTEx_R006_target_genes"] = combined.apply(exact_target_match, axis=1)
    combined["section6_exact_GTEx_R006_target_count"] = combined[
        "section6_exact_GTEx_R006_target_genes"
    ].map(lambda value: len(tokens(value)))
    combined["target_data_status"] = (
        "R006_consumed" if R006_SUMMARY.exists() and R006_EVIDENCE.exists() else "R006_missing"
    )
    return combined, target_evidence


def model_context(row: Mapping[str, object]) -> str:
    enhancer = as_bool(row.get("predicted_causal_enhancer"))
    silencer = as_bool(row.get("predicted_causal_silencer"))
    if enhancer and silencer:
        return "both"
    if enhancer:
        return "enhancer_only"
    if silencer:
        return "silencer_only"
    return "neither"


def annotate_selection(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["candidate_record_id"] = result["candidate_record_id"].astype(str)
    result["predicted_causal_priority_rank"] = result["predicted_causal_priority_rank"].astype(int)
    result["model_context"] = result.apply(lambda row: model_context(row), axis=1)

    eligibility_reasons: List[str] = []
    eligible_values: List[bool] = []
    for row in result.to_dict(orient="records"):
        reasons = []
        if clean(row.get("sequence_status")) != "ok":
            reasons.append(f"sequence_status={clean(row.get('sequence_status')) or 'missing'}")
        if as_bool(row.get("blacklisted")):
            reasons.append("ENCODE_blacklist_overlap")
        if not as_bool(row.get("predicted_causal_regulatory")):
            reasons.append("not_in_Section4_predicted_causal_union")
        eligible_values.append(not reasons)
        eligibility_reasons.append("eligible" if not reasons else ";".join(reasons))
    result["experimental_eligible"] = eligible_values
    result["experimental_eligibility_reason"] = eligibility_reasons

    ids = result["candidate_record_id"].tolist()
    union_find = UnionFind(ids)
    tag_owner: Dict[str, str] = {}
    for row in result.to_dict(orient="records"):
        candidate_id = str(row["candidate_record_id"])
        linked_tags = set(tokens(row.get("linked_gwas_tag_ids")))
        focal_tags = set(tokens(row.get("source_focal_tags")))
        for tag in sorted(linked_tags | focal_tags):
            if tag in tag_owner:
                union_find.union(candidate_id, tag_owner[tag])
            else:
                tag_owner[tag] = candidate_id

    components: Dict[str, List[str]] = defaultdict(list)
    for candidate_id in ids:
        components[union_find.find(candidate_id)].append(candidate_id)
    rank_map = result.set_index("candidate_record_id")["predicted_causal_priority_rank"].to_dict()
    eligible_map = result.set_index("candidate_record_id")["experimental_eligible"].to_dict()

    component_rows: Dict[str, Dict[str, object]] = {}
    for members in components.values():
        ordered = sorted(members, key=lambda item: (int(rank_map[item]), item))
        eligible = [member for member in ordered if bool(eligible_map[member])]
        representative = eligible[0] if eligible else ""
        component_id = f"COPD-LDC-{int(rank_map[ordered[0]]):04d}"
        for index, member in enumerate(ordered, start=1):
            component_rows[member] = {
                "ld_component_id": component_id,
                "ld_component_size": len(ordered),
                "ld_component_members": ";".join(ordered),
                "within_ld_component_priority_rank": index,
                "ld_component_representative": representative,
                "is_ld_component_representative": member == representative,
                "ld_component_definition": (
                    "transitive connection through an exact shared linked_gwas_tag_id "
                    "or source_focal_tag among Section 4 predicted-causal candidates"
                ),
            }
    component_frame = pd.DataFrame.from_dict(component_rows, orient="index")
    component_frame.index.name = "candidate_record_id"
    result = result.merge(component_frame.reset_index(), on="candidate_record_id", validate="one_to_one")

    def motif_change(row: Mapping[str, object]) -> bool:
        return any(
            as_bool(row.get(column))
            for column in [
                "enhancer_tf_disrupts_any_motif_compatible_site",
                "enhancer_tf_creates_any_motif_compatible_site",
                "silencer_tf_disrupts_any_motif_compatible_site",
                "silencer_tf_creates_any_motif_compatible_site",
            ]
        )

    criteria = [
        ("both_model_prediction", lambda row: clean(row.get("model_context")) == "both"),
        ("enhancer_only_prediction", lambda row: clean(row.get("model_context")) == "enhancer_only"),
        ("silencer_only_prediction", lambda row: clean(row.get("model_context")) == "silencer_only"),
        ("indel_or_complex", lambda row: clean(row.get("variant_class")) != "SNV"),
        ("severe_emphysema_donor_refined_overlap", lambda row: as_bool(row.get("donor_refined_any"))),
        ("exact_GTEx_v10_Lung_eQTL", lambda row: as_bool(row.get("gtex_lung_exact_significant_eqtl"))),
        ("public_MPRAbase_element_overlap", lambda row: as_bool(row.get("mprabase_any_element_evidence"))),
        ("predicted_motif_creation_or_disruption", motif_change),
    ]

    eligible_candidates = result[result["experimental_eligible"]].sort_values(
        ["predicted_causal_priority_rank", "candidate_record_id"]
    )
    representatives = eligible_candidates[
        eligible_candidates["is_ld_component_representative"]
    ]
    selected: List[str] = []
    anchor: Dict[str, str] = {}

    def satisfies(candidate_id: str, criterion) -> bool:
        row = result[result["candidate_record_id"] == candidate_id].iloc[0].to_dict()
        return bool(criterion(row))

    for criterion_name, criterion in criteria:
        if any(satisfies(candidate_id, criterion) for candidate_id in selected):
            continue
        represented_components = {
            clean(
                result.loc[
                    result["candidate_record_id"] == candidate_id, "ld_component_id"
                ].iloc[0]
            )
            for candidate_id in selected
        }
        # A diversity anchor may replace the default highest-priority member
        # inside a component. This preserves one selected candidate per
        # component while allowing a unique experimental evidence stratum to
        # be represented explicitly.
        for row in eligible_candidates.to_dict(orient="records"):
            candidate_id = str(row["candidate_record_id"])
            component_id = clean(row.get("ld_component_id"))
            if (
                candidate_id not in selected
                and component_id not in represented_components
                and criterion(row)
            ):
                selected.append(candidate_id)
                anchor[candidate_id] = criterion_name
                break

    represented_components = {
        clean(
            result.loc[
                result["candidate_record_id"] == candidate_id, "ld_component_id"
            ].iloc[0]
        )
        for candidate_id in selected
    }
    for candidate_id in representatives["candidate_record_id"].astype(str):
        if len(selected) >= SHORTLIST_SIZE:
            break
        component_id = clean(
            result.loc[
                result["candidate_record_id"] == candidate_id, "ld_component_id"
            ].iloc[0]
        )
        if candidate_id not in selected and component_id not in represented_components:
            selected.append(candidate_id)
            anchor[candidate_id] = "priority_fill_after_diversity_anchors"
            represented_components.add(component_id)

    if len(selected) != SHORTLIST_SIZE:
        raise ValueError(
            f"Selection produced {len(selected)} candidates; expected {SHORTLIST_SIZE}. "
            "Check eligibility and component definitions."
        )

    selected_order = sorted(selected, key=lambda item: (int(rank_map[item]), item))
    shortlist_rank = {candidate_id: index for index, candidate_id in enumerate(selected_order, start=1)}
    result["selected_for_experimental_shortlist"] = result["candidate_record_id"].isin(selected)
    result["experimental_shortlist_rank"] = result["candidate_record_id"].map(shortlist_rank).astype("Int64")
    result["selection_anchor_criterion"] = result["candidate_record_id"].map(anchor).fillna("")
    component_selected = {
        clean(
            result.loc[
                result["candidate_record_id"] == candidate_id, "ld_component_id"
            ].iloc[0]
        ): candidate_id
        for candidate_id in selected
    }
    result["ld_component_selected_candidate"] = result["ld_component_id"].map(
        component_selected
    ).fillna("")

    coverage_values: List[str] = []
    outcomes: List[str] = []
    for row in result.to_dict(orient="records"):
        candidate_id = str(row["candidate_record_id"])
        coverage = [name for name, criterion in criteria if criterion(row)]
        coverage_values.append(";".join(coverage))
        if not as_bool(row.get("experimental_eligible")):
            outcomes.append(f"ineligible:{clean(row.get('experimental_eligibility_reason'))}")
        elif candidate_id in selected:
            outcomes.append(f"selected:{anchor[candidate_id]}")
        elif clean(row.get("ld_component_selected_candidate")):
            outcomes.append(
                "LD_redundant_with_selected:"
                + clean(row.get("ld_component_selected_candidate"))
            )
        elif not as_bool(row.get("is_ld_component_representative")):
            outcomes.append(
                "LD_redundant_with_default_representative:"
                + clean(row.get("ld_component_representative"))
            )
        else:
            outcomes.append("eligible_nonredundant_not_selected_after_12_slots")
    result["diversity_coverage_criteria"] = coverage_values
    result["selection_outcome"] = outcomes
    result["selection_policy"] = SELECTION_POLICY
    return result


def shortlist_table(annotated: pd.DataFrame) -> pd.DataFrame:
    selected = annotated[annotated["selected_for_experimental_shortlist"]].copy()
    selected = selected.sort_values("experimental_shortlist_rank")
    selected["experimental_shortlist_rank"] = selected["experimental_shortlist_rank"].astype(int)
    selected["evidence_interpretation"] = (
        "ranked computational candidate selected for proposed validation; not causal proof"
    )
    fields = [
        "experimental_shortlist_rank",
        "predicted_causal_priority_rank",
        "canonical_section4_list",
        "canonical_rank_preserved",
        "rank_basis",
        "interpretation_boundary",
        "selection_anchor_criterion",
        "diversity_coverage_criteria",
        "ld_component_id",
        "ld_component_size",
        "ld_component_members",
        "ld_component_representative",
        "ld_component_selected_candidate",
        "candidate_record_id",
        "panel_variant_id",
        "selected_ensembl_variation_id",
        "ensembl_variation_ids",
        "candidate_rsids",
        "chromosome_grch38",
        "position_grch38",
        "ref",
        "alt",
        "variant_class",
        "model_context",
        "candidate_origin",
        "source_focal_tags",
        "linked_gwas_tag_ids",
        "ld_panels",
        "max_r2_across_links",
        "linked_gwas_genes",
        "linked_study_accessions",
        "assigned_selected_loci",
        "predicted_causal_enhancer",
        "enhancer_ref_score",
        "enhancer_alt_score",
        "enhancer_delta_alt_minus_ref",
        "enhancer_abs_delta_percentile_within_class",
        "enhancer_allelic_direction",
        "enhancer_threshold_ratio",
        "predicted_causal_silencer",
        "silencer_ref_score",
        "silencer_alt_score",
        "silencer_delta_alt_minus_ref",
        "silencer_abs_delta_percentile_within_class",
        "silencer_allelic_direction",
        "silencer_threshold_ratio",
        "coding",
        "coding_CDS_gene_names",
        "donor_refined_enhancer",
        "donor_refined_silencer",
        "donor_refined_any",
        "known_regulatory_any",
        "exclusive_class_comprehensive",
        "repeatmasker",
        "frequency_match_status",
        "global_alt_af",
        "AFR_alt_af",
        "AMR_alt_af",
        "EAS_alt_af",
        "EUR_alt_af",
        "SAS_alt_af",
        "frequency_class",
        "strict_population_specific",
        "ancestral_allele",
        "ancestral_call_status",
        "derived_candidate_allele",
        "global_derived_allele_frequency",
        "nearest_gene",
        "nearest_gene_tss_distance_bp",
        "nearest_protein_coding_gene",
        "nearest_protein_coding_tss_distance_bp",
        "genes_with_tss_within_100kb",
        "R006_provisional_target_count",
        "R006_provisional_target_genes",
        "R006_target_mapping_methods",
        "R006_target_method_gene_pairs",
        "target_mapping_limit",
        "gtex_lung_exact_significant_eqtl",
        "gtex_lung_significant_pair_count",
        "gtex_lung_significant_gene_count",
        "gtex_lung_significant_genes",
        "gtex_lung_minimum_nominal_p",
        "section6_exact_GTEx_R006_target_count",
        "section6_exact_GTEx_R006_target_genes",
        "mprabase_any_element_evidence",
        "mprabase_coordinate_element_evidence",
        "mprabase_exact_rsid_metadata_evidence",
        "mprabase_matched_element_count",
        "mprabase_matched_datasets",
        "open_targets_COPD_associated_gene_count",
        "open_targets_COPD_associated_genes",
        "integrated_validation_status",
        "enhancer_tf_any_motif_compatible_site",
        "enhancer_tf_disrupts_any_motif_compatible_site",
        "enhancer_tf_creates_any_motif_compatible_site",
        "enhancer_tf_largest_effect_motif_id",
        "enhancer_tf_largest_effect_tf_name",
        "enhancer_tf_largest_signed_delta_relative_score",
        "enhancer_tf_largest_absolute_delta_relative_score",
        "silencer_tf_any_motif_compatible_site",
        "silencer_tf_disrupts_any_motif_compatible_site",
        "silencer_tf_creates_any_motif_compatible_site",
        "silencer_tf_largest_effect_motif_id",
        "silencer_tf_largest_effect_tf_name",
        "silencer_tf_largest_signed_delta_relative_score",
        "silencer_tf_largest_absolute_delta_relative_score",
        "target_data_status",
        "evidence_interpretation",
    ]
    for field in fields:
        if field not in selected.columns:
            selected[field] = ""
    return selected[fields]


def selection_audit_table(annotated: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "predicted_causal_priority_rank",
        "canonical_section4_list",
        "canonical_rank_preserved",
        "rank_basis",
        "candidate_record_id",
        "chromosome_grch38",
        "position_grch38",
        "ref",
        "alt",
        "variant_class",
        "model_context",
        "source_focal_tags",
        "linked_gwas_tag_ids",
        "linked_gwas_genes",
        "assigned_selected_loci",
        "sequence_status",
        "blacklisted",
        "experimental_eligible",
        "experimental_eligibility_reason",
        "ld_component_id",
        "ld_component_size",
        "ld_component_members",
        "within_ld_component_priority_rank",
        "ld_component_representative",
        "ld_component_selected_candidate",
        "is_ld_component_representative",
        "ld_component_definition",
        "diversity_coverage_criteria",
        "selected_for_experimental_shortlist",
        "experimental_shortlist_rank",
        "selection_anchor_criterion",
        "selection_outcome",
        "selection_policy",
        "gtex_lung_exact_significant_eqtl",
        "mprabase_any_element_evidence",
        "integrated_validation_status",
    ]
    for field in fields:
        if field not in annotated.columns:
            annotated[field] = ""
    return annotated.sort_values("predicted_causal_priority_rank")[fields]


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def longest_homopolymer(sequence: str) -> int:
    if not sequence:
        return 0
    longest = current = 1
    for left, right in zip(sequence, sequence[1:]):
        current = current + 1 if left == right else 1
        longest = max(longest, current)
    return longest


def construct_table(shortlist: pd.DataFrame) -> pd.DataFrame:
    refs = read_fasta(REF_FASTA)
    alts = read_fasta(ALT_FASTA)
    rows: List[Dict[str, object]] = []
    anchor = 1000
    for record in shortlist.to_dict(orient="records"):
        candidate = clean(record["candidate_record_id"])
        ref = clean(record["ref"]).upper()
        alt = clean(record["alt"]).upper()
        if candidate not in refs or candidate not in alts:
            raise ValueError(f"Shortlisted candidate missing from paired FASTAs: {candidate}")
        source_sequences = {"REF": refs[candidate], "ALT": alts[candidate]}
        alleles = {"REF": ref, "ALT": alt}
        flanks: Dict[str, Tuple[str, str]] = {}
        for allele_label, allele in alleles.items():
            sequence = source_sequences[allele_label]
            if len(sequence) != 2001:
                raise ValueError(f"Unexpected source sequence length for {candidate} {allele_label}")
            observed = sequence[anchor : anchor + len(allele)]
            if observed != allele:
                raise ValueError(
                    f"FASTA allele mismatch for {candidate} {allele_label}: expected {allele}, got {observed}"
                )
            left = sequence[anchor - MPRA_FLANK_BP : anchor]
            right = sequence[
                anchor + len(allele) : anchor + len(allele) + MPRA_FLANK_BP
            ]
            if len(left) != MPRA_FLANK_BP or len(right) != MPRA_FLANK_BP:
                raise ValueError(f"Insufficient MPRA flank for {candidate} {allele_label}")
            flanks[allele_label] = (left, right)
        if flanks["REF"] != flanks["ALT"]:
            raise ValueError(f"REF and ALT genomic MPRA flanks differ for {candidate}")

        for allele_label in ("REF", "ALT"):
            allele = alleles[allele_label]
            left, right = flanks[allele_label]
            forward = left + allele + right
            for orientation, insert in (
                ("forward", forward),
                ("reverse_complement", reverse_complement(forward)),
            ):
                oriented_allele = allele if orientation == "forward" else reverse_complement(allele)
                expected_center = insert[MPRA_FLANK_BP : MPRA_FLANK_BP + len(oriented_allele)]
                if expected_center != oriented_allele:
                    raise ValueError(f"Oriented allele is not centered correctly for {candidate}")
                gc = 100.0 * (insert.count("G") + insert.count("C")) / len(insert)
                warnings: List[str] = []
                if gc < 25 or gc > 75:
                    warnings.append("extreme_GC")
                if longest_homopolymer(insert) >= 9:
                    warnings.append("homopolymer_ge_9")
                if any(site in insert for site in ("GGTCTC", "GAGACC")):
                    warnings.append("contains_BsaI_site")
                if any(site in insert for site in ("CGTCTC", "GAGACG")):
                    warnings.append("contains_BsmBI_site")
                rank = int(record["experimental_shortlist_rank"])
                label = "FWD" if orientation == "forward" else "RC"
                rows.append(
                    {
                        "construct_id": f"COPD_S6_R003_C{rank:02d}_{allele_label}_{label}",
                        "experimental_shortlist_rank": rank,
                        "candidate_record_id": candidate,
                        "chromosome_grch38": record["chromosome_grch38"],
                        "position_grch38": record["position_grch38"],
                        "genomic_ref": ref,
                        "genomic_alt": alt,
                        "variant_class": record["variant_class"],
                        "allele_label": allele_label,
                        "genomic_forward_allele": allele,
                        "orientation": orientation,
                        "oriented_allele_sequence": oriented_allele,
                        "left_flank_bp": MPRA_FLANK_BP,
                        "right_flank_bp": MPRA_FLANK_BP,
                        "allele_start_index_zero_based": MPRA_FLANK_BP,
                        "insert_length_bp": len(insert),
                        "insert_sequence": insert,
                        "GC_percent": round(gc, 6),
                        "longest_homopolymer_bp": longest_homopolymer(insert),
                        "sequence_review_flags": ";".join(warnings) if warnings else "none",
                        "genome_build": "GRCh38",
                        "source_fasta": relative(REF_FASTA if allele_label == "REF" else ALT_FASTA),
                        "source_variant_first_base_index_zero_based": anchor,
                        "sequence_background": (
                            "GRCh38 reference-sequence flanks with only the nominated allele substituted; "
                            "local phased haplotypes and nearby linked alleles are not represented"
                        ),
                        "indel_handling": (
                            "normalized VCF allele substituted at source index 1000; exactly 100 genomic "
                            "bases retained on each side; REF and ALT inserts may differ in length"
                        ),
                        "design_scope": (
                            "insert sequence only; adapters, minimal promoter, barcode architecture, "
                            "restriction-site remediation, and synthesis acceptance require platform review"
                        ),
                        "experimental_status": "proposed_not_synthesized_or_tested",
                    }
                )
    return pd.DataFrame(rows)


def editing_strategy(ref: str, alt: str) -> str:
    ref = clean(ref).upper()
    alt = clean(alt).upper()
    if len(ref) != 1 or len(alt) != 1:
        return (
            "prime-editing feasibility screen for normalized indel/complex allele; consider HDR as a "
            "secondary route; no guide or pegRNA is specified here"
        )
    if (ref, alt) in {("A", "G"), ("T", "C")}:
        return (
            "ABE feasibility screen on an appropriate strand and activity window; use prime editing if "
            "PAM, bystander, or editing-window constraints fail"
        )
    if (ref, alt) in {("C", "T"), ("G", "A")}:
        return (
            "CBE feasibility screen on an appropriate strand and activity window; use prime editing if "
            "PAM, bystander, or editing-window constraints fail"
        )
    return "prime-editing feasibility screen for transversion; no guide or pegRNA is specified here"


def primary_targets(row: Mapping[str, object]) -> str:
    exact = tokens(row.get("section6_exact_GTEx_R006_target_genes"))
    if exact:
        return ";".join(exact)
    gtex = tokens(row.get("gtex_lung_significant_genes"))
    if gtex:
        return ";".join(gtex)
    provisional = []
    provisional.extend(tokens(row.get("nearest_protein_coding_gene")))
    provisional.extend(tokens(row.get("linked_gwas_genes")))
    provisional.extend(tokens(row.get("assigned_selected_loci")))
    return ";".join(sorted(set(provisional)))


def candidate_plan_table(shortlist: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for row in shortlist.to_dict(orient="records"):
        candidate = clean(row["candidate_record_id"])
        targets = primary_targets(row)
        rows.append(
            {
                "experimental_shortlist_rank": row["experimental_shortlist_rank"],
                "candidate_record_id": candidate,
                "model_context": row["model_context"],
                "variant_class": row["variant_class"],
                "primary_target_hypotheses": targets,
                "exact_GTEx_R006_concordant_targets": row.get(
                    "section6_exact_GTEx_R006_target_genes", ""
                ),
                "all_exact_GTEx_Lung_eGenes": row.get("gtex_lung_significant_genes", ""),
                "R006_provisional_targets": row.get("R006_provisional_target_genes", ""),
                "target_interpretation": (
                    "targets are hypotheses; an exact lung eQTL is molecular association and R006 mappings "
                    "are proximity/GWAS/coding hypotheses, not enhancer-to-gene proof"
                ),
                "stage_1_MPRA": (
                    "test REF and ALT in both orientations using the R003 inserts; use at least 10 independent "
                    "barcodes per construct and at least 3 independent biological replicates per selected context"
                ),
                "stage_1_contexts": "airway_epithelium;alveolar_AT2;parenchymal_lung_fibroblast",
                "stage_1_analysis": (
                    "pre-register barcode-level quality filters and a replicate-aware allele-by-orientation-by-"
                    "context model; report effect sizes and confidence intervals with multiplicity control"
                ),
                "stage_2_CRISPRi": (
                    "after context-specific accessibility/activity screening, repress the endogenous element "
                    "with at least 3 non-overlapping sgRNAs and quantify targeted genes plus transcriptome-wide effects"
                ),
                "stage_2_controls": (
                    "non-targeting sgRNAs; safe-harbor sgRNAs; promoter-targeting positive control; multiple "
                    "independent element sgRNAs; dCas9-KRAB-only or mock control; viability and interferon-response checks"
                ),
                "stage_3_allele_editing": editing_strategy(row["ref"], row["alt"]),
                "stage_3_editing_controls": (
                    "unedited parental and mock-edited pools; independently edited pools; at least 2 independently "
                    "derived clones per genotype when cloning is used; amplicon genotype and top predicted off-target "
                    "checks; monitor bystanders, large deletions, copy number, and clone-specific expression"
                ),
                "molecular_readouts": (
                    "allele-specific reporter output; target-gene RT-qPCR and RNA-seq; local ATAC-seq; H3K27ac "
                    "CUT&Tag where informative; genotype-dose response; cell-state and viability markers"
                ),
                "contextual_COPD_readouts": (
                    "airway barrier/differentiation markers; AT2 identity and injury-repair markers; fibroblast "
                    "proliferation, alpha-SMA, collagen, and extracellular-matrix outputs, selected before unblinding"
                ),
                "progression_rule": (
                    "advance only after a reproducible allele effect or endogenous element effect in a context with "
                    "measurable element activity and target expression; require concordant independent guides or edits"
                ),
                "causal_claim_boundary": (
                    "MPRA is episomal; CRISPRi is region-level; only a genotype-verified endogenous allele perturbation "
                    "with replicated molecular consequences can support allele causality, and no such experiment is complete"
                ),
                "experimental_status": "proposed_not_performed",
            }
        )
    return pd.DataFrame(rows).sort_values("experimental_shortlist_rank")


def motif_effect_class(row: Mapping[str, object]) -> str:
    disrupts = as_bool(row.get("disrupts_motif_compatible_site"))
    creates = as_bool(row.get("creates_motif_compatible_site"))
    retains = as_bool(row.get("retains_motif_compatible_site"))
    if disrupts and creates:
        return "disrupts_and_creates_across_best_allelic_windows"
    if disrupts:
        return "ALT_disrupts_REF_compatible_site"
    if creates:
        return "ALT_creates_compatible_site"
    if retains:
        return "both_alleles_retain_compatible_site_with_score_shift"
    return "score_shift_without_0.8_compatible_site"


def tf_plan_table(shortlist: pd.DataFrame) -> pd.DataFrame:
    motifs = pd.read_csv(R008_MOTIFS, sep="\t", low_memory=False)
    motifs["candidate_record_id"] = motifs["candidate_record_id"].astype(str)
    ids = set(shortlist["candidate_record_id"].astype(str))
    motifs = motifs[motifs["candidate_record_id"].isin(ids)].copy()
    shortlist_by_id = shortlist.set_index("candidate_record_id").to_dict(orient="index")
    rows: List[Dict[str, object]] = []
    for candidate in sorted(ids, key=lambda item: int(shortlist_by_id[item]["experimental_shortlist_rank"])):
        metadata = shortlist_by_id[candidate]
        models = []
        if as_bool(metadata.get("predicted_causal_enhancer")):
            models.append("enhancer")
        if as_bool(metadata.get("predicted_causal_silencer")):
            models.append("silencer")
        for model in models:
            candidate_rows = motifs[
                (motifs["candidate_record_id"] == candidate) & (motifs["model_type"] == model)
            ].copy()
            if candidate_rows.empty:
                raise ValueError(f"No R008 motif rows for shortlisted {candidate} {model}")
            candidate_rows["effect_priority"] = candidate_rows.apply(
                lambda row: 3
                if as_bool(row.get("disrupts_motif_compatible_site"))
                or as_bool(row.get("creates_motif_compatible_site"))
                else 2
                if as_bool(row.get("retains_motif_compatible_site"))
                else 1,
                axis=1,
            )
            candidate_rows = candidate_rows.sort_values(
                ["effect_priority", "absolute_delta_relative_score", "minimum_match_q_value", "motif_id"],
                ascending=[False, False, True, True],
            )
            best = candidate_rows.iloc[0].to_dict()
            rows.append(
                {
                    "experimental_shortlist_rank": metadata["experimental_shortlist_rank"],
                    "candidate_record_id": candidate,
                    "model_type": model,
                    "motif_id": best["motif_id"],
                    "TF_name": best["tf_name"],
                    "minimum_pattern_match_q_value": best["minimum_match_q_value"],
                    "matched_TF_MoDISco_patterns": best["matched_patterns"],
                    "motif_relative_score_threshold": best["relative_score_threshold"],
                    "REF_relative_score": best["ref_relative_score"],
                    "ALT_relative_score": best["alt_relative_score"],
                    "delta_ALT_minus_REF_relative_score": best[
                        "delta_alt_minus_ref_relative_score"
                    ],
                    "absolute_delta_relative_score": best["absolute_delta_relative_score"],
                    "REF_best_window_start_zero_based": best["ref_best_window_start_zero_based"],
                    "REF_best_window_strand": best["ref_best_window_strand"],
                    "ALT_best_window_start_zero_based": best["alt_best_window_start_zero_based"],
                    "ALT_best_window_strand": best["alt_best_window_strand"],
                    "predicted_motif_effect": motif_effect_class(best),
                    "TF_gate": (
                        "confirm TF RNA/protein expression and nuclear localization in the tested cell context before perturbation"
                    ),
                    "binding_test": (
                        "allele-paired EMSA or equivalent biochemical assay, followed by CUT&RUN/CUT&Tag or ChIP in "
                        "heterozygous or isogenic cells when antibody quality and chromatin accessibility permit"
                    ),
                    "TF_perturbation": (
                        "use at least two independent TF perturbation reagents (CRISPRi, degron, or RNAi as appropriate), "
                        "measure MPRA and endogenous target response, and rescue with perturbation-resistant TF expression"
                    ),
                    "interaction_test": (
                        "test allele-by-TF-perturbation interaction rather than interpreting parallel main effects as mechanism"
                    ),
                    "controls": (
                        "IgG/no-antibody binding control; non-targeting perturbation; expression-matched rescue; motif-scrambled "
                        "reporter; positive TF-binding control; cell-viability and global stress-response checks"
                    ),
                    "interpretation_boundary": (
                        "TF-MoDISco/JASPAR similarity and PWM threshold crossing are sequence hypotheses, not TF binding evidence"
                    ),
                    "experimental_status": "proposed_not_performed",
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["experimental_shortlist_rank", "model_type"]
    )


def context_table() -> pd.DataFrame:
    rows = [
        {
            "context_id": "airway_epithelium",
            "biological_model": (
                "primary small-airway or bronchial epithelial cells from COPD donors and smoking-history-matched controls, "
                "differentiated at air-liquid interface where feasible"
            ),
            "role": "primary airway validation context",
            "pre_assay_gate": "target expression, element accessibility, epithelial identity, mycoplasma-free status, viability",
            "perturbation_conditions": (
                "baseline plus vehicle-matched, dose-calibrated cigarette-smoke or oxidative challenge only after pilot toxicity testing"
            ),
            "biological_replication": "at least 3 independent donor lines per comparison group; balance sex and smoking history where feasible",
            "core_controls": (
                "non-COPD smoking-matched donors; vehicle; non-targeting and safe-harbor guides; promoter positive control; "
                "multiple element guides; blinded sample labels during primary quantification"
            ),
            "readouts": "barrier and differentiation markers, target RNA/protein, RNA-seq, ATAC-seq, viability",
            "limit": "donor heterogeneity and epithelial-state composition require replicate-aware analysis",
            "experimental_status": "proposed_not_performed",
        },
        {
            "context_id": "alveolar_AT2",
            "biological_model": (
                "primary human AT2 cells when feasible or quality-controlled iPSC-derived AT2 cells with retained AT2 identity"
            ),
            "role": "alveolar injury-repair and emphysema-relevant validation context",
            "pre_assay_gate": "SFTPC/ABCA3 identity, target expression, element accessibility, genomic stability, viability",
            "perturbation_conditions": (
                "baseline plus a pre-specified, viability-calibrated smoke/oxidative or TGF-beta response condition when biologically justified"
            ),
            "biological_replication": "at least 3 independent donors or independently differentiated iPSC backgrounds",
            "core_controls": (
                "isogenic unedited cells; non-targeting/safe-harbor guides; differentiation-batch controls; identity-marker and viability controls"
            ),
            "readouts": "AT2 identity, injury-repair response, target RNA/protein, RNA-seq, ATAC-seq, organoid output where established",
            "limit": "iPSC-derived AT2 maturation and primary-cell expansion can alter regulatory state",
            "experimental_status": "proposed_not_performed",
        },
        {
            "context_id": "parenchymal_lung_fibroblast",
            "biological_model": (
                "low-passage primary distal/parenchymal lung fibroblasts from COPD donors and smoking-history-matched controls"
            ),
            "role": "mesenchymal and extracellular-matrix validation context",
            "pre_assay_gate": "target expression, element accessibility, passage, myofibroblast baseline, viability",
            "perturbation_conditions": "baseline plus vehicle-matched TGF-beta response when relevant to the target hypothesis",
            "biological_replication": "at least 3 independent donor lines per comparison group; passage matched",
            "core_controls": (
                "non-COPD smoking-matched donors; non-targeting/safe-harbor guides; promoter positive control; passage and density controls"
            ),
            "readouts": "target RNA/protein, proliferation, alpha-SMA, collagen/ECM output, RNA-seq, ATAC-seq",
            "limit": "culture passage, stiffness, and activation state can dominate fibroblast regulatory programs",
            "experimental_status": "proposed_not_performed",
        },
        {
            "context_id": "assay_development_cell_lines",
            "biological_model": "16HBE or BEAS-2B for optimization only; use MRC5 only for fibroblast-line pilot work",
            "role": "technical optimization and comparison with published COPD regulatory studies",
            "pre_assay_gate": "identity/authentication, mycoplasma testing, target expression, element accessibility",
            "perturbation_conditions": "baseline technical optimization",
            "biological_replication": "at least 3 independent experiments; does not replace donor replication",
            "core_controls": "the same vector, guide, viability, and positive/negative assay controls as primary models",
            "readouts": "assay dynamic range, transduction/editing efficiency, target RNA, viability",
            "limit": "cell-line findings cannot establish COPD-donor or primary-cell relevance",
            "experimental_status": "proposed_not_performed",
        },
    ]
    return pd.DataFrame(rows)


def upstream_manifest() -> pd.DataFrame:
    paths = [
        ("canonical_Section4_ranked_list", R010),
        ("ranked_predicted_causal_candidates", R004),
        ("population_genetics", R005),
        ("provisional_target_summary", R006_SUMMARY),
        ("provisional_target_evidence", R006_EVIDENCE),
        ("candidate_TFBS_summary", R008_SUMMARY),
        ("allele_specific_motif_scores", R008_MOTIFS),
        ("integrated_computational_validation", S5_INTEGRATED),
        ("paired_REF_sequences", REF_FASTA),
        ("paired_ALT_sequences", ALT_FASTA),
    ]
    rows = []
    for role, path in paths:
        if not path.exists():
            continue
        rows.append(
            {
                "input_role": role,
                "path": relative(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "consumption_scope": (
                    "full table joined by candidate_record_id"
                    if path.suffix in {".tsv", ".gz"} and "FASTA" not in role
                    else "shortlisted candidate records extracted and sequence-validated"
                ),
            }
        )
    return pd.DataFrame(rows)


def write_collaborator_readme(shortlist: pd.DataFrame, constructs: pd.DataFrame) -> None:
    content = f"""# COPD Section 6 collaborator handoff

Package status: computational design complete; all experiments are proposed and none has been performed.

## What is included

- `COPD-S6-R001_candidate_shortlist.tsv`: {len(shortlist)} nonredundant candidates with exact upstream evidence fields and explicit selection reasons.
- `COPD-S6-R002_selection_audit.tsv.gz`: every Section 4 predicted-causal record, LD-component assignment, eligibility, and selection outcome.
- `COPD-S6-R003_MPRA_constructs.tsv`: {len(constructs)} validated insert designs, representing each shortlisted variant's REF and ALT alleles in both orientations.
- `COPD-S6-R004_candidate_validation_plan.tsv`: staged MPRA, endogenous CRISPRi, and allele-editing proposals.
- `COPD-S6-R004_TF_first_plan.tsv`: one leading sequence-motif hypothesis for each candidate/model combination and a binding/perturbation test plan.
- `COPD-S6-R004_cell_context_controls.tsv`: airway, alveolar, fibroblast, and technical-optimization contexts with controls and limitations.
- `COPD-S6-R005_upstream_inputs.tsv`: immutable paths, sizes, and SHA-256 checksums for every consumed upstream artifact.
- `section_6_experimental_validation.md`: integrated interpretation and design summary.
- `COPD-S6_validation_checks.tsv`: machine-readable reconciliation checks.

## Sequence-design boundary

The R003 sequences are genomic inserts, not complete synthesis oligonucleotides. Each forward insert contains 100 genomic bases, the normalized VCF allele, and 100 genomic bases. Reverse-complement designs are exact reverse complements. REF and ALT indel constructs therefore can differ in length. Each construct isolates one nominated allele in GRCh38 reference-sequence flanks; local phased donor haplotypes and nearby linked alleles are not represented. Effects can differ on native haplotypes, so relevant phased haplotypes should be tested when biological or published evidence warrants it. The receiving laboratory must choose the MPRA backbone, promoter, adapters, cloning sites, barcode architecture, minimum barcode count, and any restriction-site remediation before ordering.

## Evidence boundary

The shortlist preserves the canonical Section 4 R010 computational priority and adds auditable diversity constraints. A TREDNet call is a model prediction. TF-MoDISco/JASPAR and PWM results are sequence hypotheses, not TF binding evidence. GTEx is an exact variant-gene association in bulk lung. MPRAbase coordinate containment is regional assay coverage unless allelic testing is independently established. R006 target mappings are provisional. Neither selection nor database overlap constitutes causal proof.

## Recommended staged handoff

1. Platform review of insert architecture and sequence flags.
2. Allele-specific MPRA in pre-screened airway, alveolar, and fibroblast contexts.
3. Endogenous region-level CRISPRi only in contexts with accessible chromatin, target expression, and interpretable reporter or model evidence.
4. Genotype-verified allele editing for candidates with reproducible reporter and region-level effects.
5. TF binding and perturbation tests using the R004 TF-first hypotheses, including an allele-by-TF interaction test.

All statistical models, exclusion rules, biological replicates, and progression criteria should be preregistered with the receiving laboratory before unblinding primary readouts.
"""
    OUT_COLLAB_README.write_text(content)


def report_text(
    annotated: pd.DataFrame,
    shortlist: pd.DataFrame,
    constructs: pd.DataFrame,
    candidate_plans: pd.DataFrame,
    tf_plans: pd.DataFrame,
) -> str:
    n_candidates = len(annotated)
    n_components = annotated["ld_component_id"].nunique()
    n_redundant = int((~annotated["is_ld_component_representative"]).sum())
    n_gtex = int(shortlist["gtex_lung_exact_significant_eqtl"].map(as_bool).sum())
    n_mpra = int(shortlist["mprabase_any_element_evidence"].map(as_bool).sum())
    n_donor = int(shortlist["donor_refined_any"].map(as_bool).sum())
    n_indel = int((shortlist["variant_class"] != "SNV").sum())
    n_target_concordant = int(
        (pd.to_numeric(shortlist["section6_exact_GTEx_R006_target_count"], errors="coerce").fillna(0) > 0).sum()
    )
    context_counts = Counter(shortlist["model_context"])
    sequence_lengths = Counter(int(value) for value in constructs["insert_length_bp"])
    flagged = int((constructs["sequence_review_flags"] != "none").sum())

    table_rows = []
    for row in shortlist.to_dict(orient="records"):
        target = clean(row.get("section6_exact_GTEx_R006_target_genes"))
        if not target:
            target = clean(row.get("nearest_protein_coding_gene"))
        tf_names = joined(
            [row.get("enhancer_tf_largest_effect_tf_name"), row.get("silencer_tf_largest_effect_tf_name")]
        )
        table_rows.append(
            f"| {int(row['experimental_shortlist_rank'])} | `{row['candidate_record_id']}` | "
            f"{row['model_context']} | {clean(row['variant_class'])} | {target or 'none'} | "
            f"{'yes' if as_bool(row.get('gtex_lung_exact_significant_eqtl')) else 'no'} | "
            f"{tf_names or 'no named leading match'} | {clean(row['selection_anchor_criterion'])} |"
        )
    table = "\n".join(table_rows)

    return f"""# Section 6: Experimental and Biological Validation

Result: COPD-S6-R006. Status: design and collaborator handoff complete; wet-lab execution is pending. No experimental result or causal proof is claimed. The computational shortlist is not causal proof.

## Scope

Section 6 consumed the canonical 337-row Section 4 R010 list and preserved its `predicted_causal_priority_rank`. R004-R008 supplied sequence eligibility, interval, population, target, and allele-specific motif detail; the analysis also consumed all 4,779 R006 target-evidence rows and the final, validated Section 5 integrated matrix. All Section 6 tables distinguish observed or computed upstream evidence from proposed future experiments.

## Transparent candidate selection

The starting set contained {n_candidates:,} predicted-causal records. Exact shared linked GWAS tags or source focal tags were joined transitively into {n_components:,} LD-redundancy components; {n_redundant:,} records were lower-priority members of a component. At most one candidate was selected per component. The policy first ensured coverage of model context, variant class, severe-emphysema-donor overlap, exact GTEx Lung eQTL evidence, public MPRAbase element evidence, and predicted motif creation/disruption when available. A unique diversity anchor could replace the default highest-priority member of its component; remaining slots were filled with default component representatives by the unchanged Section 4 rank. It did not create or optimize a new composite score.

The final panel has {context_counts.get('both', 0)} both-model, {context_counts.get('enhancer_only', 0)} enhancer-only, and {context_counts.get('silencer_only', 0)} silencer-only candidates. It includes {n_indel} indel/complex candidate(s), {n_donor} candidate(s) overlapping a refined regulatory interval from the severe-emphysema donor, {n_gtex} exact GTEx v10 Lung significant eQTL candidate(s), {n_mpra} candidate(s) with MPRAbase element evidence, and {n_target_concordant} candidate(s) for which an exact GTEx Lung eGene also occurs in the R006 provisional target set.

| S6 rank | Candidate | Model context | Class | Exact-concordant or nearest protein-coding target | Exact lung eQTL | Leading motif TF names | Selection rule |
|---:|---|---|---|---|---|---|---|
{table}

`COPD-S6-R002_selection_audit.tsv.gz` retains every candidate, component membership, eligibility, coverage strata, and exclusion or inclusion reason. Shared-tag components are an LD-redundancy heuristic tied to the available workflow links; they are not claims that every member is mutually correlated in every ancestry.

## Allele-specific MPRA design

R003 contains {len(constructs)} insert designs: REF and ALT alleles in forward and reverse-complement orientations for each of {len(shortlist)} variants. Each forward insert has exactly {MPRA_FLANK_BP} genomic bases on each side of the normalized allele. Insert lengths are {', '.join(f'{length} bp (n={count})' for length, count in sorted(sequence_lengths.items()))}. Paired FASTA allele identity, shared genomic flanks, insert length, centered allele, and reverse-complement identity are validated computationally. {flagged} constructs carry a sequence-review flag for GC, homopolymers, or common Type IIS recognition sites. Each design isolates the nominated allele in GRCh38 reference-sequence flanks and does not encode local phased donor haplotypes or nearby linked alleles. Because reporter effects can depend on haplotype background, relevant phased haplotypes should be tested where biological or prior functional evidence warrants it.

These are inserts, not ordering-ready oligonucleotides. Vector, promoter, adapters, barcodes, restriction-site remediation, randomization, and assay-specific positive/negative controls must be chosen with the receiving laboratory. The proposed minimum is 10 independent barcodes per construct and 3 independent biological replicates per selected context, subject to prospective power analysis.

## Endogenous validation plan

The R004 candidate plan is staged so that evidence classes remain separate.

1. Allelic MPRA tests sequence-dependent reporter activity in pre-screened airway epithelium, alveolar AT2, and parenchymal lung fibroblast contexts.
2. CRISPRi tests the endogenous region with at least three non-overlapping guides. A region effect does not identify the causal nucleotide.
3. Base editing is considered only for compatible transitions after PAM, activity-window, and bystander review. Prime editing is the primary proposal for transversions and indels. Genotype-verified edited pools and independent clones, off-target review, and clone-aware controls are required.
4. Target-gene readouts combine exact GTEx Lung eGenes and R006 hypotheses with targeted and transcriptome-wide measurements. Target assignments remain provisional until endogenous perturbation supports them.

Primary COPD-relevant models are donor-derived airway epithelial cultures, primary or iPSC-derived AT2 cells, and low-passage parenchymal lung fibroblasts. 16HBE, BEAS-2B, or MRC5 may support technical optimization but cannot substitute for donor replication. COPD and smoking-history-matched control donors, cell identity, viability, passage or differentiation state, and exposure conditions are explicit covariates and controls.

## TF-first perturbation plan

R004 provides {len(tf_plans)} candidate-by-model TF hypotheses. For each model-predicted candidate, it selects the strongest compatible-site creation/disruption first, then a retained site or score-only shift. The proposed sequence is expression gating, allele-paired biochemical binding, chromatin binding in a relevant heterozygous or edited context, TF perturbation with independent reagents and rescue, and a formal allele-by-TF interaction test. Motif similarity and PWM changes are not TF binding evidence or TF occupancy evidence.

## Interpretation and status

All {len(candidate_plans)} candidate plans and all {len(tf_plans)} TF plans are labeled `proposed_not_performed`. MPRA can support allele-sensitive reporter activity; CRISPRi can support endogenous region function; target expression after region repression can support a region-to-gene link. Strong allele-level causal support requires a genotype-verified endogenous allele change with replicated molecular consequences and appropriate controls. No such Section 6 experiment has been conducted in this workflow.

## Reproducible outputs

- `COPD-S6-R001_candidate_shortlist.tsv`: exact candidate evidence and selection fields.
- `COPD-S6-R002_selection_audit.tsv.gz`: all-candidate LD redundancy and selection audit.
- `COPD-S6-R003_MPRA_constructs.tsv`: paired-allele, paired-orientation inserts.
- `COPD-S6-R004_candidate_validation_plan.tsv`: staged endogenous-validation plan.
- `COPD-S6-R004_TF_first_plan.tsv`: candidate-specific TF hypotheses and tests.
- `COPD-S6-R004_cell_context_controls.tsv`: cell models, controls, and limitations.
- `COPD-S6-R005_collaborator_README.md` and package manifest/archive: collaborator handoff.
- `COPD-S6-R005_upstream_inputs.tsv`: upstream checksums.
- `COPD-S6_validation_checks.tsv`: executable consistency checks.
"""


def manifest_record(path: Path) -> Dict[str, object]:
    return {"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-missing-targets",
        action="store_true",
        help="development-only mode; production Section 6 requires Section 4 R006",
    )
    args = parser.parse_args()

    DATA.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    combined, target_evidence = load_inputs(require_targets=not args.allow_missing_targets)
    annotated = annotate_selection(combined)
    shortlist = shortlist_table(annotated)
    audit = selection_audit_table(annotated)
    constructs = construct_table(shortlist)
    candidate_plans = candidate_plan_table(shortlist)
    tf_plans = tf_plan_table(shortlist)
    contexts = context_table()
    upstream = upstream_manifest()

    write_tsv(shortlist, OUT_SHORTLIST)
    write_tsv(audit, OUT_AUDIT, gzip_output=True)
    write_tsv(constructs, OUT_CONSTRUCTS)
    write_tsv(candidate_plans, OUT_CANDIDATE_PLAN)
    write_tsv(tf_plans, OUT_TF_PLAN)
    write_tsv(contexts, OUT_CONTEXTS)
    write_tsv(upstream, OUT_UPSTREAM)
    write_collaborator_readme(shortlist, constructs)
    OUT_REPORT.write_text(
        report_text(annotated, shortlist, constructs, candidate_plans, tf_plans)
    )

    outputs = [
        OUT_SHORTLIST,
        OUT_AUDIT,
        OUT_CONSTRUCTS,
        OUT_CANDIDATE_PLAN,
        OUT_TF_PLAN,
        OUT_CONTEXTS,
        OUT_UPSTREAM,
        OUT_COLLAB_README,
        OUT_REPORT,
    ]
    input_paths = [
        R010,
        R004,
        R005,
        R006_SUMMARY,
        R006_EVIDENCE,
        R008_SUMMARY,
        R008_MOTIFS,
        S5_INTEGRATED,
        REF_FASTA,
        ALT_FASTA,
    ]
    manifest = {
        "result_id": "COPD-S6-R006",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "experimental_design_complete_wet_lab_not_performed",
        "selection_policy": SELECTION_POLICY,
        "parameters": {
            "shortlist_size": SHORTLIST_SIZE,
            "MPRA_flank_bp_each_side": MPRA_FLANK_BP,
            "LD_redundancy_definition": (
                "transitive shared linked GWAS tag or source focal tag among predicted-causal candidates; "
                "at most one shortlisted candidate per component"
            ),
            "selection_order": [
                "diversity constraints already covered by selected candidates are skipped",
                "both-model prediction",
                "enhancer-only prediction",
                "silencer-only prediction",
                "indel or complex allele",
                "refined severe-emphysema-donor regulatory overlap",
                "exact GTEx v10 Lung significant eQTL",
                "public MPRAbase element overlap",
                "motif-compatible-site creation or disruption",
                "fill by canonical Section 4 R010 predicted-causal rank",
            ],
        },
        "counts": {
            "input_predicted_causal_candidates": len(annotated),
            "LD_redundancy_components": int(annotated["ld_component_id"].nunique()),
            "shortlisted_candidates": len(shortlist),
            "MPRA_constructs": len(constructs),
            "candidate_validation_plans": len(candidate_plans),
            "TF_first_plans": len(tf_plans),
            "R006_target_evidence_rows": len(target_evidence),
        },
        "software": {"python": platform.python_version(), "pandas": pd.__version__},
        "inputs": [manifest_record(path) for path in input_paths if path.exists()],
        "outputs": [manifest_record(path) for path in outputs],
        "interpretation_guardrails": [
            "all experiments are proposed and none is represented as performed",
            "model predictions and motif matches are not causal or binding proof",
            "GTEx is molecular association evidence",
            "MPRAbase coordinate containment is regional unless allelic testing is established",
            "R006 mappings are provisional target hypotheses",
            "MPRA inserts isolate one nominated allele in a GRCh38 reference background and do not represent local phased haplotypes",
            "CRISPRi tests a region, whereas allele-level claims require endogenous allele perturbation",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(f"Input predicted-causal candidates: {len(annotated)}")
    print(f"LD-redundancy components: {annotated['ld_component_id'].nunique()}")
    print(f"Shortlisted candidates: {len(shortlist)}")
    print(f"MPRA insert designs: {len(constructs)}")
    print(f"Candidate plans: {len(candidate_plans)}")
    print(f"TF-first plans: {len(tf_plans)}")
    print(f"Report: {OUT_REPORT}")


if __name__ == "__main__":
    main()
