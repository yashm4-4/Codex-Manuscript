#!/usr/bin/env python3
"""Validate COPD Section 6 design, sequence, and evidence reconciliation."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases/COPD/06_experimental_validation"
RESULTS = SECTION / "results"
DATA = SECTION / "data"
S4 = ROOT / "diseases/COPD/04_modeling/results"

SHORTLIST = RESULTS / "COPD-S6-R001_candidate_shortlist.tsv"
AUDIT = RESULTS / "COPD-S6-R002_selection_audit.tsv.gz"
CONSTRUCTS = RESULTS / "COPD-S6-R003_MPRA_constructs.tsv"
CANDIDATE_PLAN = RESULTS / "COPD-S6-R004_candidate_validation_plan.tsv"
TF_PLAN = RESULTS / "COPD-S6-R004_TF_first_plan.tsv"
CONTEXTS = RESULTS / "COPD-S6-R004_cell_context_controls.tsv"
UPSTREAM = DATA / "COPD-S6-R005_upstream_inputs.tsv"
COLLAB_README = RESULTS / "COPD-S6-R005_collaborator_README.md"
REPORT = RESULTS / "section_6_experimental_validation.md"
MANIFEST = RESULTS / "COPD-S6-R006_analysis_manifest.json"
OUT = RESULTS / "COPD-S6_validation_checks.tsv"

R004 = S4 / "COPD-S4-R004_predicted_causal_regulatory_variants.tsv"
R006_EVIDENCE = S4 / "COPD-S4-R006_candidate_target_evidence.tsv"
R010 = S4 / "COPD-S4-R010_THE_LIST.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main() -> None:
    checks: List[Dict[str, object]] = []

    def check(name: str, passed: bool, observed: object, expected: object, note: str = "") -> None:
        checks.append(
            {
                "check_id": f"COPD-S6-C{len(checks) + 1:03d}",
                "check_name": name,
                "status": "PASS" if passed else "FAIL",
                "observed": observed,
                "expected": expected,
                "note": note,
            }
        )

    required = [
        SHORTLIST,
        AUDIT,
        CONSTRUCTS,
        CANDIDATE_PLAN,
        TF_PLAN,
        CONTEXTS,
        UPSTREAM,
        COLLAB_README,
        REPORT,
        MANIFEST,
    ]
    absent = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    check("required_design_outputs_exist", not absent, ";".join(absent) or "all_present", "all_present")
    if absent:
        pd.DataFrame(checks).to_csv(OUT, sep="\t", index=False)
        raise SystemExit(1)

    shortlist = pd.read_csv(SHORTLIST, sep="\t", low_memory=False)
    audit = pd.read_csv(AUDIT, sep="\t", low_memory=False)
    constructs = pd.read_csv(CONSTRUCTS, sep="\t", low_memory=False)
    candidate_plan = pd.read_csv(CANDIDATE_PLAN, sep="\t", low_memory=False)
    tf_plan = pd.read_csv(TF_PLAN, sep="\t", low_memory=False)
    contexts = pd.read_csv(CONTEXTS, sep="\t", low_memory=False)
    upstream = pd.read_csv(UPSTREAM, sep="\t", low_memory=False)
    r004 = pd.read_csv(R004, sep="\t", low_memory=False)
    r006_evidence = pd.read_csv(R006_EVIDENCE, sep="\t", low_memory=False)
    r010 = pd.read_csv(R010, sep="\t", low_memory=False)

    check("shortlist_has_12_rows", len(shortlist) == 12, len(shortlist), 12)
    check(
        "shortlist_candidate_ids_unique",
        shortlist["candidate_record_id"].nunique() == len(shortlist),
        shortlist["candidate_record_id"].nunique(),
        len(shortlist),
    )
    ranks = sorted(pd.to_numeric(shortlist["experimental_shortlist_rank"]).astype(int).tolist())
    check("shortlist_ranks_contiguous", ranks == list(range(1, 13)), ranks, list(range(1, 13)))
    check(
        "one_candidate_per_LD_component",
        shortlist["ld_component_id"].nunique() == len(shortlist),
        shortlist["ld_component_id"].nunique(),
        len(shortlist),
    )
    check(
        "audit_reconciles_Section4_candidates",
        len(audit) == len(r004)
        and set(audit["candidate_record_id"].astype(str))
        == set(r004["candidate_record_id"].astype(str)),
        len(audit),
        len(r004),
    )
    r010_rank = r010.set_index("candidate_record_id")["predicted_causal_priority_rank"].astype(int)
    shortlist_rank_map = shortlist.set_index("candidate_record_id")[
        "predicted_causal_priority_rank"
    ].astype(int)
    canonical_rank_ok = all(
        int(shortlist_rank_map[candidate_id]) == int(r010_rank[candidate_id])
        for candidate_id in shortlist_rank_map.index
    )
    check(
        "shortlist_preserves_canonical_R010_rank",
        canonical_rank_ok,
        canonical_rank_ok,
        True,
    )
    check(
        "shortlist_names_canonical_R010_source",
        set(shortlist["canonical_section4_list"].astype(str))
        == {"diseases/COPD/04_modeling/results/COPD-S4-R010_THE_LIST.tsv"}
        and shortlist["canonical_rank_preserved"].map(as_bool).all(),
        ";".join(sorted(set(shortlist["canonical_section4_list"].astype(str)))),
        "diseases/COPD/04_modeling/results/COPD-S4-R010_THE_LIST.tsv",
    )
    selected_audit = audit[audit["selected_for_experimental_shortlist"].map(as_bool)]
    check(
        "audit_selected_set_reconciles_shortlist",
        set(selected_audit["candidate_record_id"].astype(str))
        == set(shortlist["candidate_record_id"].astype(str)),
        len(selected_audit),
        len(shortlist),
    )
    check(
        "all_shortlist_sequences_eligible",
        all(shortlist["candidate_record_id"].astype(str).isin(
            audit.loc[audit["experimental_eligible"].map(as_bool), "candidate_record_id"].astype(str)
        )),
        int(
            shortlist["candidate_record_id"].astype(str).isin(
                audit.loc[
                    audit["experimental_eligible"].map(as_bool), "candidate_record_id"
                ].astype(str)
            ).sum()
        ),
        len(shortlist),
    )

    criteria = [
        "both_model_prediction",
        "enhancer_only_prediction",
        "silencer_only_prediction",
        "indel_or_complex",
        "severe_emphysema_donor_refined_overlap",
        "exact_GTEx_v10_Lung_eQTL",
        "public_MPRAbase_element_overlap",
        "predicted_motif_creation_or_disruption",
    ]
    selected_coverage = ";".join(shortlist["diversity_coverage_criteria"].fillna("").astype(str))
    eligible_coverage = ";".join(
        audit.loc[
            audit["experimental_eligible"].map(as_bool), "diversity_coverage_criteria"
        ].fillna("").astype(str)
    )
    for criterion in criteria:
        available = criterion in eligible_coverage
        covered = criterion in selected_coverage
        check(
            f"diversity_criterion_{criterion}",
            (not available) or covered,
            "covered" if covered else "not_covered",
            "covered_when_eligible",
        )

    check(
        "R006_targets_consumed",
        set(shortlist["target_data_status"].astype(str)) == {"R006_consumed"},
        ";".join(sorted(set(shortlist["target_data_status"].astype(str)))),
        "R006_consumed",
    )
    check("R006_evidence_row_count", len(r006_evidence) == 4779, len(r006_evidence), 4779)

    ids = set(shortlist["candidate_record_id"].astype(str))
    check("construct_count", len(constructs) == 4 * len(ids), len(constructs), 4 * len(ids))
    check(
        "construct_candidate_set_reconciles",
        set(constructs["candidate_record_id"].astype(str)) == ids,
        constructs["candidate_record_id"].nunique(),
        len(ids),
    )
    per_candidate = constructs.groupby("candidate_record_id")["construct_id"].count()
    check(
        "four_constructs_per_candidate",
        len(per_candidate) == len(ids) and (per_candidate == 4).all(),
        ";".join(f"{key}:{value}" for key, value in per_candidate.items()),
        "4_each",
    )
    check(
        "construct_ids_unique",
        constructs["construct_id"].nunique() == len(constructs),
        constructs["construct_id"].nunique(),
        len(constructs),
    )

    alphabet_ok = constructs["insert_sequence"].astype(str).map(
        lambda sequence: set(sequence) <= set("ACGT")
    )
    check("construct_sequence_alphabet", alphabet_ok.all(), int(alphabet_ok.sum()), len(constructs))
    lengths_ok = constructs.apply(
        lambda row: len(str(row["insert_sequence"]))
        == int(row["left_flank_bp"])
        + len(str(row["genomic_forward_allele"]))
        + int(row["right_flank_bp"])
        == int(row["insert_length_bp"]),
        axis=1,
    )
    check("construct_lengths_reconcile", lengths_ok.all(), int(lengths_ok.sum()), len(constructs))

    center_ok = constructs.apply(
        lambda row: str(row["insert_sequence"])[
            int(row["allele_start_index_zero_based"]) : int(row["allele_start_index_zero_based"])
            + len(str(row["oriented_allele_sequence"]))
        ]
        == str(row["oriented_allele_sequence"]),
        axis=1,
    )
    check("construct_center_alleles_reconcile", center_ok.all(), int(center_ok.sum()), len(constructs))

    orientation_ok = True
    flank_ok = True
    for candidate_id, candidate_frame in constructs.groupby("candidate_record_id"):
        for allele in ("REF", "ALT"):
            allele_rows = candidate_frame[candidate_frame["allele_label"] == allele]
            if set(allele_rows["orientation"]) != {"forward", "reverse_complement"}:
                orientation_ok = False
                continue
            forward = str(
                allele_rows.loc[allele_rows["orientation"] == "forward", "insert_sequence"].iloc[0]
            )
            reverse = str(
                allele_rows.loc[
                    allele_rows["orientation"] == "reverse_complement", "insert_sequence"
                ].iloc[0]
            )
            if reverse_complement(forward) != reverse:
                orientation_ok = False
        forward_rows = candidate_frame[candidate_frame["orientation"] == "forward"].set_index(
            "allele_label"
        )
        ref_seq = str(forward_rows.loc["REF", "insert_sequence"])
        alt_seq = str(forward_rows.loc["ALT", "insert_sequence"])
        ref_allele = str(forward_rows.loc["REF", "genomic_forward_allele"])
        alt_allele = str(forward_rows.loc["ALT", "genomic_forward_allele"])
        flank = int(forward_rows.loc["REF", "left_flank_bp"])
        if ref_seq[:flank] != alt_seq[:flank] or ref_seq[flank + len(ref_allele) :] != alt_seq[flank + len(alt_allele) :]:
            flank_ok = False
    check("reverse_complements_exact", orientation_ok, orientation_ok, True)
    check("REF_ALT_genomic_flanks_exact", flank_ok, flank_ok, True)

    check(
        "candidate_plans_reconcile",
        len(candidate_plan) == len(ids)
        and set(candidate_plan["candidate_record_id"].astype(str)) == ids,
        len(candidate_plan),
        len(ids),
    )
    expected_tf_plans = int(
        shortlist["predicted_causal_enhancer"].map(as_bool).sum()
        + shortlist["predicted_causal_silencer"].map(as_bool).sum()
    )
    check("TF_plan_count", len(tf_plan) == expected_tf_plans, len(tf_plan), expected_tf_plans)
    check(
        "TF_plan_candidate_set",
        set(tf_plan["candidate_record_id"].astype(str)) == ids,
        tf_plan["candidate_record_id"].nunique(),
        len(ids),
    )
    check(
        "required_cell_contexts_present",
        {
            "airway_epithelium",
            "alveolar_AT2",
            "parenchymal_lung_fibroblast",
            "assay_development_cell_lines",
        }.issubset(set(contexts["context_id"].astype(str))),
        ";".join(sorted(contexts["context_id"].astype(str))),
        "airway_epithelium;alveolar_AT2;parenchymal_lung_fibroblast;assay_development_cell_lines",
    )
    plan_statuses = set(candidate_plan["experimental_status"].astype(str)) | set(
        tf_plan["experimental_status"].astype(str)
    ) | set(contexts["experimental_status"].astype(str))
    check(
        "all_experimental_plans_labeled_proposed",
        plan_statuses == {"proposed_not_performed"},
        ";".join(sorted(plan_statuses)),
        "proposed_not_performed",
    )

    upstream_ok = True
    bad_upstream = []
    for row in upstream.to_dict(orient="records"):
        path = ROOT / str(row["path"])
        if not path.exists() or sha256(path) != str(row["sha256"]):
            upstream_ok = False
            bad_upstream.append(str(row["path"]))
    check(
        "upstream_checksums_current",
        upstream_ok,
        ";".join(bad_upstream) or "all_match",
        "all_match",
    )

    report_text = REPORT.read_text()
    readme_text = COLLAB_README.read_text()
    guardrail_ok = all(
        phrase in (report_text + "\n" + readme_text).lower()
        for phrase in ["proposed", "not causal proof", "not tf binding evidence"]
    )
    check("interpretation_guardrails_present", guardrail_ok, guardrail_ok, True)
    placeholders = [token for token in ("TODO", "TBD", "PLACEHOLDER") if token in report_text or token in readme_text]
    check("no_unresolved_text_placeholders", not placeholders, ";".join(placeholders) or "none", "none")

    output = pd.DataFrame(checks)
    output.to_csv(OUT, sep="\t", index=False)
    failures = output[output["status"] != "PASS"]
    print(f"Validation checks: {len(output)}")
    print(f"PASS: {(output['status'] == 'PASS').sum()}")
    print(f"FAIL: {len(failures)}")
    print(f"Output: {OUT}")
    if not failures.empty:
        print(failures[["check_id", "check_name", "observed", "expected"]].to_string(index=False))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
