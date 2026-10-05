#!/usr/bin/env python3
"""Cross-file consistency checks for COPD Section 5 outputs."""

from __future__ import annotations

import json
from pathlib import Path

from common import (
    LOGS,
    RESULTS,
    clean_chromosome,
    ensure_directories,
    integer,
    is_true,
    load_candidates,
    read_tsv,
    sha256,
    write_tsv,
)


OUT = RESULTS / "COPD-S5_validation_checks.tsv"
LOG = LOGS / "99_validate_section5.log"


def main() -> None:
    ensure_directories()
    candidate_input, candidates = load_candidates()
    candidate_ids = {row["candidate_record_id"] for row in candidates}
    paths = {
        "gtex_audit": RESULTS / "COPD-S5-R001_GTEx_v10_Lung_candidate_audit.tsv",
        "gtex_pairs": RESULTS
        / "COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz",
        "mprabase_liftover": RESULTS
        / "COPD-S5-R002_GRCh38_to_hg19_liftover_audit.tsv",
        "mprabase_audit": RESULTS / "COPD-S5-R002_MPRAbase_candidate_audit.tsv",
        "mprabase_elements": RESULTS
        / "COPD-S5-R002_MPRAbase_v4_9_3_candidate_elements.tsv.gz",
        "gene_context": RESULTS / "COPD-S5-R003_candidate_gene_context.tsv.gz",
        "gene_associations": RESULTS
        / "COPD-S5-R003_OpenTargets_COPD_gene_associations.tsv",
        "resource_audit": RESULTS / "COPD-S5-R004_resource_access_audit.tsv",
        "biobank_audit": RESULTS / "COPD-S5-R004_biobank_access_audit.tsv",
        "integrated": RESULTS / "COPD-S5-R005_integrated_candidate_validation.tsv",
        "integrated_summary": RESULTS / "COPD-S5-R005_integrated_summary.tsv",
        "report": RESULTS / "section_5_computational_validation.md",
    }
    checks: list[dict[str, object]] = []

    def add(name: str, passed: bool, details: str) -> None:
        checks.append(
            {"check": name, "status": "PASS" if passed else "FAIL", "details": details}
        )

    for label, path in paths.items():
        add(f"output_exists:{label}", path.exists(), str(path))
    if any(not path.exists() for path in paths.values()):
        write_tsv(OUT, checks, ["check", "status", "details"])
        raise FileNotFoundError("one or more required Section 5 outputs are absent")

    gtex_audit = read_tsv(paths["gtex_audit"])
    gtex_pairs = read_tsv(paths["gtex_pairs"])
    lift = read_tsv(paths["mprabase_liftover"])
    mpra_audit = read_tsv(paths["mprabase_audit"])
    mpra_elements = read_tsv(paths["mprabase_elements"])
    context = read_tsv(paths["gene_context"])
    associations = read_tsv(paths["gene_associations"])
    integrated = read_tsv(paths["integrated"])
    integrated_summary = read_tsv(paths["integrated_summary"])

    for label, rows in (
        ("GTEx audit", gtex_audit),
        ("MPRAbase liftover audit", lift),
        ("MPRAbase candidate audit", mpra_audit),
        ("integrated matrix", integrated),
    ):
        ids = [row["candidate_record_id"] for row in rows]
        add(
            f"candidate_universe:{label}",
            len(ids) == len(set(ids)) and set(ids) == candidate_ids,
            f"rows={len(ids)} unique={len(set(ids))} expected={len(candidate_ids)}",
        )

    malformed_gtex = []
    for row in gtex_pairs:
        expected = (
            f"chr{clean_chromosome(row['chromosome_grch38'])}_"
            f"{integer(row['position_grch38'])}_{row['ref']}_{row['alt']}_b38"
        )
        if row["gtex_variant_id"] != expected:
            malformed_gtex.append(row["candidate_record_id"])
    add(
        "GTEx_exact_variant_keys",
        not malformed_gtex,
        f"pair_rows={len(gtex_pairs)} malformed={len(malformed_gtex)}",
    )

    invalid_lift = [
        row["candidate_record_id"]
        for row in lift
        if is_true(row["coordinate_match_eligible"])
        != (
            row["liftover_status"] == "unique_one_base_mapping"
            and row["roundtrip_status"] == "roundtrip_exact"
        )
    ]
    add(
        "liftover_eligibility_rule",
        not invalid_lift,
        f"audit_rows={len(lift)} invalid={len(invalid_lift)}",
    )
    eligible_ids = {
        row["candidate_record_id"]
        for row in lift
        if is_true(row["coordinate_match_eligible"])
    }
    invalid_coordinate_evidence = [
        row["candidate_record_id"]
        for row in mpra_elements
        if "roundtrip_exact_hg19_coordinate_containment" in row["match_channels"]
        and row["candidate_record_id"] not in eligible_ids
    ]
    add(
        "MPRAbase_coordinate_matches_require_roundtrip",
        not invalid_coordinate_evidence,
        f"element_rows={len(mpra_elements)} invalid={len(invalid_coordinate_evidence)}",
    )

    foreign_context = {row["candidate_record_id"] for row in context} - candidate_ids
    add(
        "gene_context_candidate_subset",
        not foreign_context,
        f"context_rows={len(context)} foreign_candidates={len(foreign_context)}",
    )
    wrong_disease = [
        row
        for row in associations
        if row["open_targets_disease_id"] != "MONDO_0005002"
        or not is_true(row["exact_disease_node"])
    ]
    add(
        "OpenTargets_exact_COPD_node",
        not wrong_disease,
        f"gene_rows={len(associations)} wrong_scope={len(wrong_disease)}",
    )

    ranks = [
        integer(row["predicted_causal_priority_rank"]) or 10**12
        for row in integrated
    ]
    add(
        "integrated_preserves_priority_order",
        ranks == sorted(ranks),
        f"rows={len(ranks)}",
    )
    summary_map = {row["metric"]: row for row in integrated_summary}
    exact_eqtl_count = sum(
        is_true(row["gtex_lung_exact_significant_eqtl"]) for row in integrated
    )
    summary_eqtl = integer(
        summary_map.get("exact_significant_GTEx_v10_Lung_cis_eqtl", {}).get("n", "")
    )
    add(
        "integrated_GTEx_count_reconciles",
        exact_eqtl_count == summary_eqtl,
        f"matrix={exact_eqtl_count} summary={summary_eqtl}",
    )
    mpra_count = sum(
        is_true(row["mprabase_any_element_evidence"]) for row in integrated
    )
    summary_mpra = integer(
        summary_map.get("any_MPRAbase_element_evidence", {}).get("n", "")
    )
    add(
        "integrated_MPRAbase_count_reconciles",
        mpra_count == summary_mpra,
        f"matrix={mpra_count} summary={summary_mpra}",
    )
    report_text = paths["report"].read_text()
    markers = ("TODO", "TBD")
    add(
        "report_has_no_placeholder_markers",
        not any(marker in report_text for marker in markers),
        "checked TODO and TBD",
    )

    current_input_sha256 = sha256(candidate_input)
    for result_id in (
        "COPD-S5-R001",
        "COPD-S5-R002",
        "COPD-S5-R003",
        "COPD-S5-R004",
        "COPD-S5-R005",
    ):
        manifest_path = RESULTS / f"{result_id}_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        mismatches = []
        for label, record in manifest.get("outputs", {}).items():
            if not isinstance(record, dict) or "path" not in record:
                continue
            output_path = Path(str(record["path"]))
            if (
                not output_path.exists()
                or output_path.stat().st_size != int(record["bytes"])
                or sha256(output_path) != record["sha256"]
            ):
                mismatches.append(label)
        add(
            f"manifest_output_hashes:{result_id}",
            not mismatches,
            f"outputs={len(manifest.get('outputs', {}))} mismatches={';'.join(mismatches)}",
        )
        section4_record = manifest.get("inputs", {}).get("section4_candidates")
        if isinstance(section4_record, dict):
            add(
                f"manifest_uses_current_Section4_input:{result_id}",
                section4_record.get("sha256") == current_input_sha256,
                str(section4_record.get("path", "")),
            )

    write_tsv(OUT, checks, ["check", "status", "details"])
    failed = [row for row in checks if row["status"] == "FAIL"]
    LOG.write_text(
        f"candidate_input={candidate_input}\n"
        f"checks={len(checks)}\n"
        f"passed={len(checks) - len(failed)}\n"
        f"failed={len(failed)}\n"
        + "".join(f"FAIL\t{row['check']}\t{row['details']}\n" for row in failed)
    )
    if failed:
        raise RuntimeError(f"{len(failed)} Section 5 validation checks failed")
    print(f"Section 5 validation: {len(checks)}/{len(checks)} checks passed.")


if __name__ == "__main__":
    main()
