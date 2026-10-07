#!/usr/bin/env python3
"""Record inspected sources, append V2 registers, validate, and freeze this preflight."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from validate_preflight import REPO, GAP, STAGE, digest, run

DATE = "2026-10-07"  # investigator-facing America/New_York date
PROV = STAGE / "provenance"
PROV.mkdir(exist_ok=True)

local = [
    ("LOCAL_FMR_FREEZE", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/provenance/freeze.json", "Prior stage freeze and NOT CLEARED status"),
    ("LOCAL_FMR_REPORT", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/FINE_MAPPING_INPUT_RESOLUTION_REPORT.md", "A/B/C inputs, signed direction and method limitations"),
    ("LOCAL_R006", "diseases/COPD/04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv", "4,779 typed provisional candidate-gene rows"),
    ("LOCAL_R010", "diseases/COPD/04_modeling/results/COPD-S4-R010_THE_LIST.tsv", "Frozen 337 candidate identity/order and model fields"),
    ("LOCAL_S5_GTEX", "diseases/COPD/05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz", "Exact significant eQTL pairs and ALT slopes"),
    ("LOCAL_S5_SUMMARY", "diseases/COPD/05_computational_validation/results/COPD-S5-R005_integrated_summary.tsv", "GTEx, MPRA and Open Targets coverage"),
    ("LOCAL_S3_FUNCTIONAL", "diseases/COPD/03_regulatory_landscape/results/COPD-S3-R004_functional_variant_evidence.tsv", "Published functional variant examples"),
    ("LOCAL_PHENO", "diseases/COPD/07_gap_closure/results/COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv", "Frozen phenotype support for 337 candidates"),
    ("LOCAL_GWAS_RISK", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/tables/gwas_only_risk_direction.tsv.gz", "Previously verified signed GWAS contrasts; no new effect calculation"),
    ("LOCAL_METHOD_AUDIT", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/tables/effect_scale_summary.tsv", "A/B/C effect, SE and SPA source contracts"),
    ("LOCAL_EXT_AUDIT", "diseases/COPD/07_gap_closure/provenance/external_resource_audit.tsv", "Earlier public-resource metadata verification"),
    ("LOCAL_FRAMEWORK", "frameworks/DiseaseFramework_SL091626.md", "Workflow framework and evidence-class context"),
    ("LOCAL_GAP_PLAN", "diseases/COPD/07_gap_closure/GAP_CLOSURE_PLAN.md", "Prior target-resource feasibility and freeze boundaries"),
    ("LOCAL_MANUSCRIPT", "diseases/COPD/manuscript/COPD_regulatory_genomics_manuscript.md", "Existing exact-eQTL and provisional-target claims; read only"),
    ("LOCAL_R006_SCRIPT", "diseases/COPD/04_modeling/scripts/07_map_candidate_targets.py", "Exact GENCODE proximity, Catalog propagation and locus-mapping implementation"),
    ("LOCAL_S5_SCRIPT", "diseases/COPD/05_computational_validation/scripts/01_gtex_lung_eqtl.py", "Exact GTEx REF/ALT matching and slope field implementation"),
    ("LOCAL_S3_REPORT", "diseases/COPD/03_regulatory_landscape/results/section_3_regulatory_landscape.md", "Regulatory annotations and gene-locus overlap scope"),
    ("LOCAL_S4_REPORT", "diseases/COPD/04_modeling/results/section_4_modeling.md", "Provisional R006 target claims and model limits"),
    ("LOCAL_S5_REPORT", "diseases/COPD/05_computational_validation/results/section_5_computational_validation.md", "GTEx, MPRA and Open Targets interpretation limits"),
    ("LOCAL_S6_REPORT", "diseases/COPD/06_experimental_validation/results/section_6_experimental_validation.md", "Proposed-only target experiments"),
    ("LOCAL_PHENO_REPORT", "diseases/COPD/07_gap_closure/results/COPD-V2-PHENO_phenotype_robustness_report.md", "Direct/EHR/ML phenotype retention and cohort overlap"),
    ("LOCAL_FMR_LOCUS", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/tables/prospective_loci.tsv", "GWAS-only locus definition, not candidate derived"),
    ("LOCAL_FMR_READY", "diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/tables/track_readiness.tsv", "Prior all-NOT-CLEARED track readiness"),
]
web = [
    ("WEB_COLOC", "https://chr1swallace.github.io/coloc/", "Official coloc method overview; dense regional data and multi-signal SuSiE"),
    ("WEB_COLOC_DATA", "https://rdrr.io/cran/coloc/f/inst/doc/a02_data.Rmd", "Package data requirements; ABF LD-free under single-signal assumption; SuSiE LD"),
    ("WEB_EQTL", "https://www.ebi.ac.uk/eqtl/Data_access/", "Official complete cis-QTL, ALT effect, rsID duplicate and targeted access documentation"),
    ("WEB_EQTL_RELEASE", "https://www.ebi.ac.uk/eqtl/Release_notes/", "Release 7 and v8 prerelease status; GTEx version distinctions"),
    ("WEB_GTEX", "https://gtexportal.org/home/downloads/adult-gtex/", "GTEx significant-pair versus full-association distribution"),
    ("WEB_SAFERALI", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12481885/", "Published COPD Moloc methods, tissue Ns and 33/32-locus discrepancy"),
    ("WEB_RE2G", "https://www.encodeproject.org/annotations/ENCSR528UQX/", "ENCODE lung rE2G annotation and donor/assay context"),
    ("WEB_PCHIC", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86189", "Lung promoter capture Hi-C processed-file availability"),
    ("WEB_CELL_EQTL", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227136", "Natri processed cell-type lung QTL availability and normal/ILD scope"),
    ("WEB_NATRI", "https://doi.org/10.1038/s41588-024-01702-0", "Published genotyped sample counts and genotype access"),
]


def write_receipts() -> None:
    path = PROV / "acquisition_inspection_receipts.tsv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["receipt_id", "kind", "locator", "inspected_client_date", "retrieved_in_stage", "bytes", "sha256", "observation"])
        for rid, loc, obs in local:
            source = REPO / loc
            writer.writerow([rid, "local_read_only", loc, DATE, "no", source.stat().st_size, digest(source), obs])
        for rid, loc, obs in web:
            writer.writerow([rid, "public_metadata_page", loc, DATE, "no", "", "", obs])


def append_register(name: str, values: list[str], id_column: str) -> None:
    path = GAP / name
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        existing = {row[id_column] for row in reader}
        fieldnames = reader.fieldnames
    assert fieldnames and len(values) == len(fieldnames), (name, len(values), fieldnames)
    if values[0] in existing:
        return
    with path.open("a", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(values)


def record_baseline() -> None:
    path = PROV / "register_baseline.json"
    if path.exists():
        return
    baseline = {}
    for name in ("activity_log.tsv", "gap_closure_decision_register.tsv", "gap_closure_result_register.tsv"):
        data = (GAP / name).read_bytes()
        baseline[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    path.write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n")


def record_registers() -> None:
    rel = "target-gene-evidence-preflight-1.0/"
    append_register("activity_log.tsv", ["COPD-V2-ACT-059", DATE, "complete_preflight_audit", "Audit frozen V1 target links, GTEx exact pairs, V2 phenotype and signed GWAS contracts and public method/resource metadata", "V1 R006/R010/S5;V2 PHENO/FMR;public official pages", rel + "TARGET_GENE_EVIDENCE_PREFLIGHT_REPORT.md;" + rel + "provenance/acquisition_inspection_receipts.tsv", "No prior frozen stage altered; no downstream target or coloc run."], "activity_id")
    append_register("activity_log.tsv", ["COPD-V2-ACT-060", DATE, "frozen_stop_for_review", "Freeze allele contracts, method-specific coloc readiness, resource triage and future execution design", rel + "inventories;contracts;readiness", rel + "provenance/validation.json;" + rel + "provenance/freeze.json", "Zero execution-cleared coloc combinations; no candidate rerank, fine-mapping, manuscript edit, large download, commit or push."], "activity_id")
    append_register("gap_closure_decision_register.tsv", ["COPD-V2-DEC-063", DATE, "frozen_evidence_contract", "Keep exact eQTL, risk-expression direction, statistical coloc, contact, prediction, proximity and external biology as separate evidence classes", "Exact overlap and model/proximity consensus do not identify a causal target; allele/build/phenotype scope must be retained.", rel + "eqtl_harmonization_contract.md;" + rel + "proposed_target_evidence_hierarchy.md"], "decision_id")
    append_register("gap_closure_decision_register.tsv", ["COPD-V2-DEC-064", DATE, "no_coloc_execution_clearance", "No A/B/C × current lung eQTL package is execution-cleared for de novo colocalization", "Significant-only GTEx lacks full regions; full GTEx v8 slices and method inputs unaudited; ABF and SuSiE have distinct requirements; prior fine-map NOT CLEARED remains intact.", rel + "gwas_eqtl_readiness.tsv;" + rel + "colocalization_input_requirement_matrix.tsv"], "decision_id")
    append_register("gap_closure_result_register.tsv", ["COPD-V2-TGE-PREFLIGHT-001", DATE, "complete_preflight_stop", "target_gene_evidence_readiness", "Typed current target inventory, allele/direction rules, method-specific GWAS×QTL readiness, priority public resources and next-stage design; no biological direction results or de novo coloc.", rel + "TARGET_GENE_EVIDENCE_PREFLIGHT_REPORT.md", "no"], "result_id")
    append_register("gap_closure_result_register.tsv", ["COPD-V2-TGE-FREEZE-001", DATE, "frozen_stop_for_review", "preflight_freeze", "Validated checksum-bound preflight; all prior frozen stages and V1 337 candidates unchanged; zero execution-cleared coloc combinations.", rel + "provenance/freeze.json", "no"], "result_id")


def freeze_stage() -> None:
    result = run()
    (PROV / "validation.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if result["status"] != "PASS":
        raise RuntimeError("validation failed: " + ", ".join(t["test"] for t in result["tests"] if not t["pass"]))
    ledger = PROV / "artifact_checksums.tsv"
    files = sorted(p for p in STAGE.rglob("*") if p.is_file() and p not in {ledger, PROV / "freeze.json"} and "__pycache__" not in p.parts)
    with ledger.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for path in files:
            writer.writerow([str(path.relative_to(STAGE)), path.stat().st_size, digest(path)])
    frozen = {
        "stage": "target-gene-evidence-preflight-1.0",
        "status": "FROZEN_PREFLIGHT_STOP_FOR_INVESTIGATOR_REVIEW",
        "client_date": DATE,
        "baseline_git_commit": "07b1554e38fd8053f88fd3ad4d8abcfb6e424a7e",
        "previous_frozen_stage_sha256": digest(GAP / "fine-mapping-input-resolution-1.0/provenance/freeze.json"),
        "checksum_ledger": "provenance/artifact_checksums.tsv",
        "checksum_ledger_sha256": digest(ledger),
        "payload_file_count": len(files),
        "validation_status": result["status"],
        "coloc_execution_cleared_combinations": 0,
        "fine_mapping_execution_packages_unchanged": 0,
        "register_sha256": {name: digest(GAP / name) for name in ("activity_log.tsv", "gap_closure_decision_register.tsv", "gap_closure_result_register.tsv")},
        "boundary": "No previous frozen stage or V1 candidate order edited; no direction concordance, de novo coloc, downstream target analysis, manuscript revision, large acquisition, commit or push.",
    }
    (PROV / "freeze.json").write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {result['n_tests']} validation checks; {len(files)} frozen payload files")


if __name__ == "__main__":
    write_receipts()
    record_baseline()
    record_registers()
    freeze_stage()
