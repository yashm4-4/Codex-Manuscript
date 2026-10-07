#!/usr/bin/env python3
"""Validate the bounded COPD target-gene preflight without rerunning analyses."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path


STAGE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[5]
GAP = STAGE.parent


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        result = list(reader)
        assert reader.fieldnames and all(None not in r for r in result), path
        return result


def run() -> dict:
    tests = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        tests.append({"test": name, "pass": bool(condition), "detail": detail})

    required = [
        "TARGET_GENE_EVIDENCE_PREFLIGHT_REPORT.md",
        "current_target_evidence_inventory.tsv",
        "eqtl_harmonization_contract.md",
        "risk_expression_direction_contract.md",
        "colocalization_input_requirement_matrix.tsv",
        "gwas_eqtl_readiness.tsv",
        "candidate_external_resource_inventory.tsv",
        "proposed_target_evidence_hierarchy.md",
        "future_execution_specification.md",
        "provenance/acquisition_inspection_receipts.tsv",
    ]
    check("all_required_files", all((STAGE / p).is_file() for p in required))
    inv = rows(STAGE / "current_target_evidence_inventory.tsv")
    ext = rows(STAGE / "candidate_external_resource_inventory.tsv")
    ready = rows(STAGE / "gwas_eqtl_readiness.tsv")
    matrix = rows(STAGE / "colocalization_input_requirement_matrix.tsv")
    receipts = rows(STAGE / "provenance/acquisition_inspection_receipts.tsv")
    check("typed_inventories", len(inv) >= 12 and len(ext) >= 5 and len(matrix) >= 10)
    check("readiness_coverage", len(ready) == 21 and {r["gwas_track"] for r in ready} == {"A", "B", "C"})
    check("readiness_statuses", all(r["status"] in {"NOT CLEARED", "INCONCLUSIVE", "EXECUTION-CLEARED"} for r in ready))
    check("no_coloc_execution_clearance", all(r["status"] != "EXECUTION-CLEARED" for r in ready))
    check("resource_receipts", len(receipts) >= 15 and all(r["locator"] for r in receipts))
    check("local_receipt_hashes", all(digest(REPO / r["locator"]) == r["sha256"] for r in receipts if r["kind"] == "local_read_only"))

    r006 = rows(REPO / "diseases/COPD/04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv")
    methods = Counter(r["mapping_method"] for r in r006)
    check("V1_R006_unchanged_counts", len(r006) == 4779 and len({r["candidate_record_id"] for r in r006}) == 337 and methods["coding_CDS_overlap"] == 15)
    r010 = rows(REPO / "diseases/COPD/04_modeling/results/COPD-S4-R010_THE_LIST.tsv")
    check("frozen_candidate_order", len(r010) == 337 and [int(r["predicted_causal_priority_rank"]) for r in r010] == list(range(1, 338)))
    with gzip.open(REPO / "diseases/COPD/05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz", "rt") as fh:
        gtex = list(csv.DictReader(fh, delimiter="\t"))
    check("V1_GTEx_unchanged_counts", len(gtex) == 915 and len({r["candidate_record_id"] for r in gtex}) == 185 and len({r["gene_id_versionless"] for r in gtex}) == 224)
    phen = rows(GAP / "results/COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv")
    check("phenotype_strata_preserved", len(phen) == 337 and sum(r["primary_retained"] == "True" for r in phen) == 184 and sum(r["support_state"] == "ml_surrogate_only" for r in phen) == 124)
    old = GAP / "fine-mapping-input-resolution-1.0"
    freeze = json.loads((old / "provenance/freeze.json").read_text())
    check("prior_frozen_not_cleared", freeze["status"] == "FROZEN_INPUT_RESOLUTION_STOP_FOR_INVESTIGATOR_REVIEW" and len(rows(old / "tables/track_readiness.tsv")) == 3)
    check("prior_freeze_checksum", digest(old / "provenance/freeze.json") == next(r["sha256"] for r in receipts if r["receipt_id"] == "LOCAL_FMR_FREEZE"))

    base = json.loads((STAGE / "provenance/register_baseline.json").read_text())
    for key, info in base.items():
        payload = (GAP / key).read_bytes()
        check(f"append_only_{key}", hashlib.sha256(payload[: info["bytes"]]).hexdigest() == info["sha256"])
    changed = subprocess.run(["git", "diff", "--name-only", "HEAD", "--", "diseases/COPD"], cwd=REPO, text=True, capture_output=True, check=True).stdout.splitlines()
    allowed = {f"diseases/COPD/07_gap_closure/{name}" for name in base}
    stage_prefix = "diseases/COPD/07_gap_closure/target-gene-evidence-preflight-1.0/"
    check("no_tracked_frozen_stage_edits", all(path in allowed or path.startswith(stage_prefix) for path in changed), ",".join(changed))
    check("contracts_before_results", "No full candidate join or final concordance" in (STAGE / "TARGET_GENE_EVIDENCE_PREFLIGHT_REPORT.md").read_text())
    return {"status": "PASS" if all(t["pass"] for t in tests) else "FAIL", "tests": tests, "n_tests": len(tests)}


if __name__ == "__main__":
    result = run()
    output = STAGE / "provenance/validation.json"
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"], result["n_tests"], "checks")
    if result["status"] != "PASS":
        raise SystemExit(1)
