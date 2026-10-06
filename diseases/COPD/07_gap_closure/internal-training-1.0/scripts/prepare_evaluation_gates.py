#!/usr/bin/env python3
"""Derive hash-bound chr7 evaluation/calibration gates; never run either stage.

Use --phase evaluation only after all18 selected checkpoints and actual-network
invariance are frozen. Use --phase calibration only after BOTH C contexts PASS.
No TensorFlow, checkpoint loading, model scores, test data or external data are
opened. Existing derived gates are never overwritten.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
from pathlib import Path
import sys

from evaluate_chr7 import (CONFIGURATIONS, MODELS, SEEDS, read_json, read_rows,
                           sha256, validate_prerequisites, write_json)


def fail(message):
    raise ValueError(message)


def absolute_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}


def verified_record(base, record):
    base = Path(base).resolve()
    path = (base / record["path"]).resolve()
    if not path.is_relative_to(base):
        fail(f"Record outside permitted source root: {path}")
    if sha256(path) != record["sha256"] or ("bytes" in record and path.stat().st_size != int(record["bytes"])):
        fail(f"Source hash/size binding failed: {path}")
    return absolute_record(path)


def record_map(base, records):
    result = {}
    for record in records:
        key = str((Path(base) / record["path"]).resolve())
        if key in result:
            fail(f"Duplicate source reference: {key}")
        result[key] = record
    return result


def required_record(base, records, path):
    key = str(Path(path).resolve())
    if key not in records:
        fail(f"Missing original source binding: {key}")
    return verified_record(base, records[key])


def same_reference(left, right, description):
    if left["path"] != right["path"] or left["sha256"] != right["sha256"] or left["bytes"] != right["bytes"]:
        fail(f"Source reference disagreement: {description}")


def frozen_construction_sources(stage, repo):
    frozen = stage.parent / "pretraining-1.1"
    initial_path = stage / "provenance/input_integrity_verified_before.json"
    initial = read_json(initial_path)
    if initial.get("status") != "PASS" or initial.get("failures") != 0:
        fail("Initial independent immutable-input audit not PASS")
    checks_ref = verified_record(repo, {"path": initial["checks_path"], "sha256": initial["checks_sha256"]})
    checks = {row["check"]: row for row in read_rows(checks_ref["path"])}
    freeze_path = frozen / "provenance/freeze.json"
    ledger_path = frozen / "provenance/artifact_checksums.tsv"
    for key, path in (("pretraining-1.1:freeze:sha256", freeze_path), ("pretraining-1.1:ledger:sha256", ledger_path)):
        row = checks.get(key, {})
        if row.get("status") != "PASS" or row.get("observed") != row.get("expected") or row.get("expected") != sha256(path):
            fail(f"Original pre-execution hash anchor failed: {key}")
    freeze = read_json(freeze_path)
    if (freeze.get("version") != "pretraining-1.1" or freeze.get("status") != "READY_FOR_INVESTIGATOR_TRAINING_REVIEW"
            or any(freeze.get(name) != 0 for name in ("B_C_balance_failures", "artifact_failures", "readiness_failures", "packaging_failures", "pair_metadata_failures"))):
        fail("Frozen pretraining-1.1 construction is not ready")
    if sha256(ledger_path) != freeze["artifact_checksum_ledger_sha256"]:
        fail("Frozen ledger reference mismatch")
    ledger = record_map(repo, list(read_rows(ledger_path)))
    if len(ledger) != freeze["payload_count"]:
        fail("Frozen payload count differs")
    spec = required_record(repo, ledger, frozen / "specification/selection_adequacy.json")
    validation = required_record(repo, ledger, frozen / "provenance/validation_corrected_final_manifest.json")
    if validation["sha256"] != freeze["validation_manifest_sha256"]:
        fail("Construction validation not bound to original freeze")
    validation_obj = read_json(validation["path"])
    if validation_obj.get("artifact_failures") != 0 or validation_obj.get("readiness_failures") != 0:
        fail("Frozen independent construction validation failed")
    for record in validation_obj["outputs"]:
        current = verified_record(repo, record)
        same_reference(current, required_record(repo, ledger, current["path"]), "validation output versus frozen payload ledger")
    readiness = required_record(repo, ledger, frozen / "provenance/training_readiness.tsv")
    rows = [row for row in read_rows(readiness["path"]) if row["scope"] == "construction_review"]
    required = {"B_prespecified_balance_gates", "C_prespecified_balance_gates", "B_positive_retention", "C_positive_retention", "C_common_selection_panels", "C_calibration_control_panels", "independent_validation"}
    if not required <= {row["condition"] for row in rows} or any(not row["status"].startswith("SATISFIED") for row in rows):
        fail("Frozen B/C construction readiness conditions not satisfied")
    return spec, validation, ledger, {"initial_integrity_audit": absolute_record(initial_path), "initial_integrity_checks": checks_ref,
            "pretraining_1_1_freeze": absolute_record(freeze_path), "pretraining_1_1_checksum_ledger": absolute_record(ledger_path), "training_readiness": readiness}


def validate_audit_rows(audit, checkpoints, validation_counts):
    if (audit.get("status") != "PASS" or audit.get("real_network") is not True
            or audit.get("all_configurations_all_seeds_and_ensembles") is not True
            or audit.get("n_seed_audits") != 18 or audit.get("n_ensemble_audits") != 6):
        fail("Actual-network invariance summary not PASS for18seeds+6ensembles")
    expected = {(config, model, "seed", seed) for config in CONFIGURATIONS for model in MODELS for seed in SEEDS}
    expected |= {(config, model, "ensemble", "all3") for config in CONFIGURATIONS for model in MODELS}
    seen = set()
    for row in audit["audit_rows"]:
        key = (row["configuration"], row["model"], row["unit"], row["seed"])
        if key not in expected or key in seen:
            fail("Unexpected/duplicate actual-network invariance audit row")
        seen.add(key)
        if (row.get("status") != "PASS" or row.get("real_network") is not True or row.get("n_failed") != 0
                or row.get("n_sequences") != validation_counts[row["model"]]
                or row.get("atol") != 1e-6 or row.get("rtol") != 1e-6):
            fail("Real-network invariance row/count/tolerance gate failed")
        residual = float(row["max_tolerance_normalized_difference"])
        absolute = float(row["max_absolute_difference"])
        if not math.isfinite(residual) or not math.isfinite(absolute) or residual < 0 or residual > 1 or absolute < 0:
            fail("Real-network numeric residual is invalid/failing")
        if row["unit"] == "seed":
            expected_checkpoint = checkpoints[(row["configuration"], row["model"], int(row["seed"]))]
            if row["checkpoint"]["sha256"] != expected_checkpoint["sha256"]:
                fail("Invariance used a different selected checkpoint")
        elif tuple(row.get("fixed_seed_order", ())) != SEEDS or row.get("reduction_dtype") != "float64":
            fail("Real-network ensemble rule drift")
    if seen != expected:
        fail("Missing actual-network invariance audit units")


def build_evaluation_manifest(stage, repo):
    """Read metadata/provenance and opaque hashes only; do not read score tables."""
    spec, construction, frozen_ledger, original = frozen_construction_sources(stage, repo)
    gate_path = stage / "provenance/execution_input_gate.json"
    gate = read_json(gate_path)
    if gate.get("status") != "PASS" or gate.get("input_hashes_verified") is not True or gate.get("training_authorized") is not True:
        fail("Original authorized execution input gate not PASS")
    gate_records = record_map(stage, gate["files"])
    input_ref = required_record(stage, gate_records, stage / "inputs/input_manifest.json")
    input_manifest = read_json(input_ref["path"])
    if input_manifest.get("status") != "PASS" or input_manifest.get("completed") is not True:
        fail("Original prepared input manifest not PASS")
    for record in input_manifest["frozen_sources"]:
        current = verified_record(repo, record)
        same_reference(current, required_record(repo, frozen_ledger, current["path"]), "prepared source versus original frozen payload")
    prepared_files = record_map(stage, input_manifest["files"])
    common = {}
    for model in MODELS:
        for suffix in ("selection", "calibration", "train_C_positive_covariates"):
            relative = f"inputs/common/{model}_{suffix}.tsv.gz"
            common[relative] = required_record(stage, prepared_files, stage / relative)
    freeze_path = stage / "provenance/checkpoint_freeze.json"
    freeze = read_json(freeze_path)
    if (freeze.get("status") != "PASS" or freeze.get("all_18_selected_checkpoints_frozen") is not True
            or freeze.get("selection_prediction_authorized") is not True):
        fail("All18 selected checkpoint freeze not PASS")
    same_reference(verified_record(stage, freeze["execution_input_gate"]), absolute_record(gate_path), "checkpointfreeze execution gate")
    same_reference(verified_record(stage, freeze["input_manifest"]), input_ref, "checkpointfreeze prepared inputs")
    checkpoints, checkpoint_sources = {}, []
    expected = {(config, model, seed) for config in CONFIGURATIONS for model in MODELS for seed in SEEDS}
    for record in freeze["checkpoints"]:
        key = (record["configuration"], record["model"], int(record["seed"]))
        if key not in expected or key in checkpoints or record.get("completed") is not True:
            fail("Unexpected/duplicated/incomplete frozen checkpoint")
        checkpoint = verified_record(stage, record)
        source = verified_record(stage, record["selected_checkpoint_record"])
        selected = read_json(source["path"])
        if (selected.get("configuration"), selected.get("model"), selected.get("seed")) != key or selected.get("completed") is not True:
            fail("Checkpoint selection record identity/completion mismatch")
        same_reference(checkpoint, verified_record(stage, selected["checkpoint"]), "original selected checkpoint")
        source_bindings = {"selected_checkpoint_record": source}
        for name in ("completed_record", "history", "training_contract"):
            source_bindings[name] = verified_record(stage, record[name])
        completed = read_json(source_bindings["completed_record"]["path"])
        if completed.get("status") != "COMPLETED":
            fail("Original training completion record not COMPLETED")
        contract = read_json(source_bindings["training_contract"]["path"])
        same_reference(verified_record(stage, contract["input_gate"]), absolute_record(gate_path), "run original execution gate")
        checkpoints[key] = {"configuration": key[0], "model": key[1], "seed": key[2], "completed": True, **checkpoint}
        checkpoint_sources.append({"configuration": key[0], "model": key[1], "seed": key[2], **source_bindings})
    if set(checkpoints) != expected:
        fail("Not all18 original selected checkpoints accounted for")
    audit_path = stage / "predictions/real_network_invariance.json"
    audit = read_json(audit_path)
    completed_path = stage / "predictions/completed.json"
    completed = read_json(completed_path)
    if completed.get("status") != "COMPLETED" or completed.get("invariance_status") != "PASS":
        fail("Chr7 prediction stage did not complete with invariance PASS")
    same_reference(verified_record(stage, completed["invariance_record"]), absolute_record(audit_path), "prediction completion invariant audit")
    same_reference(verified_record(stage, audit["checkpoint_freeze"]), absolute_record(freeze_path), "invariance checkpointfreeze")
    cache_ref = required_record(stage, gate_records, stage / "cache/cache_manifest.json")
    same_reference(verified_record(stage, audit["cache_manifest"]), cache_ref, "invariance cache manifest")
    audit_scripts = record_map(stage, audit["script_files"])
    implementation_files = []
    for name in ("predict_chr7.py", "phase_two_contract.py"):
        path = stage / "scripts" / name
        current = required_record(stage, audit_scripts, path)
        same_reference(current, required_record(stage, gate_records, path), "actual inference implementation versus original execution gate")
        implementation_files.append(current)
    validation_counts = {}
    for record in audit["frozen_input_tables"]:
        current = verified_record(stage, record)
        same_reference(current, required_record(stage, prepared_files, current["path"]), "invariance validation table versus prepared inputs")
    validation_refs = record_map(stage, audit["frozen_input_tables"])
    for model in MODELS:
        current = required_record(stage, validation_refs, stage / f"inputs/validation/{model}_chr7.tsv.gz")
        count = 0
        for row in read_rows(current["path"]):
            if row["chrom"] != "chr7":
                fail("Nonchr7 metadata in invariance validation scope")
            count += 1
        if not count:
            fail("Empty actual-network validation scope")
        validation_counts[model] = count
    validate_audit_rows(audit, checkpoints, validation_counts)
    for row in audit["audit_rows"]:
        if row["unit"] == "seed":
            expected_checkpoint = checkpoints[(row["configuration"], row["model"], int(row["seed"]))]
            same_reference(verified_record(stage, row["checkpoint"]), {k: expected_checkpoint[k] for k in ("path", "sha256", "bytes")}, "invariance row checkpoint")
    if any(audit.get(key) is not False for key in ("selection_metrics_computed", "calibration_performed", "chr8_chr9_access", "external_benchmark_access")):
        fail("Prediction-stage scientific firewall record failed")
    prediction_records = record_map(stage, audit["files"])
    completion_records = record_map(stage, completed["files"])
    if set(prediction_records) != set(completion_records):
        fail("Prediction/audit output membership differs")
    predictions = {}
    for model in MODELS:
        for suffix in ("chr7", "independent_rc_pass"):
            path = stage / f"predictions/{model}_{suffix}.tsv.gz"
            current = required_record(stage, prediction_records, path)
            same_reference(current, required_record(stage, completion_records, path), "prediction completion file binding")
            if suffix == "chr7":
                predictions[model] = current
    scripts = [required_record(stage, gate_records, stage / "scripts" / name) for name in ("evaluate_chr7.py", "calibrate_chr7.py", "prepare_evaluation_gates.py")]
    implementation_specification = required_record(stage, gate_records, stage / "specification/internal_evaluation_implementation.json")
    manifest = {"version": "internal-training-1.0", "gate_phase": "evaluation", "specification": spec,
        "prerequisites": {"construction_gates_pass": True, "construction_evidence": construction, "checkpoint_freeze": absolute_record(freeze_path),
            "invariance": {"status": "PASS", "real_network": True, "all_configurations_all_seeds_and_ensembles": True, "evidence": absolute_record(audit_path)},
            "checkpoints": [checkpoints[key] for key in sorted(checkpoints)]},
        "models": {model: {"selection": common[f"inputs/common/{model}_selection.tsv.gz"], "train_covariates": common[f"inputs/common/{model}_train_C_positive_covariates.tsv.gz"], "predictions": predictions[model]} for model in MODELS},
        "source_bindings": {**original, "execution_input_gate": absolute_record(gate_path), "prepared_input_manifest": input_ref,
            "prepared_common_files": common, "checkpoint_source_records": checkpoint_sources, "prediction_completion": absolute_record(completed_path),
            "actual_inference_implementation_files": implementation_files, "evaluation_implementation_files": scripts,
            "evaluation_implementation_specification": implementation_specification},
        "planned_output_directory": str(stage / "results/chr7_evaluation"),
        "no_evaluation_calibration_or_model_execution_by_gate_helper": True}
    validate_prerequisites(manifest)
    return manifest


def validate_passed_adequacy(decision, gates, bootstrap):
    if decision.get("both_C_contexts_pass") is not True or set(decision.get("model_contexts", {})) != set(MODELS):
        fail("STOP: BOTH C contexts must PASS before any calibration gate is created")
    simple = {"frozen_construction", "real_network_invariance", "finite_valid_probabilities"}
    numeric = {"all_three_C_seeds_completed": lambda v: v == 3,
               "selection_positives": lambda v: v >= 100, "selection_controls": lambda v: v >= 200,
               "selection_components": lambda v: v >= 30, "bootstrap_valid_replicates": lambda v: v == 2000,
               "bootstrap_invalid_fraction": lambda v: 0 <= v < .1,
               "lower95_AUROC": lambda v: v > .5, "lower95_AP_gain": lambda v: v > 0,
               "lower95_BrierSkill": lambda v: v > 0, "seed_AP_sample_SD": lambda v: 0 <= v <= .03,
               "seed_AP_range": lambda v: 0 <= v <= .10}
    expected = {(model, name) for model in MODELS for name in simple | set(numeric)}
    seen = set()
    for row in gates:
        key = row["model"], row["gate"]
        if key not in expected or key in seen or row["status"] != "PASS":
            fail("Absolute-adequacy table has missing, duplicate, failed or unexpected gate")
        seen.add(key)
        if row["gate"] in simple:
            if row["observed"] != "PASS":
                fail("Absolute-adequacy prerequisite observed value not PASS")
        else:
            value = float(row["observed"])
            if not math.isfinite(value) or not numeric[row["gate"]](value):
                fail("Absolute-adequacy observed numeric value fails frozen condition")
    if seen != expected:
        fail("Absolute-adequacy table lacks required conditions")
    for model in MODELS:
        if decision["model_contexts"][model].get("status") != "PASS" or decision["model_contexts"][model].get("failed_or_inconclusive_gates") != []:
            fail("C decision contains failed/inconclusive conditions")
        b = bootstrap[model]
        if (b["status"] != "PASS" or b["seed"] != 314159 or b["valid"] != 2000 or b["required_valid"] != 2000
                or b["max_attempts"] != 20000 or b["attempted"] != b["valid"] + b["invalid"] or b["attempted"] > 20000
                or not 0 <= b["invalid"] / b["attempted"] < .1 or b["invalid_fraction"] != b["invalid"] / b["attempted"]):
            fail("Underlying bootstrap audit fails frozen adequacy contract")


def prepare_evaluation(stage, repo):
    output = stage / "provenance/evaluation_input.json"
    if output.exists():
        fail("Refusing to replace existing evaluation gate")
    manifest = build_evaluation_manifest(stage, repo)
    manifest["created_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(output, manifest)
    return {"status": "PASS", "phase": "evaluation", "output": absolute_record(output), "evaluation_executed": False}


def prepare_calibration(stage, repo):
    release_path = stage / "provenance/C_release_precalibration_freeze.json"
    calibration_path = stage / "provenance/calibration_input.json"
    if release_path.exists() or calibration_path.exists():
        fail("Refusing to replace existing C release/calibration gate")
    evaluation_input_path = stage / "provenance/evaluation_input.json"
    evaluation_input = read_json(evaluation_input_path)
    refreshed = build_evaluation_manifest(stage, repo)
    for key in ("specification", "prerequisites", "models", "source_bindings"):
        if refreshed[key] != evaluation_input[key]:
            fail(f"Original evaluation source binding drift: {key}")
    directory = stage / "results/chr7_evaluation"
    evaluation_manifest_path = directory / "evaluation_manifest.json"
    evaluation_manifest = read_json(evaluation_manifest_path)
    if evaluation_manifest["input_manifest_sha256"] != sha256(evaluation_input_path):
        fail("Evaluation output is not bound to original evaluation input")
    if evaluation_manifest["evaluation_script_sha256"] != sha256(stage / "scripts/evaluate_chr7.py"):
        fail("Actual evaluation implementation hash changed")
    artifacts = record_map(stage, evaluation_manifest["artifacts"])
    for record in evaluation_manifest["artifacts"]:
        verified_record(stage, record)
    decision_ref = required_record(stage, artifacts, directory / "adequacy_decision.json")
    decision = read_json(decision_ref["path"])
    if decision["input_manifest_sha256"] != sha256(evaluation_input_path):
        fail("Adequacy decision not bound to original input")
    gate_ref = required_record(stage, artifacts, directory / "C_absolute_adequacy.tsv")
    bootstrap_ref = required_record(stage, artifacts, directory / "bootstrap_audit.json")
    validate_passed_adequacy(decision, list(read_rows(gate_ref["path"])), read_json(bootstrap_ref["path"]))
    inference_files = evaluation_input["source_bindings"]["actual_inference_implementation_files"]
    predict_ref = next(record for record in inference_files if Path(record["path"]).name == "predict_chr7.py")
    supporting = [record for record in inference_files if Path(record["path"]).name != "predict_chr7.py"]
    release = {"status": "C_READY_FOR_CALIBRATION_FROZEN", "created_utc": datetime.now(timezone.utc).isoformat(),
        "adequacy_decision": decision_ref, "evaluation_input": absolute_record(evaluation_input_path), "evaluation_manifest": absolute_record(evaluation_manifest_path),
        "checkpoints": [record for record in evaluation_input["prerequisites"]["checkpoints"] if record["configuration"] == "V2-C"],
        "ensemble": {"seeds": list(SEEDS), "reduction_dtype": "float64", "seed_aggregation": "equal-weight", "symmetric_inference": "(p_forward+p_reverse_complement)/2"},
        "symmetric_inference_implementation": predict_ref, "symmetric_inference_supporting_files": supporting,
        "adequacy_source_bindings": [gate_ref, bootstrap_ref], "calibration_performed": False, "test_external_access_authorized": False}
    write_json(release_path, release)
    calibration = {"version": "internal-training-1.0", "created_utc": datetime.now(timezone.utc).isoformat(),
        "specification": evaluation_input["specification"], "adequacy_decision": decision_ref,
        "evaluation_input": absolute_record(evaluation_input_path), "release_freeze": absolute_record(release_path),
        "models": {model: {"calibration": evaluation_input["source_bindings"]["prepared_common_files"][f"inputs/common/{model}_calibration.tsv.gz"], "predictions": evaluation_input["models"][model]["predictions"]} for model in MODELS},
        "planned_output_directory": str(stage / "results/chr7_calibration"), "no_calibration_or_model_execution_by_gate_helper": True}
    write_json(calibration_path, calibration)
    return {"status": "PASS", "phase": "calibration", "outputs": [absolute_record(release_path), absolute_record(calibration_path)], "calibration_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=Path)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--phase", required=True, choices=("evaluation", "calibration"))
    args = parser.parse_args()
    stage = args.stage.resolve()
    repo = args.repo.resolve() if args.repo else stage.parents[3]
    if not stage.is_relative_to(repo):
        fail("Stage must be within repository")
    try:
        result = prepare_evaluation(stage, repo) if args.phase == "evaluation" else prepare_calibration(stage, repo)
    except Exception as error:
        print(f"HARD_STOP: {error}", file=sys.stderr)
        return 2
    import json
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
