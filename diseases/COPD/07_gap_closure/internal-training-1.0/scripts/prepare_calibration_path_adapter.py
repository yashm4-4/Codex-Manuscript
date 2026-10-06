#!/usr/bin/env python3
"""Append-only metadata-path adapter for the already completed chr7 evaluation.

The frozen scientific gate helpers are imported unchanged. The sole repaired
behavior is resolving producer artifact references relative to the repository
(also accepting equivalent absolute/stage-relative spellings), not blindly
relative to the internal-stage directory. No producer manifest is rewritten;
all bytes, hashes, scientific conditions and frozen ensemble rules are retained.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

from prepare_evaluation_gates import (
    MODELS, SEEDS, absolute_record, build_evaluation_manifest, fail, read_json,
    read_rows, record_map, required_record, sha256, validate_passed_adequacy,
    verified_record, write_json,
)


def normalize_evaluation_artifacts(stage, repo, records):
    """Verify original bytes/SHA and return in-memory absolute references only."""
    stage, repo = Path(stage).resolve(), Path(repo).resolve()
    directory = (stage / "results/chr7_evaluation").resolve()
    if not stage.is_relative_to(repo) or not directory.is_relative_to(stage):
        fail("Evaluation directory/stage is outside the permitted repository")
    normalized, mappings, seen = [], [], set()
    for original in records:
        if not isinstance(original.get("path"), str) or not original["path"]:
            fail("Empty/nonstring evaluation artifact reference")
        raw = Path(original["path"])
        if ".." in raw.parts:
            fail("Traversal is not permitted in evaluation artifact metadata")
        if "sha256" not in original or "bytes" not in original:
            fail("Original evaluation artifact must bind both SHA-256 and byte count")
        candidates = [raw.resolve()] if raw.is_absolute() else [(stage / raw).resolve(), (repo / raw).resolve()]
        permitted = {candidate for candidate in candidates if candidate.parent == directory}
        if len(permitted) != 1:
            fail(f"Evaluation artifact does not resolve uniquely into exact evaluation directory: {original['path']}")
        path = permitted.pop()
        if path in seen:
            fail(f"Duplicate evaluation artifact after path normalization: {path}")
        seen.add(path)
        if not path.is_file():
            fail(f"Missing evaluation artifact: {path}")
        if path.stat().st_size != int(original["bytes"]) or sha256(path) != original["sha256"]:
            fail(f"Original evaluation artifact bytes/hash mismatch: {path}")
        current = absolute_record(path)
        normalized.append(current)
        mappings.append({"original": dict(original), "resolved": current})
    if not normalized:
        fail("Evaluation artifact manifest is empty")
    return normalized, mappings


def option_value(argv, option):
    values = []
    for index, value in enumerate(argv):
        if value == option:
            if index + 1 == len(argv):
                fail(f"Command missing value for {option}")
            values.append(argv[index + 1])
        elif value.startswith(option + "="):
            values.append(value[len(option) + 1:])
    if len(values) != 1:
        fail(f"Command must have exactly one {option} value")
    return values[0]


def command_path(cwd, value):
    return (Path(cwd) / value).resolve()


def verify_command_evidence(stage, repo):
    evaluation_path = stage / "provenance/commands/evaluate_common_chr7.completed.json"
    failure_path = stage / "provenance/commands/freeze_C_release_for_calibration.completed.json"
    evaluation, failed = read_json(evaluation_path), read_json(failure_path)
    for command in (evaluation, failed):
        if Path(command["cwd"]).resolve() != repo:
            fail("Recorded command cwd differs from the repository root")
    if evaluation.get("label") != "evaluate_common_chr7" or evaluation.get("returncode") != 0:
        fail("Original common-chr7 evaluation command did not succeed")
    if command_path(repo, option_value(evaluation["argv"], "--input-manifest")) != stage / "provenance/evaluation_input.json":
        fail("Original evaluation command used another input manifest")
    if command_path(repo, option_value(evaluation["argv"], "--output")) != stage / "results/chr7_evaluation":
        fail("Original evaluation command used another output directory")
    if str(stage / "scripts/evaluate_chr7.py") not in [str(command_path(repo, item)) for item in evaluation["argv"] if not item.startswith("-")]:
        fail("Recorded successful command did not invoke the frozen evaluation script")
    if failed.get("label") != "freeze_C_release_for_calibration" or failed.get("returncode") != 2:
        fail("Expected original pre-output calibration-gate failure (exit 2) is not recorded")
    if option_value(failed["argv"], "--phase") != "calibration" or command_path(repo, option_value(failed["argv"], "--stage")) != stage:
        fail("Failed command is not the original calibration-gate invocation")
    if str(stage / "scripts/prepare_evaluation_gates.py") not in [str(command_path(repo, item)) for item in failed["argv"] if not item.startswith("-")]:
        fail("Failed command did not invoke the frozen gate helper")
    return {"successful_evaluation_command": absolute_record(evaluation_path),
            "original_failed_gate_command": absolute_record(failure_path),
            "original_failed_gate_label": failed["label"], "original_failed_gate_returncode": failed["returncode"],
            "successful_evaluation_cwd": str(repo)}


def verify_passed_evaluation_artifacts(stage, repo, evaluation_manifest, evaluation_input_path):
    """Only reference normalization differs from frozen prepare_calibration."""
    directory = stage / "results/chr7_evaluation"
    normalized, mappings = normalize_evaluation_artifacts(stage, repo, evaluation_manifest["artifacts"])
    artifacts = record_map(stage, normalized)
    for record in normalized:
        verified_record(stage, record)
    decision_ref = required_record(stage, artifacts, directory / "adequacy_decision.json")
    decision = read_json(decision_ref["path"])
    if decision["input_manifest_sha256"] != sha256(evaluation_input_path):
        fail("Adequacy decision not bound to original input")
    gate_ref = required_record(stage, artifacts, directory / "C_absolute_adequacy.tsv")
    bootstrap_ref = required_record(stage, artifacts, directory / "bootstrap_audit.json")
    # The exact original gate validator is called without replacement or patching.
    validate_passed_adequacy(decision, list(read_rows(gate_ref["path"])), read_json(bootstrap_ref["path"]))
    return decision_ref, gate_ref, bootstrap_ref, mappings


def serialized_record(path, obj):
    # Matches the unchanged write_json serializer imported from the frozen helper.
    payload = (json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}


def prepare(stage):
    stage = Path(stage).resolve()
    repo = stage.parents[3]
    if not stage.is_relative_to(repo):
        fail("Stage must remain within repository")
    release_path = stage / "provenance/C_release_precalibration_freeze.json"
    calibration_path = stage / "provenance/calibration_input.json"
    note_path = stage / "provenance/calibration_gate_path_resolution.json"
    if release_path.exists() or calibration_path.exists():
        fail("Refusing to replace existing C release/calibration gate")
    if note_path.exists():
        fail("Refusing to replace an existing metadata-path adapter note")
    helper_path = stage / "scripts/prepare_evaluation_gates.py"
    evaluation_input_path = stage / "provenance/evaluation_input.json"
    directory = stage / "results/chr7_evaluation"
    evaluation_manifest_path = directory / "evaluation_manifest.json"
    unchanged_before = {"frozen_gate_helper": absolute_record(helper_path),
                        "evaluation_manifest": absolute_record(evaluation_manifest_path),
                        "evaluation_input": absolute_record(evaluation_input_path)}
    command_evidence = verify_command_evidence(stage, repo)
    evaluation_input = read_json(evaluation_input_path)
    # Preserve every original construction/checkpoint/inference/prerequisite check.
    refreshed = build_evaluation_manifest(stage, repo)
    for key in ("specification", "prerequisites", "models", "source_bindings"):
        if refreshed[key] != evaluation_input[key]:
            fail(f"Original evaluation source binding drift: {key}")
    evaluation_manifest = read_json(evaluation_manifest_path)
    if evaluation_manifest["input_manifest_sha256"] != sha256(evaluation_input_path):
        fail("Evaluation output is not bound to original evaluation input")
    if evaluation_manifest["evaluation_script_sha256"] != sha256(stage / "scripts/evaluate_chr7.py"):
        fail("Actual evaluation implementation hash changed")
    decision_ref, gate_ref, bootstrap_ref, mappings = verify_passed_evaluation_artifacts(
        stage, repo, evaluation_manifest, evaluation_input_path)
    inference_files = evaluation_input["source_bindings"]["actual_inference_implementation_files"]
    predict_ref = next(record for record in inference_files if Path(record["path"]).name == "predict_chr7.py")
    supporting = [record for record in inference_files if Path(record["path"]).name != "predict_chr7.py"]
    c_checkpoints = [record for record in evaluation_input["prerequisites"]["checkpoints"] if record["configuration"] == "V2-C"]
    if (len(c_checkpoints) != 6 or {(r["model"], int(r["seed"])) for r in c_checkpoints} != {(m, s) for m in MODELS for s in SEEDS}
            or any(r.get("completed") is not True for r in c_checkpoints)):
        fail("The exact six completed C checkpoints are required")
    unchanged_after = {"frozen_gate_helper": absolute_record(helper_path),
                       "evaluation_manifest": absolute_record(evaluation_manifest_path),
                       "evaluation_input": absolute_record(evaluation_input_path)}
    if unchanged_before != unchanged_after:
        fail("Frozen helper or original producer metadata changed during verification")
    created = datetime.now(timezone.utc).isoformat()
    note = {"status": "PASS_METADATA_PATH_ADAPTER", "created_utc": created,
            "adapter": absolute_record(Path(__file__)), "original_frozen_sources": unchanged_before,
            "original_frozen_sources_unchanged_after_verification": True,
            **command_evidence, "n_gate_outputs_before_repair": 0,
            "artifact_reference_mappings": mappings,
            "scope": "Normalize artifact-reference namespaces only; accept equivalent absolute/stage-relative/repo-relative references into the exact existing chr7_evaluation directory, verify original bytes/SHA, reject traversal/escape/duplicates/missing files.",
            "frozen_functions_reused_unchanged": ["build_evaluation_manifest", "validate_passed_adequacy", "record_map", "required_record", "verified_record", "write_json"],
            "cause": "Original evaluation producer recorded repository-relative artifact paths, but frozen calibration-gate helper resolved those paths against the stage directory; the original invocation stopped with exit 2 before creating either gate.",
            "scientific_payload_changed": False, "evaluation_or_model_or_calibration_rerun": False,
            "evaluation_manifest_rewritten": False, "scientific_gate_conditions_changed": False,
            "checkpoints_or_predictions_changed": False, "calibration_performed": False,
            "test_external_access_authorized": False}
    note_ref = serialized_record(note_path, note)
    release = {"status": "C_READY_FOR_CALIBRATION_FROZEN", "created_utc": created,
        "adequacy_decision": decision_ref, "evaluation_input": absolute_record(evaluation_input_path), "evaluation_manifest": absolute_record(evaluation_manifest_path),
        "checkpoints": c_checkpoints,
        "ensemble": {"seeds": list(SEEDS), "reduction_dtype": "float64", "seed_aggregation": "equal-weight", "symmetric_inference": "(p_forward+p_reverse_complement)/2"},
        "symmetric_inference_implementation": predict_ref, "symmetric_inference_supporting_files": supporting,
        "adequacy_source_bindings": [gate_ref, bootstrap_ref], "calibration_performed": False, "test_external_access_authorized": False,
        "metadata_path_resolution": note_ref}
    release_ref = serialized_record(release_path, release)
    calibration = {"version": "internal-training-1.0", "created_utc": created,
        "specification": evaluation_input["specification"], "adequacy_decision": decision_ref,
        "evaluation_input": absolute_record(evaluation_input_path), "release_freeze": release_ref,
        "models": {model: {"calibration": evaluation_input["source_bindings"]["prepared_common_files"][f"inputs/common/{model}_calibration.tsv.gz"], "predictions": evaluation_input["models"][model]["predictions"]} for model in MODELS},
        "planned_output_directory": str(stage / "results/chr7_calibration"), "no_calibration_or_model_execution_by_gate_helper": True,
        "metadata_path_resolution": note_ref}
    # All checks and payload construction precede these exclusive writes. Neither
    # the original helper nor either producer manifest is changed in place.
    write_json(note_path, note)
    write_json(release_path, release)
    write_json(calibration_path, calibration)
    if absolute_record(note_path) != note_ref or absolute_record(release_path) != release_ref:
        fail("Unexpected JSON serialization binding mismatch")
    return {"status": "PASS", "phase": "calibration", "metadata_path_adapter": absolute_record(note_path),
            "outputs": [absolute_record(release_path), absolute_record(calibration_path)], "calibration_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = prepare(args.stage)
    except Exception as error:
        print(f"HARD_STOP: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
