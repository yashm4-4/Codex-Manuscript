#!/usr/bin/env python3
"""Read-only ZIP/HDF5 packaging audit; never initialize a model or optimizer.

Documents optimizer metadata serialized by Keras model_from_json. These saved
archives are inference-only by permitted use, not optimizer-state-free ZIPs.
No checkpoint, training code, scientific outcome or selection rule is changed.
"""
import argparse
import ast
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect(stage, run_id):
    attempts = sorted((stage / "runs" / run_id).glob("attempt-*/completed.json"))
    if len(attempts) != 1:
        raise ValueError("Container audit requires exactly one completed attempt: " + run_id)
    completed_path = attempts[0]
    selected_path = completed_path.with_name("selected_checkpoint.json")
    selected = json.loads(selected_path.read_text())
    record = selected["checkpoint"]
    checkpoint = stage / record["path"]
    original_sha = sha(checkpoint)
    if original_sha != record["sha256"] or checkpoint.stat().st_size != record["bytes"]:
        raise ValueError("Completed checkpoint bytes changed before container audit")
    datasets = []
    with zipfile.ZipFile(checkpoint) as archive:
        configuration = json.loads(archive.read("config.json"))
        compile_config = configuration.get("compile_config") or {}
        member_sizes = {entry.filename: entry.file_size for entry in archive.infolist()}
        with h5py.File(io.BytesIO(archive.read("model.weights.h5")), "r") as weights:
            def collect(name, obj):
                if not isinstance(obj, h5py.Dataset):
                    return
                array = np.asarray(obj[()])
                is_optimizer = name.startswith("optimizer/")
                is_slot = is_optimizer and name not in ("optimizer/vars/0", "optimizer/vars/1")
                datasets.append({"run_id": run_id, "dataset": name,
                                 "kind": "optimizer_slot" if is_slot else ("optimizer_scalar" if is_optimizer else "model_weight"),
                                 "shape": json.dumps(list(array.shape)), "dtype": str(array.dtype),
                                 "elements": int(array.size), "bytes": int(array.nbytes),
                                 "finite": bool(np.isfinite(array).all()),
                                 "slot_all_zero": bool(np.count_nonzero(array) == 0) if is_slot else "",
                                 "sha256": hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()})
            weights.visititems(collect)
            iteration = int(weights["optimizer/vars/0"][()]) if "optimizer/vars/0" in weights else None
            learning_rate = float(weights["optimizer/vars/1"][()]) if "optimizer/vars/1" in weights else None
    model_rows = [row for row in datasets if row["kind"] == "model_weight"]
    slot_rows = [row for row in datasets if row["kind"] == "optimizer_slot"]
    optimizer_rows = [row for row in datasets if row["kind"].startswith("optimizer_")]
    checks = {
        "selected_checkpoint_hash_unchanged": sha(checkpoint) == original_sha,
        "completed_selected_weight_bitwise_roundtrip_reported_PASS": selected.get("serialized_selected_weights_bitwise_roundtrip") is True,
        "checkpoint_permitted_use_inference_only": selected.get("checkpoint_is_inference_only") is True,
        "compile_configuration_present": bool(compile_config),
        "serialized_optimizer_Adadelta": compile_config.get("optimizer", {}).get("class_name") == "Adadelta",
        "optimizer_iteration_zero": iteration == 0,
        "twenty_four_optimizer_slots": len(slot_rows) == 24,
        "all_optimizer_slots_zero": bool(slot_rows) and all(row["slot_all_zero"] for row in slot_rows),
        "all_archive_tensors_finite": all(row["finite"] for row in datasets),
        "fourteen_model_tensors": len(model_rows) == 14,
        "model_parameter_count_matches_frozen_architecture": sum(row["elements"] for row in model_rows) == 29167889,
    }
    summary = {"run_id": run_id, "checkpoint": record,
               "selected_checkpoint_record": {"path": str(selected_path.relative_to(stage)), "sha256": sha(selected_path)},
               "archive_bytes": checkpoint.stat().st_size, "archive_member_bytes": member_sizes,
               "model_tensor_count": len(model_rows), "model_tensor_bytes": sum(row["bytes"] for row in model_rows),
               "optimizer_tensor_count": len(optimizer_rows), "optimizer_tensor_bytes": sum(row["bytes"] for row in optimizer_rows),
               "optimizer_slot_count": len(slot_rows), "optimizer_slot_bytes": sum(row["bytes"] for row in slot_rows),
               "optimizer_iterations": iteration, "optimizer_learning_rate": learning_rate,
               "checks": checks, "packaging_checks_pass": all(checks.values())}
    return summary, datasets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--run-id", action="append", help="Explicit completed run; repeat as needed")
    parser.add_argument("--all-completed", action="store_true")
    parser.add_argument("--final", action="store_true", help="Require all18 and emit standard independent packaging QC schema")
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    if bool(args.run_id) == bool(args.all_completed):
        raise ValueError("Choose explicit completed runs OR --all-completed")
    stage = args.stage.resolve()
    runs = args.run_id if args.run_id else sorted({path.parents[1].name for path in (stage / "runs").glob("*/attempt-*/completed.json")})
    expected_runs = {f"{config}_{model}_seed{seed}" for config in ("V2-A","V2-B","V2-C") for model in ("enhancer","h3k27me3") for seed in (104729,130363,155921)}
    if not set(runs) <= expected_runs or len(runs) != len(set(runs)):
        raise ValueError("Unexpected or duplicate prescribed run identity")
    if args.final and not args.all_completed:
        raise ValueError("Final packaging audit requires --all-completed")
    if not runs:
        raise ValueError("No completed runs to audit")
    prefix = args.output_prefix.resolve()
    paths = [prefix.with_suffix(".json"), prefix.with_name(prefix.name + "_datasets.tsv")]
    if any(path.exists() for path in paths):
        raise FileExistsError("Packaging audit is append-only; choose a new prefix")
    summaries, datasets = [], []
    for run in runs:
        summary, rows = inspect(stage, run)
        summaries.append(summary)
        datasets.extend(rows)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with paths[1].open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(datasets[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(datasets)
    result = {"status": "PASS_PACKAGING_CLARIFIED" if all(r["packaging_checks_pass"] for r in summaries) else "PACKAGING_REVIEW_REQUIRED",
              "created_utc": datetime.now(timezone.utc).isoformat(), "audit_scope": "ZIP/HDF5 container metadata and tensor integrity only",
              "model_network_executed": False, "training_or_checkpoints_modified": False,
              "selection_calibration_test_or_benchmark_predictions_read": False,
              "checkpoints_audited": len(summaries), "checkpoints": summaries,
              "dataset_table": {"path": str(paths[1].relative_to(stage)), "sha256": sha(paths[1])},
              "script_sha256": sha(__file__), "h5py": h5py.__version__, "numpy": np.__version__,
              "clarification": "Keras3.14.1 model.to_json includes compile_config; model_from_json restores compilation and initializes optimizer slots. The frozen source comment claiming architecture-only/uncompiled serialization is inaccurate. Existing archives retain the selected network weights and additional freshly initialized optimizer state.",
              "inference_assessment": "The frozen training procedure restores selected trainable/BatchNorm weights, verifies identical native checkpoint BCE, sets those weights on the serialized model, and performs compile=False reload with bitwise network-weight equality. The deployed prediction wrapper also loads compile=False. Extra optimizer tensors therefore do not enter forward inference; mandatory actual-network RC and downstream output QC still apply.",
              "risks": ["Do not resume optimization from these archives: iteration0/zero slots are not selected-epoch optimizer state.", "Do not treat checkpoint_is_inference_only as a claim that the ZIP lacks compile metadata or optimizer datasets.", "Extra optimizer tensors increase storage/checksum/I/O cost; preserve original archives and checksums, do not reserialize during the frozen stage."],
              "action": "Append-only packaging clarification; no retraining, checkpoint replacement, seed replacement, inference-policy change or scientific redesign."}
    if args.final:
        from validate_internal_stage import Audit
        audit = Audit()
        audit.check("packaging:all18_prescribed_containers", set(runs) == expected_runs and len(summaries) == 18, len(summaries), 18)
        for summary in summaries:
            for name, passed in summary["checks"].items():
                audit.check(summary["run_id"] + ":" + name, passed)
        gate = json.loads((stage / "provenance/execution_input_gate.json").read_text())
        gate_files = {row["path"]:row for row in gate["files"]}
        source_evidence = []
        for name in ("train_phase_two.py", "predict_chr7.py"):
            path = stage / "scripts" / name
            relative = str(path.relative_to(stage))
            reference = gate_files[relative]
            audit.check("packaging:immutable_consumer:" + name, sha(path) == reference["sha256"])
            tree = ast.parse(path.read_text())
            loads = [node for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and node.func.attr == "load_model"]
            compile_false = bool(loads) and all(any(keyword.arg == "compile" and isinstance(keyword.value,ast.Constant) and keyword.value.value is False for keyword in node.keywords) for node in loads)
            audit.check("packaging:compileFalse_consumer:" + name, compile_false)
            source_evidence.append({"path":relative,"sha256":sha(path),"load_model_calls":len(loads),"all_load_model_compile_false":compile_false})
        freeze = json.loads((stage / "provenance/checkpoint_freeze.json").read_text())
        frozen = {(row["configuration"],row["model"],int(row["seed"])):row for row in freeze["checkpoints"]}
        audit.check("packaging:checkpoint_freeze_complete", freeze.get("status") == "PASS" and len(frozen) == 18 and freeze.get("all_18_selected_checkpoints_frozen") is True)
        for summary in summaries:
            key = next(key for key in frozen if f"{key[0]}_{key[1]}_seed{key[2]}" == summary["run_id"])
            audit.check("packaging:frozen_archive:" + summary["run_id"], summary["checkpoint"]["sha256"] == frozen[key]["sha256"] and summary["checkpoint"]["path"] == frozen[key]["path"])
        original_status = result.pop("status")
        result.update({"phase":"packaging", "packaging_status":original_status,
                       "source_evidence":source_evidence, "final_all18_packaging_validation":True})
        failures = audit.write(prefix,result)
        return int(bool(failures))
    with paths[0].open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "checkpoints_audited": len(summaries), "output": str(paths[0]),
                      "model_tensor_bytes_each": [r["model_tensor_bytes"] for r in summaries],
                      "optimizer_tensor_bytes_each": [r["optimizer_tensor_bytes"] for r in summaries]}, sort_keys=True))
    return int(result["status"] != "PASS_PACKAGING_CLARIFIED")


if __name__ == "__main__":
    raise SystemExit(main())
