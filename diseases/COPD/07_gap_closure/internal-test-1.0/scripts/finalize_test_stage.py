#!/usr/bin/env python3
"""Administrative provenance and exclusive final freeze; never execute models.

Run provenance after scientific execution. Run freeze DIRECTLY, not through the
logged launcher, after all independent validation, report review and register
updates. Freeze writes its own finalization record before enumerating payloads;
artifact_checksums.tsv and freeze.json are the last two writes.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

sys.dont_write_bytecode = True
MODES = ("inputs", "cache", "predictions", "evaluation", "preservation")
REGISTERS = {"diseases/COPD/07_gap_closure/activity_log.tsv",
             "diseases/COPD/07_gap_closure/gap_closure_decision_register.tsv",
             "diseases/COPD/07_gap_closure/gap_closure_result_register.tsv"}
RESULTS = "results/chr8_9_evaluation"
WARNING = re.compile(r"\b(?:warning|error|exception|traceback|retry|retried|failed)\b|(?:^|\s)[WE]\d{4}\s|(?:^|\s)[WE]\s", re.IGNORECASE)


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path):
    with Path(path).open() as stream:
        return json.load(stream)


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def tsv(path, rows, fields):
    with Path(path).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def local_path(stage, value):
    path = (stage / value).resolve()
    if not path.is_relative_to(stage):
        raise RuntimeError(f"Non-stage payload reference rejected: {value}")
    return path


def record(stage, path):
    path = local_path(stage, path)
    return {"path": str(path.relative_to(stage)), "bytes": path.stat().st_size, "sha256": sha(path)}


def verify(path, expected):
    path = Path(path)
    if "bytes" in expected and path.stat().st_size != int(expected["bytes"]):
        raise RuntimeError(f"Hash-bound size mismatch: {path}")
    if sha(path) != expected["sha256"]:
        raise RuntimeError(f"Hash-bound SHA-256 mismatch: {path}")
    return path


def stage_record(stage, expected):
    return verify(local_path(stage, expected["path"]), expected)


def collect_commands(stage, require_complete=False):
    rows, issues, warnings = [], [], []
    command_dir = stage / "provenance/commands"
    started_files = sorted(command_dir.glob("*.started.json"))
    labels = {path.name[:-len(".started.json")] for path in started_files}
    if any(path.name[:-len(".completed.json")] not in labels for path in command_dir.glob("*.completed.json")):
        raise RuntimeError("Orphan completed-command record")
    if any(path.stem not in labels for path in command_dir.glob("*.log")):
        raise RuntimeError("Orphan command log")
    for started_path in started_files:
        label = started_path.name[:-len(".started.json")]
        started = read(started_path)
        completed_path = started_path.with_name(label + ".completed.json")
        log_path = started_path.with_name(label + ".log")
        done = read(completed_path) if completed_path.exists() else None
        if started["label"] != label or (done and (done["argv"] != started["argv"] or done["started_utc"] != started["started_utc"])):
            raise RuntimeError("Command ledger identity changed")
        status = "RUNNING_AT_PROVENANCE_SNAPSHOT" if done is None else ("COMPLETED" if done["returncode"] == 0 else "NONZERO_EXIT")
        if require_complete and done is None:
            raise RuntimeError(f"Cannot freeze with unfinished logged command: {label}; invoke freeze directly")
        code = done["returncode"] if done else None
        if done is None or code:
            issues.append({"label": label, "status": status, "returncode": code,
                           "interpretation": "Nonzero evaluation exit 3 can represent a prospective C gate stop; do not silently classify it as an infrastructure failure."})
        rows.append({"label": label, "argv_json": json.dumps(started["argv"], separators=(",", ":")),
            "cwd": started["cwd"], "started_utc": started["started_utc"],
            "completed_utc": done["completed_utc"] if done else "", "returncode": code if done else "",
            "elapsed_seconds": done["elapsed_seconds"] if done else "", "status": status,
            "SLURM_JOB_ID": started.get("environment", {}).get("SLURM_JOB_ID") or "",
            "SLURM_ARRAY_TASK_ID": started.get("environment", {}).get("SLURM_ARRAY_TASK_ID") or "",
            "CUDA_VISIBLE_DEVICES": started.get("environment", {}).get("CUDA_VISIBLE_DEVICES") or "",
            "PYTHONHASHSEED": started.get("environment", {}).get("PYTHONHASHSEED") or "",
            "started_path": str(started_path.relative_to(stage)), "started_sha256": sha(started_path),
            "completed_path": str(completed_path.relative_to(stage)) if done else "",
            "completed_sha256": sha(completed_path) if done else "",
            "log_path": str(log_path.relative_to(stage)) if log_path.exists() else "",
            "log_sha256": sha(log_path) if log_path.exists() else ""})
        if log_path.exists():
            with log_path.open(errors="replace") as stream:
                for number, line in enumerate(stream, 1):
                    if WARNING.search(line):
                        warnings.append({"label": label, "log_path": str(log_path.relative_to(stage)),
                                         "line_number": number, "line": line.rstrip("\n")})
    failure_records = []
    for path in sorted(stage.rglob("failure.json")):
        if "__pycache__" not in path.parts:
            failure_records.append({"file": record(stage, path), "record": read(path)})
    versions_path = stage / "provenance/independent_validator_code_versions.json"
    versions = {"file": record(stage, versions_path), "record": read(versions_path)} if versions_path.exists() else None
    model_attempts = []
    for path in sorted(stage.rglob("started.json")):
        relative = path.relative_to(stage)
        if relative.parts[0].startswith(("cache", "predictions")):
            model_attempts.append({"file": record(stage, path), "record": read(path)})
    ledger = {"created_utc": utc(), "complete_command_snapshot": require_complete,
              "command_count": len(rows), "nonzero_command_exits": sum(r["returncode"] not in (0, "") for r in rows),
              "unfinished_commands": sum(r["returncode"] == "" for r in rows),
              "command_issues": issues, "failure_records": failure_records,
              "model_execution_attempt_records": model_attempts,
              "n_recorded_phase_I_attempts": sum(Path(r["file"]["path"]).parts[0].startswith("cache") for r in model_attempts),
              "n_recorded_checkpoint_inference_attempts": sum(Path(r["file"]["path"]).parts[0].startswith("predictions") for r in model_attempts),
              "warning_or_error_keyword_lines": warnings, "warning_line_count": len(warnings),
              "warning_interpretation": "Exact diagnostic lines retained. A warning/error keyword, including error:null or nonfatal runtime diagnostics, is not itself proof of a failed command. Full logs and exit statuses are authoritative.",
              "validator_version_history": versions,
              "retry_inference": "Every distinct launched command/attempt is preserved; no failed attempt or code-version change is removed. No absence-of-retry claim is inferred solely from warning text."}
    return rows, ledger


def check_source_chain(stage, repo, require_results=True):
    """Rehash only explicitly approved source paths, never an external benchmark."""
    prospective = read(stage / "provenance/prospective_specification_freeze.json")
    if prospective["status"] != "PASS" or prospective["before_any_test_model_inference"] is not True:
        raise RuntimeError("Prospective freeze is not PASS")
    for key in ("specification", "markdown", "authorization"):
        stage_record(stage, prospective[key])
    approved = {(repo / r["path"]).resolve(): r for r in prospective["sources"]}
    verified_sources = []
    for path, expected in approved.items():
        verify(path, expected)
        verified_sources.append({**expected, "resolved_managed_storage_path": str(path)})
    checkpoints = read(stage / "provenance/checkpoint_verification.json")
    if checkpoints["status"] != "PASS" or checkpoints["n_selected_original_checkpoints"] != 18:
        raise RuntimeError("Original checkpoint gate is not PASS")
    stage_record(stage, checkpoints["prospective_specification_freeze"])
    keys = set()
    for expected in checkpoints["checkpoints"]:
        key = (expected["configuration"], expected["model"], int(expected["seed"]))
        path = (repo / expected["path"]).resolve()
        exact = stage.parent / f"internal-training-1.0/runs/{key[0]}_{key[1]}_seed{key[2]}/attempt-001/selected_checkpoint.keras"
        if path != exact or key in keys or expected["status"] != "PASS" or expected["preserved_original_archive"] is not True:
            raise RuntimeError("Original checkpoint identity/hash-chain mismatch")
        keys.add(key)
        verify(path, expected)
    if keys != {(c, m, s) for c in ("V2-A", "V2-B", "V2-C") for m in ("enhancer", "h3k27me3") for s in (104729, 130363, 155921)}:
        raise RuntimeError("Selected checkpoint population is incomplete")
    consumed = {}
    for name in ("phase_I", "inference", "evaluation"):
        path = stage / f"provenance/{name}_execution_gate.json"
        gate = read(path)
        if gate["status"] != "PASS":
            raise RuntimeError("Execution gate is not PASS")
        for expected in gate["files"]:
            checked = stage_record(stage, expected)
            consumed[str(checked.relative_to(stage))] = dict(expected)
    for relative in ("inputs/input_manifest.json", "cache/cache_manifest.json", "predictions/real_network_invariance.json"):
        manifest = read(stage / relative)
        if manifest["status"] != "PASS":
            raise RuntimeError(f"Scientific manifest is not PASS: {relative}")
        for expected in manifest["files"]:
            checked = stage_record(stage, expected)
            consumed[str(checked.relative_to(stage))] = dict(expected)
    if require_results:
        manifest = read(stage / RESULTS / "evaluation_manifest.json")
        for expected in manifest["artifacts"]:
            checked = stage_record(stage, expected)
            consumed[str(checked.relative_to(stage))] = record(stage, checked)
        verify(stage / "scripts/evaluate_test.py", {"sha256": manifest["evaluation_script_sha256"]})
        stage_record(stage, {"path": manifest["input_manifest"], "sha256": manifest["input_manifest_sha256"]})
    return {"status": "PASS", "verified_utc": utc(), "approved_sources": verified_sources,
            "original_checkpoints": checkpoints["checkpoints"], "n_original_checkpoints_unchanged": len(keys),
            "consumed_stage_files": [consumed[key] for key in sorted(consumed)],
            "no_external_benchmark_or_candidate_file_opened": True,
            "scientific_extraction_inference_evaluation_code_unchanged": True}


def resource_summaries(stage, command_rows):
    sources = ["inputs/input_manifest.json", "cache/cache_manifest.json", "predictions/completed.json"]
    sources += [f"provenance/{mode}_independent_validation.json" for mode in MODES]
    compute = []
    fields = ("elapsed_seconds", "wall_seconds", "user_cpu_seconds", "system_cpu_seconds", "peak_RSS_KiB", "started_utc", "completed_utc", "finished_utc")
    for relative in sources:
        path = stage / relative
        if not path.exists():
            continue
        source = read(path)
        values = {key: source[key] for key in fields if key in source}
        values.update(source.get("resources", {}))
        compute.append({"source": record(stage, path), "recorded_resources": values})
    environments = []
    for relative in ("cache/environment.json", "predictions/environment.json"):
        path = stage / relative
        environments.append({"source": record(stage, path), "record": read(path)})
    input_manifest = read(stage / "inputs/input_manifest.json")
    return {"created_utc": utc(), "compute_records": compute, "environment_records": environments,
            "input_preparation_environment": input_manifest.get("environment"),
            "command_wall_seconds_completed_sum": sum(float(r["elapsed_seconds"]) for r in command_rows if r["elapsed_seconds"] != ""),
            "wall_sum_interpretation": "Sum of logged subprocess durations, not elapsed calendar time or exclusively GPU time; validation and preparation can overlap. CPU/RSS fields are recorded per process and are not summed as independent incremental resources.",
            "GPU_use": "Recorded phase-I extraction and frozen-checkpoint inference only; no training",
            "scientific_analysis_performed_by_this_script": False}


def provenance(stage, repo):
    chain = check_source_chain(stage, repo)
    rows, failures = collect_commands(stage)
    if not rows:
        raise RuntimeError("Command ledger empty")
    write(stage / "provenance/execution_source_manifest.json", chain)
    tsv(stage / "provenance/command_job_ledger_provenance.tsv", rows, list(rows[0]))
    write(stage / "provenance/retry_failure_ledger.json", failures)
    write(stage / "provenance/compute_environment_summary.json", resource_summaries(stage, rows))
    return {"status": "PASS", "phase": "provenance", "commands_in_snapshot": len(rows),
            "note": "Snapshot can include its own unfinished logged wrapper; complete final command ledger is constructed during direct freeze."}


def check_validation_and_registers(stage, repo):
    checks = []
    lineage = read(stage / "provenance/independent_validator_code_versions.json")
    for mode in MODES:
        path = stage / f"provenance/{mode}_independent_validation.json"
        audit = read(path)
        if audit["status"] != "PASS" or audit["failures"] != 0 or audit["mode"] != mode:
            raise RuntimeError(f"Independent {mode} QC not PASS")
        for field in ("external_benchmark_opened", "candidate_payloads_opened", "historical_sources_written", "model_execution_performed"):
            if audit[field] is not False:
                raise RuntimeError("Independent validation scope violation")
        script = audit["script"]
        current = local_path(stage, script["path"])
        if sha(current) != script["sha256"]:
            original = lineage["initial_script"]
            if mode != "inputs" or original["sha256"] != script["sha256"]:
                raise RuntimeError("Independent validator version is unaccounted for")
            stage_record(stage, original)
        checks.append({"check": f"independent_{mode}_QC", "status": "PASS", "detail": f"{audit['checks']} checks; zero failures", "evidence": record(stage, path)})
    report_review_path = stage / "provenance/report_independent_review.json"
    review = read(report_review_path)
    if review["status"] != "PASS" or review["report"]["path"] != "INTERNAL_TEST_REPORT.md":
        raise RuntimeError("Independent final report review not PASS")
    stage_record(stage, review["report"])
    checks.append({"check": "independent_report_review", "status": "PASS", "detail": "Report hash unchanged after independent review", "evidence": record(stage, report_review_path)})
    register_path = stage / "provenance/shared_register_updates.json"
    registers = read(register_path)
    if registers["status"] != "PASS" or {r["path"] for r in registers["registers"]} != REGISTERS or len(registers["registers"]) != 3:
        raise RuntimeError("Three explicit append-only register records required")
    for expected in registers["registers"]:
        if expected["append_only"] is not True:
            raise RuntimeError("Register append-only gate failed")
        path = repo / expected["path"]
        verify(path, expected)
        content = path.read_bytes()
        before = int(expected["before_bytes"])
        if hashlib.sha256(content[:before]).hexdigest() != expected["before_sha256"]:
            raise RuntimeError("Historical register prefix changed")
        appended = content[before:]
        if len(appended) != int(expected["appended_bytes"]) or len(appended.splitlines()) != int(expected["appended_rows"]):
            raise RuntimeError("Register append accounting mismatch")
    checks.append({"check": "three_append_only_registers", "status": "PASS", "detail": "Full current hashes and historical byte prefixes verified", "evidence": record(stage, register_path)})
    return checks


def payload_files(stage):
    omitted = {"provenance/artifact_checksums.tsv", "provenance/freeze.json"}
    files = []
    for path in sorted(stage.rglob("*")):
        if "__pycache__" in path.relative_to(stage).parts:
            continue
        if path.is_symlink():
            raise RuntimeError("Frozen payload cannot contain an unaccounted symlink")
        if path.is_file() and str(path.relative_to(stage)) not in omitted:
            files.append(path)
    return files


def freeze(stage, repo):
    # Exclusive destinations checked before any finalization-side writes.
    destinations = ["command_job_ledger.tsv", "retry_failure_ledger_final.json", "compute_environment_summary_final.json",
                    "source_chain_final_validation.json", "final_validation.json", "final_validation.tsv",
                    "finalization_record.json", "artifact_checksums.tsv", "freeze.json"]
    if any((stage / "provenance" / name).exists() for name in destinations):
        raise RuntimeError("Finalization output already exists; refuse overwrite/re-freeze")
    rows, failures = collect_commands(stage, require_complete=True)
    chain = check_source_chain(stage, repo)
    checks = check_validation_and_registers(stage, repo)
    if not rows:
        raise RuntimeError("No completed command ledger")
    head = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    prospective = read(stage / "provenance/prospective_specification_freeze.json")
    if head != prospective["baseline_commit"]:
        raise RuntimeError("Repository commit changed during no-publication holdout stage")
    decision = read(stage / RESULTS / "holdout_decision.json")
    if decision["version"] != "internal-test-1.0" or decision["external_access_authorized"] is not False or decision["candidate_scoring_authorized"] is not False:
        raise RuntimeError("Final decision scope mismatch")
    tsv(stage / "provenance/command_job_ledger.tsv", rows, list(rows[0]))
    write(stage / "provenance/retry_failure_ledger_final.json", failures)
    write(stage / "provenance/compute_environment_summary_final.json", resource_summaries(stage, rows))
    write(stage / "provenance/source_chain_final_validation.json", chain)
    write(stage / "provenance/final_validation.json", {"status": "PASS", "created_utc": utc(), "checks": checks,
        "independent_validation_modes": list(MODES), "source_chain": record(stage, stage / "provenance/source_chain_final_validation.json"),
        "scientific_gate_status_not_overridden": decision["model_contexts"],
        "validation_is_QC_not_new_scientific_evaluation": True})
    tsv(stage / "provenance/final_validation.tsv", [{k: row[k] for k in ("check", "status", "detail")} for row in checks], ["check", "status", "detail"])
    write(stage / "provenance/finalization_record.json", {"status": "ALL_PRESEAL_CHECKS_PASS_READY_FOR_PAYLOAD_HASHING",
        "created_utc": utc(), "argv": sys.argv, "cwd": os.getcwd(), "script": record(stage, Path(__file__)),
        "python": platform.python_version(), "executable": sys.executable, "git_HEAD_unchanged": head,
        "direct_invocation_required": True, "all_logged_commands_completed": True,
        "final_writes": ["provenance/artifact_checksums.tsv", "provenance/freeze.json"],
        "after_freeze_no_bundle_writes_permitted": True,
        "no_models_evaluations_retraining_external_benchmark_or_candidates_executed_by_finalizer": True})
    files = payload_files(stage)
    payloads = [record(stage, path) for path in files]
    ledger = stage / "provenance/artifact_checksums.tsv"
    tsv(ledger, payloads, ["path", "bytes", "sha256"])
    frozen = {"status": "FROZEN_INTERNAL_HOLDOUT_COMPLETE", "version": "internal-test-1.0", "created_utc": utc(),
        "artifact_checksums": record(stage, ledger), "n_payloads": len(payloads),
        "total_payload_bytes": sum(r["bytes"] for r in payloads),
        "checksum_scope": "All regular stage payloads except __pycache__ and the checksum ledger/freeze self-references; original checkpoints are external immutable sources, not copied into this stage.",
        "prospective_specification_freeze": record(stage, stage / "provenance/prospective_specification_freeze.json"),
        "independent_final_validation": record(stage, stage / "provenance/final_validation.json"),
        "report": record(stage, stage / "INTERNAL_TEST_REPORT.md"),
        "model_contexts": decision["model_contexts"], "both_C_contexts_pass": decision["both_C_contexts_pass"],
        "next_action": decision["next_action"], "git_HEAD_unchanged": head,
        "commit_or_push_performed": False, "external_benchmark_opened": False,
        "candidate_scoring_performed": False, "no_further_analysis_authorized": True}
    write(stage / "provenance/freeze.json", frozen)
    # Absolutely no additional filesystem writes after the seal.
    return {"status": frozen["status"], "n_payloads": len(payloads), "total_payload_bytes": frozen["total_payload_bytes"],
            "freeze_sha256": sha(stage / "provenance/freeze.json"), "artifact_checksums_sha256": frozen["artifact_checksums"]["sha256"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--phase", choices=("provenance", "freeze"), required=True)
    args = parser.parse_args()
    stage, repo = args.stage.resolve(), args.repo.resolve()
    if stage.name != "internal-test-1.0" or not stage.is_relative_to(repo):
        raise RuntimeError("Finalizer is restricted to the authorized new stage")
    if (stage / "provenance/freeze.json").exists():
        raise RuntimeError("Completed stage is immutable")
    result = provenance(stage, repo) if args.phase == "provenance" else freeze(stage, repo)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
