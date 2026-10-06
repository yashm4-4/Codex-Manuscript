#!/usr/bin/env python3
"""Collect read-only Slurm accounting after the frozen all-18 checkpoint gate.

Does not execute a model, recalculate any scientific result, submit jobs, or
alter existing artifacts. Three exclusively created provenance outputs retain
all allocation/job-step rows, allocation-only resource totals and failures.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import time


FIELDS = ("JobID", "State", "ExitCode", "Start", "End", "ElapsedRaw", "TotalCPU", "AllocCPUS", "ReqMem", "MaxRSS", "AllocTRES", "NodeList")
CONFIGURATIONS = ("V2-A", "V2-B", "V2-C")
MODELS = ("enhancer", "h3k27me3")
SEEDS = (104729, 130363, 155921)
FAILURE_FIELDS = ("scope", "run_id", "failed_attempt", "classification", "status", "reason",
                  "failure_record", "failure_record_sha256", "retry_authorization", "retry_authorization_sha256",
                  "retry_authorized", "same_seed_identical_procedure")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def local_path(stage, relative):
    path = (stage / relative).resolve()
    if not path.is_relative_to(stage):
        raise RuntimeError(f"Non-stage provenance reference rejected: {relative}")
    return path


def file_record(stage, path):
    path = Path(path).resolve()
    return {"path": str(path.relative_to(stage)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def verify_record(stage, record):
    path = local_path(stage, record["path"])
    if sha256(path) != record["sha256"] or ("bytes" in record and path.stat().st_size != int(record["bytes"])):
        raise RuntimeError(f"Frozen provenance/checkpoint integrity failure: {record['path']}")
    return path


def tsv_bytes(rows, fields):
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return text.getvalue().encode("utf-8")


def parse_sacct(raw):
    rows = []
    for values in csv.reader(io.StringIO(raw.decode("utf-8")), delimiter="|"):
        if not values or (len(values) == 1 and not values[0].strip()):
            continue
        if len(values) != len(FIELDS):
            raise RuntimeError(f"Unexpected sacct field count {len(values)}; expected {len(FIELDS)}")
        rows.append(dict(zip(FIELDS, (value.strip() for value in values))))
    if not rows:
        raise RuntimeError("sacct returned no completed accounting records")
    return rows


def gpu_count(tres):
    parsed = {}
    for item in tres.split(","):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if not separator:
            raise RuntimeError(f"Malformed AllocTRES field: {item}")
        parsed[key.strip()] = value.strip()
    typed = {key: int(value) for key, value in parsed.items() if key.startswith("gres/gpu:")}
    if "gres/gpu" in parsed:
        count = int(parsed["gres/gpu"])
        if typed and sum(typed.values()) != count:
            raise RuntimeError("Typed and generic GPU allocation counts disagree")
    elif typed:
        count = sum(typed.values())
    else:
        raise RuntimeError("No allocated GPU count found in AllocTRES")
    if count != 1:
        raise RuntimeError(f"Expected frozen one-GPU allocation; observed {count}")
    return count


def allocation_summary(rows, array_id):
    expected = {f"{array_id}_{index}" for index in range(18)}
    allocations, steps, parent_rows = {}, [], []
    all_ids = []
    for row in rows:
        job_id = row["JobID"]
        all_ids.append(job_id)
        if job_id in expected:
            if job_id in allocations:
                raise RuntimeError(f"Duplicate allocation accounting: {job_id}")
            allocations[job_id] = row
        elif "." in job_id and job_id.split(".", 1)[0] in expected:
            steps.append(row)
        elif job_id == array_id:
            parent_rows.append(row)
        else:
            raise RuntimeError(f"Unexpected job outside the requested frozen array: {job_id}")
    if len(set(all_ids)) != len(all_ids):
        raise RuntimeError("Repeated job/step IDs make accounting ambiguous")
    if set(allocations) != expected:
        raise RuntimeError(f"Expected exactly 18 allocations; missing: {sorted(expected - set(allocations))}")
    task_records = []
    total_gpu_seconds = total_cpu_core_seconds = total_allocation_wall_seconds = 0
    for index in range(18):
        row = allocations[f"{array_id}_{index}"]
        if row["State"] != "COMPLETED" or row["ExitCode"] != "0:0":
            raise RuntimeError(f"Array task not successfully complete: {row['JobID']} {row['State']} {row['ExitCode']}")
        elapsed, cpus = int(row["ElapsedRaw"]), int(row["AllocCPUS"])
        if elapsed < 0 or cpus <= 0 or row["Start"] in ("", "Unknown", "None") or row["End"] in ("", "Unknown", "None"):
            raise RuntimeError(f"Incomplete allocation resource record: {row['JobID']}")
        gpus = gpu_count(row["AllocTRES"])
        configuration, model, seed = CONFIGURATIONS[index // 6], MODELS[(index % 6) // 3], SEEDS[index % 3]
        task_records.append({"array_task_id": index, "job_id": row["JobID"],
                             "run_id": f"{configuration}_{model}_seed{seed}", "configuration": configuration,
                             "model": model, "seed": seed, "state": row["State"], "exit_code": row["ExitCode"],
                             "elapsed_wall_seconds": elapsed, "allocated_gpus": gpus, "allocated_CPUs": cpus,
                             "allocated_GPU_seconds": elapsed * gpus, "allocated_CPU_core_seconds": elapsed * cpus,
                             "AllocTRES": row["AllocTRES"], "start": row["Start"], "end": row["End"]})
        total_gpu_seconds += elapsed * gpus
        total_cpu_core_seconds += elapsed * cpus
        total_allocation_wall_seconds += elapsed
    return {"status": "PASS", "n_task_allocations": 18, "n_job_step_rows": len(steps),
            "n_batch_step_rows": sum(r["JobID"].endswith(".batch") for r in steps),
            "n_extern_step_rows": sum(r["JobID"].endswith(".extern") for r in steps),
            "n_array_parent_rows_excluded_from_totals": len(parent_rows), "task_allocations": task_records,
            "total_allocation_wall_seconds": total_allocation_wall_seconds,
            "allocated_GPU_seconds": total_gpu_seconds, "allocated_GPU_hours": total_gpu_seconds / 3600,
            "allocated_CPU_core_seconds": total_cpu_core_seconds,
            "allocated_CPU_core_hours": total_cpu_core_seconds / 3600,
            "measured_GPU_utilization": None, "measured_GPU_kernel_active_seconds": None,
            "measurement_limit": "Totals use task allocation ElapsedRaw times allocated resources only. They are reserved allocation time, not measured GPU utilization or kernel-active time. .batch/.extern/other job steps and parent summaries are preserved but never added again.",
            "scope": "18 phase-II training array allocations only; excludes phase-I extraction, post-checkpoint inference, CPU evaluation, queue wait, and any jobs outside this array."}


def frozen_checkpoints(stage):
    path = stage / "provenance/checkpoint_freeze.json"
    freeze = read_json(path)
    if freeze.get("status") != "PASS" or freeze.get("all_18_selected_checkpoints_frozen") is not True:
        raise RuntimeError("All-18 checkpoint freeze is not PASS")
    records = freeze.get("checkpoints", [])
    expected = {(c, m, s) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}
    actual = {(r["configuration"], r["model"], int(r["seed"])) for r in records}
    if len(records) != 18 or actual != expected or any(r.get("completed") is not True for r in records):
        raise RuntimeError("Expected exactly the 18 completed prescribed frozen checkpoints")
    for record in records:
        verify_record(stage, record)  # Opaque integrity hashing, never model loading.
        for name in ("selected_checkpoint_record", "completed_record", "command_record"):
            verify_record(stage, record[name])
    return freeze, file_record(stage, path)


def failure_rows(stage, freeze):
    rows = []
    for failure in freeze.get("infrastructure_failures_and_retries", []):
        retry = failure["retry_authorization"]
        authorization = retry["authorization_file"]
        verify_record(stage, authorization)
        original = failure.get("failure_record")
        if original is not None:
            verify_record(stage, original)
        rows.append({"scope": "phase_II_same_seed_infrastructure_retry", "run_id": failure["run_id"],
                     "failed_attempt": failure["attempt"], "classification": retry["reason_class"],
                     "status": "FAILED_ATTEMPT_PRESERVED_SAME_SEED_RETRIED", "reason": retry["reason"],
                     "failure_record": original["path"] if original else "", "failure_record_sha256": original["sha256"] if original else "",
                     "retry_authorization": authorization["path"], "retry_authorization_sha256": authorization["sha256"],
                     "retry_authorized": retry["authorized"], "same_seed_identical_procedure": retry["same_seed_identical_procedure"]})
    for path in sorted(stage.rglob("failure.json")):
        relative = path.relative_to(stage)
        if relative.parts[0] == "runs":
            continue  # Fit failures/retries are already bound by the checkpoint freeze.
        failure = read_json(path)
        reference = file_record(stage, path)
        rows.append({"scope": str(relative.parent), "run_id": failure.get("run_id", ""),
                     "failed_attempt": failure.get("attempt", ""), "classification": "NON_FIT_EXECUTION_FAILURE_NOT_RECLASSIFIED",
                     "status": failure.get("status", ""), "reason": failure.get("error", ""),
                     "failure_record": reference["path"], "failure_record_sha256": reference["sha256"],
                     "retry_authorization": "", "retry_authorization_sha256": "", "retry_authorized": "",
                     "same_seed_identical_procedure": ""})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--array-job-id", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[1-9][0-9]*", args.array_job_id):
        raise ValueError("Array job ID must be a positive numeric Slurm job ID")
    stage = args.stage.resolve()
    provenance = stage / "provenance"
    paths = {name: provenance / name for name in ("slurm_accounting.tsv", "compute_usage.json", "failure_retry_ledger.tsv")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("An execution-provenance output already exists; refusing overwrite")
    freeze, freeze_reference = frozen_checkpoints(stage)
    # Explicit widths prevent a truncated State or JobID, while canonical TSV
    # headers remain exactly the requested field names.
    widths = {"JobID": 64, "State": 32, "ExitCode": 32, "Start": 40, "End": 40,
              "ElapsedRaw": 32, "TotalCPU": 64, "AllocCPUS": 32, "ReqMem": 64,
              "MaxRSS": 64, "AllocTRES": 2048, "NodeList": 2048}
    formatting = [f"{name}%{widths[name]}" for name in FIELDS]
    command = ["sacct", "--array", "--parsable2", "--noheader", f"--jobs={args.array_job_id}",
               "--format=" + ",".join(formatting)]
    collection_start = now()
    timer = time.monotonic()
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60, check=False)
    if process.returncode != 0:
        raise RuntimeError(f"Read-only sacct failed ({process.returncode}): {process.stderr.decode('utf-8', errors='replace')}")
    rows = parse_sacct(process.stdout)
    summary = allocation_summary(rows, args.array_job_id)
    by_run = {r["run_id"]: r for r in freeze["checkpoints"]}
    for task in summary["task_allocations"]:
        checkpoint = by_run[task["run_id"]]
        command_record = read_json(local_path(stage, checkpoint["command_record"]["path"]))
        if (int(command_record["environment"]["SLURM_ARRAY_TASK_ID"]) != task["array_task_id"]
                or command_record["environment"]["PYTHONHASHSEED"] != str(task["seed"])
                or command_record["returncode"] != 0):
            raise RuntimeError(f"Task-to-frozen-run command provenance mismatch: {task['run_id']}")
        task["checkpoint_sha256"] = checkpoint["sha256"]
        task["checkpoint_path"] = checkpoint["path"]
    failures = failure_rows(stage, freeze)
    accounting_bytes = tsv_bytes(rows, FIELDS)
    failure_bytes = tsv_bytes(failures, FAILURE_FIELDS)
    record = {"module": "COPD-V2-INTERNAL", "version": "internal-training-1.0", "status": "PASS",
              "array_job_id": args.array_job_id, "created_utc": now(), "checkpoint_freeze": freeze_reference,
              "checkpoint_count": 18, "all_checkpoint_hashes_verified": True, "array_accounting": summary,
              "collection": {"argv": command, "started_utc": collection_start, "completed_utc": now(),
                             "elapsed_seconds": time.monotonic() - timer, "returncode": process.returncode,
                             "raw_stdout_sha256": hashlib.sha256(process.stdout).hexdigest(),
                             "raw_stdout_bytes": len(process.stdout), "raw_stdout": process.stdout.decode("utf-8"),
                             "raw_stderr_sha256": hashlib.sha256(process.stderr).hexdigest(),
                             "raw_stderr": process.stderr.decode("utf-8", errors="replace")},
              "accounting_table": {"path": str(paths["slurm_accounting.tsv"].relative_to(stage)),
                                   "bytes": len(accounting_bytes), "sha256": hashlib.sha256(accounting_bytes).hexdigest()},
              "failure_retry_ledger": {"path": str(paths["failure_retry_ledger.tsv"].relative_to(stage)),
                                       "rows": len(failures), "bytes": len(failure_bytes),
                                       "sha256": hashlib.sha256(failure_bytes).hexdigest()},
              "collector": file_record(stage, Path(__file__)), "scientific_computation_performed": False,
              "model_execution": False, "test_external_or_candidate_data_access": False}
    compute_bytes = (json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    provenance.mkdir(parents=True, exist_ok=True)
    for name, payload in (("slurm_accounting.tsv", accounting_bytes), ("failure_retry_ledger.tsv", failure_bytes), ("compute_usage.json", compute_bytes)):
        with paths[name].open("xb") as stream:
            stream.write(payload)
    print(json.dumps({"status": "PASS", "array_job_id": args.array_job_id, "completed_allocations": 18,
                      "allocated_GPU_hours": summary["allocated_GPU_hours"], "failure_retry_rows": len(failures)}))


if __name__ == "__main__":
    main()
