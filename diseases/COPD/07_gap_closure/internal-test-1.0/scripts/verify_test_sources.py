#!/usr/bin/env python3
"""Independent pre-execution source/preservation gate; never executes a model.

Only the authorized pretraining, internal-training and prospective test-stage
files and explicitly allowlisted frozen reference inputs are opened. Historical
Git trees are inventoried as Git object metadata, without reading benchmark or
candidate payloads. All evidence is exclusively created in the new test stage.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


BASELINE_COMMIT = "07d72dad0f57d605f5f50bbc7fa6d2c7843a94ba"
GAP = "diseases/COPD/07_gap_closure"
TRAINING = GAP + "/internal-training-1.0"
PRETRAINING = GAP + "/pretraining-1.1"
TEST = GAP + "/internal-test-1.0"
TRAINING_FREEZE_SHA = "c151e61dba04a9e5500534e507a392d51e1af0601a6e62565288bc6b276e6de4"
PRETRAINING_FREEZE_SHA = "6933f98bd6a1a7b7afe5adddc8a6bf9d7b1201df4b25f06d0cd8c35e183941a0"
CHECKPOINT_FREEZE_SHA = "d2330e4881c477cb7d8e5b2b621940674761133b779c89d0db1ad88635eee379"
PHASE_I_SHA = "483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6"
REFERENCES = {
    "diseases/COPD/04_modeling/trednet/fasta/hg38.fa",
    "diseases/COPD/04_modeling/trednet/fasta/hg38.fa.fai",
    "diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5",
}
REGISTERS = [GAP + "/" + name for name in (
    "activity_log.tsv", "gap_closure_decision_register.tsv", "gap_closure_result_register.tsv")]


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def repo_file(repo, relative):
    relative = str(relative)
    if Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise RuntimeError("Absolute/traversing source path rejected: " + relative)
    allowed = any(relative.startswith(prefix + "/") for prefix in (TRAINING, PRETRAINING, TEST))
    if not allowed and relative not in REFERENCES and relative not in REGISTERS:
        raise RuntimeError("Source path outside explicit authorization: " + relative)
    # Preserve lexical paths because the authorized hg38/phase-I references are
    # existing symlinks to the original managed model/reference storage.
    path = repo / relative
    if not path.is_file():
        raise RuntimeError("Required frozen source is unavailable: " + relative)
    return path


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--prospective-freeze", type=Path)
    args = parser.parse_args()
    repo, stage = args.repo.resolve(), args.stage.resolve()
    if stage != repo / TEST:
        raise RuntimeError("Only the authorized new internal-test-1.0 stage is writable")
    provenance = stage / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    output = provenance / "input_integrity_verified_before.json"
    if output.exists() or (provenance / "input_integrity_verified_before.tsv").exists():
        raise RuntimeError("Pre-execution evidence already exists; refusing replacement")
    for name in ("predictions", "results"):
        path = stage / name
        if path.exists() and any(p.is_file() for p in path.rglob("*")):
            raise RuntimeError("Prospective source gate cannot first run after test outcomes")
    started, started_utc = time.monotonic(), utcnow()
    checks, checked = [], {}

    def check(relative, expected_sha, expected_bytes=None, scope="source"):
        path = repo_file(repo, relative)
        key = str(relative)
        if key not in checked:
            checked[key] = {"path": key, "bytes": path.stat().st_size, "sha256": sha256(path)}
        actual = checked[key]
        passed = actual["sha256"] == expected_sha and (
            expected_bytes is None or actual["bytes"] == int(expected_bytes))
        checks.append({"scope": scope, "path": key, "expected_bytes": expected_bytes,
                       "observed_bytes": actual["bytes"], "expected_sha256": expected_sha,
                       "observed_sha256": actual["sha256"], "status": "PASS" if passed else "FAIL"})
        if not passed:
            raise RuntimeError("Frozen source SHA-256/size mismatch: " + key)
        return actual

    try:
        head = git(repo, "rev-parse", "HEAD").decode().strip()
        if head != BASELINE_COMMIT:
            raise RuntimeError("Published historical baseline commit differs from authorized HEAD")
        freeze_path = args.prospective_freeze or provenance / "prospective_specification_freeze.json"
        if freeze_path.resolve() != provenance / "prospective_specification_freeze.json":
            raise RuntimeError("Unexpected prospective specification freeze path")
        prospective = read_json(freeze_path)
        if prospective.get("status") != "PASS" or prospective.get("before_any_test_model_inference") is not True:
            raise RuntimeError("Prospective specification is not explicitly frozen before inference")
        for key in ("specification", "markdown"):
            record = prospective[key]
            check(TEST + "/" + record["path"], record["sha256"], record["bytes"], "prospective_" + key)
        for record in prospective["sources"]:
            check(record["path"], record["sha256"], record.get("bytes"), "prospective_source")
        prospective_record = {"path": str(freeze_path.relative_to(stage)), "bytes": freeze_path.stat().st_size,
                              "sha256": sha256(freeze_path)}
        training_freeze = check(TRAINING + "/provenance/freeze.json", TRAINING_FREEZE_SHA, scope="training_freeze")
        training = read_json(repo / training_freeze["path"])
        ledger = training["artifact_checksums"]
        check(TRAINING + "/" + ledger["path"], ledger["sha256"], scope="training_ledger")
        with (repo / TRAINING / ledger["path"]).open(newline="", encoding="utf-8") as stream:
            training_rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(training_rows) != 470 or len({r["path"] for r in training_rows}) != 470:
            raise RuntimeError("Training ledger does not contain exactly 470 unique frozen payloads")
        for record in training_rows:
            check(TRAINING + "/" + record["path"], record["sha256"], record["bytes"], "training_payload")
        if sum(int(r["bytes"]) for r in training_rows) != int(training["payload_bytes"]):
            raise RuntimeError("Training ledger total byte count mismatch")
        pretraining_record = check(PRETRAINING + "/provenance/freeze.json", PRETRAINING_FREEZE_SHA,
                                   scope="pretraining_freeze")
        pretraining = read_json(repo / pretraining_record["path"])
        pretraining_ledger = pretraining["artifact_checksum_ledger"]
        check(pretraining_ledger, pretraining["artifact_checksum_ledger_sha256"], scope="pretraining_ledger")
        with (repo / pretraining_ledger).open(newline="", encoding="utf-8") as stream:
            pretraining_rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(pretraining_rows) != 124 or len({r["path"] for r in pretraining_rows}) != 124:
            raise RuntimeError("Pretraining ledger does not contain exactly 124 unique frozen payloads")
        for record in pretraining_rows:
            if not record["path"].startswith(PRETRAINING + "/"):
                raise RuntimeError("Pretraining ledger source escaped authorized stage")
            check(record["path"], record["sha256"], record["bytes"], "pretraining_payload")
        checkpoint_record = check(TRAINING + "/provenance/checkpoint_freeze.json", CHECKPOINT_FREEZE_SHA,
                                  scope="checkpoint_freeze")
        checkpoint_freeze = read_json(repo / checkpoint_record["path"])
        if checkpoint_freeze.get("status") != "PASS" or checkpoint_freeze.get("all_18_selected_checkpoints_frozen") is not True:
            raise RuntimeError("Original all-18 checkpoint freeze is not PASS")
        expected = {(c, m, seed) for c in ("V2-A", "V2-B", "V2-C")
                    for m in ("enhancer", "h3k27me3") for seed in (104729, 130363, 155921)}
        checkpoint_rows, seen = [], set()
        for record in checkpoint_freeze["checkpoints"]:
            key = record["configuration"], record["model"], int(record["seed"])
            if key not in expected or key in seen or record.get("completed") is not True:
                raise RuntimeError("Unexpected, duplicate or incomplete frozen checkpoint")
            seen.add(key)
            actual = check(TRAINING + "/" + record["path"], record["sha256"], record["bytes"], "selected_checkpoint")
            checkpoint_rows.append({"configuration": key[0], "model": key[1], "seed": key[2],
                                    **actual, "status": "PASS", "preserved_original_archive": True})
        if seen != expected:
            raise RuntimeError("Not all 18 original checkpoints were verified")
        check("diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5", PHASE_I_SHA, scope="frozen_phase_I")
        # Git object metadata and status filenames only: no benchmark or candidate
        # files are opened/parsed by this auditor, including for opaque hashing.
        tree = []
        for entry in git(repo, "ls-tree", "-r", "-l", "-z", "HEAD").split(b"\0"):
            if not entry:
                continue
            metadata, raw_path = entry.split(b"\t", 1)
            mode, kind, oid, size = metadata.decode().split()
            tree.append({"path": raw_path.decode(), "mode": mode, "type": kind,
                         "git_object": oid, "git_blob_bytes": size})
        raw_status = git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=all", "--ignore-submodules=all")
        status_entries = [x.decode() for x in raw_status.split(b"\0") if x]
        for entry in status_entries:
            if len(entry) < 4 or not entry[3:].startswith(TEST + "/"):
                raise RuntimeError("Pre-execution worktree contains unrelated historical modifications: " + entry)
        registers = []
        for relative in REGISTERS:
            path = repo_file(repo, relative)
            registers.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path),
                              "preservation_rule": "Original byte prefix must remain exact; append-only stage records permitted later"})
        write_tsv(provenance / "historical_git_tree_before.tsv", tree,
                  ["path", "mode", "type", "git_object", "git_blob_bytes"])
        write_json(provenance / "historical_preservation_baseline.json", {
            "status": "PASS", "created_utc": utcnow(), "baseline_commit": head,
            "tracked_entries": len(tree), "working_status_filenames_only": status_entries,
            "tree_metadata_file": "provenance/historical_git_tree_before.tsv", "registers": registers,
            "external_benchmark_files_opened": False, "external_benchmark_payloads_hashed": False,
            "candidate_payloads_opened": False,
            "historical_tree_interpretation": "Git object identities only; external payload preservation is not claimed from newly reading those forbidden files"})
        write_json(provenance / "checkpoint_verification.json", {
            "status": "PASS", "created_utc": utcnow(), "n_selected_original_checkpoints": 18,
            "checkpoint_freeze": checkpoint_record, "prospective_specification_freeze": prospective_record,
            "checkpoints": checkpoint_rows, "model_loads": 0, "inference_performed": False})
        write_tsv(provenance / "input_integrity_verified_before.tsv", checks,
                  ["scope", "path", "expected_bytes", "observed_bytes", "expected_sha256", "observed_sha256", "status"])
        result = {"status": "PASS", "started_utc": started_utc, "completed_utc": utcnow(),
                  "argv": sys.argv, "wall_seconds": time.monotonic() - started,
                  "prospective_specification_freeze": prospective_record,
                  "published_baseline_commit": head, "checks": len(checks), "failures": 0,
                  "training_frozen_payloads_verified": len(training_rows),
                  "pretraining_frozen_payloads_verified": len(pretraining_rows),
                  "selected_original_checkpoints_verified": len(checkpoint_rows),
                  "input_hashes_verified": True, "all_historical_authorized_sources_unchanged": True,
                  "files": list(checked.values()), "external_benchmark_access": False,
                  "candidate_access": False, "model_inference_performed": False,
                  "test_preparation_authorized": True}
        write_json(output, result)
        print(json.dumps({k: result[k] for k in ("status", "checks", "failures", "training_frozen_payloads_verified",
              "pretraining_frozen_payloads_verified", "selected_original_checkpoints_verified", "wall_seconds")}))
    except BaseException as exc:
        failure_path = provenance / "input_integrity_before_failure.json"
        if not failure_path.exists():
            write_json(failure_path, {"status": "FAIL", "started_utc": started_utc, "failed_utc": utcnow(),
                       "error": str(exc), "exception_type": type(exc).__name__, "checks_completed": len(checks),
                       "test_preparation_authorized": False, "inference_authorized": False})
        raise


if __name__ == "__main__":
    main()
