#!/usr/bin/env python3
"""Read-only publication verification; never loads a model or runs analysis.

Run from the repository root after fetching all 18 checkpoint objects into a
fresh, isolated LFS storage directory. This script does not fetch or write files.
It prints a JSON receipt to stdout, or exits unsuccessfully on any mismatch.
"""

import argparse
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import subprocess


STAGE = "diseases/COPD/07_gap_closure/internal-training-1.0"
PUBLICATION = (
    "diseases/COPD/07_gap_closure/provenance/"
    "COPD-V2-INTERNAL_internal-training-1.0_publication"
)
BASELINE = "7c0caf4af5f8f346089c92f43b88161331ca420b"
FREEZE_HASHES = {
    "provenance/artifact_checksums.tsv":
        "661d2e95956068819713ed5d2ba01591f0019f5cb98f61e5b61e221aec8fd691",
    "provenance/freeze.json":
        "c151e61dba04a9e5500534e507a392d51e1af0601a6e62565288bc6b276e6de4",
}
EXTERNAL = {
    STAGE + "/cache/features_canonical.npy",
    STAGE + "/cache/features_rc.npy",
    STAGE + "/inputs/canonical_sequences.npy",
}


def require(condition, description):
    if not condition:
        raise RuntimeError(description)


def git(*args):
    return subprocess.check_output(["git", *args])


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def blob(commit, path):
    return git("show", commit + ":" + path)


def tree(commit):
    entries = {}
    for item in git("ls-tree", "-r", "-z", commit).split(b"\0"):
        if item:
            metadata, path = item.split(b"\t", 1)
            entries[path.decode()] = metadata.decode()
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--lfs-storage", required=True, type=Path)
    parser.add_argument("--remote", default="codex-manuscript")
    args = parser.parse_args()
    root = Path(git("rev-parse", "--show-toplevel").decode().strip())
    require(Path.cwd().resolve() == root.resolve(), "Run from repository root")
    commit = git("rev-parse", args.commit + "^{commit}").decode().strip()
    remote = git("ls-remote", args.remote, "refs/heads/main").decode().split()[0]
    require(remote == commit, "Remote main does not resolve to requested commit")
    published = tree(commit)
    inventory = list(csv.DictReader(io.StringIO(blob(
        commit, PUBLICATION + "/storage_inventory.tsv"
    ).decode()), delimiter="\t"))
    require(len(inventory) == 477, "Inventory must contain 477 records")
    require(len({r["path"] for r in inventory}) == 477, "Duplicate inventory path")
    for path, expected in FREEZE_HASHES.items():
        require(hashlib.sha256(blob(commit, STAGE + "/" + path)).hexdigest()
                == expected, "Original freeze metadata changed: " + path)
    ledger = list(csv.DictReader(io.StringIO(blob(
        commit, STAGE + "/provenance/artifact_checksums.tsv"
    ).decode()), delimiter="\t"))
    by_path = {r["path"]: r for r in inventory}
    require(len(ledger) == 470, "Original scientific ledger count changed")
    for row in ledger:
        inv = by_path[STAGE + "/" + row["path"]]
        require(inv["sha256"] == row["sha256"]
                and int(inv["bytes"]) == int(row["bytes"]),
                "Inventory differs from frozen ledger: " + row["path"])

    counts, byte_totals, checkpoints, excluded = {}, {}, [], []
    hosted_stage = set()
    for row in inventory:
        path, category = row["path"], row["storage_class"]
        size, expected = int(row["bytes"]), row["sha256"]
        local = root / path
        require(local.stat().st_size == size, "Local size mismatch: " + path)
        local_hash = digest(local)
        require(local_hash == expected, "Local SHA-256 mismatch: " + path)
        counts[category] = counts.get(category, 0) + 1
        byte_totals[category] = byte_totals.get(category, 0) + size
        if category == "external_reconstructable_intermediate":
            require(path in EXTERNAL and path not in published,
                    "External array unexpectedly Git tracked: " + path)
            excluded.append({"path": path, "bytes": size, "sha256": expected,
                             "local_bytes_preserved": True, "git_hosted": False})
            continue
        require(path in published, "Missing published path: " + path)
        if path.startswith(STAGE + "/"):
            hosted_stage.add(path)
        content = blob(commit, path)
        if category == "git_lfs_original_checkpoint":
            pointer = ("version https://git-lfs.github.com/spec/v1\n"
                       f"oid sha256:{expected}\nsize {size}\n").encode()
            require(content == pointer, "Incorrect LFS pointer: " + path)
            downloaded = (args.lfs_storage / "objects" / expected[:2]
                          / expected[2:4] / expected)
            require(downloaded.stat().st_size == size,
                    "Downloaded object size mismatch: " + path)
            remote_hash = digest(downloaded)
            require(remote_hash == expected, "Downloaded SHA-256 mismatch: " + path)
            checkpoints.append({"path": path, "bytes": size,
                                "frozen_sha256": expected,
                                "local_sha256": local_hash,
                                "downloaded_sha256": remote_hash,
                                "downloaded_object_path": str(downloaded),
                                "pointer_resolves_and_hash_matches": True})
        else:
            require(category == "git_ordinary", "Unknown storage class")
            require(len(content) == size and hashlib.sha256(content).hexdigest()
                    == expected, "Published ordinary bytes differ: " + path)
    require(counts == {"git_ordinary": 456,
                       "git_lfs_original_checkpoint": 18,
                       "external_reconstructable_intermediate": 3},
            "Storage class counts differ")
    require(hosted_stage == {p for p in published if p.startswith(STAGE + "/")}
            and len(hosted_stage) == 469, "Unexpected stage tree paths")

    registers = json.loads(blob(commit, STAGE +
                                "/provenance/shared_register_updates.json"))
    permitted_updates = {r["path"] for r in registers["files"]}
    permitted_updates.update({".gitattributes", ".gitignore"})
    baseline = tree(BASELINE)
    for path, entry in baseline.items():
        require(path in published, "Historical path deleted: " + path)
        if path in permitted_updates:
            require(blob(commit, path).startswith(blob(BASELINE, path)),
                    "Historical prefix changed: " + path)
        else:
            require(published[path] == entry, "Historical artifact changed: " + path)
    changed = git("diff", "--name-only", "-z", BASELINE, commit).split(b"\0")
    changed_paths = [p.decode() for p in changed if p]
    require(all(p in permitted_updates or p.startswith(STAGE + "/")
                or p.startswith(PUBLICATION + "/") for p in changed_paths),
            "Unrelated paths included in publication")
    ordinary_paths = [p for p in changed_paths if not p.endswith(".keras")]
    ordinary_bytes = sum(len(blob(commit, p)) for p in ordinary_paths)
    pointer_bytes = sum(len(blob(commit, r["path"])) for r in checkpoints)
    print(json.dumps({
        "status": "PASS", "verified_utc": datetime.datetime.now(
            datetime.timezone.utc).isoformat(),
        "artifact_commit": commit, "remote_main_verified": remote,
        "baseline_commit": BASELINE,
        "verification_kind": "Read-only byte/hash/tree preservation checks only",
        "isolated_lfs_storage": str(args.lfs_storage),
        "inventory_class_counts": counts, "inventory_class_bytes": byte_totals,
        "local_frozen_payloads_rehashed_unchanged": len(ledger),
        "original_freeze_and_ledger_hashes_unchanged": FREEZE_HASHES,
        "ordinary_stage_files_verified": 451,
        "remote_lfs_objects_verified": len(checkpoints),
        "checkpoints": checkpoints, "external_intermediates": excluded,
        "historical_tree_entries_checked": len(baseline),
        "historical_unchanged_except_authorized_append_only_paths":
            sorted(permitted_updates),
        "V1_benchmark_pretraining_1_0_pretraining_1_1_unchanged": True,
        "publication_payload_at_artifact_commit": {
            "ordinary_files": len(ordinary_paths), "ordinary_bytes": ordinary_bytes,
            "lfs_pointer_bytes": pointer_bytes,
            "lfs_original_checkpoint_bytes": byte_totals["git_lfs_original_checkpoint"],
            "combined_original_content_bytes": ordinary_bytes +
                byte_totals["git_lfs_original_checkpoint"],
            "accounting_note": "Changed paths counted at full content size; excludes "
                "Git history/pack overhead and the later verification receipt. "
                "LFS pointer bytes reported separately."
        }
    }, indent=2))


if __name__ == "__main__":
    main()
