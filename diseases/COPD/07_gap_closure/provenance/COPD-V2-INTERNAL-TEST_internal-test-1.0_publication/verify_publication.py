#!/usr/bin/env python3
"""Read-only Git-object publication audit; no scientific modules are imported.

Run against a staged tree or a freshly fetched remote revision. Historical
scientific artifacts are compared by Git tree entries, never opened. Only
internal-test payloads, its publication metadata, the three authorized registers,
.gitignore, and the existing training checkpoint LFS pointers are read.
"""

import argparse
import csv
import hashlib
import io
import json
import subprocess
from datetime import datetime, timezone


BASELINE = "07d72dad0f57d605f5f50bbc7fa6d2c7843a94ba"
STAGE = "diseases/COPD/07_gap_closure/internal-test-1.0"
PUB = "diseases/COPD/07_gap_closure/provenance/COPD-V2-INTERNAL-TEST_internal-test-1.0_publication"
TRAINING = "diseases/COPD/07_gap_closure/internal-training-1.0"
EXTERNAL = {"cache/features_canonical.npy", "cache/features_rc.npy"}
LEDGER_SHA = "c354375030411b346145709da0ec6ef1512ca981648813789e38675750a7f1f2"
FREEZE_SHA = "db514e7877ce948718be3632fb8e8626350da03b8fce3878ffa2acb2698a3d55"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="Working or bare Git repository")
    parser.add_argument("--revision", required=True, help="Commit or staged tree object")
    args = parser.parse_args()

    def git(*arguments):
        return subprocess.check_output(["git", "-C", args.repo, *arguments])

    def blob(revision, path):
        return git("show", f"{revision}:{path}")

    def digest(data):
        return hashlib.sha256(data).hexdigest()

    def require(condition, message):
        if not condition:
            raise RuntimeError(message)

    def tree(revision):
        entries = {}
        for record in git("ls-tree", "-rz", revision).split(b"\0"):
            if record:
                metadata, path = record.split(b"\t", 1)
                entries[path.decode()] = metadata.decode()
        return entries

    revision = git("rev-parse", args.revision).decode().strip()
    previous, published = tree(BASELINE), tree(revision)
    ledger = blob(revision, f"{STAGE}/provenance/artifact_checksums.tsv")
    freeze = blob(revision, f"{STAGE}/provenance/freeze.json")
    require(digest(ledger) == LEDGER_SHA, "Original ledger changed")
    require(digest(freeze) == FREEZE_SHA, "Original freeze changed")
    rows = list(csv.DictReader(io.StringIO(ledger.decode()), delimiter="\t"))
    require(len(rows) == 172, "Unexpected frozen payload count")
    require(len({r['path'] for r in rows}) == 172, "Duplicate frozen path")
    require(sum(int(r['bytes']) for r in rows) == 1093208401, "Frozen total mismatch")
    require(EXTERNAL <= {r['path'] for r in rows}, "External caches absent from ledger")
    expected = {f"{STAGE}/{r['path']}" for r in rows if r['path'] not in EXTERNAL}
    expected |= {f"{STAGE}/provenance/artifact_checksums.tsv", f"{STAGE}/provenance/freeze.json"}
    actual = {p for p in published if p.startswith(STAGE + "/")}
    require(actual == expected, f"Stage membership mismatch: {sorted(actual ^ expected)}")
    hosted_bytes = 0
    for row in rows:
        if row['path'] in EXTERNAL:
            continue
        path = f"{STAGE}/{row['path']}"
        data = blob(revision, path)
        require(len(data) == int(row['bytes']), f"Size mismatch: {path}")
        require(digest(data) == row['sha256'], f"Hash mismatch: {path}")
        hosted_bytes += len(data)

    updates = json.loads(blob(revision, f"{STAGE}/provenance/shared_register_updates.json"))
    require(updates['baseline_commit'] == BASELINE, "Register baseline changed")
    register_paths = {r['path'] for r in updates['registers']}
    require(len(register_paths) == 3, "Unexpected register count")
    append_rows = {}
    for row in updates['registers']:
        old = blob(BASELINE, row['path'])
        new = blob(revision, row['path'])
        require(len(old) == row['before_bytes'] and digest(old) == row['before_sha256'], "Register historical bytes mismatch")
        require(new.startswith(old), "Register is not append-only")
        require(len(new) == row['bytes'] and digest(new) == row['sha256'], "Frozen register update changed")
        require(len(new) - len(old) == row['appended_bytes'], "Register append size mismatch")
        require(new[len(old):].count(b'\n') == row['appended_rows'], "Register append row count mismatch")
        append_rows[row['path']] = row['appended_rows']

    mutable = register_paths | {'.gitignore'}
    unchanged = 0
    for path, entry in previous.items():
        if path not in mutable:
            require(published.get(path) == entry, f"Historical artifact changed: {path}")
            unchanged += 1
    new_paths = set(published) - set(previous)
    require(all(p.startswith(STAGE + '/') or p.startswith(PUB + '/') for p in new_paths), "Unrelated new artifact")
    old_ignore = blob(BASELINE, '.gitignore')
    new_ignore = blob(revision, '.gitignore')
    require(new_ignore.startswith(old_ignore), ".gitignore historical prefix changed")
    for path in EXTERNAL:
        require(f"/{STAGE}/{path}\n".encode() in new_ignore, "External cache not explicitly ignored")
    require(f"/{STAGE}/inputs/canonical_sequences.npy\n".encode() not in new_ignore, "Hosted sequence array incorrectly ignored")
    for name in ('README.md', 'storage_manifest.json', 'storage_inventory.tsv'):
        metadata = blob(revision, f"{PUB}/{name}")
        for path in EXTERNAL:
            require(path.encode() in metadata, f"Cache undocumented in {name}")

    checkpoints = sorted(p for p in previous if p.startswith(TRAINING + '/') and p.endswith('.keras'))
    require(len(checkpoints) == 18, "Unexpected original checkpoint count")
    checkpoint_bytes = 0
    for path in checkpoints:
        pointer = blob(revision, path).decode()
        require(pointer.startswith('version https://git-lfs.github.com/spec/v1\n'), "Checkpoint is not original LFS pointer")
        require('oid sha256:' in pointer, "Missing LFS SHA-256")
        checkpoint_bytes += int(next(line[5:] for line in pointer.splitlines() if line.startswith('size ')))
    require(checkpoint_bytes == 6301157652, "Checkpoint aggregate size changed")
    require(not any(p.endswith('.keras') for p in new_paths), "Duplicate/new checkpoint archive")

    result = {
        'status': 'PASS',
        'verified_utc': datetime.now(timezone.utc).isoformat(),
        'revision': revision,
        'baseline_commit': BASELINE,
        'frozen_scientific_payload_count': len(rows),
        'frozen_scientific_payload_bytes': 1093208401,
        'git_hosted_scientific_payload_count': 170,
        'git_hosted_stage_files_including_ledger_and_freeze': len(actual),
        'git_hosted_stage_bytes': hosted_bytes + len(ledger) + len(freeze),
        'hosted_payload_sha256_verification': '170/170 PASS',
        'original_ledger_and_freeze_sha256': 'PASS',
        'external_intermediates': [r for r in rows if r['path'] in EXTERNAL],
        'external_intermediates_absent_from_git_tree': True,
        'external_intermediates_documented_and_ignored': True,
        'append_only_register_rows': append_rows,
        'historical_git_tree_entries_unchanged': unchanged,
        'historical_content_opened': False,
        'historical_content_read_exceptions': ['three authorized register prefixes', '.gitignore', '18 training LFS pointers (not model archives)'],
        'checkpoint_lfs_pointers_unchanged': 18,
        'checkpoint_lfs_payload_bytes_unchanged': checkpoint_bytes,
        'new_checkpoint_archives': 0,
        'scientific_computation_performed': False,
        'external_benchmark_or_candidate_content_opened': False,
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
