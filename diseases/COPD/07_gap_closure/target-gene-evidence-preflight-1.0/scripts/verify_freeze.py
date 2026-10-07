#!/usr/bin/env python3
"""Recheck the target-gene preflight payload checksum ledger."""

import csv
import json
from pathlib import Path
from validate_preflight import STAGE, GAP, digest, run

freeze = json.loads((STAGE / "provenance/freeze.json").read_text())
ledger = STAGE / freeze["checksum_ledger"]
assert digest(ledger) == freeze["checksum_ledger_sha256"]
with ledger.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
assert len(rows) == freeze["payload_file_count"]
for row in rows:
    path = STAGE / row["relative_path"]
    assert path.is_file() and path.stat().st_size == int(row["bytes"]) and digest(path) == row["sha256"], path
for name, sha in freeze["register_sha256"].items():
    assert digest(GAP / name) == sha, name
assert run()["status"] == "PASS"
print(f"PASS: {len(rows)} payload checksums, 3 register hashes, validation")
