#!/usr/bin/env python3
"""Prepare only the frozen pretraining-1.1 common C-task chr8-9 panels.

This is input preparation, not prediction or rematching. No frozen-stage module
is imported or executed. All source strings and their order are retained in the
common-panel outputs; only cache-row and nucleotide-orientation columns are added.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import resource
import sys
import time
from datetime import datetime, timezone

import numpy as np

MODELS = ("enhancer", "h3k27me3")
ALLOWED = frozenset(("chr8", "chr9"))
COLS = ("interval_id", "chrom", "core_start", "core_end", "input_start", "input_end",
        "partition", "validation_role", "component_id", "sequence_available", "sequence_length",
        "sequence_sha256", "canonical_rc_sequence_sha256", "gc_fraction", "non_acgt_fraction")
NORMALIZE = bytes(value if value in b"ACGT" else ord("N") for value in range(256))
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")
EXPECTED = {
    "enhancer": {"rows": 25390, "positive": 12695, "negative": 12695, "components": 8754,
                 "encoded_identities": 25390, "chromosome_counts": {"chr8": 11963, "chr9": 13427}},
    "h3k27me3": {"rows": 3508, "positive": 1754, "negative": 1754, "components": 1809,
                 "encoded_identities": 3508, "chromosome_counts": {"chr8": 1724, "chr9": 1784}},
}
PRETRAINING = "diseases/COPD/07_gap_closure/pretraining-1.1"
SOURCE_RECORDS = {
    f"{PRETRAINING}/data/evaluation/enhancer_common_challenge_panel.tsv.gz":
        (4962271, "4f374ebf993b1e3b6a5fff5aba51ac6e0efaf6859300bf752f56245a7db65c20"),
    f"{PRETRAINING}/data/evaluation/h3k27me3_common_challenge_panel.tsv.gz":
        (719591, "895979a83b12d0bb372dd93a1baaa1a83d174e96a57e17acc4f5f3bea8775c26"),
    f"{PRETRAINING}/data/configurations/V2-C_enhancer_interval_manifest.tsv.gz":
        (44914883, "5165080dbeba1a43decb57286c9ed44e820eaa73617f252cecda6ebe4321e72f"),
    f"{PRETRAINING}/data/configurations/V2-C_h3k27me3_interval_manifest.tsv.gz":
        (6293534, "1d325e61ae5c41bfabe955c12e0b84f683a80472fb855dbac9a19fabd1cbc7ec"),
    "diseases/COPD/04_modeling/trednet/fasta/hg38.fa":
        (3273481150, "5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51"),
    "diseases/COPD/04_modeling/trednet/fasta/hg38.fa.fai":
        (19381, "3b425de206296a5c8053023fa5ca61da43cfe78c1737c12e58c83367c7e83c21"),
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path, base):
    path, base = Path(path), Path(base)
    return {"path": str(path.relative_to(base)), "bytes": path.stat().st_size, "sha256": sha(path)}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def read_table(path, test_only=False):
    with gzip.open(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames
        if not columns or len(columns) != len(set(columns)):
            raise RuntimeError("Missing/duplicate source columns")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise RuntimeError("Ragged source table")
            # Include any test indicator, then require all indicators to agree.
            if not test_only or row.get("partition") == "test" or row.get("validation_role") == "test" or row.get("chrom") in ALLOWED:
                rows.append(row)
        return columns, rows


def write_table(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as binary:
        with gzip.GzipFile(filename="", fileobj=binary, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def table_digest(columns, rows):
    """Hash unambiguous JSON strings in frozen row/column order, not float casts."""
    digest = hashlib.sha256()
    for item in [list(columns), *([row[key] for key in columns] for row in rows)]:
        digest.update(json.dumps(item, ensure_ascii=True, separators=(",", ":")).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def verify_geometry(row):
    try:
        start, end, core_start, core_end = (int(row[key]) for key in ("input_start", "input_end", "core_start", "core_end"))
    except (KeyError, ValueError) as exc:
        raise RuntimeError("Invalid frozen geometry") from exc
    if (row["chrom"] not in ALLOWED or row["partition"] != "test" or row["validation_role"] != "test" or
            end - start != 2001 or core_end - core_start != 1000 or start != core_start - 501 or end != core_end + 500 or
            start < 0 or row["sequence_available"] != "1" or row["sequence_length"] != "2001" or
            row["interval_id"] != f"{row['chrom']}:{core_start}-{core_end}" or
            not row["component_id"].startswith(row["chrom"] + ":")):
        raise RuntimeError("Frozen geometry/role/eligibility failure: " + row["interval_id"])
    for column in ("sequence_sha256", "canonical_rc_sequence_sha256"):
        value = row[column]
        if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
            raise RuntimeError("Malformed sequence hash")


def validate_panel(model, columns, rows, configuration_columns, configuration_rows, expected=None):
    expected = EXPECTED[model] if expected is None else expected
    if list(columns) != list(configuration_columns):
        raise RuntimeError("Panel/configuration column disagreement")
    if "cache_row" in columns or "forward_is_canonical" in columns:
        raise RuntimeError("Frozen source unexpectedly contains cache columns")
    if len(rows) != expected["rows"] or len(configuration_rows) != expected["rows"]:
        raise RuntimeError("Frozen test population count mismatch")
    by_id = {row["interval_id"]: row for row in configuration_rows}
    if len(by_id) != len(configuration_rows) or len({row["interval_id"] for row in rows}) != len(rows):
        raise RuntimeError("Duplicate frozen test interval")
    if set(by_id) != {row["interval_id"] for row in rows}:
        raise RuntimeError("Frozen test membership disagreement")
    for row in rows:
        if row != by_id[row["interval_id"]]:
            raise RuntimeError("Frozen test row metadata/label disagreement: " + row["interval_id"])
        verify_geometry(row)
        if row["configuration"] != "V2-C" or row["model"] != model or row["label"] not in ("0", "1"):
            raise RuntimeError("Common C-task label/configuration failure")
        if row["label"] == "1" and int(row[model + "_same_lobe_peak_support"]) < 1:
            raise RuntimeError("Same-lobe positive-label eligibility failure")
        if row["label"] == "0" and row[model + "_accessible_control_eligible"] != "1":
            raise RuntimeError("Accessible-control eligibility failure")
        if row["non_acgt_any"] != str(int(float(row["non_acgt_fraction"]) > 0)):
            raise RuntimeError("Ambiguous-base flag disagreement")
        rank = float(row["atac_signal_percentile_max_train_only"])
        if not np.isfinite(rank) or not 0 <= rank <= 1:
            raise RuntimeError("Frozen train-derived ATAC rank is invalid")
    observed = {"rows": len(rows), "positive": sum(row["label"] == "1" for row in rows),
                "negative": sum(row["label"] == "0" for row in rows),
                "components": len({row["component_id"] for row in rows}),
                "encoded_identities": len({row["canonical_rc_sequence_sha256"] for row in rows}),
                "chromosome_counts": {chrom: sum(row["chrom"] == chrom for row in rows) for chrom in sorted(ALLOWED)}}
    if observed != expected:
        raise RuntimeError("Frozen population/component/identity counts differ: " + json.dumps(observed))
    return {"status": "PASS", "model": model, **observed,
            "frozen_columns": list(columns), "all_frozen_cells_identical_to_configuration": True,
            "ordered_frozen_rows_sha256": table_digest(columns, rows),
            "labels_and_roles_unchanged": True, "matching_rebuilt": False,
            "train_derived_covariates_recomputed": False,
            "ambiguous_base_rows": sum(float(row["non_acgt_fraction"]) > 0 for row in rows)}


def unique_intervals(panels):
    by_id = {}
    for rows in panels.values():
        for row in rows:
            projected = {key: row[key] for key in COLS}
            if row["interval_id"] in by_id and by_id[row["interval_id"]] != projected:
                raise RuntimeError("Cross-panel frozen interval metadata conflict")
            by_id[row["interval_id"]] = projected
    return sorted(by_id.values(), key=lambda row: (int(row["chrom"][3:]), int(row["core_start"]), row["interval_id"]))


class TestIndexedFasta:
    """FAI reader with an explicit two-chromosome, exact-interval access allowlist."""
    def __init__(self, path, permitted_intervals):
        self.allowed_intervals = frozenset(permitted_intervals)
        if not self.allowed_intervals or any(chrom not in ALLOWED for chrom, _, _ in self.allowed_intervals):
            raise RuntimeError("FASTA allowlist contains forbidden/non-test input")
        self.handle = path.open("rb")
        self.index = {}
        self.fetches = 0
        for line in Path(str(path) + ".fai").read_text().splitlines():
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[name] = tuple(map(int, (length, offset, bases, width)))

    def fetch(self, chrom, start, end):
        if chrom not in ALLOWED or (chrom, start, end) not in self.allowed_intervals:
            raise RuntimeError("Forbidden chromosome or non-panel sequence fetch")
        key = chrom if chrom in self.index else chrom.removeprefix("chr")
        length, offset, bases, width = self.index[key]
        if not 0 <= start < end <= length:
            raise RuntimeError("Sequence bounds failure")
        first = offset + (start // bases) * width + start % bases
        last = offset + ((end - 1) // bases) * width + (end - 1) % bases
        self.handle.seek(first)
        self.fetches += 1
        return self.handle.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"").upper()

    def close(self):
        self.handle.close()


def verify_sequence(row, raw):
    normalized = raw.translate(NORMALIZE)
    reverse = normalized.translate(COMPLEMENT)[::-1]
    canonical = min(normalized, reverse)
    if (len(raw) != 2001 or hashlib.sha256(raw).hexdigest() != row["sequence_sha256"] or
            hashlib.sha256(canonical).hexdigest() != row["canonical_rc_sequence_sha256"] or
            abs((raw.count(b"G") + raw.count(b"C")) / 2001 - float(row["gc_fraction"])) > 1e-9 or
            abs(normalized.count(b"N") / 2001 - float(row["non_acgt_fraction"])) > 1e-9):
        raise RuntimeError("Frozen sequence/hash/covariate failure: " + row["interval_id"])
    return canonical, int(normalized == canonical)


def verify_preparation_gate(stage):
    integrity_path = stage / "provenance/input_integrity_verified_before.json"
    integrity = json.loads(integrity_path.read_text())
    if integrity.get("status") != "PASS":
        raise RuntimeError("Independent immutable-input gate is not PASS")
    freeze_path = stage / "provenance/prospective_specification_freeze.json"
    freeze_record = integrity.get("prospective_specification_freeze", {})
    if freeze_record.get("path") != "provenance/prospective_specification_freeze.json":
        raise RuntimeError("Independent input gate does not bind the prospective freeze")
    if freeze_path.stat().st_size != freeze_record.get("bytes") or sha(freeze_path) != freeze_record.get("sha256"):
        raise RuntimeError("Prospective freeze hash/size mismatch against independent input gate")
    freeze = json.loads(freeze_path.read_text())
    if freeze.get("status") != "PASS" or freeze.get("before_any_test_model_inference") is not True:
        raise RuntimeError("Prospective test specification freeze is not PASS")
    spec_record = freeze["specification"]
    if spec_record["path"] != "specification/test_specification.json":
        raise RuntimeError("Unexpected prospective specification path")
    path = stage / spec_record["path"]
    if path.stat().st_size != spec_record["bytes"] or sha(path) != spec_record["sha256"]:
        raise RuntimeError("Prospective specification hash/size mismatch")
    return [record(path, stage), record(freeze_path, stage), record(integrity_path, stage)]


def prepare(stage, repo):
    stage, repo = stage.resolve(), repo.resolve()
    if stage != repo / "diseases/COPD/07_gap_closure/internal-test-1.0":
        raise RuntimeError("Only the new internal-test-1.0 destination is allowed")
    if (stage / "inputs").exists():
        raise RuntimeError("Refusing to replace any existing test input preparation")
    started = time.monotonic()
    gates = verify_preparation_gate(stage)
    sources = []
    for relative, (size, digest) in SOURCE_RECORDS.items():
        path = repo / relative
        if path.stat().st_size != size or sha(path) != digest:
            raise RuntimeError("Frozen source hash/size mismatch: " + relative)
        sources.append({"path": relative, "bytes": size, "sha256": digest})
    panels, columns, audits = {}, {}, []
    for model in MODELS:
        columns[model], panels[model] = read_table(repo / PRETRAINING / f"data/evaluation/{model}_common_challenge_panel.tsv.gz", test_only=True)
        config_columns, config_rows = read_table(repo / PRETRAINING / f"data/configurations/V2-C_{model}_interval_manifest.tsv.gz", test_only=True)
        audits.append(validate_panel(model, columns[model], panels[model], config_columns, config_rows))
    unique = unique_intervals(panels)
    encoded_count = len({row["canonical_rc_sequence_sha256"] for row in unique})
    if len(unique) != 26225 or encoded_count != 26225 or len({row["component_id"] for row in unique}) != 8926:
        raise RuntimeError("Frozen union interval/identity/component count disagreement")
    input_dir = stage / "inputs"
    input_dir.mkdir()
    sequence_path = input_dir / "canonical_sequences.npy"
    sequences = np.lib.format.open_memmap(sequence_path, mode="w+", dtype=np.uint8, shape=(encoded_count, 2001))
    permitted = [(row["chrom"], int(row["input_start"]), int(row["input_end"])) for row in unique]
    fasta = TestIndexedFasta(repo / "diseases/COPD/04_modeling/trednet/fasta/hg38.fa", permitted)
    encoded_map = {}
    try:
        for row in unique:
            raw = fasta.fetch(row["chrom"], int(row["input_start"]), int(row["input_end"]))
            canonical, flag = verify_sequence(row, raw)
            identity = row["canonical_rc_sequence_sha256"]
            if identity not in encoded_map:
                encoded_map[identity] = len(encoded_map)
                sequences[encoded_map[identity]] = np.frombuffer(canonical, dtype=np.uint8)
            elif not np.array_equal(sequences[encoded_map[identity]], np.frombuffer(canonical, dtype=np.uint8)):
                raise RuntimeError("Encoded-identity hash collision")
            row["cache_row"] = str(encoded_map[identity])
            row["forward_is_canonical"] = str(flag)
        sequences.flush()
    finally:
        fasta.close()
    if len(encoded_map) != encoded_count or fasta.fetches != len(unique):
        raise RuntimeError("Sequence extraction count disagreement")
    write_table(input_dir / "cache_index.tsv.gz", [*COLS, "cache_row", "forward_is_canonical"], unique)
    mapping = {row["interval_id"]: row for row in unique}
    for model in MODELS:
        mapped_rows = [{**row, "cache_row": mapping[row["interval_id"]]["cache_row"],
                        "forward_is_canonical": mapping[row["interval_id"]]["forward_is_canonical"]} for row in panels[model]]
        output = input_dir / f"common/{model}_test.tsv.gz"
        write_table(output, [*columns[model], "cache_row", "forward_is_canonical"], mapped_rows)
        reloaded_columns, reloaded = read_table(output)
        if reloaded_columns != [*columns[model], "cache_row", "forward_is_canonical"] or table_digest(columns[model], reloaded) != table_digest(columns[model], panels[model]):
            raise RuntimeError("Written frozen panel string/order roundtrip failed")
    audit = {"status": "PASS", "panels": audits, "sequence_intervals": len(unique),
             "encoded_identity_count": encoded_count, "union_components": 8926,
             "sequence_fetches": fasta.fetches, "forbidden_sequence_fetches": 0,
             "sequence_hash_failures": 0, "sequence_eligibility_failures": 0,
             "source_panel_membership_and_order_preserved": True,
             "cross_model_shared_intervals": sum(len(rows) for rows in panels.values()) - len(unique),
             "all_existing_component_ids_preserved": True,
             "arbitrary_or_candidate_sequences_accessed": False,
             "external_benchmark_opened": False, "model_inference_performed": False}
    write_json(input_dir / "test_panel_integrity_audit.json", audit)
    outputs = [record(path, stage) for path in sorted(input_dir.rglob("*")) if path.is_file()]
    manifest = {"status": "PASS", "completed": True, "created_utc": datetime.now(timezone.utc).isoformat(),
                "sequence_intervals": len(unique), "encoded_identity_count": encoded_count,
                "sequence_shape": [encoded_count, 2001], "sequence_dtype": "uint8_ASCII",
                "normalization": "nonACGT to N; canonical=min(normalized_forward,nucleotide_RC)",
                "chromosome_counts": {chrom: sum(row["chrom"] == chrom for row in unique) for chrom in sorted(ALLOWED)},
                "forbidden_sequence_fetches": 0, "missing_sequences": 0, "hash_failures": 0,
                "class_role_counts": [{"configuration": "V2-C", "model": model, "role": "test", **EXPECTED[model]} for model in MODELS],
                "files": outputs, "frozen_sources": sources, "preparation_gates": gates,
                "environment": {"python": platform.python_version(), "executable": sys.executable, "numpy": np.__version__, "platform": platform.platform()},
                "elapsed_seconds": time.monotonic() - started,
                "peak_rss_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "model_inference_performed": False, "historical_stages_modified": False}
    write_json(input_dir / "input_manifest.json", manifest)
    print(json.dumps({key: manifest[key] for key in ("status", "sequence_intervals", "encoded_identity_count", "elapsed_seconds")}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.stage, args.repo)


if __name__ == "__main__":
    main()
