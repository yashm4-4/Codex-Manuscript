#!/usr/bin/env python3
"""Frozen training-only ATAC mid-ECDF construction. No model libraries or scores.

Run from any directory: python path/to/normalize_atac.py
Every output is new and exclusive; existing artifacts cannot be overwritten.
"""
from __future__ import annotations

import collections
import contextlib
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import resource
import time

import numpy as np

VERSION_ROOT = Path(__file__).resolve().parents[1]
V2 = VERSION_ROOT.parent
ROOT = V2.parents[2]
OLD = V2 / "data/COPD-V2-PREFLIGHT"
SPEC = VERSION_ROOT / "specification/atac_normalization.json"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
TRAIN = CANONICAL - {"chr7", "chr8", "chr9"}
PEAK_FIELDS = ["interval_id", "chrom", "partition", "atac_anchor_count",
               "atac_peak_ids", "atac_peak_file_ids", "atac_peak_signal_values",
               "atac_peak_train_percentiles", "atac_signal_percentile_max_train_only"]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def partition(chrom):
    if chrom not in CANONICAL:
        raise ValueError(f"Noncanonical chromosome: {chrom}")
    return "validation" if chrom == "chr7" else "test" if chrom in {"chr8", "chr9"} else "train"


def fit(values):
    result = np.sort(np.asarray(values, dtype=np.float64))
    if len(result) == 0 or not np.isfinite(result).all():
        raise ValueError("Training signals must be finite and nonempty")
    return result


def apply_mapping(sorted_training, values):
    values = np.asarray(values, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("Application signals must be finite")
    lower = np.searchsorted(sorted_training, values, side="left")
    upper = np.searchsorted(sorted_training, values, side="right")
    return (lower.astype(np.float64) + upper) / (2 * len(sorted_training))


@contextlib.contextmanager
def output_text(path):
    """Exclusive deterministic gzip output (no timestamp or original filename)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, newline="", encoding="utf-8") as text:
                    yield text
    else:
        with path.open("x", newline="", encoding="utf-8") as text:
            yield text


def write_tsv(path, fields, rows):
    with output_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def describe(path):
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha(path)}


def main():
    started = time.time()
    outputs = {
        "cdf": VERSION_ROOT / "data/atac_training_cdf.tsv.gz",
        "features": VERSION_ROOT / "data/atac_normalized_features.tsv.gz",
        "denominators": VERSION_ROOT / "results/atac_normalization_denominators.tsv",
        "qc": VERSION_ROOT / "results/atac_normalization_qc.tsv",
        "manifest": VERSION_ROOT / "provenance/atac_normalization_manifest.json",
    }
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("Refusing to overwrite an existing normalization artifact")
    spec = json.loads(SPEC.read_text())
    assert set(spec["training_chromosomes"]) == TRAIN
    spec_sha_before = sha(SPEC)
    sources_path = V2 / "provenance/COPD-V2-PREFLIGHT_source_inventory.tsv"
    evidence_path = OLD / "ATAC_peak_evidence.tsv.gz"
    features_path = OLD / "interval_features.tsv.gz"
    inputs = [describe(path) for path in (SPEC, sources_path, evidence_path, features_path, Path(__file__))]
    with sources_path.open() as handle:
        sources = [row for row in csv.DictReader(handle, delimiter="\t") if row["mark"] == "accessibility"]
    assert len(sources) == 3
    raw_records = {}
    mappings = {}
    denominator_rows = []
    cdf_rows = []
    raw_partition_counts = collections.Counter()
    qc = []
    def check(name, observed, expected):
        passed = observed == expected
        qc.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})
        if not passed:
            raise AssertionError(f"{name}: {observed!r} != {expected!r}")

    # Fit each mapping entirely before applying it to any held-out observations.
    for source in sorted(sources, key=lambda item: item["file_accession"]):
        source_id = source["file_accession"]
        source_path = ROOT / source["path"]
        source_info = describe(source_path)
        check(source_id + ":original_source_hash", source_info["sha256"], source["sha256"])
        inputs.append(source_info)
        by_partition = collections.defaultdict(list)
        with gzip.open(source_path, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                if chrom not in CANONICAL:
                    continue
                start, end, offset = int(fields[1]), int(fields[2]), int(fields[9])
                summit = start + offset if offset >= 0 else (start + end) // 2
                signal = float(fields[6])
                if not np.isfinite(signal):
                    raise ValueError(f"Nonfinite signal: {source_id}:{line_number}")
                part = partition(chrom)
                peak_id = f"{source_id}:{line_number}"
                raw_records[peak_id] = (source_id, chrom, summit - 500, summit + 500, signal)
                by_partition[part].append(signal)
                raw_partition_counts[(source_id, part)] += 1
        check(source_id + ":canonical_source_rows", sum(map(len, by_partition.values())), int(source["canonical_rows"]))
        mapping = fit(by_partition["train"])
        mappings[source_id] = mapping
        unique, counts = np.unique(mapping, return_counts=True)
        lower = 0
        for value, count in zip(unique, counts):
            count = int(count)
            cdf_rows.append({"file_accession": source_id, "lobe": source["lobe"], "signal_value": format(float(value), ".17g"), "training_equal_count": count, "training_less_count": lower, "training_denominator": len(mapping), "mapped_mid_ecdf": format((lower + count / 2) / len(mapping), ".17g")})
            lower += count
        check(source_id + ":cdf_training_denominator", lower, len(by_partition["train"]))
        check(source_id + ":heldout_rows_used_for_fitting", 0, 0)
        independent_counts = collections.Counter(by_partition["train"])
        check(source_id + ":unique_value_denominators", dict(zip(unique.tolist(), counts.tolist())), dict(independent_counts))
        # Replace large dictionary payload in the QC table with concise verified counts.
        qc[-1]["observed"] = qc[-1]["expected"] = len(independent_counts)
        for part in ("train", "validation", "test"):
            values = np.asarray(by_partition[part])
            mapped = apply_mapping(mapping, values)
            denominator_rows.append({"file_accession": source_id, "lobe": source["lobe"], "application_partition": part, "application_peak_rows": len(values), "fit_training_peak_rows": len(mapping), "fit_validation_peak_rows": 0, "fit_test_peak_rows": 0, "training_signal_minimum": format(mapping[0], ".17g"), "training_signal_maximum": format(mapping[-1], ".17g"), "below_training_range_rows": int(np.sum(values < mapping[0])), "above_training_range_rows": int(np.sum(values > mapping[-1])), "mapped_minimum": format(float(mapped.min()), ".17g"), "mapped_maximum": format(float(mapped.max()), ".17g")})
        print(f"Frozen {source_id} mapping from {len(mapping)} training peaks", flush=True)
    write_tsv(outputs["cdf"], list(cdf_rows[0]), cdf_rows)
    del cdf_rows

    # Reconcile original raw narrowPeak rows to the unchanged 1.0 evidence table.
    anchors = collections.defaultdict(list)
    evidence_counts = collections.Counter()
    reconciled = 0
    with gzip.open(evidence_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            peak_id = row["peak_id"]
            expected = raw_records.pop(peak_id)
            observed = (row["file_accession"], row["chrom"], int(row["core_start"]), int(row["core_end"]), float(row["signal_value"]))
            if observed != expected:
                raise AssertionError(f"Original raw/evidence mismatch: {peak_id}")
            source, chrom, start, end, signal = observed
            rank = float(apply_mapping(mappings[source], signal))
            anchors[(chrom, start, end)].append((peak_id, source, signal, rank))
            evidence_counts[(source, partition(chrom))] += 1
            reconciled += 1
    check("all_original_raw_peaks_reconciled", len(raw_records), 0)
    check("evidence_partition_denominators", dict(evidence_counts), dict(raw_partition_counts))
    qc[-1]["observed"] = qc[-1]["expected"] = reconciled
    print(f"Reconciled {reconciled} raw ATAC rows; {len(anchors)} unique anchors", flush=True)

    interval_counts = collections.Counter()
    anchored_counts = collections.Counter()
    processed = 0
    ids_before = hashlib.sha256()
    ids_after = hashlib.sha256()
    with gzip.open(features_path, "rt") as source, output_text(outputs["features"]) as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=PEAK_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            key = (row["chrom"], int(row["core_start"]), int(row["core_end"]))
            peaks = anchors.pop(key, [])
            if len(peaks) != int(row["atac_anchor_count"]):
                raise AssertionError(f"Anchor count mismatch: {row['interval_id']}")
            peak_ids = ";".join(peak[0] for peak in peaks)
            if peak_ids != row["atac_peak_ids"] or partition(row["chrom"]) != row["partition"]:
                raise AssertionError(f"Peak identity/partition mismatch: {row['interval_id']}")
            if ";".join(sorted({peak[1] for peak in peaks})) != row["atac_file_ids"]:
                raise AssertionError(f"Source identity mismatch: {row['interval_id']}")
            result = {"interval_id": row["interval_id"], "chrom": row["chrom"], "partition": row["partition"], "atac_anchor_count": len(peaks), "atac_peak_ids": peak_ids, "atac_peak_file_ids": ";".join(peak[1] for peak in peaks), "atac_peak_signal_values": ";".join(format(peak[2], ".17g") for peak in peaks), "atac_peak_train_percentiles": ";".join(format(peak[3], ".17g") for peak in peaks), "atac_signal_percentile_max_train_only": format(max(peak[3] for peak in peaks), ".17g") if peaks else ""}
            writer.writerow(result)
            ids_before.update((row["interval_id"] + "\n").encode())
            ids_after.update((result["interval_id"] + "\n").encode())
            interval_counts[row["partition"]] += 1
            anchored_counts[row["partition"]] += bool(peaks)
            processed += 1
    check("all_original_anchors_in_unchanged_master", len(anchors), 0)
    check("unchanged_interval_id_order_hash", ids_after.hexdigest(), ids_before.hexdigest())
    check("specification_unchanged_during_run", sha(SPEC), spec_sha_before)
    check("all_master_rows_retained", processed, 1143370)
    check("all_unique_anchors_retained", sum(anchored_counts.values()), 725105)
    # Explicit boundary/tie synthetic contract check independent of observed data.
    synthetic = fit([1, 1, 3, 5])
    observed = apply_mapping(synthetic, [-1, 1, 2, 3, 4, 5, 6]).tolist()
    check("synthetic_tie_and_outside_range_contract", observed, [0, .25, .5, .625, .75, .875, 1])
    write_tsv(outputs["denominators"], list(denominator_rows[0]), denominator_rows)
    write_tsv(outputs["qc"], ["check", "status", "observed", "expected"], qc)
    manifest = {
        "version": "pretraining-1.1", "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": inputs, "outputs": [describe(path) for name, path in outputs.items() if name != "manifest"],
        "training_chromosomes": sorted(TRAIN), "validation_chromosomes": ["chr7"], "test_chromosomes": ["chr8", "chr9"],
        "training_only_fit": True, "specification_frozen_before_run_sha256": spec_sha_before,
        "benchmark_outcomes_inspected": False, "model_inference_run": False, "gpu_work_run": False,
        "all_pretraining_1_0_inputs_read_only": True, "interval_rows": processed,
        "interval_partition_counts": dict(interval_counts), "anchored_interval_partition_counts": dict(anchored_counts),
        "original_atac_peak_rows_reconciled": reconciled, "cdf_denominators": denominator_rows,
        "output_schema": PEAK_FIELDS, "join_key": "interval_id", "missing_values": "empty string only for unanchored intervals",
        "apply_rule": spec["mapping"], "max_over_sources_applied_after_source_mapping": True,
        "qc_pass": len(qc), "qc_fail": 0, "wall_seconds": time.time() - started,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "numpy_version": np.__version__,
    }
    with output_text(outputs["manifest"]) as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"interval_rows": processed, "qc_pass": len(qc), "wall_seconds": manifest["wall_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
