#!/usr/bin/env python3
"""Independent full-table validation via Counter/bisect, not numpy searchsorted.

Reads only donor ATAC sources and frozen construction artifacts. Writes one new
validation record; never refits or changes the frozen production mapping.
"""
from __future__ import annotations

import bisect
import collections
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

VERSION = Path(__file__).resolve().parents[1]
V2 = VERSION.parent
ROOT = V2.parents[2]


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8388608), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    start_time = time.time()
    destination = VERSION / "results/atac_independent_validation.json"
    if destination.exists():
        raise FileExistsError(destination)
    spec = json.loads((VERSION / "specification/atac_normalization.json").read_text())
    train_chromosomes = set(spec["training_chromosomes"])
    assert not train_chromosomes & {"chr7", "chr8", "chr9"}
    canonical = train_chromosomes | {"chr7", "chr8", "chr9"}
    with (V2 / "provenance/COPD-V2-PREFLIGHT_source_inventory.tsv").open() as handle:
        sources = [row for row in csv.DictReader(handle, delimiter="\t") if row["mark"] == "accessibility"]
    counts = {}
    raw_peaks = {}
    fit_chromosome_counts = collections.Counter()
    for source in sources:
        file_id = source["file_accession"]
        counts[file_id] = collections.Counter()
        raw_source = ROOT / source["path"]
        assert sha(raw_source) == source["sha256"]
        with gzip.open(raw_source, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                if chrom not in canonical:
                    continue
                signal = float(fields[6])
                summit = int(fields[1]) + int(fields[9]) if int(fields[9]) >= 0 else (int(fields[1]) + int(fields[2])) // 2
                raw_peaks[f"{file_id}:{line_number}"] = (file_id, chrom, summit - 500, summit + 500, signal)
                if chrom in train_chromosomes:
                    counts[file_id][signal] += 1
                    fit_chromosome_counts[chrom] += 1
    assert fit_chromosome_counts["chr7"] == fit_chromosome_counts["chr8"] == fit_chromosome_counts["chr9"] == 0
    cdf = {}
    for file_id, counter in counts.items():
        values = sorted(counter)
        cumulative = [0]
        for value in values:
            cumulative.append(cumulative[-1] + counter[value])
        cdf[file_id] = (values, cumulative)
    cdf_records = 0
    loaded_signals = collections.defaultdict(set)
    cdf_path = VERSION / "data/atac_training_cdf.tsv.gz"
    with gzip.open(cdf_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            file_id = row["file_accession"]
            signal = float(row["signal_value"])
            values, cumulative = cdf[file_id]
            index = bisect.bisect_left(values, signal)
            assert values[index] == signal
            assert int(row["training_equal_count"]) == counts[file_id][signal]
            assert int(row["training_less_count"]) == cumulative[index]
            assert int(row["training_denominator"]) == cumulative[-1]
            rank = (cumulative[index] + counts[file_id][signal] / 2) / cumulative[-1]
            assert math.isclose(float(row["mapped_mid_ecdf"]), rank, rel_tol=0, abs_tol=1e-15)
            assert signal not in loaded_signals[file_id]
            loaded_signals[file_id].add(signal)
            cdf_records += 1
    assert all(loaded_signals[file_id] == set(counter) for file_id, counter in counts.items())
    sidecar_path = VERSION / "data/atac_normalized_features.tsv.gz"
    original_path = V2 / "data/COPD-V2-PREFLIGHT/interval_features.tsv.gz"
    row_count = 0
    peak_count = 0
    source_partition_counts = collections.Counter()
    with gzip.open(original_path, "rt") as original, gzip.open(sidecar_path, "rt") as normalized:
        left = csv.DictReader(original, delimiter="\t")
        right = csv.DictReader(normalized, delimiter="\t")
        for old, new in itertools.zip_longest(left, right):
            assert old is not None and new is not None
            for field in ("interval_id", "chrom", "partition", "atac_anchor_count", "atac_peak_ids"):
                assert old[field] == new[field]
            row_count += 1
            ids = new["atac_peak_ids"].split(";") if new["atac_peak_ids"] else []
            files = new["atac_peak_file_ids"].split(";") if ids else []
            signals = [float(item) for item in new["atac_peak_signal_values"].split(";")] if ids else []
            ranks = [float(item) for item in new["atac_peak_train_percentiles"].split(";")] if ids else []
            assert len(ids) == len(files) == len(signals) == len(ranks) == int(old["atac_anchor_count"])
            expected_ranks = []
            for peak_id, file_id, signal, actual in zip(ids, files, signals, ranks):
                expected_raw = raw_peaks.pop(peak_id)
                assert expected_raw == (file_id, old["chrom"], int(old["core_start"]), int(old["core_end"]), signal)
                values, cumulative = cdf[file_id]
                lower = bisect.bisect_left(values, signal)
                less = cumulative[lower]
                equal = counts[file_id].get(signal, 0)
                expected = (less + equal / 2) / cumulative[-1]
                assert math.isclose(actual, expected, rel_tol=0, abs_tol=1e-15)
                assert 0 <= actual <= 1
                expected_ranks.append(expected)
                peak_count += 1
                source_partition_counts[(file_id, old["partition"])] += 1
            if ids:
                assert math.isclose(float(new["atac_signal_percentile_max_train_only"]), max(expected_ranks), rel_tol=0, abs_tol=1e-15)
            else:
                assert new["atac_signal_percentile_max_train_only"] == ""
    assert not raw_peaks
    result = {
        "status": "PASS", "version": "pretraining-1.1", "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "method": "Independent Python Counter/bisect reconstruction from original donor ATAC narrowPeak rows; full frozen-CDF and full normalized-sidecar audit",
        "interval_rows_validated": row_count, "source_peak_rows_validated": peak_count,
        "unique_training_signal_rows_validated": cdf_records,
        "fit_chromosome_peak_counts": dict(sorted(fit_chromosome_counts.items())),
        "heldout_rows_in_fit": 0, "mapping_denominators": {key: sum(value.values()) for key, value in counts.items()},
        "application_partition_counts": [{"file_accession": key[0], "partition": key[1], "peak_rows": value} for key, value in sorted(source_partition_counts.items())],
        "original_interval_ids_order_and_partitions_preserved": True,
        "original_raw_source_hashes_verified": True, "every_normalized_peak_reconciled_to_raw_source": True,
        "every_interval_max_aggregation_verified": True, "no_duplicate_or_omitted_source_peaks": True,
        "all_unanchored_missing_values_preserved": True, "cdf_sha256": sha(cdf_path),
        "sidecar_sha256": sha(sidecar_path), "validator_sha256": sha(Path(__file__)),
        "wall_seconds": time.time() - start_time, "model_or_benchmark_information_read": False,
    }
    with destination.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "interval_rows": row_count, "source_peak_rows": peak_count, "wall_seconds": result["wall_seconds"]}))


if __name__ == "__main__":
    main()
