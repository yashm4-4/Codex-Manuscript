#!/usr/bin/env python3
"""Independent CPU-only audit of the frozen COPD model-redesign preflight.

This validator never imports a training implementation, TensorFlow, or Keras.
External benchmark artifacts may only be read as opaque bytes for SHA-256.
Artifact correctness and training readiness are distinct: a correctly reported
data limitation can pass artifact QC while blocking authorization to train.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


ROOT = Path(__file__).resolve().parents[4]
V2 = ROOT / "diseases/COPD/07_gap_closure"
DATA = V2 / "data/COPD-V2-PREFLIGHT"
PROV = V2 / "provenance"
RESULTS = V2 / "results"
PREFIX = "COPD-V2-PREFLIGHT"
SEEDS = [104729, 130363, 155921]
MODELS = ("enhancer", "h3k27me3")
PARTITIONS = ("train", "validation", "test")
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
CHECKS = []
READINESS = []


class AggregatedChecks:
    def __init__(self, prefix):
        self.prefix = prefix
        self.n = collections.Counter()
        self.failures = collections.Counter()
        self.examples = collections.defaultdict(list)

    def add(self, name, condition, example=""):
        self.n[name] += 1
        if not condition:
            self.failures[name] += 1
            if len(self.examples[name]) < 5:
                self.examples[name].append(str(example))

    def finish(self):
        for name in sorted(self.n):
            bad = self.failures[name]
            check(f"{self.prefix}:{name}", bad == 0,
                  json.dumps({"tested": self.n[name], "failures": bad,
                              "examples": self.examples[name]}), "zero failures")


class IndependentIntervals:
    """Union intervals independently; query first interval ending after start."""
    def __init__(self, rows):
        groups = collections.defaultdict(list)
        for chrom, start, end in rows:
            if end > start:
                groups[chrom].append((int(start), int(end)))
        self.data = {}
        for chrom, intervals in groups.items():
            merged = []
            for start, end in sorted(intervals):
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(end, merged[-1][1])
            self.data[chrom] = ([r[0] for r in merged], [r[1] for r in merged])

    def overlap(self, chrom, start, end):
        starts, ends = self.data.get(chrom, ([], []))
        candidate = bisect.bisect_right(ends, start)
        return candidate < len(starts) and starts[candidate] < end

    def covered(self, chrom, start, end):
        starts, ends = self.data.get(chrom, ([], []))
        index = bisect.bisect_right(ends, start)
        amount = 0
        while index < len(starts) and starts[index] < end:
            amount += max(0, min(end, ends[index]) - max(start, starts[index]))
            index += 1
        return amount


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check(name, condition, observed="", expected="", domain="artifact"):
    record = {"check": name, "status": "PASS" if condition else "FAIL",
              "observed": str(observed), "expected": str(expected), "domain": domain}
    CHECKS.append(record)
    return bool(condition)


def readiness(name, condition, observed="", required=""):
    READINESS.append({"condition": name, "status": "SATISFIED" if condition else "BLOCKED",
                      "observed": str(observed), "required": str(required)})


def open_table(path):
    if "COPD-V2-BENCH" in str(path):
        raise RuntimeError("Benchmark content parsing is forbidden")
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()


def table(path):
    with open_table(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_json(path):
    if "COPD-V2-BENCH" in str(path):
        raise RuntimeError("Benchmark JSON parsing is forbidden")
    with Path(path).open() as handle:
        return json.load(handle)


def partition(chrom):
    if chrom not in CANONICAL:
        return "excluded"
    if chrom == "chr7":
        return "validation"
    if chrom in {"chr8", "chr9"}:
        return "test"
    return "train"


def ids_from_bed(path):
    result = []
    with open_table(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                chrom, start, end, *_ = line.rstrip().split("\t")
                result.append(f"{chrom}:{int(start)}-{int(end)}")
    return result


def verify_ledger(path, label):
    count = 0
    for entry in table(path):
        target = ROOT / entry["path"]
        count += 1
        exists = target.is_file()
        check(f"{label}:exists:{entry['path']}", exists)
        if not exists:
            continue
        check(f"{label}:size:{entry['path']}", target.stat().st_size == int(entry["bytes"]),
              target.stat().st_size, entry["bytes"])
        digest = sha256(target)
        check(f"{label}:sha256:{entry['path']}", digest == entry["sha256"], digest, entry["sha256"])
        if "COPD-V2-BENCH" in str(target):
            check(f"{label}:benchmark_opaque_access:{entry['path']}",
                  entry.get("access_mode") == "opaque_bytes_only", entry.get("access_mode"))
    return count


def verify_manifest_files(manifest, label):
    count = 0
    for entry in manifest.get("files", []):
        target = ROOT / entry["path"]
        count += 1
        exists = target.is_file()
        check(f"{label}:exists:{entry['path']}", exists)
        if not exists:
            continue
        check(f"{label}:size:{entry['path']}", target.stat().st_size == int(entry["bytes"]),
              target.stat().st_size, entry["bytes"])
        digest = sha256(target)
        check(f"{label}:sha256:{entry['path']}", digest == entry["sha256"], digest, entry["sha256"])
    return count


def verify_orientation_contract():
    contract = read_json(PROV / f"{PREFIX}_orientation_contract_qc.json")
    tests = contract.get("checks", [])
    check("orientation:synthetic_checks_present", len(tests) >= 30, len(tests), ">=30")
    for test in tests:
        check("orientation:" + test["check"], test["status"] == "PASS", test["status"], "PASS")
    actual = contract.get("actual_trained_network_invariance_validation", "")
    check("orientation:real_network_test_truthfully_not_run", actual.startswith("NOT_RUN"), actual, "NOT_RUN")
    check("orientation:relative_tolerance", contract.get("relative_tolerance") == 1e-6,
          contract.get("relative_tolerance"), 1e-6)
    check("orientation:absolute_tolerance", contract.get("absolute_tolerance") == 1e-6,
          contract.get("absolute_tolerance"), 1e-6)
    check("orientation:three_fixed_seeds", contract.get("training_seeds") == SEEDS)
    check("orientation:no_training_model_loading", contract.get("training_or_model_loading_performed") is False)
    path = V2 / "scripts/preflight_orientation_contract.py"
    spec = importlib.util.spec_from_file_location("preflight_cpu_contract_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    repeated = module.run_synthetic_qc()
    check("orientation:CPU_contract_reexecution_passes", all(r["status"] == "PASS" for r in repeated["checks"]))
    check("orientation:CPU_contract_check_set_reproduced", {r["check"] for r in repeated["checks"]} == {r["check"] for r in tests})
    return contract


def require_complete_features():
    path = PROV / f"{PREFIX}_feature_manifest.json"
    if not path.is_file():
        raise RuntimeError("Feature construction has not completed; no validator outputs written")
    manifest = read_json(path)
    check("features:CPU_only_status", manifest.get("status") == "COMPLETE_CPU_ONLY_NO_TRAINING",
          manifest.get("status"))
    check("features:no_training_or_inference", manifest.get("training_or_inference") is False)
    check("features:no_benchmark_content_access", manifest.get("benchmark_contents_inspected") is False)
    verify_manifest_files(manifest, "features")
    return manifest


def original_bed_sets():
    membership = {}
    for model, word in (("enhancer", "Enhancer"), ("h3k27me3", "Silencer")):
        for label in ("positive", "control"):
            path = ROOT / ("diseases/COPD/04_modeling/trednet/input_training_data/"
                           f"COPD_SevereEmphysema_Lung_{word}_DHS_x2_{label}_1kb.bed")
            rows = ids_from_bed(path)
            membership[(model, label)] = set(rows)
            check(f"original_BED:unique:{model}:{label}", len(rows) == len(set(rows)),
                  len(rows), len(set(rows)))
    return membership


def id_coordinate(identifier):
    chrom, position = identifier.split(":")
    start, end = position.split("-")
    return chrom, int(start), int(end)


def original_source_evidence():
    """Reconstruct anatomical co-occurrence from original nine peak files."""
    sources = list(table(ROOT / "diseases/COPD/03_regulatory_landscape/data/encode_peaks/manifest.tsv"))
    check("sources:nine_files", len(sources) == 9, len(sources), 9)
    check("sources:one_donor", {r["donor_accessions"] for r in sources} == {"ENCDO520EJG"})
    check("sources:GRCh38", {r["assembly"] for r in sources} == {"GRCh38"})
    lobes = ("lower_left", "lower_right", "upper_right")
    marks = {"enhancer": "H3K27ac", "h3k27me3": "H3K27me3"}
    by_key = {(r["lobe"], r["mark"]): r for r in sources}
    check("sources:complete_three_lobe_three_assay_design",
          set(by_key) == {(l, m) for l in lobes for m in ("accessibility", "H3K27ac", "H3K27me3")})
    raw = {}
    for key, source in by_key.items():
        peaks = []
        with gzip.open(source["path"], "rt") as handle:
            for lineno, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip().split("\t")
                if fields[0] not in CANONICAL:
                    continue
                start, end = int(fields[1]), int(fields[2])
                offset = int(fields[9])
                summit = start + offset if offset >= 0 else (start + end) // 2
                peaks.append((fields[0], start, end, summit, float(fields[6]),
                              f"{source['file_accession']}:{lineno}"))
        raw[key] = peaks
    mark_indices = {(l, m): IndependentIntervals((r[0], r[1], r[2]) for r in raw[(l, mark)])
                    for l in lobes for m, mark in marks.items()}
    union_mark = {m: IndependentIntervals((r[0], r[1], r[2])
                                         for l in lobes for r in raw[(l, mark)])
                  for m, mark in marks.items()}
    union_atac = IndependentIntervals((r[0], r[1], r[2])
                                     for l in lobes for r in raw[(l, "accessibility")])
    anchors = {}
    peak_rules = AggregatedChecks("ATAC_source_reconstruction")
    # Compact expected anchor structure: peak IDs, files, lobes, signal/rank maxima,
    # and per-model same-lobe/all-lobe support sets.
    for lobe in lobes:
        source = by_key[(lobe, "accessibility")]
        peaks = raw[(lobe, "accessibility")]
        signals = sorted(r[4] for r in peaks)
        for chrom, start, end, summit, signal, pid in peaks:
            identifier = f"{chrom}:{summit-500}-{summit+500}"
            anchor = anchors.setdefault(identifier, {"peaks": set(), "files": set(), "lobes": set(),
                                                     "signal": float("-inf"), "rank": 0.0,
                                                     "enhancer_all": set(), "enhancer_same": set(),
                                                     "h3k27me3_all": set(), "h3k27me3_same": set()})
            anchor["peaks"].add(pid)
            anchor["files"].add(source["file_accession"])
            anchor["lobes"].add(lobe)
            anchor["signal"] = max(signal, anchor["signal"])
            rank = (bisect.bisect_left(signals, signal) + bisect.bisect_right(signals, signal) + 1) / (2 * len(signals))
            anchor["rank"] = max(rank, anchor["rank"])
            for model in MODELS:
                supported = {other for other in lobes if mark_indices[(other, model)].overlap(chrom, start, end)}
                anchor[f"{model}_all"].update(supported)
                if lobe in supported:
                    anchor[f"{model}_same"].add(lobe)
            peak_rules.add("summit_anchor_geometry", summit - 500 < summit + 500, pid)
    peak_rules.finish()
    return anchors, union_mark, union_atac, by_key


def independent_annotation_indices():
    blacklist = []
    with gzip.open(ROOT / "data/blacklist/hg38-blacklist.v2.bed.gz", "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                fields = line.split("\t")
                blacklist.append((fields[0], int(fields[1]), int(fields[2])))
    promoters = []
    with gzip.open(ROOT / "data/gencode/gencode.v50.annotation.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] in CANONICAL and fields[2] == "gene":
                tss = int(fields[3] if fields[6] == "+" else fields[4]) - 1
                promoters.append((fields[0], max(0, tss-2000), tss+2000))
    def repeats():
        with gzip.open(ROOT / "data/repeats/rmsk.txt.gz", "rt") as handle:
            for line in handle:
                fields = line.split("\t")
                if len(fields) >= 8 and fields[5] in CANONICAL:
                    yield fields[5], int(fields[6]), int(fields[7])
    return IndependentIntervals(blacklist), IndependentIntervals(promoters), IndependentIntervals(repeats())


def verify_features(manifest, membership, anchors, mark_indices, atac_index, source_by_key):
    """Validate every source/provenance row; retain compact records for set QC."""
    import pysam

    rules = AggregatedChecks("feature_rows")
    blacklist_index, promoter_index, repeat_index = independent_annotation_indices()
    features = {}
    feature_counts = collections.Counter()
    duplicate_ids = {"exact_sequence": collections.defaultdict(list),
                     "canonical_forward_RC": collections.defaultdict(list)}
    positive_index = {m: IndependentIntervals((c, s-501, e+500)
                                             for c, s, e in map(id_coordinate, membership[(m, "positive")]))
                      for m in MODELS}
    fasta = ROOT / "diseases/COPD/04_modeling/trednet/fasta/hg38.fa"
    reference = pysam.FastaFile(str(fasta))
    lengths = dict(zip(reference.references, reference.lengths))
    complement = str.maketrans("ACGTN", "TGCAN")
    marks = {"enhancer": "H3K27ac", "h3k27me3": "H3K27me3"}
    observed_membership = {(m, label): set() for m in MODELS for label in ("positive", "control")}
    for index, row in enumerate(table(DATA / "interval_features.tsv.gz"), 1):
        identifier = row["interval_id"]
        chrom, start, end = row["chrom"], int(row["core_start"]), int(row["core_end"])
        left, right = int(row["input_start"]), int(row["input_end"])
        part = partition(chrom)
        rules.add("unique_interval_identifier", identifier not in features, identifier)
        rules.add("identifier_coordinates", identifier == f"{chrom}:{start}-{end}", identifier)
        rules.add("partition", row["partition"] == part and part != "excluded", identifier)
        rules.add("window_geometry", end-start == 1000 and left == start-501 and right == end+500 and right-left == 2001, identifier)
        for name, annotation in (("blacklist", blacklist_index), ("promoter", promoter_index)):
            for suffix, a, b in (("1kb", start, end), ("2001", left, right)):
                field = f"{name}_bp_{suffix}"
                rules.add(field, int(row[field]) == annotation.covered(chrom, a, b), identifier)
        rules.add("RepeatMasker_full_input_fraction", abs(float(row["repeat_fraction_2001"]) - repeat_index.covered(chrom, left, right)/2001) < 1e-9, identifier)
        anchor = anchors.get(identifier)
        rules.add("ATAC_anchor_count", int(row["atac_anchor_count"]) == (len(anchor["peaks"]) if anchor else 0), identifier)
        tokens = lambda value: set(value.split(";")) if value else set()
        for source_field, anchor_field in (("atac_peak_ids", "peaks"), ("atac_file_ids", "files"), ("atac_lobes", "lobes")):
            rules.add(source_field, tokens(row[source_field]) == (anchor[anchor_field] if anchor else set()), identifier)
        if anchor:
            rules.add("ATAC_signal", abs(float(row["atac_signal_raw_max"]) - anchor["signal"]) < 1e-8, identifier)
            rules.add("ATAC_source_normalized_rank", abs(float(row["atac_signal_percentile_max"]) - anchor["rank"]) < 1e-9, identifier)
        for field, query_start, query_end in (("atac_overlap_1kb", start, end), ("atac_overlap_2001", left, right), ("atac_overlap_summit", start+500, start+501)):
            rules.add(field, int(row[field]) == int(atac_index.overlap(chrom, query_start, query_end)), identifier)
        for model in MODELS:
            for label in ("positive", "control"):
                expected = identifier in membership[(model, label)]
                rules.add(f"V1_{model}_{label}_membership", int(row[f"v1_{model}_{label}"]) == int(expected), identifier)
                if expected:
                    observed_membership[(model, label)].add(identifier)
                    feature_counts[(model, part, f"V1_{label}")] += 1
            for support, key in (("all_lobe_peak_support", "all"), ("same_lobe_peak_support", "same")):
                support_lobes = anchor[f"{model}_{key}"] if anchor else set()
                rules.add(f"{model}_{support}", int(row[f"{model}_{support}"]) == int(bool(support_lobes)), identifier)
                field = f"{model}_{'all_lobe' if key == 'all' else 'same_lobe'}_histone_file_ids"
                expected_files = {source_by_key[(l, marks[model])]["file_accession"] for l in support_lobes}
                rules.add(field, tokens(row[field]) == expected_files, identifier)
            same_lobes = anchor[f"{model}_same"] if anchor else set()
            rules.add(f"{model}_same_lobe_support_lobes", tokens(row[f"{model}_same_lobe_support_lobes"]) == same_lobes, identifier)
            for field, a, b in (("mark_overlap_1kb", start, end), ("mark_overlap_2001", left, right)):
                rules.add(f"{model}_{field}", int(row[f"{model}_{field}"]) == int(mark_indices[model].overlap(chrom, a, b)), identifier)
            rules.add(f"{model}_positive_input_overlap", int(row[f"{model}_v1_positive_input_overlap"]) == int(positive_index[model].overlap(chrom, left, right)), identifier)
            if identifier in membership[(model, "positive")]:
                feature_counts[(model, part, "V1_positive_same_lobe" if same_lobes else "V1_positive_cross_lobe_only")] += 1
        reference_chrom = chrom if chrom in lengths else chrom.removeprefix("chr")
        expected_available = reference_chrom in lengths and left >= 0 and right <= lengths.get(reference_chrom, 0)
        rules.add("sequence_available", int(row["sequence_available"]) == int(expected_available), identifier)
        if expected_available:
            sequence = reference.fetch(reference_chrom, left, right).upper()
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            acgt_bases = sum(sequence.count(base) for base in "ACGT")
            normalized = sequence if acgt_bases == 2001 else "".join(base if base in "ACGT" else "N" for base in sequence)
            canonical = min(normalized, normalized.translate(complement)[::-1])
            canonical_digest = hashlib.sha256(canonical.encode()).hexdigest()
            rules.add("sequence_length", int(row["sequence_length"]) == len(sequence) == 2001, identifier)
            rules.add("sequence_SHA256", row["sequence_sha256"] == digest, identifier)
            rules.add("encoded_RC_canonical_SHA256", row["canonical_rc_sequence_sha256"] == canonical_digest, identifier)
            rules.add("GC_fraction", abs(float(row["gc_fraction"]) - (sequence.count("G")+sequence.count("C"))/2001) < 1e-9, identifier)
            rules.add("ambiguous_base_fraction", abs(float(row["non_acgt_fraction"]) - normalized.count("N")/2001) < 1e-9, identifier)
            duplicate_ids["exact_sequence"][digest].append(identifier)
            duplicate_ids["canonical_forward_RC"][canonical_digest].append(identifier)
        else:
            rules.add("unavailable_sequence_no_fake_hash", not row["sequence_sha256"] and not row["canonical_rc_sequence_sha256"], identifier)
        for model in MODELS:
            eligible = bool(anchor) and expected_available and start >= 501 and int(row["blacklist_bp_1kb"]) == 0 and int(row["promoter_bp_1kb"]) == 0
            eligible = eligible and not anchor[f"{model}_all"] and not mark_indices[model].overlap(chrom, left, right) and not positive_index[model].overlap(chrom, left, right) if anchor else False
            rules.add(f"{model}_control_eligibility", int(row[f"{model}_accessible_control_eligible"]) == int(eligible), identifier)
            if eligible:
                feature_counts[(model, part, "accessible_control_eligible")] += 1
            if not expected_available and (identifier in membership[(model, "positive")] or identifier in membership[(model, "control")]):
                feature_counts[(model, part, "V1_sequence_unavailable")] += 1
        feature_counts[("all", part, "union_intervals")] += 1
        feature_counts[("all", part, "sequence_available" if expected_available else "sequence_unavailable")] += 1
        features[identifier] = row
        if index % 100000 == 0:
            print(f"Independently verified {index} source/sequence feature rows", flush=True)
    reference.close()
    rules.finish()
    check("features:manifest_row_count", len(features) == int(manifest["interval_count"]), len(features), manifest["interval_count"])
    check("features:all_anchors_and_V1_inputs_exact_union", set(features) == set(anchors).union(*membership.values()))
    for key, expected in membership.items():
        check("features:complete_original_membership:" + ":".join(key), observed_membership[key] == expected,
              len(observed_membership[key]), len(expected))
    reported_counts = {(r["model"], r["partition"], r["metric"]): int(r["count"])
                       for r in table(RESULTS / f"{PREFIX}_feature_construction_counts.tsv")}
    check("features:all_construction_counts_reconciled", dict(feature_counts) == reported_counts)
    expected_duplicates = {(kind, digest): ids for kind, groups in duplicate_ids.items()
                           for digest, ids in groups.items() if len(ids) > 1}
    reported_duplicates = {(r["identity_type"], r["sha256"]): r
                           for r in table(RESULTS / f"{PREFIX}_sequence_duplicate_audit.tsv")}
    check("duplicates:complete_identity_group_keys", set(expected_duplicates) == set(reported_duplicates),
          len(reported_duplicates), len(expected_duplicates))
    duplicate_rules = AggregatedChecks("duplicates")
    for key, identifiers in expected_duplicates.items():
        if key not in reported_duplicates:
            continue
        report = reported_duplicates[key]
        parts = {features[i]["partition"] for i in identifiers}
        duplicate_rules.add("member_identity", set(report["interval_ids"].split(";")) == set(identifiers), key)
        duplicate_rules.add("member_count", int(report["n_intervals"]) == len(identifiers), key)
        duplicate_rules.add("partition_count", int(report["n_partitions"]) == len(parts), key)
        duplicate_rules.add("cross_partition_status", int(report["cross_partition"]) == int(len(parts) > 1), key)
    duplicate_rules.finish()
    return features, expected_duplicates


def verify_roles(features, duplicate_groups):
    """Rebuild overlap components and sequence-equivalence unions independently."""
    ordered = sorted(features, key=lambda i: (features[i]["chrom"], int(features[i]["input_start"]), int(features[i]["input_end"])))
    initial = {}
    parent = {}
    rank = {}
    last_chrom, last_end, component = None, -1, None
    for identifier in ordered:
        row = features[identifier]
        chrom, start, end = row["chrom"], int(row["input_start"]), int(row["input_end"])
        if chrom != last_chrom or start >= last_end:
            component = f"{chrom}:{start}"
            parent[component] = component
            rank[component] = len(rank)
            last_end = end
        else:
            last_end = max(last_end, end)
        initial[identifier] = component
        last_chrom = chrom

    def find(component_id):
        trail = []
        while component_id != parent[component_id]:
            trail.append(component_id)
            component_id = parent[component_id]
        for old in trail:
            parent[old] = component_id
        return component_id

    for (kind, _), identifiers in duplicate_groups.items():
        if kind != "canonical_forward_RC":
            continue
        representative = find(initial[identifiers[0]])
        for identifier in identifiers[1:]:
            other = find(initial[identifier])
            if other != representative:
                keep, remove = sorted((representative, other), key=lambda c: rank[c])
                parent[remove] = keep
                representative = keep
    expected = {identifier: find(group) for identifier, group in initial.items()}
    roles = {}
    checks = AggregatedChecks("roles")
    for row in table(DATA / "interval_role_assignment.tsv.gz"):
        identifier = row["interval_id"]
        checks.add("known_unique_interval", identifier in features and identifier not in roles, identifier)
        if identifier not in features:
            continue
        feature = features[identifier]
        component_id = expected[identifier]
        if feature["chrom"] == "chr7":
            bucket = int(hashlib.sha256(("chr7-role-v1|" + component_id).encode()).hexdigest()[:16], 16) % 4
            role = "checkpoint" if bucket < 2 else "selection" if bucket == 2 else "calibration"
        else:
            role = feature["partition"]
        checks.add("independent_component_reconstruction", row["component_id"] == component_id, identifier)
        checks.add("independent_role_hash_rule", row["validation_role"] == role, identifier)
        for field in ("chrom", "partition", "canonical_rc_sequence_sha256"):
            checks.add(field, row[field] == feature[field], identifier)
        roles[identifier] = row
    checks.finish()
    check("roles:complete_master_interval_coverage", set(roles) == set(features), len(roles), len(features))
    component_roles = collections.defaultdict(set)
    for identifier, role in roles.items():
        if role["chrom"] == "chr7":
            component_roles[role["component_id"]].add(role["validation_role"])
    check("roles:no_chr7_component_cross_role", all(len(values) == 1 for values in component_roles.values()))
    return roles


def close_number(left, right, tolerance=1e-8):
    if left in ("", None) or right in ("", None):
        return left in ("", None) and right in ("", None)
    return abs(float(left) - float(right)) <= tolerance * max(1.0, abs(float(left)), abs(float(right)))


def verify_configurations(features, membership, roles, configuration_manifest):
    registry = list(table(DATA / "configuration_registry.tsv"))
    expected_registry = {(c, m) for c in ("V2-A", "V2-B", "V2-C") for m in MODELS}
    check("configurations:exact_six_registry_entries", len(registry) == 6 and {(r["configuration"], r["model"]) for r in registry} == expected_registry)
    configuration_sets = {}
    records = {}
    checks = AggregatedChecks("configuration_rows")
    numeric_fields = {"core_start", "core_end", "input_start", "input_end", "atac_anchor_count", "atac_signal_percentile_max", "gc_fraction", "repeat_fraction_2001", "blacklist_bp_1kb", "blacklist_bp_2001", "promoter_bp_1kb", "promoter_bp_2001", "sequence_available", "sequence_length", "non_acgt_fraction"}
    for registration in registry:
        config, model = registration["configuration"], registration["model"]
        key = (config, model)
        path = ROOT / registration["manifest"]
        check(f"configuration:{config}:{model}:manifest_hash", sha256(path) == registration["sha256"])
        sets = {"positive": set(), "control": set()}
        config_rows = {}
        for row in table(path):
            identifier = row["interval_id"]
            checks.add("unique_known_id", identifier in features and identifier not in config_rows, identifier)
            if identifier not in features:
                continue
            original = features[identifier]
            checks.add("binary_label", row["label"] in {"0", "1"}, identifier)
            label = "positive" if row["label"] == "1" else "control"
            sets[label].add(identifier)
            for field, value in row.items():
                if field in original:
                    equal = close_number(value, original[field]) if field in numeric_fields else value == original[field]
                    checks.add("unaltered_feature:" + field, equal, identifier)
            for field in ("component_id", "validation_role"):
                checks.add("role_consistency:" + field, row[field] == roles[identifier][field], identifier)
            checks.add("sequence_available_and_correct_length", row["sequence_available"] == "1" and row["sequence_length"] == "2001", identifier)
            if config in {"V2-B", "V2-C"} and label == "control":
                checks.add("new_controls_donor_ATAC_anchored", int(row["atac_anchor_count"]) > 0 and bool(row["atac_file_ids"]), identifier)
                checks.add("new_controls_no_positive_mark", row[f"{model}_all_lobe_peak_support"] == "0" and row[f"{model}_mark_overlap_2001"] == "0", identifier)
                checks.add("new_controls_no_positive_input_overlap", row[f"{model}_v1_positive_input_overlap"] == "0", identifier)
                checks.add("new_controls_eligible", row[f"{model}_accessible_control_eligible"] == "1", identifier)
                checks.add("new_controls_core_blacklist_promoter_free", row["blacklist_bp_1kb"] == "0" and row["promoter_bp_1kb"] == "0", identifier)
            if config == "V2-C" and label == "positive":
                checks.add("C_positive_same_lobe", row[f"{model}_same_lobe_peak_support"] == "1" and bool(row[f"{model}_same_lobe_histone_file_ids"]) and bool(row[f"{model}_same_lobe_support_lobes"]), identifier)
            config_rows[identifier] = row
        records[key] = config_rows
        configuration_sets[key] = sets
        check(f"configuration:{config}:{model}:registry_counts", len(config_rows) == int(registration["n_intervals"]) and len(sets["positive"]) == int(registration["positive"]) and len(sets["control"]) == int(registration["control"]))
        check(f"configuration:{config}:{model}:no_conflicting_class_membership", not sets["positive"] & sets["control"])
        if config in {"V2-A", "V2-B"}:
            check(f"configuration:{config}:{model}:exact_V1_positives", sets["positive"] == membership[(model, "positive")])
        if config == "V2-A":
            check(f"configuration:{config}:{model}:exact_V1_controls", sets["control"] == membership[(model, "control")])
        if config == "V2-C":
            same_lobe = {i for i in membership[(model, "positive")] if features[i][f"{model}_same_lobe_peak_support"] == "1"}
            check(f"configuration:{config}:{model}:exact_same_lobe_subset", sets["positive"] == same_lobe)
    checks.finish()
    for model in MODELS:
        check(f"configuration:{model}:B_C_identical_controls", configuration_sets[("V2-B", model)]["control"] == configuration_sets[("V2-C", model)]["control"])
    for flag in ("training_started", "model_inference_performed", "benchmark_outcomes_read"):
        check("configuration:no_execution:" + flag, configuration_manifest.get(flag) is False)
    return registry, configuration_sets, records


def matching_stratum(identifier, features, roles):
    row = features[identifier]
    return (row["chrom"], roles[identifier]["validation_role"], row["atac_lobes"],
            int(int(row["blacklist_bp_2001"]) > 0), int(float(row["non_acgt_fraction"] or 0) > 0))


def verify_matching(features, membership, roles, sets, manifest):
    import numpy as np
    from scipy.spatial import cKDTree

    spec = read_json(PROV / f"{PREFIX}_matching_specification.json")
    check("matching:three_fixed_training_seeds", spec.get("training_seeds") == SEEDS)
    check("matching:fixed_matching_seed", spec.get("matching_seed") == 271828)
    check("matching:fixed_bootstrap_seed", spec.get("bootstrap_seed") == 314159)
    check("matching:one_to_one_target", spec.get("requested_controls_per_B_positive") == 1)
    check("matching:no_replacement_or_positive_subsampling", spec.get("without_replacement") is True and spec.get("positive_subsampling") is False)
    check("matching:recorded_rule_hash", sha256(PROV / f"{PREFIX}_matching_specification.json") == manifest["matching_rules_sha256"])
    check("matching:feature_input_hash", sha256(DATA / "interval_features.tsv.gz") == spec["input_feature_sha256"])
    check("matching:implementation_hash", sha256(V2 / "scripts/construct_preflight_configurations.py") == spec["implementation_sha256"])
    for flag in ("no_training", "no_model_inference", "no_benchmark_outcome_access"):
        check("matching:" + flag, spec.get(flag) is True)
    pairs = {model: {} for model in MODELS}
    used_controls = {model: set() for model in MODELS}
    checks = AggregatedChecks("matching_pairs")
    variables = ("gc_fraction", "atac_signal_percentile_max", "repeat_fraction_2001")
    calipers = np.array([0.05, 0.20, 0.20])
    for row in table(DATA / "control_matching_pairs.tsv.gz"):
        model, positive, control = row["model"], row["positive_id"], row["control_id"]
        checks.add("known_model", model in MODELS, model)
        if model not in MODELS or positive not in features or control not in features:
            checks.add("known_intervals", False, positive + "|" + control)
            continue
        checks.add("one_match_per_positive", positive not in pairs[model], positive)
        checks.add("without_control_replacement", control not in used_controls[model], control)
        checks.add("positive_is_V1_positive", positive in membership[(model, "positive")], positive)
        checks.add("control_is_eligible", features[control][f"{model}_accessible_control_eligible"] == "1", control)
        checks.add("all_exact_strata", matching_stratum(positive, features, roles) == matching_stratum(control, features, roles), positive)
        differences = np.array([abs(float(features[positive][v]) - float(features[control][v])) for v in variables])
        checks.add("all_fixed_calipers", bool(np.all(differences <= calipers + 1e-10)), positive)
        for variable, difference in zip(variables, differences):
            checks.add(variable + "_reported_difference", close_number(row[variable + "_absolute_difference"], difference), positive)
        checks.add("reported_Chebyshev_distance", close_number(row["normalized_chebyshev_distance"], float(np.max(differences/calipers))), positive)
        for field in ("chrom", "atac_lobes"):
            checks.add("reported_" + field, row[field] == features[control][field], control)
        checks.add("reported_validation_role", row["validation_role"] == roles[control]["validation_role"], control)
        pairs[model][positive] = row
        used_controls[model].add(control)
    unmatched = {model: {} for model in MODELS}
    tie_scope_exceptions = []
    for row in table(DATA / "unmatched_positive_targets.tsv.gz"):
        model, identifier = row["model"], row["positive_id"]
        checks.add("unique_unmatched_positive", identifier not in unmatched[model], identifier)
        unmatched[model][identifier] = row["reason"]
    for model in MODELS:
        matched_ids = set(pairs[model])
        unmatched_ids = set(unmatched[model])
        check(f"matching:{model}:complete_positive_disposition", not matched_ids & unmatched_ids and matched_ids | unmatched_ids == membership[(model, "positive")])
        check(f"matching:{model}:selected_controls_equal_pairs", used_controls[model] == sets[("V2-B", model)]["control"])
        check(f"matching:{model}:manifest_matched_count", len(matched_ids) == int(manifest["matched_controls"].get(model, 0)))
        check(f"matching:{model}:manifest_unmatched_count", len(unmatched_ids) == int(manifest["unmatched_positive_targets"].get(model, 0)))
        fraction = len(matched_ids) / len(membership[(model, "positive")])
        readiness(f"{model}:B_matching_fraction_at_least_0.90", fraction >= 0.90, fraction, 0.90)
        # Independently test the greedy assignment against all available controls
        # within the observed distance, and exhaustive radius-1 neighborhoods for
        # unmatched targets. Never perform rematching or modify a constructed set.
        positive_strata = collections.defaultdict(list)
        control_strata = collections.defaultdict(list)
        for identifier in membership[(model, "positive")]:
            positive_strata[matching_stratum(identifier, features, roles)].append(identifier)
        for identifier, feature in features.items():
            if feature[f"{model}_accessible_control_eligible"] == "1":
                control_strata[matching_stratum(identifier, features, roles)].append(identifier)
        for stratum, positives in positive_strata.items():
            controls = sorted(control_strata.get(stratum, []))
            if not controls:
                for identifier in positives:
                    checks.add("unmatched_exact_stratum_reason", unmatched[model].get(identifier) == "no_control_in_exact_stratum", identifier)
                continue
            values = np.array([[float(features[i][v]) for v in variables] for i in controls]) / calipers
            tree = cKDTree(values)
            used = set()
            ordered = sorted(positives, key=lambda identifier: hashlib.sha256(f"271828|{model}|{identifier}".encode()).hexdigest())
            for identifier in ordered:
                point = np.array([float(features[identifier][v]) for v in variables]) / calipers
                record = pairs[model].get(identifier)
                radius = float(record["normalized_chebyshev_distance"]) + 1e-8 if record else 1.0 + 1e-12
                nearby = [int(i) for i in tree.query_ball_point(point, radius, p=np.inf) if int(i) not in used]
                if record:
                    distances = np.max(np.abs(values[nearby] - point), axis=1) if nearby else np.array([])
                    best = min(zip(distances, nearby), default=(None, None))
                    observed_control = record["control_id"]
                    expected_control = controls[best[1]] if best[1] is not None else ""
                    observed_index = bisect.bisect_left(controls, observed_control)
                    observed_distance = float(np.max(np.abs(values[observed_index] - point)))
                    checks.add("independent_minimum_available_distance", best[0] is not None and abs(observed_distance-float(best[0])) <= 1e-12,
                               identifier + "|" + repr(observed_distance) + "|" + repr(best[0]))
                    if observed_control != expected_control:
                        # The original code guarantees nearest available distance,
                        # but breaks coordinate ties within a bounded returned
                        # batch. Reproduce that documented implementation only for
                        # observed global-tie exceptions; never rematch data.
                        k = min(64, len(controls))
                        bounded_index = None
                        query_mode = "bounded_nearest_query"
                        while True:
                            kd_distances, kd_indices = tree.query(point, k=k, p=np.inf, distance_upper_bound=1.0+1e-12)
                            kd_distances = np.atleast_1d(kd_distances)
                            kd_indices = np.atleast_1d(kd_indices)
                            valid = np.isfinite(kd_distances) & (kd_indices < len(controls))
                            available_pairs = [(float(d), int(i)) for d, i in zip(kd_distances[valid], kd_indices[valid]) if int(i) not in used]
                            if available_pairs:
                                bounded_index = min(available_pairs)[1]
                                break
                            if k == len(controls) or not valid.all():
                                break
                            if k >= 1024:
                                query_mode = "radius_one_fallback"
                                candidates = [int(i) for i in tree.query_ball_point(point, r=1.0+1e-12, p=np.inf) if int(i) not in used]
                                if candidates:
                                    candidate_distances = np.max(np.abs(values[candidates]-point), axis=1)
                                    bounded_index = min(zip(candidate_distances, candidates))[1]
                                break
                            k = min(k*4, len(controls))
                        bounded_control = controls[bounded_index] if bounded_index is not None else ""
                        checks.add("documented_bounded_query_tie_assignment", bounded_control == observed_control,
                                   identifier + "|" + observed_control + "|" + bounded_control)
                        tie_scope_exceptions.append({
                            "model": model, "positive_id": identifier, "observed_control_id": observed_control,
                            "global_coordinate_first_at_minimum_distance": expected_control,
                            "documented_bounded_query_control_id": bounded_control,
                            "query_mode": query_mode, "query_k": k,
                            "global_minimum_available_distance": float(best[0]),
                            "observed_distance": observed_distance,
                            "unused_candidates_at_exact_minimum_distance": int(np.sum(distances == best[0])),
                            "calipers_unchanged": True, "frozen_control_unchanged": True,
                            "interpretation": "Equal-distance coordinate tie outside returned-batch scope; not a caliper or nearest-distance failure."
                        })
                    if observed_index < len(controls) and controls[observed_index] == observed_control:
                        used.add(observed_index)
                else:
                    checks.add("unmatched_has_no_unused_in_caliper_control", not nearby, identifier)
                    checks.add("unmatched_caliper_reason", unmatched[model].get(identifier) == "no_unused_control_within_all_calipers", identifier)
        print(f"Independently verified greedy matching for {model}", flush=True)
    checks.finish()
    notes = read_json(PROV / f"{PREFIX}_matching_implementation_notes.json")
    check("matching:bounded_tie_exception_explicitly_documented", bool(notes.get("nearest_neighbor_ties", {}).get("initial_wording_exception")))
    check("matching:original_strict_tie_validator_preserved", sha256(V2 / "scripts/validate_model_redesign_preflight_strict_tie_attempt.py") == read_json(PROV / f"{PREFIX}_preliminary_strict_tie_validation_manifest.json")["validator"]["sha256"])
    return tie_scope_exceptions


def verify_summary_tables(features, roles, sets, records, configuration_manifest):
    import numpy as np

    count_checks = AggregatedChecks("class_and_role_summaries")
    counts = list(table(RESULTS / f"{PREFIX}_class_counts.tsv"))
    expected_count_keys = {(c, m, p, label) for c in ("V2-A", "V2-B", "V2-C") for m in MODELS
                           for p in PARTITIONS for label in ("positive", "control")}
    count_keys = [(r["configuration"], r["model"], r["partition"], r["class"]) for r in counts]
    check("class_counts:complete_unique_strata", len(count_keys) == len(set(count_keys)) and set(count_keys) == expected_count_keys)
    for row in counts:
        config, model, part, label = row["configuration"], row["model"], row["partition"], row["class"]
        expected = {i for i in sets[(config, model)][label] if features[i]["partition"] == part}
        observed = ids_from_bed(ROOT / row["bed"])
        count_checks.add("exact_BED_membership", len(observed) == len(set(observed)) and set(observed) == expected, row["bed"])
        count_checks.add("BED_hash", sha256(ROOT / row["bed"]) == row["bed_sha256"], row["bed"])
        count_checks.add("interval_count", len(expected) == int(row["n_intervals"]), row["bed"])
        count_checks.add("component_count", len({roles[i]["component_id"] for i in expected}) == int(row["n_components"]), row["bed"])
    common = {}
    common_checks = AggregatedChecks("common_challenge")
    for model in MODELS:
        rows = {}
        for row in table(DATA / "evaluation" / f"{model}_common_challenge_panel.tsv.gz"):
            identifier = row["interval_id"]
            common_checks.add("unique_known_interval", identifier not in rows and identifier in records[("V2-C", model)], identifier)
            common_checks.add("validation_or_test_only", row["partition"] in {"validation", "test"}, identifier)
            if identifier in records[("V2-C", model)]:
                common_checks.add("exact_C_labels_and_metadata", row == records[("V2-C", model)][identifier], identifier)
            rows[identifier] = row
        expected = {i for ids in sets[("V2-C", model)].values() for i in ids if features[i]["partition"] in {"validation", "test"}}
        check(f"common_challenge:{model}:exact_shared_C_task", set(rows) == expected)
        common[model] = rows
    common_checks.finish()
    all_records = {**records, **{("COMMON", m): r for m, r in common.items()}}
    expected_role_counts = {}
    for key, rows in all_records.items():
        grouped = collections.defaultdict(list)
        for identifier, row in rows.items():
            grouped[(row["partition"], row["validation_role"])].append(identifier)
        for (part, role), identifiers in grouped.items():
            positive = sum(int(rows[i]["label"]) for i in identifiers)
            expected_role_counts[(*key, part, role)] = (positive, len(identifiers)-positive, len({roles[i]["component_id"] for i in identifiers}))
    actual_role_rows = list(table(RESULTS / f"{PREFIX}_role_class_counts.tsv"))
    actual_role_counts = {(r["configuration"], r["model"], r["partition"], r["role"]): r for r in actual_role_rows}
    check("role_class_counts:complete_unique_keys", len(actual_role_counts) == len(actual_role_rows) and set(actual_role_counts) == set(expected_role_counts))
    for key, (positive, control, components) in expected_role_counts.items():
        row = actual_role_counts.get(key, {})
        count_checks.add("role_positive", row.get("positive") == str(positive), key)
        count_checks.add("role_control", row.get("control") == str(control), key)
        count_checks.add("role_components", row.get("components") == str(components), key)
        count_checks.add("role_ratio", close_number(row.get("control_to_positive_ratio"), control/positive if positive else None), key)
        if key[0] == "COMMON" and key[3] == "selection":
            readiness(f"{key[1]}:common_selection_class_and_component_minima", positive >= 100 and control >= 200 and components >= 30,
                      {"positive": positive, "control": control, "components": components}, "100+/200-/30 components")
        if key[0] == "COMMON" and key[3] == "calibration":
            readiness(f"{key[1]}:common_calibration_at_least_200_controls", control >= 200, control, 200)
    count_checks.finish()
    overlap_checks = AggregatedChecks("configuration_overlap")
    overlap_rows = list(table(RESULTS / f"{PREFIX}_configuration_overlap.tsv"))
    check("configuration_overlap:six_rows", len(overlap_rows) == 6)
    for row in overlap_rows:
        model, first, second = row["model"], row["configuration_1"], row["configuration_2"]
        a = set(records[(first, model)])
        b = set(records[(second, model)])
        for field, value in (("n_1", len(a)), ("n_2", len(b)), ("shared", len(a & b)), ("only_1", len(a-b)), ("only_2", len(b-a))):
            overlap_checks.add(field, int(row[field]) == value, (model, first, second))
        overlap_checks.add("Jaccard", close_number(row["Jaccard"], len(a & b)/len(a | b)), (model, first, second))
    overlap_checks.finish()
    sequence_checks = AggregatedChecks("selected_sequence_context")
    sequence_rows = list(table(RESULTS / f"{PREFIX}_sequence_context_qc.tsv"))
    check("selected_sequence_context:six_rows", len(sequence_rows) == 6)
    for row in sequence_rows:
        key = (row["configuration"], row["model"])
        selected = records[key]
        counts_by_hash = collections.Counter(features[i]["canonical_rc_sequence_sha256"] for i in selected)
        expected = {
            "n_intervals": len(selected),
            "unavailable": sum(features[i]["sequence_available"] != "1" for i in selected),
            "wrong_length": sum(features[i]["sequence_length"] != "2001" for i in selected),
            "full_input_blacklist_overlap": sum(int(features[i]["blacklist_bp_2001"]) > 0 for i in selected),
            "full_input_promoter_overlap": sum(int(features[i]["promoter_bp_2001"]) > 0 for i in selected),
            "non_ACGT_sequences": sum(float(features[i]["non_acgt_fraction"] or 0) > 0 for i in selected),
            "within_set_exact_or_RC_duplicate_rows": sum(n for n in counts_by_hash.values() if n > 1)
        }
        for field, value in expected.items():
            sequence_checks.add(field, int(row[field]) == value, key)
        readiness(":".join(key) + ":sequence_available_2001bp", expected["unavailable"] == expected["wrong_length"] == 0,
                  {"unavailable": expected["unavailable"], "wrong_length": expected["wrong_length"]})
    sequence_checks.finish()
    active = set().union(*(set(r) for r in records.values()))
    groups = collections.defaultdict(list)
    for identifier in active:
        groups[features[identifier]["canonical_rc_sequence_sha256"]].append(identifier)
    cross = {key: ids for key, ids in groups.items() if len({features[i]["partition"] for i in ids}) > 1}
    reported_duplicate_rows = list(table(DATA / "selected_sequence_duplicate_audit.tsv.gz"))
    expected_ids = {i for ids in groups.values() if len(ids) > 1 for i in ids}
    check("selected_duplicates:complete_unique_rows", len(reported_duplicate_rows) == len(expected_ids) and {r["interval_id"] for r in reported_duplicate_rows} == expected_ids)
    duplicate_checks = AggregatedChecks("selected_duplicates")
    for row in reported_duplicate_rows:
        identifier = row["interval_id"]
        for field in ("partition", "component_id", "validation_role", "canonical_rc_sequence_sha256"):
            duplicate_checks.add(field, row[field] == roles[identifier][field], identifier)
    duplicate_checks.finish()
    check("selected_duplicates:cross_partition_manifest_count", len(cross) == int(configuration_manifest["cross_partition_selected_duplicate_groups"]), len(cross))
    check("selected:unique_interval_manifest_count", len(active) == int(configuration_manifest["selected_unique_intervals"]), len(active))
    check("selected:encoded_sequence_manifest_count", len(groups) == int(configuration_manifest["unique_selected_encoded_sequences"]), len(groups))
    check("selected:cache_bytes_manifest_arithmetic", len(groups)*36480 == int(configuration_manifest["future_two_orientation_feature_cache_bytes"]))
    readiness("no_selected_cross_partition_encoded_sequence_duplicates", len(cross) == 0, len(cross), 0)
    component_partitions = collections.defaultdict(set)
    for identifier in active:
        component_partitions[roles[identifier]["component_id"]].add(features[identifier]["partition"])
    cross_components = sum(len(p) > 1 for p in component_partitions.values())
    readiness("no_selected_component_spans_chromosome_partitions", cross_components == 0, cross_components, 0)
    return common, active


def verify_matching_diagnostics(features, records, configuration_manifest):
    import numpy as np

    checks = AggregatedChecks("matching_diagnostics")
    expected = {}
    quantitative = ("gc_fraction", "atac_signal_percentile_max", "repeat_fraction_2001", "blacklist_input_any", "non_acgt_any")
    def number(identifier, variable):
        row = features[identifier]
        if variable == "blacklist_input_any":
            return float(int(row["blacklist_bp_2001"]) > 0)
        if variable == "non_acgt_any":
            return float(float(row["non_acgt_fraction"] or 0) > 0)
        return float(row[variable]) if row[variable] else float("nan")
    for (config, model), rows in records.items():
        for part in PARTITIONS:
            positive = [i for i, r in rows.items() if r["partition"] == part and r["label"] == "1"]
            control = [i for i, r in rows.items() if r["partition"] == part and r["label"] == "0"]
            for variable in quantitative:
                a = np.array([number(i, variable) for i in positive])
                b = np.array([number(i, variable) for i in control])
                mean_a, mean_b = np.nanmean(a), np.nanmean(b)
                scale = np.sqrt((np.nanvar(a, ddof=1) + np.nanvar(b, ddof=1))/2)
                difference = mean_a - mean_b
                smd = difference/scale if scale > 0 else 0.0 if difference == 0 else float("nan")
                expected[(config, model, part, variable)] = {"positive_n": len(a), "control_n": len(b), "positive_mean": mean_a, "control_mean": mean_b, "SMD": smd, "status": "PASS" if np.isfinite(smd) and abs(smd) <= 0.10 else "FAIL"}
            for variable, bound in (("chrom", 0.02), ("atac_lobes", 0.05)):
                a = collections.Counter(features[i][variable] for i in positive)
                b = collections.Counter(features[i][variable] for i in control)
                for level in sorted(set(a) | set(b)):
                    ap, bp = a[level]/len(positive), b[level]/len(control)
                    gap = ap - bp
                    expected[(config, model, part, variable + "=" + level)] = {"positive_n": len(positive), "control_n": len(control), "positive_mean": ap, "control_mean": bp, "proportion_gap": gap, "status": "PASS" if abs(gap) <= bound else "FAIL"}
    rows = list(table(RESULTS / f"{PREFIX}_matching_diagnostics.tsv"))
    reported = {(r["configuration"], r["model"], r["partition"], r["variable"]): r for r in rows}
    check("matching_diagnostics:complete_unique_keys", len(reported) == len(rows) and set(reported) == set(expected))
    failures = []
    for key, values in expected.items():
        row = reported.get(key, {})
        for field, value in values.items():
            if field == "status":
                equal = row.get(field) == value
            elif not np.isfinite(value):
                equal = row.get(field) in {"", "nan", None}
            else:
                equal = close_number(row.get(field), value)
            checks.add(field, equal, key)
        if key[0] != "V2-A" and values["status"] == "FAIL":
            failures.append({"configuration": key[0], "model": key[1], "partition": key[2], "variable": key[3], **values})
    checks.finish()
    check("matching_diagnostics:manifest_B_C_failure_count", len(failures) == int(configuration_manifest["B_C_balance_failures"]), len(failures), configuration_manifest["B_C_balance_failures"])
    readiness("all_B_C_prespecified_balance_gates", not failures, len(failures), 0)
    return failures


def verify_run_matrix(registry, records):
    rows = list(table(DATA / "run_matrix.tsv"))
    expected_keys = {(c, m, seed) for c in ("V2-A", "V2-B", "V2-C") for m in MODELS for seed in SEEDS}
    observed_keys = [(r["configuration"], r["model"], int(r["seed"])) for r in rows]
    check("run_matrix:exact_18_prespecified_runs", len(rows) == 18 and len(set(observed_keys)) == 18 and set(observed_keys) == expected_keys)
    by_key = {(r["configuration"], r["model"]): r for r in registry}
    checks = AggregatedChecks("run_matrix")
    expected = {"epochs_max": "50", "batch_size": "256", "patience": "15", "optimizer": "Adadelta", "gpu_type": "a100", "gpus": "1", "cpus": "8", "memory_GiB": "48", "walltime_cap_hours": "4", "status": "NOT_AUTHORIZED_NOT_STARTED"}
    for row in rows:
        key = (row["configuration"], row["model"])
        for field, value in expected.items():
            checks.add(field, row[field] == value, row["run_id"])
        checks.add("learning_rate", close_number(row["learning_rate"], 0.001), row["run_id"])
        checks.add("manifest", row["manifest"] == by_key[key]["manifest"] and row["manifest_sha256"] == by_key[key]["sha256"], row["run_id"])
        checks.add("training_intervals", int(row["training_intervals"]) == sum(r["partition"] == "train" for r in records[key].values()), row["run_id"])
    checks.finish()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-skeleton", action="store_true",
                        help="Check static validator setup only; write nothing and do not claim full validation.")
    args = parser.parse_args()
    if args.check_skeleton:
        print(json.dumps({"status": "STATIC_ONLY_NO_DATA_VALIDATION", "no_model_execution": True,
                          "implemented": "full source, sequence, matching, role and configuration audit"}))
        return
    started = time.time()
    output = PROV / f"{PREFIX}_final_validation.tsv"
    readiness_output = PROV / f"{PREFIX}_training_readiness.tsv"
    manifest_output = PROV / f"{PREFIX}_final_validation_manifest.json"
    tie_output = PROV / f"{PREFIX}_matching_tie_contract_audit.json"
    if any(p.exists() for p in (output, readiness_output, manifest_output, tie_output)):
        raise RuntimeError("Refusing to overwrite an existing validation; preserve/version explicitly")
    configuration_path = PROV / f"{PREFIX}_configuration_manifest.json"
    if not configuration_path.is_file():
        raise RuntimeError("Configuration construction is incomplete; no validation outputs written")
    configuration_manifest = read_json(configuration_path)
    feature_manifest = require_complete_features()
    print("Verifying original source labels and immutable V1 membership", flush=True)
    membership = original_bed_sets()
    anchors, marks, atac, sources = original_source_evidence()
    features, duplicates = verify_features(feature_manifest, membership, anchors, marks, atac, sources)
    del anchors, marks, atac, sources
    print("Verifying global overlap components and fixed validation roles", flush=True)
    roles = verify_roles(features, duplicates)
    registry, sets, records = verify_configurations(features, membership, roles, configuration_manifest)
    print("Verifying actual matched-control assignments, without rematching", flush=True)
    tie_exceptions = verify_matching(features, membership, roles, sets, configuration_manifest)
    common, active = verify_summary_tables(features, roles, sets, records, configuration_manifest)
    balance_failures = verify_matching_diagnostics(features, records, configuration_manifest)
    verify_run_matrix(registry, records)
    verify_orientation_contract()
    print("Independently rehashing protected artifacts and consumed inputs as opaque bytes", flush=True)
    protected_count = verify_ledger(PROV / f"{PREFIX}_protected_tracked_hashes_before.tsv", "protected")
    consumed_count = verify_ledger(PROV / f"{PREFIX}_consumed_input_hashes.tsv", "consumed")
    check("runtime:no_TensorFlow_or_Keras_import", not any(name == "tensorflow" or name.startswith("tensorflow.") or name == "keras" or name.startswith("keras.") for name in sys.modules))
    readiness("investigator_training_authorization", False, "Not requested or granted by this preflight", "Explicit investigator approval after review")
    failures = [r for r in CHECKS if r["status"] == "FAIL"]
    blocked = [r for r in READINESS if r["status"] == "BLOCKED"]
    preserved_outputs = [
        (f"{PREFIX}_final_validation.tsv", f"{PREFIX}_preliminary_strict_tie_validation.tsv"),
        (f"{PREFIX}_training_readiness.tsv", f"{PREFIX}_preliminary_training_readiness.tsv"),
        (f"{PREFIX}_final_validation_manifest.json", f"{PREFIX}_preliminary_strict_tie_validation_manifest.json")
    ]
    tie_audit = {
        "module": PREFIX, "strict_global_coordinate_tie_property": "DOES_NOT_HOLD_FOR_ALL_OBSERVED_PAIRS" if tie_exceptions else "OBSERVED_FOR_ALL_PAIRS",
        "strict_minimum_available_distance_tested_for_all_pairs": True,
        "exception_count": len(tie_exceptions), "exceptions": tie_exceptions,
        "final_contract": "Minimum unused distance; coordinate-ID tie resolution within actual returned nearest-k batch, as the unchanged executed constructor implements.",
        "no_rematching_or_dataset_edits": True,
        "preliminary_evidence_relocations": [{"original_historical_path": str((PROV / original).relative_to(ROOT)),
                                              "preserved_path": str((PROV / preserved).relative_to(ROOT)),
                                              "sha256": sha256(PROV / preserved)} for original, preserved in preserved_outputs],
        "initial_validator_preserved_path": str((V2 / "scripts/validate_model_redesign_preflight_strict_tie_attempt.py").relative_to(ROOT)),
        "implementation_notes_sha256": sha256(PROV / f"{PREFIX}_matching_implementation_notes.json")
    }
    with tie_output.open("x") as handle:
        json.dump(tie_audit, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected", "domain"], delimiter="\t")
        writer.writeheader()
        writer.writerows(CHECKS)
    with readiness_output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition", "status", "observed", "required"], delimiter="\t")
        writer.writeheader()
        writer.writerows(READINESS)
    import resource
    result = {
        "module": PREFIX, "completed_utc": utc(), "artifact_validation_status": "PASS" if not failures else "FAIL",
        "artifact_checks": len(CHECKS), "artifact_failures": len(failures),
        "training_readiness": "BLOCKED_NOT_AUTHORIZED" if blocked else "READY_PENDING_EXPLICIT_AUTHORIZATION",
        "training_readiness_conditions": READINESS, "B_C_balance_failures": balance_failures,
        "documented_equal_distance_tie_scope_exceptions": len(tie_exceptions),
        "master_intervals_independently_source_and_sequence_checked": len(features),
        "selected_unique_intervals": len(active), "configuration_count": len(registry), "planned_fit_count": 18,
        "protected_files_opaque_hash_verified": protected_count, "consumed_files_hash_verified": consumed_count,
        "real_network_invariance": "NOT_RUN; no training/inference authorized", "GPU_model_execution": False,
        "benchmark_content_parsing": False, "elapsed_seconds": time.time()-started,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "validator": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": sha256(__file__)},
        "outputs": [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in (output, readiness_output, tie_output)],
        "failure_details": failures
    }
    with manifest_output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: result[k] for k in ("artifact_validation_status", "artifact_checks", "artifact_failures", "training_readiness", "elapsed_seconds", "peak_rss_kib")}), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
