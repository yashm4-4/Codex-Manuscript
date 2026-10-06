#!/usr/bin/env python3
"""Independent, CPU-only validation of the COPD pretraining-1.1 construction.

No training implementation is imported. Historical benchmark files are accessed
only as opaque bytes when verifying immutability. Validation outputs are created
exclusively, never used to overwrite a previous validation attempt.
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
import time

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[5]
OLD = ROOT / "diseases/COPD/07_gap_closure"
BASE = OLD / "pretraining-1.1"
PROV = BASE / "provenance"
DATA = BASE / "data"
MODELS = ("enhancer", "h3k27me3")
SEEDS = [104729, 130363, 155921]
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
CHECKS, READINESS, BALANCE, TIES = [], [], [], []
HISTORICAL_FREEZE_SHA = "2d022d19b4fbdd485d0762c0946519334e654b1de0581c5075758b7550d15fcb"
HISTORICAL_LEDGER_SHA = "7a6440c17be3374e5b4732c0bd6a5450fa61b74e4849e99a19b856a02e9e5055"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def guard(path):
    if "COPD-V2-BENCH" in str(path):
        raise RuntimeError("Benchmark parsing is forbidden")
    return Path(path)


def read_json(path):
    with guard(path).open() as f:
        return json.load(f)


def table(path, **kwargs):
    return pd.read_csv(guard(path), sep="\t", keep_default_na=False, **kwargs)


def records(path):
    p = guard(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as f:
        yield from csv.DictReader(f, delimiter="\t")


def check(name, condition, observed="", expected="", domain="artifact"):
    CHECKS.append({"check": name, "status": "PASS" if condition else "FAIL",
                   "observed": str(observed), "expected": str(expected), "domain": domain})
    return bool(condition)


def ready(name, condition, observed="", expected=""):
    READINESS.append({"condition": name, "status": "SATISFIED" if condition else "BLOCKED",
                      "observed": str(observed), "required": str(expected)})


class Aggregate:
    def __init__(self, prefix):
        self.prefix = prefix
        self.n = collections.Counter()
        self.bad = collections.Counter()
        self.examples = collections.defaultdict(list)

    def add(self, name, condition, example=""):
        self.n[name] += 1
        if not condition:
            self.bad[name] += 1
            if len(self.examples[name]) < 5:
                self.examples[name].append(str(example))

    def finish(self):
        for name in sorted(self.n):
            check(self.prefix + ":" + name, self.bad[name] == 0,
                  json.dumps({"n": self.n[name], "failures": self.bad[name],
                              "examples": self.examples[name]}), "zero failures")


class UnionIntervals:
    """Independent half-open union, with overlap and covered-base queries."""
    def __init__(self, rows):
        groups = collections.defaultdict(list)
        for chrom, start, end in rows:
            if end > start:
                groups[chrom].append((int(start), int(end)))
        self.data = {}
        for chrom, rows in groups.items():
            merged = []
            for start, end in sorted(rows):
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            self.data[chrom] = ([r[0] for r in merged], [r[1] for r in merged])

    def covered(self, chrom, start, end):
        starts, ends = self.data.get(chrom, ([], []))
        i = bisect.bisect_right(ends, start)
        n = 0
        while i < len(starts) and starts[i] < end:
            n += max(0, min(end, ends[i]) - max(start, starts[i]))
            i += 1
        return n

    def overlap(self, chrom, start, end):
        starts, ends = self.data.get(chrom, ([], []))
        i = bisect.bisect_right(ends, start)
        return i < len(starts) and starts[i] < end


def partition(chrom):
    return "validation" if chrom == "chr7" else "test" if chrom in {"chr8", "chr9"} else "train" if chrom in CANONICAL else "excluded"


def ledger(path, prefix):
    n = 0
    for r in records(path):
        p = ROOT / r["path"]
        n += 1
        ok = p.is_file()
        check(prefix + ":exists:" + r["path"], ok)
        if ok:
            check(prefix + ":bytes:" + r["path"], p.stat().st_size == int(r["bytes"]), p.stat().st_size, r["bytes"])
            check(prefix + ":sha256:" + r["path"], sha(p) == r["sha256"], expected=r["sha256"])
    return n


def historical_integrity():
    freeze = OLD / "provenance/COPD-V2-PREFLIGHT_freeze.json"
    hashes = OLD / "provenance/COPD-V2-PREFLIGHT_artifact_checksums.tsv"
    check("historical_1.0:freeze_unchanged", sha(freeze) == HISTORICAL_FREEZE_SHA)
    check("historical_1.0:ledger_unchanged", sha(hashes) == HISTORICAL_LEDGER_SHA)
    check("historical_1.0:all_100_payloads_present", ledger(hashes, "historical_1.0") == 100)
    check("historical_1.0:all_28_consumed_inputs_unchanged", ledger(OLD / "provenance/COPD-V2-PREFLIGHT_consumed_input_hashes.tsv", "consumed_input") == 28)
    prior = read_json(freeze)
    check("historical_1.0:failed_preflight_status_preserved", prior["status"] == "FROZEN_REVIEW_ONLY_TRAINING_BLOCKED" and prior["matching_balance_failures"] == 18)
    benchmark_files = 0
    for r in records(OLD / "provenance/COPD-V2-PREFLIGHT_protected_tracked_hashes_before.tsv"):
        if "BENCH" in r["path"] or "benchmark" in r["path"].lower():
            target = ROOT / r["path"]
            benchmark_files += 1
            check("benchmark_opaque_integrity:" + r["path"], target.is_file() and target.stat().st_size == int(r["bytes"]) and sha(target) == r["sha256"])
    check("benchmark_opaque_integrity:files_checked", benchmark_files > 0, benchmark_files)


def load_master():
    names = ["interval_id", "chrom", "core_start", "core_end", "input_start", "input_end", "partition", "atac_lobes", "atac_file_ids", "atac_peak_ids", "atac_anchor_count", "atac_signal_percentile_max", "gc_fraction", "repeat_fraction_2001", "blacklist_bp_1kb", "blacklist_bp_2001", "promoter_bp_1kb", "promoter_bp_2001", "sequence_available", "sequence_length", "non_acgt_fraction", "sequence_sha256", "canonical_rc_sequence_sha256"]
    for m in MODELS:
        names += [f"v1_{m}_positive", f"v1_{m}_control", f"{m}_same_lobe_peak_support", f"{m}_all_lobe_peak_support", f"{m}_mark_overlap_2001", f"{m}_v1_positive_input_overlap", f"{m}_accessible_control_eligible"]
    frame = table(OLD / "data/COPD-V2-PREFLIGHT/interval_features.tsv.gz", usecols=names)
    check("master:unique_interval_ids", frame.interval_id.is_unique, len(frame))
    frame = frame.set_index("interval_id", drop=False)
    roles = table(OLD / "data/COPD-V2-PREFLIGHT/interval_role_assignment.tsv.gz").set_index("interval_id")
    check("roles:complete_master_coverage", set(frame.index) == set(roles.index))
    frame["component_id"] = roles.loc[frame.index, "component_id"]
    frame["validation_role"] = roles.loc[frame.index, "validation_role"]
    frame["blacklist_input_any"] = frame.blacklist_bp_2001.gt(0).astype(int)
    frame["non_acgt_any"] = frame.non_acgt_fraction.gt(0).astype(int)
    frame["atac_signal_percentile_max"] = pd.to_numeric(frame.atac_signal_percentile_max, errors="coerce")
    return frame


def verify_global_roles(frame):
    """Rebuild all-master overlap and encoded-sequence components independently."""
    ordered = frame.sort_values(["chrom", "input_start", "input_end"])
    initial, names = {}, []
    last_chrom, last_end, g = None, -1, -1
    for identifier, chrom, start, end in ordered[["interval_id", "chrom", "input_start", "input_end"]].itertuples(index=False, name=None):
        if chrom != last_chrom or start >= last_end:
            g += 1
            names.append(f"{chrom}:{start}")
            last_end = end
        else:
            last_end = max(last_end, end)
        initial[identifier] = g
        last_chrom = chrom
    parent = list(range(g + 1))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    dup = frame[frame.canonical_rc_sequence_sha256.duplicated(keep=False)]
    for _, subset in dup.groupby("canonical_rc_sequence_sha256", sort=False):
        groups = [initial[i] for i in subset.index]
        for other in groups[1:]:
            a, b = find(groups[0]), find(other)
            parent[max(a, b)] = min(a, b)
    agg = Aggregate("global_roles")
    for r in frame[["interval_id", "chrom", "partition", "component_id", "validation_role"]].itertuples(index=False):
        component = names[find(initial[r.interval_id])]
        part = partition(r.chrom)
        if r.chrom == "chr7":
            bucket = int(hashlib.sha256(("chr7-role-v1|" + component).encode()).hexdigest()[:16], 16) % 4
            role = "checkpoint" if bucket < 2 else "selection" if bucket == 2 else "calibration"
        else:
            role = part
        agg.add("component_reconstructed", r.component_id == component, r.interval_id)
        agg.add("partition_fixed", r.partition == part, r.interval_id)
        agg.add("role_reconstructed", r.validation_role == role, r.interval_id)
    agg.finish()


def verify_normalization(frame):
    """Reconstruct all three training-only CDFs from original narrowPeak rows."""
    spec = read_json(BASE / "specification/atac_normalization.json")
    check("ATAC:training_chromosomes_exact", set(spec["training_chromosomes"]) == CANONICAL - {"chr7", "chr8", "chr9"})
    sources = table(OLD / "provenance/COPD-V2-PREFLIGHT_source_inventory.tsv")
    peaks, distributions = {}, {}
    for source in sources[sources.mark.eq("accessibility")].itertuples(index=False):
        source_peaks, training = [], []
        with gzip.open(ROOT / source.path, "rt") as f:
            for number, line in enumerate(f, 1):
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                v = line.rstrip().split("\t")
                if v[0] not in CANONICAL:
                    continue
                signal = float(v[6])
                start, end, offset = int(v[1]), int(v[2]), int(v[9])
                summit = start + offset if offset >= 0 else (start+end)//2
                interval = f"{v[0]}:{summit-500}-{summit+500}"
                source_peaks.append((f"{source.file_accession}:{number}", signal, interval))
                if partition(v[0]) == "train":
                    training.append(signal)
        ordered = sorted(training)
        distributions[source.file_accession] = (ordered, collections.Counter(training))
        for peak, signal, interval in source_peaks:
            value = (bisect.bisect_left(ordered, signal)+bisect.bisect_right(ordered, signal))/(2*len(ordered))
            peaks[peak] = (source.file_accession, signal, interval, value)
    agg = Aggregate("ATAC_training_mapping")
    observed_cdf = collections.defaultdict(set)
    for r in records(DATA / "atac_training_cdf.tsv.gz"):
        file_id, signal = r["file_accession"], float(r["signal_value"])
        ordered, counts = distributions[file_id]
        value = (bisect.bisect_left(ordered, signal)+bisect.bisect_right(ordered, signal))/(2*len(ordered))
        agg.add("CDF_unique_source_signal", signal not in observed_cdf[file_id], file_id)
        observed_cdf[file_id].add(signal)
        agg.add("CDF_equal_count", int(r["training_equal_count"]) == counts[signal], file_id)
        agg.add("CDF_less_count", int(r["training_less_count"]) == bisect.bisect_left(ordered, signal), file_id)
        agg.add("CDF_train_denominator", int(r["training_denominator"]) == len(ordered), file_id)
        agg.add("CDF_value", abs(float(r["mapped_mid_ecdf"])-value) <= 1e-15, file_id)
    for file_id, (_, counter) in distributions.items():
        check("ATAC:CDF_complete:" + file_id, observed_cdf[file_id] == set(counter), len(observed_cdf[file_id]), len(counter))
    values = []
    peak_seen = set()
    for n, r in enumerate(records(DATA / "atac_normalized_features.tsv.gz")):
        identifier = r["interval_id"]
        agg.add("master_row_order", n < len(frame) and identifier == frame.index[n], identifier)
        f = frame.loc[identifier]
        names = r["atac_peak_ids"].split(";") if r["atac_peak_ids"] else []
        files = r["atac_peak_file_ids"].split(";") if names else []
        signals = [float(x) for x in r["atac_peak_signal_values"].split(";")] if names else []
        ranks = [float(x) for x in r["atac_peak_train_percentiles"].split(";")] if names else []
        agg.add("peak_list_lengths", len(names) == len(files) == len(signals) == len(ranks) == int(r["atac_anchor_count"]) == int(f.atac_anchor_count), identifier)
        agg.add("peak_names_match_frozen_master", r["atac_peak_ids"] == f.atac_peak_ids, identifier)
        expected = []
        for peak, file_id, signal, rank in zip(names, files, signals, ranks):
            raw = peaks[peak]
            agg.add("unique_original_peak_row", peak not in peak_seen, peak)
            peak_seen.add(peak)
            agg.add("peak_source_signal_interval", raw[:3] == (file_id, signal, identifier), peak)
            agg.add("peak_train_only_percentile", abs(rank-raw[3]) <= 1e-15, peak)
            expected.append(raw[3])
        value = max(expected) if expected else np.nan
        observed = float(r["atac_signal_percentile_max_train_only"]) if r["atac_signal_percentile_max_train_only"] else np.nan
        agg.add("anchor_maximum", (np.isnan(value) and np.isnan(observed)) or abs(value-observed) <= 1e-15, identifier)
        values.append(observed)
    check("ATAC:all_master_rows", len(values) == len(frame), len(values), len(frame))
    check("ATAC:all_original_peak_rows_once", peak_seen == set(peaks), len(peak_seen), len(peaks))
    # Independent CDF/source validation above uses Python's correctly rounded
    # decimal reader. Matching must replay the constructor's actual binary64
    # operands, however: pandas' default C-reader decimal conversion can differ
    # by an ulp. Preserve that explicit numerical-ingestion contract; never
    # replace exact matching ties by approximate-isclose ties.
    numeric = pd.read_csv(DATA / "atac_normalized_features.tsv.gz", sep="\t",
                          usecols=["interval_id", "atac_signal_percentile_max_train_only"])
    numeric = numeric.set_index("interval_id").loc[frame.index,"atac_signal_percentile_max_train_only"].to_numpy(float)
    independently_reconstructed = np.asarray(values,dtype=float)
    check("ATAC:matching_numeric_ingestion_within_one_ulp_scale", bool(np.isclose(numeric,independently_reconstructed,atol=1e-15,rtol=0,equal_nan=True).all()))
    frame["atac_signal_percentile_max_train_only"] = numeric
    agg.finish()
    contract = [1, 1, 3, 5]
    observed = [(bisect.bisect_left(contract,x)+bisect.bisect_right(contract,x))/(2*len(contract)) for x in (-1,1,2,3,4,5,6)]
    check("ATAC:synthetic_ties_steps_outside_range", observed == [0,.25,.5,.625,.75,.875,1], observed)
    return frame


def raw_source_indices():
    source = table(OLD / "provenance/COPD-V2-PREFLIGHT_source_inventory.tsv")
    check("source:one_severe_emphysema_donor", set(source.donor_accessions) == {"ENCDO520EJG"})
    check("source:GRCh38", set(source.assembly) == {"GRCh38"})
    mark_rows = collections.defaultdict(list)
    atac = []
    for r in source.itertuples(index=False):
        with gzip.open(ROOT / r.path, "rt") as f:
            for line_number, line in enumerate(f, 1):
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                v = line.rstrip().split("\t")
                if v[0] not in CANONICAL:
                    continue
                chrom, start, end = v[0], int(v[1]), int(v[2])
                if r.mark == "accessibility":
                    summit = start + int(v[9]) if int(v[9]) >= 0 else (start + end) // 2
                    atac.append((f"{chrom}:{summit-500}-{summit+500}", chrom, start, end, r.lobe, f"{r.file_accession}:{line_number}", r.file_accession))
                else:
                    model = "enhancer" if r.mark == "H3K27ac" else "h3k27me3"
                    mark_rows[(model, r.lobe)].append((chrom, start, end))
    mark_lobe = {key: UnionIntervals(rows) for key, rows in mark_rows.items()}
    marks = {m: UnionIntervals(v for (model, _), rows in mark_rows.items() if model == m for v in rows) for m in MODELS}
    anchors = {}
    for identifier, chrom, start, end, lobe, peak, file_id in atac:
        row = anchors.setdefault(identifier, {"lobes": set(), "peaks": set(), "files": set(), "enhancer_same": False, "h3k27me3_same": False})
        row["lobes"].add(lobe)
        row["peaks"].add(peak)
        row["files"].add(file_id)
        for m in MODELS:
            row[m + "_same"] |= mark_lobe[(m, lobe)].overlap(chrom, start, end)
    return anchors, marks


def annotation_indices():
    rows = []
    with gzip.open(ROOT / "data/blacklist/hg38-blacklist.v2.bed.gz", "rt") as f:
        for line in f:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                x = line.split("\t")
                rows.append((x[0], int(x[1]), int(x[2])))
    black = UnionIntervals(rows)
    rows = []
    with gzip.open(ROOT / "data/gencode/gencode.v50.annotation.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t")
            if x[0] in CANONICAL and x[2] == "gene":
                tss = int(x[3] if x[6] == "+" else x[4]) - 1
                rows.append((x[0], max(0, tss-2000), tss+2000))
    promoters = UnionIntervals(rows)
    def repeats():
        with gzip.open(ROOT / "data/repeats/rmsk.txt.gz", "rt") as f:
            for line in f:
                x = line.split("\t")
                if len(x) >= 8 and x[5] in CANONICAL:
                    yield x[5], int(x[6]), int(x[7])
    return black, promoters, UnionIntervals(repeats())


def verify_selected_sources(frame, sets):
    import pysam
    selected = set().union(*(set(x.index) for x in sets.values()))
    anchors, marks = raw_source_indices()
    black, promoters, repeats = annotation_indices()
    pos_indices = {m: UnionIntervals(frame.loc[frame[f"v1_{m}_positive"].eq(1), ["chrom", "input_start", "input_end"]].itertuples(index=False, name=None)) for m in MODELS}
    new_controls = {m: set().union(*(set(x.index[x.label.eq(0)]) for (c, model), x in sets.items() if c != "V2-A" and model == m)) for m in MODELS}
    c_positives = {m: set(sets[("V2-C", m)].index[sets[("V2-C", m)].label.eq(1)]) for m in MODELS}
    fasta = pysam.FastaFile(str(OLD / "../04_modeling/trednet/fasta/hg38.fa"))
    refs = set(fasta.references)
    complement = str.maketrans("ACGTN", "TGCAN")
    agg = Aggregate("selected_source_and_sequence")
    for n, r in enumerate(frame.loc[sorted(selected)].itertuples(index=False), 1):
        reference_chrom = r.chrom if r.chrom in refs else r.chrom.removeprefix("chr")
        seq = fasta.fetch(reference_chrom, r.input_start, r.input_end).upper()
        normalized = seq if sum(seq.count(x) for x in "ACGT") == len(seq) else "".join(x if x in "ACGT" else "N" for x in seq)
        canonical = min(normalized, normalized.translate(complement)[::-1])
        agg.add("2001bp_available", r.sequence_available == 1 and len(seq) == r.sequence_length == 2001, r.interval_id)
        agg.add("raw_sequence_hash", hashlib.sha256(seq.encode()).hexdigest() == r.sequence_sha256, r.interval_id)
        agg.add("encoded_forward_RC_hash", hashlib.sha256(canonical.encode()).hexdigest() == r.canonical_rc_sequence_sha256, r.interval_id)
        agg.add("GC_recomputed", abs((seq.count("G") + seq.count("C"))/2001-r.gc_fraction) <= 1e-9, r.interval_id)
        agg.add("ambiguity_recomputed", abs(normalized.count("N")/2001-r.non_acgt_fraction) <= 1e-9, r.interval_id)
        agg.add("repeat_recomputed", abs(repeats.covered(r.chrom, r.input_start, r.input_end)/2001-r.repeat_fraction_2001) <= 1e-9, r.interval_id)
        for word, index in (("blacklist", black), ("promoter", promoters)):
            for suffix, start, end in (("1kb", r.core_start, r.core_end), ("2001", r.input_start, r.input_end)):
                agg.add(word + "_" + suffix, index.covered(r.chrom, start, end) == getattr(r, word + "_bp_" + suffix), r.interval_id)
        anchor = anchors.get(r.interval_id)
        if anchor:
            agg.add("source_anchor_provenance", set(r.atac_lobes.split(";")) == anchor["lobes"] and set(r.atac_file_ids.split(";")) == anchor["files"] and set(r.atac_peak_ids.split(";")) == anchor["peaks"], r.interval_id)
        for m in MODELS:
            if r.interval_id in new_controls[m]:
                agg.add(m + ":donor_anchor", anchor is not None, r.interval_id)
                agg.add(m + ":full_mark_absent", not marks[m].overlap(r.chrom, r.input_start, r.input_end), r.interval_id)
                agg.add(m + ":full_known_positive_absent", not pos_indices[m].overlap(r.chrom, r.input_start, r.input_end), r.interval_id)
                agg.add(m + ":blacklist_TSS_core_exclusions", r.blacklist_bp_1kb == 0 and r.promoter_bp_1kb == 0 and r.core_start >= 501, r.interval_id)
            if r.interval_id in c_positives[m]:
                agg.add(m + ":same_lobe_positive_reconstructed", bool(anchor and anchor[m + "_same"]), r.interval_id)
        if n % 100000 == 0:
            print(f"Validated source/sequence of {n}/{len(selected)} selected intervals", flush=True)
    fasta.close()
    agg.finish()
    active = frame.loc[sorted(selected)]
    for field in ("canonical_rc_sequence_sha256", "component_id"):
        bad = int(active.groupby(field).partition.nunique().gt(1).sum())
        check("leakage:" + field + ":no_cross_partition", bad == 0, bad, 0)
        ready("no_cross_partition_" + field, bad == 0, bad, 0)
        validation = active[active.chrom.eq("chr7")]
        cross_role = int(validation.groupby(field).validation_role.nunique().gt(1).sum())
        check("leakage:" + field + ":no_chr7_cross_role", cross_role == 0, cross_role, 0)
    return {"selected_unique_intervals": len(selected), "selected_unique_encoded_sequences": active.canonical_rc_sequence_sha256.nunique()}


def compare_balance(frame, configuration, model):
    for part, group in frame.groupby("partition", sort=True):
        pos, neg = group[group.label.eq(1)], group[group.label.eq(0)]
        check(f"balance:{configuration}:{model}:{part}:classes_n_at_least_two",len(pos) >= 2 and len(neg) >= 2,f"{len(pos)}/{len(neg)}",">=2 each")
        for variable in ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001", "blacklist_input_any", "non_acgt_any"]:
            a, b = pos[variable].to_numpy(float), neg[variable].to_numpy(float)
            if configuration != "V2-A":
                check(f"balance:{configuration}:{model}:{part}:{variable}:finite_covariates",bool(np.isfinite(a).all() and np.isfinite(b).all()))
            pooled = np.sqrt((np.nanvar(a, ddof=1) + np.nanvar(b, ddof=1))/2)
            delta = float(np.nanmean(a)-np.nanmean(b))
            smd = delta/pooled if pooled > 0 else 0.0 if delta == 0 else np.inf
            ok = np.isfinite(smd) and abs(smd) <= 0.10
            row = {"configuration": configuration, "model": model, "partition": part, "variable": variable, "positive_n": len(a), "control_n": len(b), "positive_mean": np.nanmean(a), "control_mean": np.nanmean(b), "SMD": smd, "proportion_gap": "", "limit": 0.10, "status": "PASS" if ok else "FAIL"}
            BALANCE.append(row)
        for variable, limit in [("chrom", .02), ("atac_lobes", .05)]:
            a, b = pos[variable].value_counts(normalize=True), neg[variable].value_counts(normalize=True)
            for level in sorted(set(a.index) | set(b.index)):
                x, y = float(a.get(level, 0)), float(b.get(level, 0))
                delta = x-y
                BALANCE.append({"configuration": configuration, "model": model, "partition": part, "variable": variable + "=" + str(level), "positive_n": len(pos), "control_n": len(neg), "positive_mean": x, "control_mean": y, "SMD": "", "proportion_gap": delta, "limit": limit, "status": "PASS" if abs(delta) <= limit else "FAIL"})
        for lobe in ("lower_left", "lower_right", "upper_right"):
            x = float(pos.atac_lobes.str.split(";").map(lambda s: lobe in s).mean())
            y = float(neg.atac_lobes.str.split(";").map(lambda s: lobe in s).mean())
            BALANCE.append({"configuration": configuration, "model": model, "partition": part, "variable": "marginal_lobe=" + lobe, "positive_n": len(pos), "control_n": len(neg), "positive_mean": x, "control_mean": y, "SMD": "", "proportion_gap": x-y, "limit": .05, "status": "PASS" if abs(x-y) <= .05 else "FAIL"})


def verify_selection_specification():
    old = read_json(OLD / "provenance/COPD-V2-PREFLIGHT_selection_specification.json")
    new = read_json(BASE / "specification/selection_adequacy.json")
    for field in ("seeds", "phase_I_weights_sha256", "training", "augmentation", "checkpoint"):
        check("selection:preserved:" + field, new[field] == old[field])
    for field, value in old["inference"].items():
        check("selection:preserved:inference:" + field, new["inference"].get(field) == value)
    for field in ("formula", "call", "target_empirical_FPR", "minimum_controls", "tie_policy", "test_recalibration", "variant_FPR_guarantee"):
        check("selection:preserved:threshold:" + field, new["threshold"].get(field) == old["threshold"].get(field))
    check("selection:explicit_external_recalibration_prohibition",new["threshold"]["external_recalibration"] is False)
    check("selection:training_not_authorized_started", new["training_authorized"] is False and new["training_started"] is False)
    check("selection:only_C_final_eligible", {c for c, v in new["configurations"].items() if v["eligible_final_configuration"]} == {"V2-C"})
    check("selection:no_C_over_B_requirement", new["selection"]["paired_ablation_report"]["comparative_margins_used_for_retention"] is False and new["selection"]["paired_ablation_report"]["superiority_equivalence_noninferiority_required"] is False)
    adequacy = new["selection"]["adequacy"]
    for field, value in (("lower95_AUROC_strictly_gt", .5), ("lower95_AP_minus_replicate_prevalence_strictly_gt", 0), ("lower95_BrierSkill_strictly_gt", 0)):
        check("selection:absolute_adequacy:" + field, adequacy[field] == value)
    check("selection:seed_instability_limits_preserved", new["selection"]["seed_stability_hard_stop"] == old["selection"]["seed_stability_hard_stop"])
    check("selection:minimum_selection_preserved", new["selection"]["minimum_selection"] == old["selection"]["minimum_selection"])
    check("selection:calibration_component_minimum", new["threshold"]["minimum_control_components"] == 30)
    check("selection:no_test_external_tiebreaker", new["selection"]["test_or_external_tiebreaker"] is False)
    check("selection:no_allele_delta_cutoff", new["candidate_plan"]["allele_delta_cutoff"] is None and new["candidate_plan"]["execute_now"] is False)
    check("selection:common_panels_roles", set(new["common_panels"]["roles"]) == {"selection", "calibration", "test"} and new["common_panels"]["checkpoint_included"] is False)
    return new


def verify_configurations(frame):
    sets = {}
    original = {}
    for model, word in (("enhancer", "Enhancer"), ("h3k27me3", "Silencer")):
        for label in ("positive", "control"):
            p = OLD / f"../04_modeling/trednet/input_training_data/COPD_SevereEmphysema_Lung_{word}_DHS_x2_{label}_1kb.bed"
            ids = []
            with p.open() as f:
                for line in f:
                    if line.strip() and not line.startswith(("#", "track", "browser")):
                        c, a, b, *_ = line.rstrip().split("\t")
                        ids.append(f"{c}:{int(a)}-{int(b)}")
            original[(model, label)] = set(ids)
            check(f"V1_membership:{model}:{label}", len(ids) == len(set(ids)) and set(ids) == set(frame.index[frame[f"v1_{model}_{label}"].eq(1)]), len(ids))
    for config in ("V2-A", "V2-B", "V2-C"):
        for model in MODELS:
            path = DATA / "configurations" / f"{config}_{model}_interval_manifest.tsv.gz"
            x = table(path).set_index("interval_id", drop=False)
            check(f"configuration:{config}:{model}:unique_ids", x.index.is_unique)
            check(f"configuration:{config}:{model}:known_ids", set(x.index) <= set(frame.index))
            check(f"configuration:{config}:{model}:binary_labels", set(x.label) <= {0,1})
            master = frame.loc[x.index]
            for column in x.columns.intersection(master.columns):
                if column == "interval_id":
                    continue
                if pd.api.types.is_numeric_dtype(master[column]):
                    a = pd.to_numeric(x[column], errors="coerce").to_numpy(float)
                    b = master[column].to_numpy(float)
                    equal = np.isclose(a, b, atol=1e-9, rtol=1e-9, equal_nan=True)
                else:
                    equal = x[column].to_numpy() == master[column].to_numpy()
                check(f"configuration:{config}:{model}:master_column:{column}", bool(np.all(equal)), int(np.sum(~equal)), 0)
            positive = set(x.index[x.label.eq(1)])
            controls = set(x.index[x.label.eq(0)])
            universe = set(original[(model, "positive")])
            if config == "V2-C":
                universe &= set(frame.index[frame[f"{model}_same_lobe_peak_support"].eq(1)])
            check(f"configuration:{config}:{model}:positive_universe", positive <= universe)
            if config == "V2-A":
                check(f"configuration:{config}:{model}:unchanged_V1_positive", positive == original[(model,"positive")])
                check(f"configuration:{config}:{model}:unchanged_V1_control", controls == original[(model,"control")])
            else:
                check(f"configuration:{config}:{model}:one_to_one", len(positive) == len(controls), f"{len(positive)}+/{len(controls)}-")
                eligible = set(frame.index[frame[f"{model}_accessible_control_eligible"].eq(1)])
                check(f"configuration:{config}:{model}:controls_eligible", controls <= eligible)
                for part in ("train", "validation", "test"):
                    total = len(universe & set(frame.index[frame.partition.eq(part)]))
                    retained = len(positive & set(frame.index[frame.partition.eq(part)]))
                    ready(f"{config}:{model}:{part}:positive_support", retained/total >= .90, f"{retained}/{total}={retained/total:.12g}", ">=0.90")
                ready(f"{config}:{model}:total_positive_support", len(positive)/len(universe) >= .90, f"{len(positive)}/{len(universe)}={len(positive)/len(universe):.12g}", ">=0.90")
            enriched = master.copy()
            enriched["label"] = x.label
            sets[(config, model)] = enriched
            compare_balance(enriched, config, model)
    for model in MODELS:
        b, c = sets[("V2-B", model)], sets[("V2-C", model)]
        check(f"configuration:{model}:C_not_unchanged_reuse_B_controls", set(b.index[b.label.eq(0)]) != set(c.index[c.label.eq(0)]))
    for config in ("V2-B", "V2-C"):
        for model in MODELS:
            bad = [r for r in BALANCE if r["configuration"] == config and r["model"] == model and r["status"] == "FAIL"]
            ready(f"{config}:{model}:all_balance_gates", not bad, len(bad), 0)
    return sets


def verify_common_panels(sets, specification):
    panel_summary = []
    for model in MODELS:
        panel = table(DATA / "evaluation" / f"{model}_common_challenge_panel.tsv.gz").set_index("interval_id", drop=False)
        expected = sets[("V2-C", model)]
        expected = expected[expected.validation_role.isin(["selection", "calibration", "test"])]
        check(f"panel:{model}:unique_complete_C_task", panel.index.is_unique and set(panel.index) == set(expected.index))
        for col in ("label", "chrom", "partition", "validation_role", "component_id", "sequence_sha256", "canonical_rc_sequence_sha256"):
            check(f"panel:{model}:C_manifest:{col}", bool((panel.loc[expected.index,col] == expected[col]).all()))
        check(f"panel:{model}:no_checkpoint_training_rows", set(panel.validation_role) == {"selection", "calibration", "test"})
        for role, group in panel.groupby("validation_role"):
            positive = int(group.label.eq(1).sum())
            control = int(group.label.eq(0).sum())
            components = int(group.component_id.nunique())
            control_components = int(group.loc[group.label.eq(0),"component_id"].nunique())
            row = {"model": model, "role": role, "positive": positive, "control": control, "components": components, "control_components": control_components}
            panel_summary.append(row)
            if role == "selection":
                minima = specification["selection"]["minimum_selection"]
                ready(f"panel:{model}:selection_minima", positive >= minima["positive"] and control >= minima["control"] and components >= minima["components"], row, minima)
            if role == "calibration":
                minima = specification["threshold"]
                ready(f"panel:{model}:calibration_minima", control >= minima["minimum_controls"] and control_components >= minima["minimum_control_components"], row, {"controls": minima["minimum_controls"], "control_components": minima["minimum_control_components"]})
    return panel_summary


def coordinate_key(identifier):
    chrom, limits = identifier.split(":")
    start, end = map(int, limits.split("-"))
    chromosome = int(chrom[3:]) if chrom[3:].isdigit() else 23 if chrom == "chrX" else 24
    return chromosome, start, end, identifier


def verify_initial_matching(frame, path):
    """Audit every retained initial assignment without reproducing bounded-k logic.

    The recorded selected distance gives a radius guaranteed to include the
    selected control; every eligible unused control within that complete radius
    is tested. Any closer or globally earlier tied control falsifies the record.
    Available-count bookkeeping independently chooses the first nonempty pool in
    the exact-chromosome/lobe -> lobe -> validation-role hierarchy.
    """
    pairs = table(path)
    if "configuration" not in pairs and "config" in pairs:
        pairs = pairs.rename(columns={"config": "configuration"})
    variables = ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001"]
    scales = np.array([.05, .20, .20])
    agg = Aggregate("initial_global_matching")
    audit = {}
    for (config, model), rows in pairs.groupby(["configuration", "model"], sort=True):
        controls_frame = frame[frame[f"{model}_accessible_control_eligible"].eq(1)]
        positive_frame = frame[frame[f"v1_{model}_positive"].eq(1)]
        if config == "V2-C":
            positive_frame = positive_frame[positive_frame[f"{model}_same_lobe_peak_support"].eq(1)]
        agg.add("known_configuration_model", config in {"V2-B","V2-C"} and model in MODELS)
        agg.add("unique_positive", rows.positive_id.is_unique)
        agg.add("unique_control", rows.control_id.is_unique)
        agg.add("full_positive_population", set(rows.positive_id) == set(positive_frame.index))
        observed = rows.set_index("positive_id")
        groups = collections.defaultdict(list)
        keys = {}
        for r in controls_frame[["interval_id", "chrom", "validation_role", "atac_lobes"]].itertuples(index=False):
            tiers = ((0,r.validation_role,r.chrom,r.atac_lobes), (1,r.validation_role,r.atac_lobes), (2,r.validation_role))
            keys[r.interval_id] = tiers
            for key in tiers:
                groups[key].append(r.interval_id)
        available = {key: len(ids) for key, ids in groups.items()}
        trees, used = {}, set()
        ordered = sorted(positive_frame.index, key=lambda i: hashlib.sha256(f"271828|{model}|{config}|{i}".encode()).hexdigest())
        matching = {}
        for order, pid in enumerate(ordered):
            p = frame.loc[pid]
            levels = ((0,p.validation_role,p.chrom,p.atac_lobes), (1,p.validation_role,p.atac_lobes), (2,p.validation_role))
            selected_key = next((key for key in levels if available.get(key,0) > 0), None)
            if selected_key is None:
                agg.add("nonempty_candidate_pool", False, pid)
                continue
            if selected_key not in trees:
                ids = sorted(groups[selected_key], key=coordinate_key)
                values = frame.loc[ids, variables].to_numpy(float)/scales
                trees[selected_key] = (ids, values, cKDTree(values), {i:j for j,i in enumerate(ids)})
            ids, values, tree, lookup = trees[selected_key]
            row = observed.loc[pid]
            control = row.control_id
            agg.add("reported_processing_order", int(row.processing_order) == order, pid)
            agg.add("reported_candidate_priority", int(row.candidate_priority_level) == selected_key[0],pid)
            point = p[variables].to_numpy(float)/scales
            agg.add("control_in_priority_pool", control in lookup, pid + "|" + control)
            agg.add("control_unused", control not in used, pid + "|" + control)
            if control not in lookup:
                continue
            ci = lookup[control]
            distance = float(np.max(np.abs(values[ci]-point)))
            radius = np.nextafter(distance,np.inf) + 1e-12
            nearby = [int(i) for i in tree.query_ball_point(point,radius,p=np.inf) if ids[int(i)] not in used]
            distances = np.max(np.abs(values[nearby]-point),axis=1)
            candidates = sorted(zip(distances, nearby))
            best_distance, best = candidates[0]
            expected = ids[best]
            agg.add("globally_nearest_available_control", distance == float(best_distance), pid)
            agg.add("global_natural_coordinate_tie", expected == control, pid + "|" + control + "|" + expected)
            agg.add("same_partition_and_validation_role", p.partition == frame.loc[control,"partition"] and p.validation_role == frame.loc[control,"validation_role"], pid)
            tied = int(np.sum(distances == best_distance))
            agg.add("reported_initial_distance",np.isclose(float(row.initial_normalized_chebyshev_distance),distance,atol=1e-14,rtol=1e-12),pid)
            agg.add("reported_global_minimum_tie_count",int(row.global_min_distance_tie_count) == tied,pid)
            if tied > 1:
                TIES.append({"stage":"initial", "configuration":config, "model":model,"positive_id":pid,"selected_control":control,"expected_control":expected,"minimum_distance":float(best_distance),"number_globally_tied_available_controls":tied,"status":"PASS" if control == expected else "FAIL"})
            used.add(control)
            matching[pid] = control
            for key in keys[control]:
                available[key] -= 1
        audit[(config,model)] = matching
        print(f"Independent complete initial tie replay: {config} {model}, {len(matching)} assignments",flush=True)
    agg.finish()
    check("initial_matching:all_four_constructions", set(audit) == {(c,m) for c in ("V2-B","V2-C") for m in MODELS})
    return audit


def sample_smd(a, b):
    if len(a) < 2 or len(b) < 2 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Undefined SMD prerequisites: both classes require >=2 finite rows")
    delta = a.mean(axis=0)-b.mean(axis=0)
    scale = np.sqrt((a.var(axis=0,ddof=1)+b.var(axis=0,ddof=1))/2)
    result = np.zeros_like(delta)
    np.divide(delta,scale,out=result,where=scale > 0)
    result[(scale == 0) & (delta != 0)] = np.inf
    return result


def verify_exchange_replay(frame, initial, attempt):
    swap_table = table(attempt / "exchange_ledger.tsv.gz")
    final_table = table(attempt / "control_matching_pairs.tsv.gz")
    round_table = table(attempt / "exchange_rounds.tsv")
    variables = ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001", "blacklist_input_any", "non_acgt_any"]
    agg = Aggregate("exchange_independent_replay")
    total = 0
    for (config,model), original_matching in initial.items():
        matching = dict(original_matching)
        inverse = {control:positive for positive,control in matching.items()}
        selected = set(inverse)
        all_eligible = frame[frame[f"{model}_accessible_control_eligible"].eq(1)]
        for role in sorted({frame.loc[p,"validation_role"] for p in matching}):
            positives = [p for p in matching if frame.loc[p,"validation_role"] == role]
            ids = sorted(all_eligible.index[all_eligible.validation_role.eq(role)],key=coordinate_key)
            lookup = {i:j for j,i in enumerate(ids)}
            values = frame.loc[ids,variables].to_numpy(float)
            pv = frame.loc[positives,variables].to_numpy(float)
            scale = np.std(pv,axis=0,ddof=1)/np.sqrt(2)
            pool_scale = np.std(values,axis=0,ddof=1)/np.sqrt(2)
            scale = np.where(scale == 0,pool_scale,scale)
            constant = scale == 0
            agg.add("constant_feature_feasibility", not np.any(constant & (pv.mean(axis=0) != values.mean(axis=0))), (config,model,role))
            scale[constant] = 1
            normalized = values/scale
            target = pv.mean(axis=0)/scale
            n = len(positives)
            groups = collections.defaultdict(list)
            for j,(chrom,lobes) in enumerate(frame.loc[ids,["chrom","atac_lobes"]].itertuples(index=False,name=None)):
                groups[(chrom,lobes)].append(j)
            local_swaps = swap_table[swap_table.configuration.eq(config) & swap_table.model.eq(model) & swap_table.validation_role.eq(role)]
            local_rounds = round_table[round_table.configuration.eq(config) & round_table.model.eq(model) & round_table.validation_role.eq(role)]
            observed_swaps = list(local_swaps.itertuples(index=False))
            accepted_all = 0
            for round_number in range(201):
                active = np.array([i in selected for i in ids])
                residual = normalized[active].mean(axis=0)-target
                loss = float(residual @ residual)
                smd = sample_smd(pv,values[active])
                reported = local_rounds[local_rounds["round"].eq(round_number)]
                agg.add("round_record_unique_present", len(reported) == 1,(config,model,role,round_number))
                if len(reported) != 1:
                    break
                row = reported.iloc[0]
                agg.add("round_loss", np.isclose(row.loss_before,loss,atol=1e-12,rtol=1e-10),(config,model,role,round_number))
                agg.add("round_maximum_actual_SMD", np.isclose(row.maximum_absolute_SMD_before,np.max(np.abs(smd)),atol=1e-12,rtol=1e-10))
                if np.all(np.abs(smd) <= .10):
                    agg.add("stop_when_role_gates_pass", row.status == "ALL_ROLE_SMD_GATES_PASS")
                    break
                if round_number == 200:
                    agg.add("prespecified_round_cap", row.status == "STOP_MAX_PRESPECIFIED_ROUNDS")
                    break
                projections = normalized @ (2*residual/n)
                proposals = []
                for key in sorted(groups):
                    old = sorted((j for j in groups[key] if active[j]),key=lambda j:(-projections[j],j))
                    new = sorted((j for j in groups[key] if not active[j]),key=lambda j:(projections[j],j))
                    proposals.extend((float(projections[b]-projections[a]),a,b) for a,b in zip(old,new))
                proposals.sort(key=lambda v:(v[0],v[1],v[2]))
                agg.add("complete_global_proposal_count", len(proposals) == int(row.proposals), (config,model,role,round_number))
                accepted, rejected = 0, 0
                for proposal_rank,(derivative,a,b) in enumerate(proposals):
                    new_residual = residual+(normalized[b]-normalized[a])/n
                    next_loss = float(new_residual @ new_residual)
                    if not next_loss < loss-1e-15:
                        rejected += 1
                        continue
                    old_id,new_id = ids[a],ids[b]
                    positive = inverse[old_id]
                    if accepted_all >= len(observed_swaps):
                        agg.add("all_expected_swaps_recorded",False,(config,model,role,round_number))
                    else:
                        record = observed_swaps[accepted_all]
                        agg.add("global_proposal_order_and_ties", record.proposal_rank == proposal_rank and record.round == round_number and record.swap_in_round == accepted and record.positive_id == positive and record.removed_control_id == old_id and record.added_control_id == new_id, f"{config}|{model}|{role}|{round_number}|{proposal_rank}")
                        agg.add("frozen_gradient_derivative", np.isclose(record.frozen_gradient_derivative,derivative,atol=1e-12,rtol=1e-10))
                        agg.add("recorded_quadratic_loss", np.isclose(record.loss_before,loss,atol=1e-12,rtol=1e-10) and np.isclose(record.loss_after,next_loss,atol=1e-12,rtol=1e-10))
                    selected.remove(old_id); selected.add(new_id)
                    del inverse[old_id]; inverse[new_id] = positive
                    matching[positive] = new_id
                    active[a] = False; active[b] = True
                    residual,loss = new_residual,next_loss
                    accepted += 1; accepted_all += 1; total += 1
                    if np.all(np.abs(sample_smd(pv,values[active])) <= .10) or accepted >= 256:
                        break
                agg.add("accepted_and_rejected_round_counts", accepted == int(row.accepted) and rejected == int(row.rejected_nonimproving),(config,model,role,round_number))
                if accepted == 0:
                    break
            agg.add("complete_accepted_swap_ledger", accepted_all == len(observed_swaps),(config,model,role))
        final = final_table[final_table.configuration.eq(config) & final_table.model.eq(model)]
        agg.add("final_mapping_exact_replay", matching == dict(zip(final.positive_id,final.control_id)),(config,model))
        agg.add("final_controls_unique", len(set(matching.values())) == len(matching),(config,model))
        print(f"Independent full exchange replay: {config} {model}",flush=True)
    agg.finish()
    check("exchange:all_accepted_swaps_replayed", total == len(swap_table),total,len(swap_table))
    return total


def verify_synthetic_nearest():
    source = BASE / "scripts/construct_matched_design.py"
    spec = importlib.util.spec_from_file_location("outcome_blind_selector_under_test",source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cases = []
    cases.append(("257_equal_distance_candidates", np.ones((257,3)),np.zeros(257,dtype=bool),np.arange(256,-1,-1),np.zeros(3)))
    p = np.ones((257,3)); used = np.zeros(257,dtype=bool); used[:200] = True
    cases.append(("unused_global_tie_beyond_bounded_k",p,used,np.arange(257),np.zeros(3)))
    p = np.array([[1+5e-13,0,0],[1.,0,0],[-1.,0,0],[0,2,0]])
    cases.append(("near_distance_not_false_tie",p,np.zeros(4,dtype=bool),np.array([3,0,2,1]),np.zeros(3)))
    cases.append(("all_candidates_used",np.zeros((4,3)),np.ones(4,dtype=bool),np.arange(4),np.zeros(3)))
    p = np.array([[0,0,0],[0,0,0],[1,0,0],[0,-1,0]],dtype=float)
    used = np.array([True,False,False,False])
    cases.append(("zero_distance_unused_duplicate",p,used,np.arange(4),np.zeros(3)))
    rng = np.random.default_rng(271828)
    for n in range(20):
        points = rng.integers(-2,3,size=(300,3)).astype(float)
        used = rng.random(300) < .4
        cases.append(("discrete_full_ties_" + str(n),points,used,rng.permutation(300),rng.integers(-1,2,size=3).astype(float)))
    for name,points,used,indices,point in cases:
        candidates = [(float(np.max(np.abs(points[i]-point))),int(i)) for i in indices if not used[i]]
        expected = min(candidates) if candidates else None
        result = module.AvailableNearest(indices,points,used).choose(point)
        passed = result is None if expected is None else result is not None and (result[1],result[0]) == expected and result[2] == sum(d == expected[0] for d,_ in candidates)
        check("synthetic_global_tie:" + name,passed,result,expected)


def verify_attempt_metadata(attempt):
    spec = read_json(BASE / "specification/matching_attempt_001.json")
    manifest = read_json(attempt / "attempt_manifest.json")
    check("matching:specification_recorded_before_execution_hash", manifest["rule_sha256"] == sha(BASE / "specification/matching_attempt_001.json"))
    check("matching:implementation_hash", manifest["implementation_sha256"] == spec["implementation_sha256"] == sha(BASE / "scripts/construct_matched_design.py"))
    expected = {"old_features": OLD / "data/COPD-V2-PREFLIGHT/interval_features.tsv.gz", "old_roles": OLD / "data/COPD-V2-PREFLIGHT/interval_role_assignment.tsv.gz", "train_only_ATAC_sidecar": DATA / "atac_normalized_features.tsv.gz"}
    for name,path in expected.items():
        check("matching:input_hash:" + name,spec["input_sha256"][name] == sha(path))
    check("matching:original_balance_gates_unrelaxed",spec["balance_gates"] == {"absolute_SMD":.10,"absolute_chromosome_proportion_gap":.02,"absolute_ATAC_lobe_signature_proportion_gap":.05,"absolute_each_lobe_support_proportion_gap":.05})
    check("matching:forbidden_execution_flags",all(value is False for value in spec["firewall"].values()))
    for flag in ("training_started","model_inference_performed","benchmark_outcomes_read"):
        check("matching:" + flag,manifest[flag] is False)
    check("matching:common_support_exclusions_empty",len(table(attempt / "positive_common_support_exclusions.tsv.gz")) == 0)
    return manifest


def verify_reported_balance(attempt):
    reference = {(r["configuration"],r["model"],r["partition"],r["variable"]):r for r in BALANCE if r["configuration"] != "V2-A"}
    reported = table(attempt / "covariate_balance.tsv")
    reported = reported[reported.stage.eq("final") & reported.scope_type.eq("partition")]
    agg = Aggregate("reported_balance")
    observed_keys = set()
    for r in reported.itertuples(index=False):
        variable = r.variable.replace("lobe_support=","marginal_lobe=")
        key = (r.configuration,r.model,r.scope,variable)
        observed_keys.add(key)
        expected = reference.get(key)
        agg.add("known_unique_key",expected is not None,key)
        if expected is None:
            continue
        value = expected["SMD"] if r.kind == "SMD" else expected["proportion_gap"]
        agg.add("independent_metric_value",np.isclose(r.value,value,atol=1e-12,rtol=1e-10),key)
        agg.add("unrelaxed_limit",r.limit == expected["limit"],key)
        agg.add("status",r.status == expected["status"],key)
        agg.add("class_counts",r.positive_n == expected["positive_n"] and r.control_n == expected["control_n"],key)
    check("reported_balance:all_partition_metrics_present",observed_keys == set(reference),len(observed_keys),len(reference))
    agg.finish()


def write_outputs(label, started, summary):
    outputs = {f"validation_{label}_checks.tsv": pd.DataFrame(CHECKS),
               f"validation_{label}_readiness.tsv": pd.DataFrame(READINESS),
               f"validation_{label}_balance.tsv": pd.DataFrame(BALANCE),
               f"validation_{label}_ties.tsv": pd.DataFrame(TIES)}
    for name in outputs:
        if (PROV / name).exists():
            raise RuntimeError("Refusing to overwrite previous validation: " + name)
    target = PROV / f"validation_{label}_manifest.json"
    if target.exists():
        raise RuntimeError("Refusing to overwrite previous validation: " + str(target))
    PROV.mkdir(parents=True, exist_ok=True)
    for name, data in outputs.items():
        data.to_csv(PROV / name, sep="\t", index=False)
    result = {"version": "pretraining-1.1", "validator_sha256": sha(__file__), "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "elapsed_seconds": time.time()-started, "artifact_checks": len(CHECKS), "artifact_failures": sum(r["status"] == "FAIL" for r in CHECKS), "readiness_failures": sum(r["status"] == "BLOCKED" for r in READINESS), "training_authorized": False, "GPU_training_or_model_inference_performed": False, "benchmark_outcomes_parsed": False, "phase_I_extraction_performed": False, "summary": summary, "outputs": [{"path": str((PROV / name).relative_to(ROOT)), "bytes": (PROV / name).stat().st_size, "sha256": sha(PROV / name)} for name in outputs]}
    with target.open("x") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps(result, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label",default="final",help="Unique output label; existing outputs cannot be replaced")
    args = parser.parse_args()
    if not args.label.replace("_","").replace("-","").isalnum():
        raise ValueError("Invalid validation output label")
    started = time.time()
    attempt = BASE / "attempts/001_full_population"
    print("Verifying immutable 1.0 artifacts, consumed inputs and opaque benchmark hashes",flush=True)
    historical_integrity()
    manifest = verify_attempt_metadata(attempt)
    specification = verify_selection_specification()
    verify_synthetic_nearest()
    print("Reading master features and independently reconstructing train-only ATAC mapping",flush=True)
    frame = verify_normalization(load_master())
    print("Reconstructing all master components and fixed validation roles",flush=True)
    verify_global_roles(frame)
    initial = verify_initial_matching(frame,attempt / "initial_pairs.tsv.gz")
    exchanges = verify_exchange_replay(frame,initial,attempt)
    print("Verifying all final manifests, support, balance and common panels",flush=True)
    sets = verify_configurations(frame)
    panels = verify_common_panels(sets,specification)
    verify_reported_balance(attempt)
    final_pairs = table(attempt / "control_matching_pairs.tsv.gz")
    for (config,model), x in sets.items():
        if config == "V2-A":
            continue
        pairs = final_pairs[final_pairs.configuration.eq(config) & final_pairs.model.eq(model)]
        check(f"final_pairs:{config}:{model}:complete_manifest_membership",set(pairs.positive_id) == set(x.index[x.label.eq(1)]) and set(pairs.control_id) == set(x.index[x.label.eq(0)]))
    print("Independently validating original source annotations and every selected FASTA sequence",flush=True)
    summary = verify_selected_sources(frame,sets)
    summary.update({"initial_assignments_global_tie_audited":sum(map(len,initial.values())),"accepted_exchanges_fully_replayed":exchanges,"initial_assignments_with_exact_global_ties":len(TIES),"common_C_panels":panels,"B_C_partition_balance_failures":sum(r["status"] == "FAIL" and r["configuration"] != "V2-A" for r in BALANCE)})
    check("matching:reported_balance_failure_total",summary["B_C_partition_balance_failures"] == manifest["balance_failures"],summary["B_C_partition_balance_failures"],manifest["balance_failures"])
    check("execution:no_accelerator_framework_imported",not any(name == "tensorflow" or name.startswith("tensorflow.") or name == "torch" or name.startswith("torch.") or name == "keras" or name.startswith("keras.") for name in __import__("sys").modules))
    ready("artifact_correctness",not any(r["status"] == "FAIL" for r in CHECKS),sum(r["status"] == "FAIL" for r in CHECKS),0)
    summary["readiness_interpretation"] = "Construction validation only. Freeze all payloads/specifications before READY FOR INVESTIGATOR TRAINING REVIEW. No training authorization; actual-network adequacy and invariance remain future mandatory gates."
    write_outputs(args.label,started,summary)
    return 1 if any(r["status"] == "FAIL" for r in CHECKS) else 0


if __name__ == "__main__":
    raise SystemExit(main())
