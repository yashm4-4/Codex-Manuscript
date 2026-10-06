#!/usr/bin/env python3
"""CPU-only, no-model pretraining feature construction from frozen V1 inputs.

This stage does not select/match controls, train, extract learned embeddings, or
inspect external benchmark contents. All new outputs are confined to PREFLIGHT.
Existing outputs are never silently replaced. Coordinates are BED half-open.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[4]
S3 = ROOT / "diseases/COPD/03_regulatory_landscape"
S4 = ROOT / "diseases/COPD/04_modeling"
V2 = ROOT / "diseases/COPD/07_gap_closure"
DATA = V2 / "data/COPD-V2-PREFLIGHT"
PROV = V2 / "provenance"
RESULTS = V2 / "results"
PREFIX = "COPD-V2-PREFLIGHT"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_SET = set(CHROMS)
LOBES = ["lower_left", "lower_right", "upper_right"]
MODELS = {"enhancer": "H3K27ac", "h3k27me3": "H3K27me3"}
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path):
    return str(Path(path).relative_to(ROOT))


def write_json(path, obj):
    with Path(path).open("x") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")


def write_tsv(path, rows, fields):
    with Path(path).open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()


def bed(path):
    with opener(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] in CHROM_SET:
                yield p[0], int(p[1]), int(p[2])


def split(chrom):
    return "validation" if chrom == "chr7" else "test" if chrom in {"chr8", "chr9"} else "train"


class MergedIndex:
    """Disjoint merged intervals, logarithmic overlap and covered-length queries."""

    def __init__(self, intervals):
        groups = collections.defaultdict(list)
        for c, s, e, *_ in intervals:
            if e > s:
                groups[c].append((s, e))
        self.data = {}
        self.n_intervals = 0
        for chrom, ivs in groups.items():
            ivs.sort()
            merged = []
            for s, e in ivs:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            starts, ends, prefix = [], [], [0]
            for s, e in merged:
                starts.append(s)
                ends.append(e)
                prefix.append(prefix[-1] + e - s)
            self.data[chrom] = (starts, ends, prefix)
            self.n_intervals += len(merged)

    def hit(self, chrom, start, end):
        if chrom not in self.data or end <= start:
            return False
        starts, ends, _ = self.data[chrom]
        i = bisect.bisect_left(starts, end)
        return i > 0 and ends[i - 1] > start

    def bp(self, chrom, start, end):
        if chrom not in self.data or end <= start:
            return 0
        starts, ends, prefix = self.data[chrom]
        left = bisect.bisect_right(ends, start)
        right = bisect.bisect_left(starts, end)
        if left >= right:
            return 0
        return prefix[right] - prefix[left] - max(0, start - starts[left]) - max(0, ends[right - 1] - end)


def promoter_intervals(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if p[0] in CHROM_SET and p[2] == "gene":
                tss = int(p[3]) - 1 if p[6] == "+" else int(p[4]) - 1
                yield p[0], max(0, tss - 2000), tss + 2000


def repeat_intervals(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.split("\t")
            if len(p) >= 8 and p[5] in CHROM_SET:
                yield p[5], int(p[6]), int(p[7])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["features"], default="features")
    args = parser.parse_args()
    started = time.time()
    output = DATA / "interval_features.tsv.gz"
    manifest_path = PROV / f"{PREFIX}_feature_manifest.json"
    rules_path = PROV / f"{PREFIX}_feature_rules.json"
    if any(p.exists() for p in (output, manifest_path, rules_path, DATA / "ATAC_peak_evidence.tsv.gz")):
        raise RuntimeError("Refusing to overwrite an existing PREFLIGHT feature run; inspect/version explicitly.")
    DATA.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    source_manifest = S3 / "data/encode_peaks/manifest.tsv"
    with source_manifest.open() as f:
        sources = list(csv.DictReader(f, delimiter="\t"))
    assert len(sources) == 9 and {r["donor_accessions"] for r in sources} == {"ENCDO520EJG"}
    fasta = S4 / "trednet/fasta/hg38.fa"
    gtf = ROOT / "data/gencode/gencode.v50.annotation.gtf.gz"
    blacklist = ROOT / "data/blacklist/hg38-blacklist.v2.bed.gz"
    repeats = ROOT / "data/repeats/rmsk.txt.gz"
    beds = {}
    for model, eid_word in (("enhancer", "Enhancer"), ("h3k27me3", "Silencer")):
        for label in ("positive", "control"):
            beds[(model, label)] = S4 / f"trednet/input_training_data/COPD_SevereEmphysema_Lung_{eid_word}_DHS_x2_{label}_1kb.bed"
    consumed = [source_manifest, S3 / "data/encode_lung_experiments.tsv", fasta, Path(str(fasta) + ".fai"), gtf, blacklist, repeats, Path(__file__)]
    consumed += [Path(r["path"]) for r in sources] + list(beds.values())
    consumed += [S4 / "trednet/model_phase_I/phase_one_weights.h5", S4 / "trednet/TREDNet_v2_seeded.py", S4 / "trednet/TREDNet_v2_seeded_v1.py"]
    consumed += [S4 / f"trednet/models_output/COPD_SevereEmphysema_Lung_{word}_DHS_x2/phase_two_model.keras" for word in ("Enhancer", "Silencer")]
    consumed += [ROOT / "models/TREDNET_v2/make_input_training_data.py", ROOT / "models/TREDNET_v2/input_training_data/control/allDHS.merge.nonPromoterExon"]
    rules = {
        "module": PREFIX, "stage": args.stage, "initialized_utc": utc(), "no_training_or_inference": True,
        "seeds": [104729, 130363, 155921], "matching_seed": 271828,
        "genome": "GRCh38", "input_rule": "[core_start-501,core_end+500),2001bp; preserve V1 expansion",
        "encoding": "Uppercase A/C/G/T one-hot; every other symbol zero; no N-based exclusion",
        "peak_overlap": ">=1bp original ATAC peak overlap; preserve exact summit-window identity",
        "same_lobe": "Union exact summit windows supported by ATAC and relevant mark in same lobe; no duplicate weighting",
        "eligible_control": "ATAC summit anchor; original V1 core filters; no relevant all-lobe mark in ANY anchor peak or full2001bp input; no overlap with full2001bp V1 positive inputs for that model; sequence available",
        "blacklist_flanks": "Audit full input; core-only filter retained to avoid an unprespecified additional exclusion",
        "promoter": "GENCODE v50 gene TSS +/-2000bp, strand aware, exactly V1 rule",
        "signal_percentile": "Within each source ATAC file, average rank/n among canonical rows; maximum rank over identical anchor windows",
        "duplicate_policy": "Audit raw exact sequence and model-encoded canonical forward/RC identity across fixed chromosome partitions; normalize all nonACGT symbols to N before canonical identity, matching V1 zero encoding; do not silently reassign/drop frozen V1 examples",
        "external_benchmark_firewall": "No benchmark contents inspected; tracked benchmark files may be opaque byte-hashed for immutability only",
        "V1_and_prior_modules_immutable": True,
    }
    write_json(rules_path, rules)
    print("Freezing consumed-input hashes before interval construction", flush=True)
    hash_rows = []
    for p in consumed:
        hash_rows.append({"path": rel(p), "resolved_path": str(p.resolve()), "bytes": p.stat().st_size, "sha256": sha(p)})
    write_tsv(PROV / f"{PREFIX}_consumed_input_hashes.tsv", hash_rows, ["path", "resolved_path", "bytes", "sha256"])
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    protected = []
    for name in tracked:
        if not name:
            continue
        prior = name.startswith("diseases/COPD/") and "/07_gap_closure/" not in name
        prior = prior or any(t in name for t in ("COPD-V2-PHENO", "COPD-V2-RC", "COPD-V2-BENCH"))
        if prior:
            p = ROOT / name
            if p.is_file():
                protected.append({"path": name, "bytes": p.stat().st_size, "sha256": sha(p), "access_mode": "opaque_bytes_only"})
    write_tsv(PROV / f"{PREFIX}_protected_tracked_hashes_before.tsv", protected, ["path", "bytes", "sha256", "access_mode"])
    print(f"Frozen {len(hash_rows)} consumed and {len(protected)} prior tracked files", flush=True)
    source_inventory = []
    peak_sets = {}
    source_by_key = {}
    for row in sources:
        key = (row["lobe"], row["mark"])
        path = Path(row["path"])
        rows = []
        with gzip.open(path, "rt") as f:
            for lineno, line in enumerate(f, 1):
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if p[0] not in CHROM_SET:
                    continue
                start, end, offset = int(p[1]), int(p[2]), int(p[9])
                summit = start + offset if offset >= 0 else (start + end) // 2
                rows.append((p[0], start, end, summit, float(p[6]), f"{row['file_accession']}:{lineno}"))
        peak_sets[key] = rows
        source_by_key[key] = row
        source_inventory.append({**{k: row[k] for k in ("lobe", "assay", "mark", "experiment", "file_accession", "biosample_accessions", "donor_accessions", "assembly", "output_type", "biological_replicates", "technical_replicates", "sha256")}, "path": rel(path), "canonical_rows": len(rows)})
    write_tsv(PROV / f"{PREFIX}_source_inventory.tsv", source_inventory, list(source_inventory[0]))
    mark_ix = {(l, m): MergedIndex(peak_sets[(l, m)]) for l in LOBES for m in MODELS.values()}
    allmark_ix = {model: MergedIndex(p for l in LOBES for p in peak_sets[(l, mark)]) for model, mark in MODELS.items()}
    atac_ix = MergedIndex(p for l in LOBES for p in peak_sets[(l, "accessibility")])
    black_ix = MergedIndex(bed(blacklist))
    promoter_ix = MergedIndex(promoter_intervals(gtf))
    print("Building merged RepeatMasker coverage index", flush=True)
    repeat_ix = MergedIndex(repeat_intervals(repeats))
    membership = {(model, label): set(bed(path)) for (model, label), path in beds.items()}
    positive_input_ix = {model: MergedIndex((c, s - 501, e + 500) for c, s, e in membership[(model, "positive")]) for model in MODELS}
    anchors = collections.defaultdict(list)
    # Record tuple: peak_id,file_id,lobe,signal,percentile,enhancer_lobe_mask,me3_lobe_mask.
    peak_fields = ["peak_id", "file_accession", "experiment", "lobe", "chrom", "peak_start", "peak_end", "summit", "core_start", "core_end", "signal_value", "signal_percentile", "enhancer_all_lobe_support", "enhancer_same_lobe_support", "enhancer_histone_lobes", "enhancer_histone_file_ids", "h3k27me3_all_lobe_support", "h3k27me3_same_lobe_support", "h3k27me3_histone_lobes", "h3k27me3_histone_file_ids"]
    with gzip.open(DATA / "ATAC_peak_evidence.tsv.gz", "xt", compresslevel=6) as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(peak_fields)
        for lobe in LOBES:
            source = source_by_key[(lobe, "accessibility")]
            rows = peak_sets[(lobe, "accessibility")]
            values = np.array([p[4] for p in rows])
            sorted_values = np.sort(values)
            percentiles = (np.searchsorted(sorted_values, values, side="left") + np.searchsorted(sorted_values, values, side="right") + 1) / (2 * len(rows))
            for peak, percentile in zip(rows, percentiles):
                c, s, e, summit, signal, pid = peak
                masks = []
                support = []
                for model, mark in MODELS.items():
                    support_lobes = [l for l in LOBES if mark_ix[(l, mark)].hit(c, s, e)]
                    mask = sum(1 << LOBES.index(l) for l in support_lobes)
                    masks.append(mask)
                    support += [int(bool(mask)), int(lobe in support_lobes), ";".join(support_lobes), ";".join(source_by_key[(l, mark)]["file_accession"] for l in support_lobes)]
                anchors[(c, summit - 500, summit + 500)].append((pid, source["file_accession"], lobe, signal, float(percentile), *masks))
                w.writerow([pid, source["file_accession"], source["experiment"], lobe, c, s, e, summit, summit - 500, summit + 500, signal, f"{percentile:.10g}", *support])
    all_windows = set(anchors)
    for vals in membership.values():
        all_windows.update(vals)
    keys = sorted(all_windows, key=lambda k: (CHROMS.index(k[0]), k[1], k[2]))
    del all_windows, peak_sets
    print(f"Master union: {len(keys)} intervals; {len(anchors)} distinct ATAC summit anchors", flush=True)
    structural_counts = collections.Counter()
    structural_eligible = {m: set() for m in MODELS}
    for key, records in anchors.items():
        c, s, e = key
        if s < 501 or black_ix.hit(c, s, e) or promoter_ix.hit(c, s, e):
            continue
        for mi, model in enumerate(MODELS):
            if any(r[5 + mi] for r in records):
                continue
            if allmark_ix[model].hit(c, s - 501, e + 500) or positive_input_ix[model].hit(c, s - 501, e + 500):
                continue
            structural_eligible[model].add(key)
            structural_counts[(model, split(c))] += 1
    structural_rows = [{"model": m, "partition": p, "structural_eligible_before_sequence_qc": n} for (m, p), n in sorted(structural_counts.items())]
    write_tsv(RESULTS / f"{PREFIX}_structural_control_pool_counts.tsv", structural_rows, ["model", "partition", "structural_eligible_before_sequence_qc"])
    print("STRUCTURAL_CONTROL_POOLS " + json.dumps(structural_rows), flush=True)
    fields = ["interval_id", "chrom", "core_start", "core_end", "input_start", "input_end", "partition"]
    fields += [f"v1_{m}_{label}" for m in MODELS for label in ("positive", "control")]
    fields += ["atac_anchor_count", "atac_peak_ids", "atac_file_ids", "atac_lobes", "atac_signal_raw_max", "atac_signal_percentile_max"]
    for m in MODELS:
        fields += [f"{m}_all_lobe_peak_support", f"{m}_same_lobe_peak_support", f"{m}_all_lobe_histone_file_ids", f"{m}_same_lobe_histone_file_ids", f"{m}_same_lobe_support_lobes"]
    fields += ["atac_overlap_1kb", "atac_overlap_2001", "atac_overlap_summit"]
    for m in MODELS:
        fields += [f"{m}_mark_overlap_1kb", f"{m}_mark_overlap_2001", f"{m}_v1_positive_input_overlap"]
    fields += ["blacklist_bp_1kb", "blacklist_bp_2001", "promoter_bp_1kb", "promoter_bp_2001", "repeat_fraction_2001", "sequence_available", "sequence_length", "gc_fraction", "non_acgt_fraction", "sequence_sha256", "canonical_rc_sequence_sha256"]
    fields += [f"{m}_accessible_control_eligible" for m in MODELS]
    write_json(PROV / f"{PREFIX}_feature_schema.json", {"columns": fields, "binary_fields": "0/1", "missing_values": "empty string", "fraction_denominator": "full2001bp including ambiguous bases", "canonical_rc_sequence_sha256": "SHA256(min(normalize(sequence),RC(normalize(sequence)))); normalize maps every nonACGT to N to match model zero encoding", "no_model_features": True})
    counts = collections.Counter()
    sequence_first = {}
    duplicate_groups = collections.defaultdict(list)
    exact_first = {}
    exact_duplicate_groups = collections.defaultdict(list)
    generated = 0
    available = 0
    current_chrom = None
    reference = pysam.FastaFile(str(fasta))
    lengths = dict(zip(reference.references, reference.lengths))
    with gzip.open(output, "xt", compresslevel=6) as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for c, s, e in keys:
            if c != current_chrom:
                if current_chrom is not None:
                    print(f"Completed {current_chrom}; rows={generated}; elapsed_seconds={time.time()-started:.1f}", flush=True)
                current_chrom = c
            key = (c, s, e); a, b = s - 501, e + 500; part = split(c)
            rid = f"{c}:{s}-{e}"
            records = anchors.get(key, [])
            row = {"interval_id": rid, "chrom": c, "core_start": s, "core_end": e, "input_start": a, "input_end": b, "partition": part}
            for model in MODELS:
                for label in ("positive", "control"):
                    row[f"v1_{model}_{label}"] = int(key in membership[(model, label)])
                    if row[f"v1_{model}_{label}"]:
                        counts[(model, part, f"V1_{label}")] += 1
            row.update(atac_anchor_count=len(records), atac_peak_ids=";".join(r[0] for r in records), atac_file_ids=";".join(sorted({r[1] for r in records})), atac_lobes=";".join(sorted({r[2] for r in records})), atac_signal_raw_max=max((r[3] for r in records), default=""), atac_signal_percentile_max=f"{max(r[4] for r in records):.10g}" if records else "")
            for mi, (model, mark) in enumerate(MODELS.items()):
                all_lobes = {l for r in records for li, l in enumerate(LOBES) if r[5 + mi] & (1 << li)}
                same_lobes = {r[2] for r in records if r[5 + mi] & (1 << LOBES.index(r[2]))}
                row.update({f"{model}_all_lobe_peak_support": int(bool(all_lobes)), f"{model}_same_lobe_peak_support": int(bool(same_lobes)), f"{model}_all_lobe_histone_file_ids": ";".join(source_by_key[(l, mark)]["file_accession"] for l in sorted(all_lobes)), f"{model}_same_lobe_histone_file_ids": ";".join(source_by_key[(l, mark)]["file_accession"] for l in sorted(same_lobes)), f"{model}_same_lobe_support_lobes": ";".join(sorted(same_lobes)), f"{model}_mark_overlap_1kb": int(allmark_ix[model].hit(c, s, e)), f"{model}_mark_overlap_2001": int(allmark_ix[model].hit(c, a, b)), f"{model}_v1_positive_input_overlap": int(positive_input_ix[model].hit(c, a, b))})
                if row[f"v1_{model}_positive"]:
                    counts[(model, part, "V1_positive_same_lobe" if same_lobes else "V1_positive_cross_lobe_only")] += 1
            row.update(atac_overlap_1kb=int(atac_ix.hit(c, s, e)), atac_overlap_2001=int(atac_ix.hit(c, a, b)), atac_overlap_summit=int(atac_ix.hit(c, s + 500, s + 501)), blacklist_bp_1kb=black_ix.bp(c, s, e), blacklist_bp_2001=black_ix.bp(c, a, b), promoter_bp_1kb=promoter_ix.bp(c, s, e), promoter_bp_2001=promoter_ix.bp(c, a, b), repeat_fraction_2001=f"{repeat_ix.bp(c, a, b)/2001:.10g}")
            rc = c if c in lengths else c.removeprefix("chr")
            seq = reference.fetch(rc, a, b).upper() if rc in lengths and a >= 0 and b <= lengths[rc] else ""
            ok = len(seq) == 2001
            row.update(sequence_available=int(ok), sequence_length=len(seq), gc_fraction="", non_acgt_fraction="", sequence_sha256="", canonical_rc_sequence_sha256="")
            if ok:
                available += 1
                gc = seq.count("G") + seq.count("C")
                acgt = gc + seq.count("A") + seq.count("T")
                sqhash = hashlib.sha256(seq.encode()).hexdigest()
                normalized = seq if acgt == 2001 else "".join(base if base in "ACGT" else "N" for base in seq)
                canonical = min(normalized, normalized.translate(COMP)[::-1])
                canhash = hashlib.sha256(canonical.encode()).hexdigest()
                row.update(gc_fraction=f"{gc/2001:.10g}", non_acgt_fraction=f"{1-acgt/2001:.10g}", sequence_sha256=sqhash, canonical_rc_sequence_sha256=canhash)
                info = (rid, part, c)
                for digest, first, groups in ((sqhash, exact_first, exact_duplicate_groups), (canhash, sequence_first, duplicate_groups)):
                    if digest in first:
                        if digest not in groups:
                            groups[digest].append(first[digest])
                        groups[digest].append(info)
                    else:
                        first[digest] = info
            for model in MODELS:
                eligible = int(ok and key in structural_eligible[model])
                row[f"{model}_accessible_control_eligible"] = eligible
                if eligible:
                    counts[(model, part, "accessible_control_eligible")] += 1
                if row[f"v1_{model}_positive"] or row[f"v1_{model}_control"]:
                    if not ok:
                        counts[(model, part, "V1_sequence_unavailable")] += 1
            counts[("all", part, "union_intervals")] += 1
            counts[("all", part, "sequence_available" if ok else "sequence_unavailable")] += 1
            writer.writerow(row)
            generated += 1
    reference.close()
    print(f"Completed {current_chrom}; rows={generated}; elapsed_seconds={time.time()-started:.1f}", flush=True)
    count_rows = [{"model": m, "partition": p, "metric": metric, "count": n} for (m, p, metric), n in sorted(counts.items())]
    write_tsv(RESULTS / f"{PREFIX}_feature_construction_counts.tsv", count_rows, ["model", "partition", "metric", "count"])
    dup_rows = []
    for kind, groups in (("exact_sequence", exact_duplicate_groups), ("canonical_forward_RC", duplicate_groups)):
        for digest, group in sorted(groups.items()):
            partitions = sorted({x[1] for x in group})
            dup_rows.append({"identity_type": kind, "sha256": digest, "n_intervals": len(group), "n_partitions": len(partitions), "cross_partition": int(len(partitions) > 1), "partitions": ";".join(partitions), "interval_ids": ";".join(x[0] for x in group)})
    write_tsv(RESULTS / f"{PREFIX}_sequence_duplicate_audit.tsv", dup_rows, ["identity_type", "sha256", "n_intervals", "n_partitions", "cross_partition", "partitions", "interval_ids"])
    files = [output, DATA / "ATAC_peak_evidence.tsv.gz", rules_path, PROV / f"{PREFIX}_consumed_input_hashes.tsv", PROV / f"{PREFIX}_protected_tracked_hashes_before.tsv", PROV / f"{PREFIX}_source_inventory.tsv", PROV / f"{PREFIX}_feature_schema.json", RESULTS / f"{PREFIX}_structural_control_pool_counts.tsv", RESULTS / f"{PREFIX}_feature_construction_counts.tsv", RESULTS / f"{PREFIX}_sequence_duplicate_audit.tsv"]
    manifest = {"module": PREFIX, "stage": "features", "status": "COMPLETE_CPU_ONLY_NO_TRAINING", "completed_utc": utc(), "elapsed_seconds": time.time() - started, "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, "cpu_workers": 1, "available_cpu_count": os.cpu_count(), "interval_count": generated, "sequence_available_count": available, "distinct_ATAC_anchor_windows": len(anchors), "duplicate_groups": len(dup_rows), "cross_partition_duplicate_groups": sum(r["cross_partition"] for r in dup_rows), "training_or_inference": False, "benchmark_contents_inspected": False, "source_files": len(sources), "consumed_inputs": len(hash_rows), "protected_tracked_files": len(protected), "files": [{"path": rel(p), "bytes": p.stat().st_size, "sha256": sha(p)} for p in files], "feature_columns": fields, "construction_counts": count_rows}
    write_json(manifest_path, manifest)
    print(json.dumps({k: manifest[k] for k in ("status", "interval_count", "sequence_available_count", "elapsed_seconds", "peak_rss_kib", "cross_partition_duplicate_groups")}), flush=True)
    print("FINAL_CONTROL_POOLS " + json.dumps([r for r in count_rows if r["metric"] == "accessible_control_eligible"]), flush=True)


if __name__ == "__main__":
    main()
