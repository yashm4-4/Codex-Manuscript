#!/usr/bin/env python3
"""Independent final verification of the frozen-model COPD V2 RC diagnostic.

This validator never imports either scoring or analysis implementation. It reads
the frozen V1 sources and raw RC scores, reconstructs the decisions and summary
statistics, and writes only final V2 QC/checksum artifacts. Run only after all
analysis, report, figures, and register updates are complete.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


V2 = Path(__file__).resolve().parents[1]
ROOT = V2.parents[2]
COPD = V2.parent
PREFIX = "COPD-V2-RC"
SPEC_SHA = "09d3656f3e06c24fa5402c28b976d6be06fff91604ad0b7b29c3150120d76ac1"
MODELS = ("enhancer", "silencer")
THRESHOLDS = {"enhancer": 0.643623, "silencer": 0.58505}
DELTA = {
    "enhancer": {"SNV": 0.05706318769999998, "indel_or_complex": 0.04936093850000001},
    "silencer": {"SNV": 0.028802613899999996, "indel_or_complex": 0.026183359500000003},
}
ID = "candidate_record_id"
CHECKS = []
HASH_CACHE = {}


def check(name, passed, detail=""):
    CHECKS.append({"check": name, "status": "PASS" if bool(passed) else "FAIL", "detail": str(detail)})
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def sha(path):
    path = Path(path)
    info = path.stat()
    key = (str(path.resolve()), info.st_size, info.st_mtime_ns, info.st_ctime_ns)
    if key not in HASH_CACHE:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1048576), b""):
                digest.update(block)
        HASH_CACHE[key] = digest.hexdigest()
    return HASH_CACHE[key]


def read(path):
    return pd.read_csv(path, sep="\t", low_memory=False, float_precision="round_trip", na_values=["not_evaluable"])


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def result(suffix):
    return V2 / "results" / f"{PREFIX}-{suffix}"


def provenance(suffix):
    return V2 / "provenance" / f"{PREFIX}_{suffix}"


def boolean(values):
    series = pd.Series(values)
    if series.dtype == bool:
        return series.to_numpy()
    converted = series.astype(str).str.lower()
    invalid = ~converted.isin(["true", "false", "1", "0"])
    if invalid.any():
        raise ValueError(f"Invalid boolean serialization: {series[invalid].head().tolist()}")
    return converted.isin(["true", "1"]).to_numpy()


def equal_numeric(actual, expected, atol=2e-12):
    return np.allclose(np.asarray(actual, float), np.asarray(expected, float), rtol=1e-10, atol=atol, equal_nan=True)


def correlation(a, b, spearman=False):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    if spearman:
        a, b = rankdata(a, method="average"), rankdata(b, method="average")
    return float(np.corrcoef(a, b)[0, 1])


def describe(values):
    values = np.asarray(values, float)
    if not len(values):
        return {key: np.nan for key in ("mean", "median", "p90", "p95", "p99", "max")}
    return {"mean": float(values.mean()), "median": float(np.median(values)),
            "p90": float(np.quantile(values, .90)), "p95": float(np.quantile(values, .95)),
            "p99": float(np.quantile(values, .99)), "max": float(values.max())}


def rate(numerator, denominator):
    return float(numerator / denominator) if denominator else np.nan


def verify_array(frame, field, expected, numeric=False):
    check(f"Comparison column exists: {field}", field in frame, len(frame))
    observed = frame[field].to_numpy()
    okay = equal_numeric(observed, expected) if numeric else np.array_equal(observed.astype(str), np.asarray(expected).astype(str))
    check(f"Independent per-record reconstruction: {field}", okay, f"n={len(frame)}")


def construct_independent_data():
    original = read(COPD / "04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz")
    metadata = read(COPD / "04_modeling/results/COPD-S4-R004_prioritized_candidates.tsv.gz")
    rc = read(result("R001_reverse_complement_scores.tsv.gz"))
    check("Authoritative frozen R003 unique complete universe", len(original) == 15303 and original[ID].is_unique)
    check("RC complete universe and authoritative source order", len(rc) == 15303 and rc[ID].is_unique and rc[ID].tolist() == original[ID].tolist())
    check("Original R004 includes all 15389 unique records", len(metadata) == 15389 and metadata[ID].is_unique)
    data = metadata.set_index(ID).loc[original[ID]].reset_index()
    eligible = ~boolean(data["encode_blacklist"])
    classes = data["variant_class_group"].to_numpy()
    check("Preserved eligibility/class denominators", int(eligible.sum()) == 15283 and int((eligible & (classes == "SNV")).sum()) == 13747 and int((eligible & (classes == "indel_or_complex")).sum()) == 1536)
    expected = {ID: original[ID].to_numpy(), "eligible": eligible, "classes": classes}
    for model in MODELS:
        cutoff = np.array([DELTA[model][v] for v in classes])
        expected[f"{model}_cutoff"] = cutoff
        epsilon = np.maximum(2e-5, .1 * cutoff)
        for orientation, source in (("forward", original), ("rc", rc)):
            ref, alt = source[f"{model}_ref_score"].to_numpy(float), source[f"{model}_alt_score"].to_numpy(float)
            check(f"{model}/{orientation}: finite scores in [0,1]", np.isfinite(ref).all() and np.isfinite(alt).all() and ((ref >= 0) & (ref <= 1) & (alt >= 0) & (alt <= 1)).all())
            delta = alt - ref
            serialized = source[f"{model}_delta_alt_minus_ref"].to_numpy(float)
            check(f"{model}/{orientation}: serialized delta finite", np.isfinite(serialized).all())
            if orientation == "rc":
                check(f"{model}: raw RC serialized delta equals ALT minus REF", equal_numeric(serialized, delta, atol=1e-15))
            expected[f"{model}_{orientation}_serialized_delta"] = serialized
            expected[f"{model}_{orientation}_serialization_delta_rounding"] = delta - serialized
            region = np.maximum(ref, alt)
            call = eligible & (region >= THRESHOLDS[model]) & (np.abs(delta) >= cutoff)
            expected.update({f"{model}_{orientation}_ref_score": ref,
                             f"{model}_{orientation}_alt_score": alt,
                             f"{model}_{orientation}_delta": delta,
                             f"{model}_{orientation}_region_score": region,
                             f"{model}_{orientation}_abs_delta": np.abs(delta),
                             f"{model}_{orientation}_call": call,
                             f"{model}_{orientation}_meaningful_sign": np.where(delta > epsilon, 1, np.where(delta < -epsilon, -1, 0))})
        fd, rd = expected[f"{model}_forward_delta"], expected[f"{model}_rc_delta"]
        region = expected[f"{model}_forward_region_score"]
        expected[f"{model}_abs_delta_disagreement"] = np.abs(rd - fd)
        expected[f"{model}_normalized_delta_disagreement"] = np.abs(rd - fd) / cutoff
        expected[f"{model}_forward_well_separated"] = eligible & (region >= THRESHOLDS[model] + .05) & (np.abs(fd) >= 1.5 * cutoff)
        expected[f"{model}_borderline"] = (np.abs(region - THRESHOLDS[model]) <= .02) | (np.abs(np.abs(fd) / cutoff - 1) <= .10)
        check(f"{model}: independent forward decisions reproduce every V1 R004 call", np.array_equal(expected[f"{model}_forward_call"], boolean(data[f"predicted_causal_{model}"])))
    for orientation in ("forward", "rc"):
        e, s = [expected[f"{m}_{orientation}_call"] for m in MODELS]
        expected[f"{orientation}_union_call"] = e | s
        expected[f"{orientation}_model_context"] = np.where(e & s, "both", np.where(e, "enhancer_only", np.where(s, "H3K27me3_only", "neither")))
    e, s = [expected[f"{m}_forward_call"] for m in MODELS]
    check("Frozen forward model/union counts", (int(e.sum()), int(s.sum()), int((e&s).sum()), int((e|s).sum())) == (175,199,37,337))
    return original, metadata, rc, data, expected



def verify_original_thresholds():
    frozen = read(COPD / "04_modeling/results/COPD-S4-R004_delta_thresholds.tsv")
    check("Four frozen model/class threshold rows", len(frozen) == 4)
    for row in frozen.itertuples(index=False):
        check(f"Frozen cutoff preserved: {row.model_type}/{row.variant_class_group}",
              abs(row.region_score_cutoff - THRESHOLDS[row.model_type]) < 1e-15
              and abs(row.abs_delta_cutoff - DELTA[row.model_type][row.variant_class_group]) < 1e-15)
    threshold_scores = read(COPD / "04_modeling/results/COPD-S4-R001_test_threshold_performance.tsv")
    for model in MODELS:
        values = threshold_scores.loc[(threshold_scores.model_type == model) & (threshold_scores.target_fpr_percent == 5), "score_threshold"]
        check(f"V1 region threshold independently cross-referenced: {model}", len(values) == 1 and abs(values.iloc[0] - THRESHOLDS[model]) < 1e-15)


def verify_forward_gate(original, data, expected):
    preflight = json.loads(provenance("preflight_manifest.json").read_text())
    scoring = json.loads(provenance("scoring_manifest.json").read_text())
    check("Prespecified specification lock retained in manifests",
          preflight["spec_sha256"] == SPEC_SHA and scoring["spec_sha256"] == SPEC_SHA)
    check("Specification lock and passing forward gate precede RC inference",
          datetime.fromisoformat(preflight["spec_lock_utc"])
          <= datetime.fromisoformat(preflight["forward_gate_completed_utc"])
          <= datetime.fromisoformat(scoring["rc_inference_started_utc"]))
    shortlist = read(COPD / "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv")
    panel = sorted(set(sorted(data.loc[data.variant_class_group == "SNV", ID])[:64])
                   | set(sorted(data.loc[data.variant_class_group == "indel_or_complex", ID])[:64])
                   | set(shortlist[ID]) | {"4:88963935:G:T", "4:88962828:C:T"})
    check("Numerical-gate panel exactly matches prespecified selection", panel == preflight["forward_verification_panel_ids"])
    verification = read(V2 / "results" / f"{PREFIX}_forward_verification.tsv")
    check("Numerical-gate pair/model coverage", len(verification) == 2 * len(panel)
          and not verification.duplicated([ID, "model_type"]).any())
    original = original.set_index(ID)
    maxima = {}
    for model in MODELS:
        sample = verification[verification.model_type == model]
        check(f"Forward verification exact IDs: {model}", sample[ID].tolist() == panel)
        source = original.loc[sample[ID]]
        ref = source[f"{model}_ref_score"].to_numpy(float)
        alt = source[f"{model}_alt_score"].to_numpy(float)
        observed_ref = sample.verification_ref_score.to_numpy(float)
        observed_alt = sample.verification_alt_score.to_numpy(float)
        check(f"Verification uses authoritative forward scores: {model}",
              equal_numeric(sample.frozen_ref_score, ref)
              and equal_numeric(sample.frozen_alt_score, alt)
              and equal_numeric(sample.frozen_delta_recomputed, alt - ref))
        check(f"Verification output finite and in range: {model}",
              np.isfinite(observed_ref).all() and np.isfinite(observed_alt).all()
              and ((observed_ref >= 0) & (observed_ref <= 1)
                   & (observed_alt >= 0) & (observed_alt <= 1)).all())
        ref_error, alt_error = abs(observed_ref - ref), abs(observed_alt - alt)
        delta_error = abs((observed_alt - observed_ref) - (alt - ref))
        check(f"Numerical-gate recorded errors independently reproduce: {model}",
              equal_numeric(sample.ref_absolute_error, ref_error)
              and equal_numeric(sample.alt_absolute_error, alt_error)
              and equal_numeric(sample.delta_absolute_error, delta_error)
              and equal_numeric(sample.verification_delta, observed_alt - observed_ref))
        maximum_score, maximum_delta = max(ref_error.max(), alt_error.max()), delta_error.max()
        check(f"Forward numerical gate genuinely passes fixed tolerances: {model}",
              maximum_score <= 1e-5 and maximum_delta <= 2e-5,
              f"max_score_error={maximum_score:.17g}; max_delta_error={maximum_delta:.17g}")
        maxima[model] = {"max_score_error": float(maximum_score), "max_delta_error": float(maximum_delta)}
    check("Preflight claims passed numerical gate", preflight["gate_status"] == "passed")
    check("Scoring run completeness and frozen batch/seed settings",
          scoring["n_rc_pairs_scored"] == 15303 and scoring["seed"] == 20261001
          and scoring["outer_batch_size"] == 256 and scoring["prediction_batch_size"] == 64)
    for package, version in {"tensorflow": "2.20.0", "keras": "3.14.1", "numpy": "2.5.0"}.items():
        check(f"Exact frozen numerical runtime package: {package}", scoring["software"][package] == version)
    return preflight, scoring, maxima


def fasta(path):
    order, sequences = [], {}
    current, pieces = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(pieces).upper()
                current = line[1:].split()[0]
                if current in sequences or current in order:
                    raise AssertionError(f"Duplicate FASTA identifier {current}")
                order.append(current)
                pieces = []
            else:
                pieces.append(line.strip())
        if current is not None:
            sequences[current] = "".join(pieces).upper()
    return order, sequences


def verify_sequence_transformations(original):
    qc = read(result("R002_sequence_transformation_qc.tsv.gz"))
    check("All sequence transformation QC rows and frozen score order", qc[ID].tolist() == original[ID].tolist() and qc[ID].is_unique)
    original_audit = read(COPD / "04_modeling/results/COPD-S4-R002_candidate_sequence_audit.tsv.gz")
    audit = original_audit.set_index(ID).loc[original[ID]]
    alphabet = dict(zip("ACGTN", "TGCAN"))
    encode = np.zeros((256, 4), dtype=np.float32)
    for index, base in enumerate("ACGT"):
        encode[ord(base), index] = 1.0
    problems = []
    for role in ("ref", "alt"):
        ids, sequences = fasta(COPD / f"04_modeling/data/COPD_candidate_variants_{role}_2001bp.fa")
        check(f"{role} FASTA ID count and authoritative order", ids == original[ID].tolist())
        for index, record in enumerate(ids):
            sequence = sequences[record]
            allele = str(audit.iloc[index][role]).upper()
            if len(sequence) != 2001 or set(sequence) - set(alphabet):
                problems.append((record, role, "invalid input sequence"))
                continue
            complement = "".join(alphabet[base] for base in reversed(sequence))
            restored = "".join(alphabet[base] for base in reversed(complement))
            span_start = len(sequence) - 1000 - len(allele)
            span_end = len(sequence) - 1000
            complement_allele = "".join(alphabet[base] for base in reversed(allele))
            known = qc.iloc[index]
            conditions = {
                "allele letters": sequence[1000:1000+len(allele)] == allele and complement[span_start:span_end] == complement_allele,
                "correct involution": restored == sequence,
                "onehot independently reverses axes": np.array_equal(encode[np.frombuffer(complement.encode("ascii"), dtype=np.uint8)], encode[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)][::-1, ::-1]),
                "original string hash": known[f"{role}_forward_sequence_sha256"] == hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "RC string hash": known[f"{role}_rc_sequence_sha256"] == hashlib.sha256(complement.encode("ascii")).hexdigest(),
                "span positions": int(known[f"{role}_forward_allele_start"]) == 1000 and int(known[f"{role}_forward_allele_end"]) == 1000+len(allele) and int(known[f"{role}_rc_allele_start"]) == span_start and int(known[f"{role}_rc_allele_end"]) == span_end,
                "allele length": int(known[f"{role}_allele_length"]) == len(allele),
                "N count": int(known[f"{role}_N_count"]) == sequence.count("N"),
                "edge distances": int(known[f"{role}_forward_min_edge_distance"]) == min(1000, 1001-len(allele)) and int(known[f"{role}_rc_min_edge_distance"]) == min(1000, 1001-len(allele)),
            }
            problems.extend((record, role, label) for label, passed in conditions.items() if not passed)
        del sequences
    check("Independent exact full-sequence RC, allele identity/spans, onehot, hashes, N and edges", not problems, f"30606 sequences; discrepancies={problems[:20]}")
    qc_flags = [column for column in qc if column.endswith("_valid")]
    check("Recorded all-pair transformation checks agree with independent verification", len(qc_flags) == 14 and all(boolean(qc[column]).all() for column in qc_flags))
    reconciliation = read(result("R003_record_reconciliation.tsv"))
    check("15389-record reconciliation preserves original order and unique IDs",
          reconciliation[ID].tolist() == original_audit[ID].tolist() and reconciliation[ID].is_unique)
    check("Original sequence-status reasons retained, including all 86 unscorable cases",
          reconciliation.v1_sequence_status.tolist() == original_audit.sequence_status.tolist()
          and int((original_audit.sequence_status != "ok").sum()) == 86)
    scorable = original_audit[ID].isin(set(original[ID])).to_numpy()
    for field in ("in_frozen_forward_scores", "in_ref_fasta", "in_alt_fasta", "rc_scoring_expected", "rc_scored"):
        check(f"Complete reconciliation identity flag: {field}", np.array_equal(boolean(reconciliation[field]), scorable))
    for field, filename in (("in_frozen_337", "04_modeling/results/COPD-S4-R010_THE_LIST.tsv"),
                            ("in_shortlist_12", "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv")):
        target_ids = set(read(COPD / filename)[ID])
        check(f"Complete reconciliation subgroup: {field}", np.array_equal(boolean(reconciliation[field]), original_audit[ID].isin(target_ids)))
    return qc


def verify_provenance(preflight, scoring):
    snapshot = rows(V2 / "provenance/v1_snapshot.tsv")
    check("Original V1 freeze snapshot hashes remain unchanged",
          all(sha(COPD / row["artifact"]) == row["sha256"] for row in snapshot), f"{len(snapshot)} freeze anchors")
    frozen_scoring = json.loads((COPD / "04_modeling/results/COPD-S4-R003_scoring_manifest.json").read_text())
    check("Every original R003 consumed-file hash still matches",
          all(sha(Path(entry["path"])) == entry["sha256"] for entry in frozen_scoring["files"].values()), f"{len(frozen_scoring['files'])} frozen source/model hashes")
    consumed = rows(provenance("consumed_file_hashes.tsv"))
    altered = []
    for row in consumed:
        path = ROOT / row["path"]
        if sha(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]) or str(path.resolve()) != row["resolved_path"]:
            altered.append(row["label"])
        if row["expected_frozen_sha256"] and row["expected_frozen_sha256"] != row["sha256"]:
            altered.append(row["label"] + ":original frozen hash mismatch")
    check("All consumed inputs and resolved model paths immutable", not altered, f"{len(consumed)} inputs; altered={altered}")
    tracked = rows(provenance("v1_tracked_hashes_before.tsv"))
    observed = {}
    for entry in subprocess.check_output(["git", "ls-files", "-s", "-z", "--", "diseases/COPD"], cwd=ROOT).decode().split("\0"):
        if not entry:
            continue
        info, path = entry.split("\t", 1)
        if path.startswith("diseases/COPD/07_gap_closure/"):
            continue
        mode, blob, stage = info.split()
        observed[path] = (mode, blob, stage)
    check("Tracked V1 path inventory and stages unchanged",
          set(observed) == {row["path"] for row in tracked} and all(v[2] == "0" for v in observed.values()))
    altered = []
    for row in tracked:
        path = ROOT / row["path"]
        if path.is_symlink():
            raw = os.readlink(path).encode()
            digest, size = hashlib.sha256(raw).hexdigest(), len(raw)
        else:
            digest, size = sha(path), path.stat().st_size
        if digest != row["sha256"] or size != int(row["bytes"]) or observed[row["path"]][:2] != (row["git_mode"], row["git_blob"]):
            altered.append(row["path"])
    check("Every tracked V1 content hash and Git index blob unchanged", not altered, f"{len(tracked)} files; altered={altered[:20]}")
    baseline_tree = rows(provenance("v1_tree_stat_before.tsv"))
    current_tree = {}
    for directory, directories, files in os.walk(COPD, followlinks=False):
        directories[:] = sorted(name for name in directories if not (Path(directory) == COPD and name == "07_gap_closure"))
        for name in directories + files:
            path = Path(directory) / name
            info = path.lstat()
            current_tree[str(path.relative_to(ROOT))] = {
                "type": "symlink" if path.is_symlink() else "directory" if path.is_dir() else "file",
                "mode": str(info.st_mode), "size": str(info.st_size), "mtime_ns": str(info.st_mtime_ns),
                "ctime_ns": str(info.st_ctime_ns), "symlink_target": os.readlink(path) if path.is_symlink() else ""}
    check("Entire V1 tree path inventory unchanged, including ignored outputs",
          set(current_tree) == {row["path"] for row in baseline_tree}, f"{len(baseline_tree)} entries")
    altered = [row["path"] for row in baseline_tree if any(current_tree[row["path"]][key] != row[key] for key in current_tree[row["path"]])]
    check("Entire V1 tree metadata unchanged, including symlinks and ignored outputs", not altered, f"altered={altered[:20]}")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    check("No unauthorized commit created during RC audit", head == preflight["git_head"] == scoring["git_head"])
    outside_v2 = [path for path in subprocess.check_output(["git", "diff", "--name-only", "875e995406ddf5d00ce77920d9e650ff128604c0", "--", "diseases/COPD"], cwd=ROOT, text=True).splitlines() if not path.startswith("diseases/COPD/07_gap_closure/")]
    check("Tracked V1 remains unchanged from original freeze commit", not outside_v2, outside_v2)
    return {"tracked_files": len(tracked), "tree_entries": len(baseline_tree), "consumed_inputs": len(consumed)}



def prepare_comparison_expectations(original, data, expected, sequence_qc):
    ids = expected[ID]
    canonical = read(COPD / "04_modeling/results/COPD-S4-R010_THE_LIST.tsv")
    shortlist = read(COPD / "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv")
    pheno = read(V2 / "results/COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv")
    pheno_short = read(V2 / "results/COPD-V2-PHENO-R007_shortlist_phenotype_support.tsv")
    check("Canonical R010 exactly equals independently derived union", set(canonical[ID]) == set(ids[expected["forward_union_call"]]))
    check("Original 337/12 phenotype records preserve exact original ID order",
          pheno[ID].tolist() == canonical[ID].tolist() and pheno_short[ID].tolist() == shortlist[ID].tolist())
    direct = set(pheno.loc[boolean(pheno.primary_retained), ID])
    sole = set(pheno.loc[boolean(pheno.GCST90244098_only_study), ID])
    check("Authoritative descriptive phenotype subgroup sizes", len(direct) == 184 and len(sole) == 124)
    masks = {"all_15303": np.ones(len(ids), dtype=bool), "SNV": expected["classes"] == "SNV",
             "indel_or_complex": expected["classes"] == "indel_or_complex", "call_eligible": expected["eligible"],
             "blacklisted": ~expected["eligible"], "frozen_337": np.isin(ids, canonical[ID]),
             "direct_COPD_184": np.isin(ids, list(direct)), "sole_GCST90244098_124": np.isin(ids, list(sole)),
             "shortlist_12": np.isin(ids, shortlist[ID]),
             "exact_v1_literature": np.isin(ids, ["4:88963935:G:T", "4:88962828:C:T"])}
    n_mask = sequence_qc.ref_N_count.to_numpy() + sequence_qc.alt_N_count.to_numpy() > 0
    if n_mask.any():
        masks["N_containing"], masks["N_free"] = n_mask, ~n_mask
    expected.update({"blacklisted": ~expected["eligible"], "causal_call_eligible": expected["eligible"],
                     "variant_class_group": expected["classes"], "in_frozen_337": masks["frozen_337"],
                     "primary_retained": masks["direct_COPD_184"], "GCST90244098_only_study": masks["sole_GCST90244098_124"],
                     "in_shortlist_12": masks["shortlist_12"], "in_exact_v1_literature": masks["exact_v1_literature"],
                     "ref_N_count": sequence_qc.ref_N_count.to_numpy(), "alt_N_count": sequence_qc.alt_N_count.to_numpy(),
                     "has_ambiguous_N": n_mask})
    p99 = read(provenance("forward_p99_cutoffs.tsv")).set_index(["model_type", "variant_class_group"])
    for model in MODELS:
        cutoff = expected[f"{model}_cutoff"]
        threshold = THRESHOLDS[model]
        f, r = expected[f"{model}_forward_call"], expected[f"{model}_rc_call"]
        expected[f"{model}_region_threshold"] = np.full(len(ids), threshold)
        expected[f"{model}_abs_delta_threshold"] = cutoff
        expected[f"{model}_meaningful_epsilon"] = np.maximum(2e-5, .1 * cutoff)
        for orientation in ("forward", "rc"):
            delta = expected[f"{model}_{orientation}_delta"]
            region = expected[f"{model}_{orientation}_region_score"]
            expected[f"{model}_{orientation}_region_margin"] = region - threshold
            expected[f"{model}_{orientation}_delta_margin"] = abs(delta) - cutoff
            expected[f"{model}_{orientation}_normalized_delta_margin"] = abs(delta) / cutoff - 1
            expected[f"{model}_{orientation}_region_gate"] = region >= threshold
            expected[f"{model}_{orientation}_delta_gate"] = abs(delta) >= cutoff
            expected[f"{model}_{orientation}_raw_sign"] = np.sign(delta).astype(int)
        fdelta = expected[f"{model}_forward_delta"]
        expected[f"{model}_forward_score_borderline"] = abs(expected[f"{model}_forward_region_score"] - threshold) <= .02
        expected[f"{model}_forward_delta_borderline"] = abs(abs(fdelta) / cutoff - 1) <= .10
        expected[f"{model}_forward_any_borderline"] = expected[f"{model}_borderline"]
        expected[f"{model}_forward_tail95"] = expected["eligible"] & (abs(fdelta) >= cutoff)
        tail99_cutoff = np.zeros(len(ids))
        for group in ("SNV", "indel_or_complex"):
            selected = expected["eligible"] & (expected["classes"] == group)
            quantile = float(np.quantile(abs(fdelta[selected]), .99))
            row = p99.loc[(model, group)]
            check(f"Frozen forward-only p99 cutoff independently recomputed: {model}/{group}",
                  abs(quantile - row.forward_abs_delta_p99) < 1e-13
                  and int(row.n_eligible_background) == int(selected.sum())
                  and abs(row.frozen_abs_delta_cutoff - DELTA[model][group]) < 1e-15
                  and abs(row.region_score_cutoff - threshold) < 1e-15)
            # Use the frozen serialized p99 cutoff in decisions, as specified.
            tail99_cutoff[expected["classes"] == group] = row.forward_abs_delta_p99
        expected[f"{model}_forward_abs_delta_p99_cutoff"] = tail99_cutoff
        expected[f"{model}_forward_tail99"] = expected["eligible"] & (abs(fdelta) >= tail99_cutoff)
        for quantity in ("ref_score", "alt_score", "region_score", "delta", "abs_delta"):
            difference = expected[f"{model}_rc_{quantity}"] - expected[f"{model}_forward_{quantity}"]
            expected[f"{model}_{quantity}_signed_disagreement"] = difference
            expected[f"{model}_{quantity}_absolute_disagreement"] = abs(difference)
        expected[f"{model}_call_changed"] = f != r
        expected[f"{model}_lost_forward_call"] = f & ~r
        expected[f"{model}_gained_rc_only_call"] = ~f & r
        expected[f"{model}_lost_well_separated_call"] = expected[f"{model}_forward_well_separated"] & ~r
        expected[f"{model}_raw_sign_changed"] = expected[f"{model}_forward_raw_sign"] != expected[f"{model}_rc_raw_sign"]
        expected[f"{model}_meaningful_sign_changed"] = expected[f"{model}_forward_meaningful_sign"] != expected[f"{model}_rc_meaningful_sign"]
    expected["forward_union_well_separated"] = expected["enhancer_forward_well_separated"] | expected["silencer_forward_well_separated"]
    expected["union_call_changed"] = expected["forward_union_call"] != expected["rc_union_call"]
    expected["lost_union_call"] = expected["forward_union_call"] & ~expected["rc_union_call"]
    expected["retained_union_call"] = expected["forward_union_call"] & expected["rc_union_call"]
    expected["lost_any_original_model_call"] = expected["enhancer_lost_forward_call"] | expected["silencer_lost_forward_call"]
    expected["gained_any_additional_rc_only_model_call"] = expected["enhancer_gained_rc_only_call"] | expected["silencer_gained_rc_only_call"]
    expected["lost_any_well_separated_model_call"] = expected["enhancer_lost_well_separated_call"] | expected["silencer_lost_well_separated_call"]
    expected["exact_model_context_stable"] = expected["forward_model_context"] == expected["rc_model_context"]
    expected["support_switching"] = expected["lost_any_original_model_call"] & expected["gained_any_additional_rc_only_model_call"]
    expected["pure_enhancer_to_H3K27me3_switch"] = (expected["forward_model_context"] == "enhancer_only") & (expected["rc_model_context"] == "H3K27me3_only")
    expected["pure_H3K27me3_to_enhancer_switch"] = (expected["forward_model_context"] == "H3K27me3_only") & (expected["rc_model_context"] == "enhancer_only")
    full = read(result("R004_all_scorable_comparison.tsv.gz"))
    check("Full comparison count, uniqueness and R003 row order",
          len(full) == 15303 and full[ID].is_unique and full[ID].tolist() == list(ids))
    internal = {"classes", "eligible", "enhancer_cutoff", "silencer_cutoff", "enhancer_borderline", "silencer_borderline"}
    for field, values in expected.items():
        if field not in internal:
            verify_array(full, field, values, np.asarray(values).dtype.kind in "fiu")
    # The two subset tables must be exact row selections, not regenerated rankings.
    for suffix, source, rank, new_rank in (
            ("R005_frozen_337_comparison.tsv", canonical, "predicted_causal_priority_rank", "v1_priority_rank"),
            ("R006_shortlist_12_comparison.tsv", shortlist, "experimental_shortlist_rank", "v1_shortlist_rank")):
        table = read(result(suffix))
        check(f"{suffix}: original row order and original rank", table[ID].tolist() == source[ID].tolist()
              and equal_numeric(table[new_rank], source[rank]))
        selected = full.set_index(ID).loc[source[ID]].reset_index()[table.columns]
        pd.testing.assert_frame_equal(table, selected, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-14)
        check(f"{suffix}: every field equals corresponding full-universe row", True, len(table))
    # Original literature register coverage is preserved even when no exact ID is evaluable.
    literature = read(COPD / "04_modeling/results/COPD-S4-R009_literature_functional_variant_recovery.tsv")
    coverage, exact = read(result("R007_literature_coverage.tsv")), read(result("R008_exact_literature_comparison.tsv"))
    pd.testing.assert_frame_equal(coverage[literature.columns], literature, check_dtype=False)
    check("Every original literature register field retained", len(coverage) == 8)
    check("Only two already-cataloged exact literature identities evaluated",
          exact[ID].tolist() == ["4:88963935:G:T", "4:88962828:C:T"]
          and coverage.rc_evaluable_exact_n.astype(int).tolist() == [1,1,0,0,0,0,0,0])
    pd.testing.assert_frame_equal(exact[full.columns].reset_index(drop=True),
                                 full.set_index(ID).loc[exact[ID]].reset_index()[full.columns],
                                 check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-14)
    check("Exact literature comparisons agree with full-universe values", True)
    check("rs2013701 retains V1 original rank 210", int(full.loc[full[ID] == "4:88963935:G:T", "v1_priority_rank"].iloc[0]) == 210)
    return full, masks


def add_distribution(record, prefix, values):
    record.update({prefix + ("maximum" if name == "max" else name): value for name, value in describe(values).items()})


def independent_score_metrics(a, b):
    difference = np.asarray(b) - np.asarray(a)
    absolute = abs(difference)
    n = len(difference)
    status = "not_evaluable_n_lt_2" if n < 2 else "not_evaluable_constant_vector" if np.ptp(a) == 0 or np.ptp(b) == 0 else "evaluable"
    result = {"n": n, "pearson_r": correlation(a,b), "spearman_rho": correlation(a,b,True),
              "correlation_status": status, "signed_mean_difference": float(difference.mean()) if n else np.nan}
    add_distribution(result, "abs_difference_", absolute)
    for label, threshold in (("1e_5", 1e-5), ("0_02", .02), ("0_10", .10)):
        result["fraction_abs_difference_gt_" + label] = rate(int((absolute > threshold).sum()), n)
    return result


def independent_delta_metrics(a, b, cutoff):
    a, b, cutoff = np.asarray(a), np.asarray(b), np.asarray(cutoff)
    n = len(a)
    result = {"n": n, "pearson_r": correlation(a,b), "spearman_rho": correlation(a,b,True),
              "correlation_status": "not_evaluable_n_lt_2" if n < 2 else "not_evaluable_constant_vector" if np.ptp(a) == 0 or np.ptp(b) == 0 else "evaluable",
              "signed_mean_difference": float((b-a).mean()) if n else np.nan}
    add_distribution(result, "delta_abs_difference_", abs(b-a))
    add_distribution(result, "absolute_effect_magnitude_difference_", abs(abs(b)-abs(a)))
    add_distribution(result, "normalized_delta_disagreement_", abs(b-a)/cutoff)
    raw_f, raw_r = np.sign(a), np.sign(b)
    nonzero = (a != 0) & (b != 0)
    raw_equal = (raw_f == raw_r) & nonzero
    epsilon = np.maximum(2e-5, .1*cutoff)
    mf = np.where(a > epsilon, 1, np.where(a < -epsilon, -1, 0))
    mr = np.where(b > epsilon, 1, np.where(b < -epsilon, -1, 0))
    meaningful = mf != 0
    retained = meaningful & (mf == mr)
    collapse = meaningful & (mr == 0)
    opposite = meaningful & (mf == -mr)
    result.update({"raw_both_nonzero_n": int(nonzero.sum()), "raw_sign_agree_n": int(raw_equal.sum()),
                   "raw_nonzero_sign_agreement": rate(int(raw_equal.sum()), int(nonzero.sum())),
                   "raw_gain_to_loss_n": int(((raw_f == 1) & (raw_r == -1)).sum()),
                   "raw_loss_to_gain_n": int(((raw_f == -1) & (raw_r == 1)).sum()),
                   "forward_meaningful_n": int(meaningful.sum()), "meaningful_direction_retained_n": int(retained.sum()),
                   "meaningful_direction_retention": rate(int(retained.sum()), int(meaningful.sum())),
                   "meaningful_near_zero_collapse_n": int(collapse.sum()),
                   "meaningful_near_zero_collapse_fraction": rate(int(collapse.sum()), int(meaningful.sum())),
                   "meaningful_opposite_direction_n": int(opposite.sum()),
                   "meaningful_opposite_direction_fraction": rate(int(opposite.sum()), int(meaningful.sum()))})
    return result


def independent_call_metrics(f, r, w):
    both, lost, gained, neither = [int(v.sum()) for v in (f&r, f&~r, ~f&r, ~f&~r)]
    wn, wl = int(w.sum()), int((w&~r).sum())
    return {"n": len(f), "forward_positive_n": int(f.sum()), "rc_positive_n": int(r.sum()),
            "both_positive_n": both, "forward_only_n": lost, "reverse_only_n": gained, "neither_n": neither,
            "transition_1_to_1": both, "transition_1_to_0": lost, "transition_0_to_1": gained, "transition_0_to_0": neither,
            "jaccard": rate(both,both+lost+gained), "forward_retention": rate(both,int(f.sum())),
            "well_separated_forward_positive_n": wn, "well_separated_lost_n": wl,
            "well_separated_retention": rate(wn-wl,wn)}


def independent_metric_tables(e, masks):
    outputs = {key: [] for key in ("score","delta","call","context","sign","tail","proximity","bins","support","gate")}
    for name, mask in masks.items():
        n = int(mask.sum())
        fctx, rctx = e["forward_model_context"][mask], e["rc_model_context"][mask]
        support = {"stratum": name, "n": n,
                   "forward_union_positive_n": int(e["forward_union_call"][mask].sum()),
                   "rc_union_positive_n": int(e["rc_union_call"][mask].sum()),
                   "either_model_call_changed_n": int((fctx != rctx).sum()),
                   "both_model_calls_changed_n": int((e["enhancer_call_changed"][mask] & e["silencer_call_changed"][mask]).sum())}
        for field in ("exact_model_context_stable","lost_union_call","retained_union_call","lost_any_original_model_call",
                      "gained_any_additional_rc_only_model_call","lost_any_well_separated_model_call","support_switching",
                      "pure_enhancer_to_H3K27me3_switch","pure_H3K27me3_to_enhancer_switch"):
            count = int(e[field][mask].sum())
            support[field+"_n"], support[field+"_fraction"] = count, rate(count,n)
        for model in MODELS:
            f, r, w = (e[f"{model}_{field}"] for field in ("forward_call","rc_call","forward_well_separated"))
            delta_f, delta_r, cutoff = (e[f"{model}_{field}"] for field in ("forward_delta","rc_delta","cutoff"))
            common = {"stratum": name, "model_type": model}
            call = independent_call_metrics(f[mask],r[mask],w[mask])
            outputs["call"].append(dict(common, **call))
            outputs["delta"].append(dict(common, **independent_delta_metrics(delta_f[mask],delta_r[mask],cutoff[mask])))
            for kind in ("ref","alt","region"):
                outputs["score"].append(dict(common, score_type=kind, **independent_score_metrics(e[f"{model}_forward_{kind}_score"][mask],e[f"{model}_rc_{kind}_score"][mask])))
            for percentile in (95,99):
                selected = mask & e[f"{model}_forward_tail{percentile}"]
                record = dict(common, forward_tail=f"matched_class_p{percentile}", **independent_delta_metrics(delta_f[selected],delta_r[selected],cutoff[selected]))
                positive_f, positive_r = int(f[selected].sum()), int(r[selected].sum())
                retained = int((f[selected]&r[selected]).sum())
                record.update({"forward_call_n":positive_f,"rc_call_n":positive_r,"retained_forward_call_n":retained,
                               "forward_call_retention":rate(retained,positive_f),
                               "rc_positive_fraction_of_whole_forward_tail":rate(positive_r,int(selected.sum())),
                               "retained_call_fraction_of_whole_forward_tail":rate(retained,int(selected.sum()))})
                outputs["tail"].append(record)
            for mode in ("raw","meaningful"):
                fs, rs = e[f"{model}_forward_{mode}_sign"][mask], e[f"{model}_rc_{mode}_sign"][mask]
                labels = {-1:"loss", 1:"gain", 0:"zero" if mode == "raw" else "near_zero"}
                for first in (-1,0,1):
                    for second in (-1,0,1):
                        count, denominator = int(((fs==first)&(rs==second)).sum()), int((fs==first).sum())
                        outputs["sign"].append(dict(common,sign_mode=mode,forward_sign=labels[first],rc_sign=labels[second],
                                                    n=count,stratum_n=n,forward_sign_n=denominator,
                                                    fraction_within_forward_sign=rate(count,denominator)))
            changed = f != r
            borderline = e[f"{model}_borderline"]
            count, near = int(changed[mask].sum()), int((changed&borderline&mask).sum())
            support[f"{model}_changed_call_n"] = count
            support[f"{model}_changed_call_forward_borderline_n"] = near
            support[f"{model}_changed_call_forward_borderline_fraction"] = rate(near,count)
            for near_score in (False,True):
                for near_delta in (False,True):
                    selected = mask & (e[f"{model}_forward_score_borderline"] == near_score) & (e[f"{model}_forward_delta_borderline"] == near_delta)
                    bn = int(selected.sum())
                    sr = e[f"{model}_forward_region_gate"] != e[f"{model}_rc_region_gate"]
                    dr = e[f"{model}_forward_delta_gate"] != e[f"{model}_rc_delta_gate"]
                    record = dict(common,forward_score_borderline=near_score,forward_delta_borderline=near_delta,
                                  **independent_call_metrics(f[selected],r[selected],w[selected]))
                    record.update({"eligible_n":int(e["eligible"][selected].sum()),
                                   "region_gate_crossing_n":int(sr[selected].sum()),"region_gate_crossing_fraction":rate(int(sr[selected].sum()),bn),
                                   "delta_gate_crossing_n":int(dr[selected].sum()),"delta_gate_crossing_fraction":rate(int(dr[selected].sum()),bn),
                                   "call_change_n":int(changed[selected].sum()),"call_change_fraction":rate(int(changed[selected].sum()),bn)})
                    add_distribution(record,"region_abs_difference_",e[f"{model}_region_score_absolute_disagreement"][selected])
                    outputs["proximity"].append(record)
            threshold = THRESHOLDS[model]
            boundaries = [0,threshold-.10,threshold-.02,threshold+.02,threshold+.10,1]
            forward_region, rc_region = e[f"{model}_forward_region_score"], e[f"{model}_rc_region_score"]
            for index in range(5):
                lo, hi = boundaries[index:index+2]
                selected = mask & (forward_region>=lo) & ((forward_region<=hi) if index==4 else (forward_region<hi))
                outputs["bins"].append(dict(common,forward_region_bin=index+1,lower_inclusive=lo,upper=hi,upper_inclusive=index==4,
                                            **independent_score_metrics(forward_region[selected],rc_region[selected]),
                                            changed_call_n=int(changed[selected].sum()),
                                            changed_call_fraction=rate(int(changed[selected].sum()),int(selected.sum()))))
            for direction, changing, opposite in (("lost_forward_call",f&~r,"rc"),("gained_rc_only_call",~f&r,"forward")):
                rg, dg = e[f"{model}_{opposite}_region_gate"],e[f"{model}_{opposite}_delta_gate"]
                transition_n = int((changing&mask).sum())
                for failure, failed in (("region_only",~rg&dg),("delta_only",rg&~dg),("both",~rg&~dg)):
                    count = int((mask&changing&failed).sum())
                    outputs["gate"].append(dict(common,transition=direction,failed_orientation=opposite,failed_gate=failure,
                                                n=count,transition_n=transition_n,fraction_of_transition=rate(count,transition_n)))
        support["changed_model_decisions_n"] = support["enhancer_changed_call_n"]+support["silencer_changed_call_n"]
        support["changed_model_decisions_forward_borderline_n"] = support["enhancer_changed_call_forward_borderline_n"]+support["silencer_changed_call_forward_borderline_n"]
        support["changed_model_decisions_forward_borderline_fraction"] = rate(support["changed_model_decisions_forward_borderline_n"],support["changed_model_decisions_n"])
        outputs["support"].append(support)
        outputs["call"].append(dict(stratum=name,model_type="union",
                                    **independent_call_metrics(e["forward_union_call"][mask],e["rc_union_call"][mask],e["forward_union_well_separated"][mask])))
        for first in ("neither","enhancer_only","H3K27me3_only","both"):
            for second in ("neither","enhancer_only","H3K27me3_only","both"):
                count, denominator = int(((fctx==first)&(rctx==second)).sum()),int((fctx==first).sum())
                outputs["context"].append({"stratum":name,"forward_context":first,"rc_context":second,"n":count,
                                           "forward_context_n":denominator,"fraction_within_forward_context":rate(count,denominator),"stratum_n":n})
    return outputs


def verify_metric_tables(computed):
    specs = {
        "score":("R009_score_metrics.tsv",["stratum","model_type","score_type"]),
        "delta":("R010_delta_metrics.tsv",["stratum","model_type"]),
        "sign":("R011_sign_transitions.tsv",["stratum","model_type","sign_mode","forward_sign","rc_sign"]),
        "tail":("R012_forward_tail_metrics.tsv",["stratum","model_type","forward_tail"]),
        "call":("R013_call_transitions.tsv",["stratum","model_type"]),
        "context":("R014_model_context_transitions.tsv",["stratum","forward_context","rc_context"]),
        "proximity":("R015_threshold_proximity.tsv",["stratum","model_type","forward_score_borderline","forward_delta_borderline"]),
        "bins":("R016_forward_region_bins.tsv",["stratum","model_type","forward_region_bin"]),
        "support":("R017_subset_support_summary.tsv",["stratum"]),
        "gate":("R019_gate_failure_transitions.tsv",["stratum","model_type","transition","failed_gate"])}
    for kind,(filename,keys) in specs.items():
        actual = read(result(filename))
        independent = pd.DataFrame(computed[kind])
        check(f"{filename}: complete unique summary keys", len(actual)==len(independent) and not actual.duplicated(keys).any())
        actual = actual.set_index(keys).sort_index()
        independent = independent.set_index(keys).sort_index()
        check(f"{filename}: every intended model/stratum is represented", actual.index.equals(independent.index))
        mismatches = []
        for field in independent:
            if field not in actual:
                mismatches.append(field+":missing")
            elif independent[field].dtype.kind in "fiub":
                if not equal_numeric(actual[field],independent[field],atol=2e-12):
                    mismatches.append(field+":numeric")
            elif not np.array_equal(actual[field].astype(str),independent[field].astype(str)):
                mismatches.append(field+":text")
        check(f"{filename}: every statistic independently recomputed", not mismatches, f"{len(independent)} rows; {len(independent.columns)} columns; mismatches={mismatches}")


def verify_severity(computed):
    scores = pd.DataFrame(computed["score"])
    deltas = pd.DataFrame(computed["delta"]).set_index(["stratum","model_type"])
    tails = pd.DataFrame(computed["tail"]).set_index(["stratum","model_type","forward_tail"])
    calls = pd.DataFrame(computed["call"]).set_index(["stratum","model_type"])
    support = pd.DataFrame(computed["support"]).set_index("stratum")
    expected = {}
    def classify(material, negligible):
        return "scientifically_material" if material else "negligible" if negligible else "modest"
    for model in MODELS:
        s = scores[(scores.stratum=="all_15303")&(scores.model_type==model)]
        expected[("scores",model)] = classify(bool((s.abs_difference_mean>=.05).any() or (s.abs_difference_p95>=.10).any()),
                                             bool((s.abs_difference_mean<=.005).all() and (s.abs_difference_p95<=.02).all()))
        d,t = deltas.loc[("all_15303",model)],tails.loc[("all_15303",model,"matched_class_p95")]
        expected[("allele_deltas",model)] = classify(d.normalized_delta_disagreement_median>=.5 or d.normalized_delta_disagreement_p95>=1 or t.meaningful_direction_retention<.9,
                                                   d.normalized_delta_disagreement_median<=.05 and d.normalized_delta_disagreement_p95<=.25 and t.meaningful_direction_retention>=.99)
    for model in (*MODELS,"union"):
        c = calls.loc[("all_15303",model)]
        expected[("candidate_calls",model)] = classify(c.jaccard<.8 or c.forward_retention<.9 or (c.well_separated_retention<.95 and c.well_separated_lost_n>=3),
                                                      c.jaccard>=.98 and c.forward_retention>=.99 and c.well_separated_lost_n==0)
    s = support.loc["frozen_337"]
    expected[("frozen_337_support_vectors","both_models")] = classify(s.lost_union_call_fraction>.05 or s.lost_any_original_model_call_fraction>.1,
                                                                     s.exact_model_context_stable_fraction>=.99 and s.lost_any_well_separated_model_call_n==0)
    ordering = {"negligible":0,"modest":1,"scientifically_material":2}
    worst = max(expected.values(),key=ordering.get)
    expected[("overall","both_models_and_union")] = worst
    actual = read(result("R018_severity.tsv")).set_index(["domain","model_type"])
    check("Severity table has all and only prespecified domain/model rows", set(actual.index)==set(expected))
    check("Severity independently follows unchanged locked rubric", all(actual.loc[key,"severity"]==value for key,value in expected.items()),f"overall={worst}")
    return worst



def verify_artifact_manifests(preflight, scoring):
    check("Final scoring-script hash preserved", scoring["script_sha256"] == sha(V2 / "scripts/run_reverse_complement_scoring.py"))
    for path, metadata in scoring["outputs"].items():
        actual = ROOT / path
        check("Scoring manifest output: " + actual.name,
              actual.is_file() and sha(actual) == metadata["sha256"] and actual.stat().st_size == metadata["bytes"])
    analyzer = json.loads(provenance("analysis_manifest.json").read_text())
    check("Analysis locked specification and source script", analyzer["spec_sha256"] == SPEC_SHA
          and sha(ROOT / analyzer["script"]["path"]) == analyzer["script"]["sha256"])
    check("Analyzer consumed-input hashes independently match", all(sha(ROOT / value["path"]) == value["sha256"] for value in analyzer["inputs"].values()))
    for metadata in analyzer["outputs"]:
        actual = V2 / metadata["path"]
        check("Analysis manifest output: " + actual.name,
              actual.is_file() and sha(actual) == metadata["sha256"] and actual.stat().st_size == metadata["size_bytes"])
    check("Analysis summary hash matches", sha(V2 / analyzer["summary"]["path"]) == analyzer["summary"]["sha256"])
    analysis_qc = read(result("R020_analysis_validation.tsv"))
    check("Analyzer internal checks complete and pass", len(analysis_qc) == analyzer["analysis_qc_passed"] and (analysis_qc.status == "PASS").all())
    figures = json.loads(provenance("figure_manifest.json").read_text())
    check("Figure specification and all consumed hashes", figures["spec_sha256"] == SPEC_SHA
          and all(sha(ROOT / value["path"]) == value["sha256"] for value in figures["inputs"]))
    expected_figures = {f"results/{PREFIX}-{stem}.{extension}" for stem in
                        ("F001_score_concordance","F002_delta_concordance","F003_threshold_proximity","F004_call_transitions")
                        for extension in ("png","pdf")}
    check("Four required figure families available in PNG/PDF",
          len(figures["outputs"]) == 8
          and {str((ROOT/value["path"]).relative_to(V2)) for value in figures["outputs"]} == expected_figures)
    check("Every graphical artifact has matching manifest hash/size",
          all(sha(ROOT/value["path"]) == value["sha256"] and (ROOT/value["path"]).stat().st_size == value["bytes"]
              and value["bytes"] > 1000 for value in figures["outputs"]))
    failed_gate = provenance("preflight_manifest_attempt1_failed.json")
    failed_table = V2 / "results" / f"{PREFIX}_forward_verification_attempt1_failed.tsv"
    if failed_gate.exists():
        historical = json.loads(failed_gate.read_text())
        check("Historical failed attempt preserved without modifying specification",
              failed_table.is_file() and historical["spec_sha256"] == SPEC_SHA
              and historical["gate_status"] == "failed_forward_numerical_gate")
        check("Numerical-panel batching correction documented; unique panel unchanged",
              "forward_verification_execution_policy" in preflight
              and preflight["forward_verification_panel_ids"] == historical["forward_verification_panel_ids"])
    return {"analysis_checks": len(analysis_qc), "graphical_files": len(figures["outputs"])}


def verify_report_and_registers(preflight, severity):
    # This report uses an underscore after the RC namespace, unlike numbered results.
    report = V2 / "results" / f"{PREFIX}_reverse_complement_robustness_report.md"
    check("Required detailed audit report exists", report.is_file())
    text = report.read_text()
    check("Detailed report states complete census and prespecified severity",
          "15303" in text.replace(",","") and severity.replace("_"," ") in text.lower().replace("_"," "))
    own_outputs = {provenance("final_validation.tsv"), provenance("final_validation_manifest.json"), provenance("artifact_checksums.tsv")}
    missing_links = []
    for target in re.findall(r"\]\(([^)]+)\)", text):
        target = target.strip("<>").split("#")[0]
        if not target or "://" in target or target.startswith("mailto:"):
            continue
        path = (report.parent / target).resolve()
        if not path.exists() and path not in own_outputs:
            missing_links.append(target)
    check("Report local artifact links resolve", not missing_links, missing_links)
    result_rows = []
    for filename in ("activity_log.tsv","gap_closure_decision_register.tsv","gap_closure_result_register.tsv"):
        path = V2 / filename
        with path.open(newline="") as handle:
            actual = list(csv.reader(handle, delimiter="\t"))
        check("Register rectangular unique-ID structure: " + filename,
              len({len(row) for row in actual}) == 1 and len({row[0] for row in actual[1:]}) == len(actual)-1)
        original_text = subprocess.check_output(["git","show",preflight["git_head"]+":"+str(path.relative_to(ROOT))],cwd=ROOT,text=True)
        original_rows = list(csv.reader(original_text.splitlines(),delimiter="\t"))
        lookup = {row[0]:row for row in actual[1:]}
        check("Previous register content preserved: " + filename,
              actual[0] == original_rows[0] and all(lookup.get(row[0]) == row for row in original_rows[1:]))
        additions = [row for row in actual[1:] if row[0] not in {r[0] for r in original_rows[1:]}]
        check("RC completion appended to register: " + filename,
              bool(additions) and all("COPD-V2" in row[0] for row in additions)
              and any("reverse" in " ".join(row).lower() or PREFIX in " ".join(row) for row in additions))
        if filename == "gap_closure_result_register.tsv":
            result_rows = rows(path)
    check("All registered results resolve, including planned final QC outputs",
          all((V2/row["primary_file"]).is_file() or (V2/row["primary_file"]).resolve() in own_outputs for row in result_rows))
    # PHENO checksums attest the prior committed register state and must not be rewritten.
    prior_ledger = V2 / "provenance/COPD-V2-PHENO_artifact_checksums.tsv"
    original_ledger = subprocess.check_output(["git","show",preflight["git_head"]+":"+str(prior_ledger.relative_to(ROOT))],cwd=ROOT)
    check("Historical phenotype checksum ledger retained unchanged", sha(prior_ledger) == hashlib.sha256(original_ledger).hexdigest())
    return report


def run_independent_validation():
    check("Locked specification SHA256", sha(provenance("analysis_specification.md")) == SPEC_SHA)
    verify_original_thresholds()
    original, metadata, rc, data, expected = construct_independent_data()
    preflight, scoring, maxima = verify_forward_gate(original, data, expected)
    sequence_qc = verify_sequence_transformations(original)
    full, masks = prepare_comparison_expectations(original,data,expected,sequence_qc)
    metrics = independent_metric_tables(expected,masks)
    verify_metric_tables(metrics)
    severity = verify_severity(metrics)
    artifact_summary = verify_artifact_manifests(preflight,scoring)
    verify_report_and_registers(preflight,severity)
    immutable = verify_provenance(preflight,scoring)
    check("Locked specification still unchanged at completion", sha(provenance("analysis_specification.md")) == SPEC_SHA)
    return {"status":"PASS","v1_unchanged":True,"overall_severity":severity,"scorable_pairs":15303,
            "original_records":15389,"frozen_candidates":337,"shortlist_entries":12,"direct_COPD_candidates":184,
            "sole_GCST90244098_candidates":124,"numerical_gate":maxima,**immutable,**artifact_summary,
            "completed_modules":["frozen_model_reverse_complement_robustness_audit"],
            "no_new_model_training_or_replacement_ranking":True}


def write_final_outputs(summary):
    qc = provenance("final_validation.tsv")
    with qc.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(CHECKS)
    summary.update({"validator": str(Path(__file__).relative_to(ROOT)), "validator_sha256": sha(Path(__file__)),
                    "specification_sha256": SPEC_SHA, "completed_utc": datetime.now(timezone.utc).isoformat(),
                    "checks": len(CHECKS), "passed": sum(r["status"] == "PASS" for r in CHECKS),
                    "failed": sum(r["status"] == "FAIL" for r in CHECKS),
                    "implementation_independence": "No scoring or analyzer import; reconstruct from original V1 scores and raw RC scores",
                    "command": sys.argv, "executable": sys.executable, "working_directory": os.getcwd(),
                    "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                    "checksum_self_exclusion": "COPD-V2-RC_artifact_checksums.tsv is excluded to avoid a recursive self-hash"})
    provenance("final_validation_manifest.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    selected = []
    for path in V2.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts or "COPD-V2-RC_runtime" in path.parts or any(part.startswith(".") for part in path.relative_to(V2).parts[:-1]):
            continue
        if path.name == f"{PREFIX}_artifact_checksums.tsv":
            continue
        if path.name.startswith(PREFIX) or path.name in {"run_reverse_complement_scoring.py", "analyze_reverse_complement_audit.py", "validate_reverse_complement_audit.py", "plot_reverse_complement_audit.py", "activity_log.tsv", "gap_closure_decision_register.tsv", "gap_closure_result_register.tsv", ".gitignore", "README.md"}:
            selected.append(path)
    with provenance("artifact_checksums.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "bytes", "sha256"])
        for path in sorted(selected):
            writer.writerow([str(path.relative_to(ROOT)), path.stat().st_size, sha(path)])
    return len(selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", action="store_true", help="Confirm report/register/figure preparation is complete")
    args = parser.parse_args()
    if not args.final:
        parser.error("Final validation requires --final after report and registers are ready")
    summary = {"result_id": "COPD-V2-RC-QC-FINAL", "status": "FAIL", "v1_unchanged": False}
    try:
        summary.update(run_independent_validation())
        count = write_final_outputs(summary)
        # Verify the completed checksum ledger after writing it, excluding its own hash.
        checksum_rows = rows(provenance("artifact_checksums.tsv"))
        if not all(sha(ROOT/row["path"]) == row["sha256"] and (ROOT/row["path"]).stat().st_size == int(row["bytes"]) for row in checksum_rows):
            raise AssertionError("Completed checksum ledger verification failed")
        print(json.dumps({"status": "PASS", "checks": len(CHECKS), "checksummed_artifacts": count,
                          "v1_unchanged": True, "overall_severity": summary["overall_severity"]}, indent=2))
        return 0
    except Exception as exc:
        if not CHECKS or CHECKS[-1]["status"] != "FAIL":
            CHECKS.append({"check": "Unhandled final-validation condition", "status": "FAIL", "detail": f"{type(exc).__name__}: {exc}"})
        summary["error"] = str(exc)
        count = write_final_outputs(summary)
        print(json.dumps({"status": "FAIL", "error": str(exc), "checks": len(CHECKS), "checksummed_artifacts": count}, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
