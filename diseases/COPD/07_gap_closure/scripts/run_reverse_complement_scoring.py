#!/usr/bin/env python3
"""Prespecified frozen-model RC preflight, numerical verification and inference.

R001 uses the original R003 column names, but every score in R001 is a reverse-
complement score. No forward score or V1 artifact is overwritten. Run with -B.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
COPD = ROOT / "diseases/COPD"
V1 = COPD / "04_modeling"
V2 = COPD / "07_gap_closure"
PROV, RESULTS, LOGS = (V2 / d for d in ("provenance", "results", "logs"))
RUNTIME = PROV / "COPD-V2-RC_runtime"
for directory in (PROV, RESULTS, LOGS, RUNTIME, RUNTIME / "tmp", RUNTIME / "cache"):
    directory.mkdir(parents=True, exist_ok=True)
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
os.environ["TMPDIR"] = str(RUNTIME / "tmp")
os.environ["XDG_CACHE_HOME"] = str(RUNTIME / "cache")
os.environ["MPLCONFIGDIR"] = str(RUNTIME / "matplotlib")
os.environ["KERAS_HOME"] = str(RUNTIME / "keras")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("TF_XLA_FLAGS", "--tf_xla_enable_xla_devices=false")
os.environ["TF_DETERMINISTIC_OPS"] = "1"
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

SPEC = PROV / "COPD-V2-RC_analysis_specification.md"
SPEC_HASH = "09d3656f3e06c24fa5402c28b976d6be06fff91604ad0b7b29c3150120d76ac1"
SEED, OUTER_BATCH, PREDICTION_BATCH = 20261001, 256, 64
MODELS = {
    "enhancer": "COPD_SevereEmphysema_Lung_Enhancer_DHS_x2",
    "silencer": "COPD_SevereEmphysema_Lung_Silencer_DHS_x2",
}
SCORE_COLUMNS = ["candidate_record_id"] + [
    f"{model}_{field}" for model in MODELS
    for field in ("ref_score", "alt_score", "delta_alt_minus_ref")
]
OUT_RAW = RESULTS / "COPD-V2-RC-R001_reverse_complement_scores.tsv.gz"
OUT_SEQ = RESULTS / "COPD-V2-RC-R002_sequence_transformation_qc.tsv.gz"
OUT_RECON = RESULTS / "COPD-V2-RC-R003_record_reconciliation.tsv"
OUT_FORWARD = RESULTS / "COPD-V2-RC_forward_verification.tsv"
OUT_P99 = PROV / "COPD-V2-RC_forward_p99_cutoffs.tsv"
PREFLIGHT = PROV / "COPD-V2-RC_preflight_manifest.json"
MANIFEST = PROV / "COPD-V2-RC_scoring_manifest.json"
TRACKED = PROV / "COPD-V2-RC_v1_tracked_hashes_before.tsv"
TREE = PROV / "COPD-V2-RC_v1_tree_stat_before.tsv"
CONSUMED = PROV / "COPD-V2-RC_consumed_file_hashes.tsv"
ORIGINAL_FREEZE = PROV / "COPD-V2-RC_original_freeze_validation.tsv"
LOG = LOGS / "COPD-V2-RC_scoring.log"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def now():
    return datetime.now(timezone.utc).isoformat()


def log(message):
    line = f"{now()} {message}"
    print(line, flush=True)
    with LOG.open("a") as handle:
        handle.write(line + "\n")


def check(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def string_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_tsv(path, rows, columns=None):
    table = pd.DataFrame(rows, columns=columns)
    compression = {"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None
    table.to_csv(path, sep="\t", index=False, compression=compression)


def read_table(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, low_memory=False)


def flags(series):
    return series.astype(str).str.lower().isin(("true", "1", "yes"))


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def load_fasta(path):
    names, sequences = [], {}
    with path.open() as handle:
        current, pieces = None, []
        for line in handle:
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(pieces).upper()
                current = line[1:].split()[0]
                check(current not in names, f"duplicate FASTA ID: {current}")
                names.append(current)
                pieces = []
            else:
                pieces.append(line.strip())
        if current is not None:
            sequences[current] = "".join(pieces).upper()
    return names, sequences


def tree_inventory():
    rows = []
    for directory, dirs, files in os.walk(COPD, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not (Path(directory) == COPD and d == "07_gap_closure"))
        for name in sorted(dirs + files):
            path = Path(directory) / name
            info = path.lstat()
            kind = "symlink" if stat.S_ISLNK(info.st_mode) else "directory" if stat.S_ISDIR(info.st_mode) else "file"
            rows.append({
                "path": relative(path), "type": kind, "mode": info.st_mode,
                "size": info.st_size, "mtime_ns": info.st_mtime_ns,
                "ctime_ns": info.st_ctime_ns,
                "symlink_target": os.readlink(path) if kind == "symlink" else "",
            })
    return rows


def tracked_inventory():
    raw = subprocess.check_output(["git", "ls-files", "-s", "-z", "--", "diseases/COPD"], cwd=ROOT).decode()
    rows = []
    for entry in raw.split("\0"):
        if not entry:
            continue
        metadata, path = entry.split("\t", 1)
        if path.startswith("diseases/COPD/07_gap_closure/"):
            continue
        mode, blob, stage = metadata.split()
        check(stage == "0", f"unmerged frozen file {path}")
        full = ROOT / path
        check(full.exists(), f"missing tracked frozen file {path}")
        if full.is_symlink():
            content = os.readlink(full).encode()
            digest, size = hashlib.sha256(content).hexdigest(), len(content)
        else:
            digest, size = sha256(full), full.stat().st_size
        rows.append({"path": path, "sha256": digest, "bytes": size, "git_blob": blob, "git_mode": mode})
    return rows


def files_to_consume():
    paths = {
        "specification": SPEC,
        "original_v1_snapshot": PROV / "v1_snapshot.tsv",
        "scoring_script": Path(__file__),
        "frozen_scoring_manifest": V1 / "results/COPD-S4-R003_scoring_manifest.json",
        "frozen_sequence_manifest": V1 / "results/COPD-S4-R002_candidate_sequence_manifest.json",
        "frozen_analysis_manifest": V1 / "results/COPD-S4-R004_analysis_manifest.json",
        "ref_fasta": V1 / "data/COPD_candidate_variants_ref_2001bp.fa",
        "alt_fasta": V1 / "data/COPD_candidate_variants_alt_2001bp.fa",
        "scores": V1 / "results/COPD-S4-R003_candidate_allele_scores.tsv.gz",
        "sequence_audit": V1 / "results/COPD-S4-R002_candidate_sequence_audit.tsv.gz",
        "prioritized_candidates": V1 / "results/COPD-S4-R004_prioritized_candidates.tsv.gz",
        "delta_thresholds": V1 / "results/COPD-S4-R004_delta_thresholds.tsv",
        "score_thresholds": V1 / "results/COPD-S4-R001_test_threshold_performance.tsv",
        "canonical_337": V1 / "results/COPD-S4-R010_THE_LIST.tsv",
        "shortlist_12": COPD / "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv",
        "phenotype_337": RESULTS / "COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv",
        "phenotype_12": RESULTS / "COPD-V2-PHENO-R007_shortlist_phenotype_support.tsv",
        "literature": V1 / "results/COPD-S4-R009_literature_functional_variant_recovery.tsv",
        "inference_source": ROOT / "models/TREDNET_v2/TREDNet_v2_inference.py",
        "phase_one_weights": V1 / "trednet/model_phase_I/phase_one_weights.h5",
        "v1_scoring_script": V1 / "scripts/03_score_candidate_alleles.py",
        "v1_sequence_script": V1 / "scripts/02_prepare_candidate_alleles.py",
        "v1_decision_script": V1 / "scripts/05_prioritize_candidate_variants.py",
    }
    for model, eid in MODELS.items():
        base = V1 / "trednet/models_output" / eid
        paths[f"{model}_best_weights"] = base / f"{eid}_phase_two_weights.weights.h5"
        paths[f"{model}_architecture"] = base / "phase_two_model.keras"
        paths[f"{model}_threshold_file"] = base / "fpr_threshold_scores.txt"
    return paths


def check_transform(sequence, allele, encode):
    check(len(sequence) == 2001, "sequence length !=2001")
    check(not (set(sequence) - set("ACGTN")), "unsupported sequence symbol")
    check(sequence[1000:1000 + len(allele)] == allele, "forward allele-span mismatch")
    transformed = rc(sequence)
    start, end = 1001 - len(allele), 1001
    check(start >= 0, "allele extends beyond frozen window")
    check(transformed[start:end] == rc(allele), "RC allele-span mismatch")
    check(transformed[1000] == rc(sequence[1000]), "RC anchor mismatch")
    check(rc(transformed) == sequence, "RC involution failed")
    check(np.array_equal(encode(transformed), encode(sequence)[::-1, ::-1]), "RC encoding failed")
    return transformed


def synthetic_checks(encode):
    output = []
    flank_left = ("ACGT" * 250)
    flank_right = ("TGCA" * 501)
    for case, ref, alt, with_n in [
        ("SNV", "G", "T", False), ("insertion", "A", "ACT", False),
        ("deletion", "AGT", "A", False), ("multibase", "ACT", "TGA", False),
        ("N_context", "G", "A", True),
    ]:
        left = flank_left if not with_n else "N" + flank_left[1:]
        for name, allele in (("ref", ref), ("alt", alt)):
            sequence = (left + allele + flank_right)[:2001]
            check_transform(sequence, allele, encode)
            output.append({"case": case, "allele_role": name, "pass": True})
    return output


def main():
    start_utc = now()
    log("Starting frozen preflight; no RC inference authorized before every gate passes")
    check(sha256(SPEC) == SPEC_HASH, "locked specification hash mismatch")
    check(not OUT_RAW.exists(), "raw RC output already exists; refuse overwrite/re-inference")
    if PREFLIGHT.exists() and json.loads(PREFLIGHT.read_text()).get("gate_status") == "failed_forward_numerical_gate":
        for source in (PREFLIGHT, OUT_FORWARD):
            destination = source.with_name(source.stem + "_attempt1_failed" + source.suffix)
            check(not destination.exists(), "failed attempt archive already exists")
            shutil.copyfile(source, destination)
        log("IMPLEMENTATION CORRECTION: preserve failed partial-batch verification; pad forward-only execution to full64 minibatches, retain same142 panel IDs and unchanged tolerances; RC batches remain original V1 shapes")
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    baseline_tree, baseline_tracked = tree_inventory(), tracked_inventory()
    write_tsv(TREE, baseline_tree)
    write_tsv(TRACKED, baseline_tracked)
    log(f"V1 baseline captured: {len(baseline_tracked)} tracked files, {len(baseline_tree)} filesystem entries")

    paths = files_to_consume()
    freeze_checks = []
    for row in read_table(paths["original_v1_snapshot"]).to_dict(orient="records"):
        frozen_path = COPD / row["artifact"]
        observed = sha256(frozen_path)
        check(observed == row["sha256"], f"original V1 freeze mismatch: {row['artifact']}")
        freeze_checks.append({"snapshot_id":row["snapshot_id"],"artifact":row["artifact"],
                              "expected_sha256":row["sha256"],"observed_sha256":observed,
                              "pass":True,"checked_utc":now()})
    check(len(freeze_checks)==3,"original freeze record count changed")
    write_tsv(ORIGINAL_FREEZE,freeze_checks)
    frozen_manifest = json.loads(paths["frozen_scoring_manifest"].read_text())
    rows = []
    for label, path in paths.items():
        check(path.is_file(), f"missing consumed input {label}: {path}")
        expected = frozen_manifest["files"].get(label, {}).get("sha256", "")
        digest = sha256(path)
        check(not expected or digest == expected, f"frozen hash mismatch: {label}")
        rows.append({"label": label, "path": relative(path), "resolved_path": str(path.resolve()),
                     "sha256": digest, "bytes": path.stat().st_size, "expected_frozen_sha256": expected})
    write_tsv(CONSUMED, rows)
    log("All frozen R003 model, FASTA, forward-score and inference-source hashes verified")
    inference_spec = importlib.util.spec_from_file_location("frozen_trednet_inference", paths["inference_source"])
    inference = importlib.util.module_from_spec(inference_spec)
    inference_spec.loader.exec_module(inference)
    synthetic = synthetic_checks(inference.seq2one_hot)

    scores = read_table(paths["scores"])
    audit = read_table(paths["sequence_audit"])
    prior = read_table(paths["prioritized_candidates"])
    canonical = read_table(paths["canonical_337"])
    shortlist = read_table(paths["shortlist_12"])
    literature = read_table(paths["literature"])
    pheno = read_table(paths["phenotype_337"])
    pheno12 = read_table(paths["phenotype_12"])
    thresholds = read_table(paths["delta_thresholds"])
    for name, frame, count in [("scores", scores, 15303), ("audit", audit, 15389), ("prior", prior, 15389),
                               ("canonical", canonical, 337), ("shortlist", shortlist, 12),
                               ("pheno", pheno, 337), ("pheno12", pheno12, 12)]:
        check(len(frame) == count, f"wrong {name} row count")
        check(frame.candidate_record_id.is_unique, f"duplicate {name} IDs")
    ids = scores.candidate_record_id.tolist()
    ref_ids, ref_seqs = load_fasta(paths["ref_fasta"])
    alt_ids, alt_seqs = load_fasta(paths["alt_fasta"])
    check(ids == ref_ids == alt_ids, "score/REF/ALT order disagreement")
    check(set(audit.candidate_record_id) == set(prior.candidate_record_id), "audit/prior universe mismatch")
    check(set(audit.loc[audit.sequence_status.eq("ok"), "candidate_record_id"]) == set(ids), "scorable status mismatch")
    check(set(canonical.candidate_record_id) == set(pheno.candidate_record_id), "337 phenotype reconciliation failed")
    check(set(shortlist.candidate_record_id) == set(pheno12.candidate_record_id), "12 phenotype reconciliation failed")
    check(set(shortlist.candidate_record_id) <= set(canonical.candidate_record_id), "shortlist outside canonical337")
    check(int(flags(pheno.primary_retained).sum()) == 184, "direct COPD stratum count mismatch")
    check(int(flags(pheno.GCST90244098_only_study).sum()) == 124, "sole GCST90244098 stratum count mismatch")

    a = audit.set_index("candidate_record_id")
    p = prior.set_index("candidate_record_id").loc[ids]
    f = scores.set_index("candidate_record_id")
    calls, p99_rows, rounding = {}, [], {}
    expected_thresholds = {
        ("enhancer", "SNV"): (0.643623, 0.05706318769999998),
        ("enhancer", "indel_or_complex"): (0.643623, 0.04936093850000001),
        ("silencer", "SNV"): (0.58505, 0.028802613899999996),
        ("silencer", "indel_or_complex"): (0.58505, 0.026183359500000003),
    }
    eligible = ~flags(p.encode_blacklist)
    check(int(eligible.sum()) == 15283, "eligible background count mismatch")
    check(flags(p.causal_call_eligible).equals(eligible), "eligibility flags mismatch")
    for model in MODELS:
        ref, alt = f[f"{model}_ref_score"], f[f"{model}_alt_score"]
        delta, region = alt - ref, np.maximum(ref, alt)
        check(np.isfinite(ref).all() and np.isfinite(alt).all(), "nonfinite authoritative scores")
        rounding[model] = float(np.abs(delta - f[f"{model}_delta_alt_minus_ref"]).max())
        call = pd.Series(False, index=ids)
        for group in ("SNV", "indel_or_complex"):
            threshold_row = thresholds.loc[thresholds.model_type.eq(model) & thresholds.variant_class_group.eq(group)]
            check(len(threshold_row) == 1, "threshold row ambiguity")
            t, d = float(threshold_row.iloc[0].region_score_cutoff), float(threshold_row.iloc[0].abs_delta_cutoff)
            et, ed = expected_thresholds[(model, group)]
            check(abs(t-et)<1e-15 and abs(d-ed)<1e-15, "frozen cutoff disagreement")
            select = eligible & p.variant_class_group.eq(group)
            check(int(select.sum()) == (13747 if group == "SNV" else 1536), "class background count mismatch")
            call |= select & region.ge(t) & delta.abs().ge(d)
            p99_rows.append({"model_type": model, "variant_class_group": group,
                             "n_eligible_background": int(select.sum()),
                             "forward_abs_delta_p99": float(np.quantile(delta.loc[select].abs(), .99, method="linear")),
                             "frozen_abs_delta_cutoff": d, "region_score_cutoff": t})
        check(np.array_equal(call.to_numpy(), flags(p[f"predicted_causal_{model}"]).to_numpy()), "frozen model call reconstruction failed")
        check(int(call.sum()) == (175 if model == "enhancer" else 199), "forward call count mismatch")
        calls[model] = call
    union = calls["enhancer"] | calls["silencer"]
    check(int((calls["enhancer"] & calls["silencer"]).sum()) == 37, "intersection mismatch")
    check(set(union.index[union]) == set(canonical.candidate_record_id), "canonical337 call reconstruction failed")
    write_tsv(OUT_P99, p99_rows)

    seq_rows = []
    for idx, record_id in enumerate(ids):
        row = {"candidate_record_id": record_id, "variant_class_group": p.loc[record_id, "variant_class_group"]}
        for role, sequences in (("ref", ref_seqs), ("alt", alt_seqs)):
            allele = str(a.loc[record_id, role]).upper()
            sequence = sequences[record_id]
            transformed = check_transform(sequence, allele, inference.seq2one_hot)
            row.update({f"{role}_N_count": sequence.count("N"), f"{role}_allele_length": len(allele),
                        f"{role}_forward_allele_start": 1000, f"{role}_forward_allele_end": 1000 + len(allele),
                        f"{role}_rc_allele_start": 1001 - len(allele), f"{role}_rc_allele_end": 1001,
                        f"{role}_forward_min_edge_distance": min(1000, 1001-len(allele)),
                        f"{role}_rc_min_edge_distance": min(1001-len(allele), 1000),
                        f"{role}_forward_sequence_sha256": string_hash(sequence),
                        f"{role}_rc_sequence_sha256": string_hash(transformed),
                        f"{role}_length_valid": True, f"{role}_symbols_valid": True,
                        f"{role}_forward_allele_span_valid": True, f"{role}_rc_allele_span_valid": True,
                        f"{role}_rc_anchor_valid": True, f"{role}_rc_involution_valid": True,
                        f"{role}_one_hot_transform_valid": True})
        seq_rows.append(row)
        if (idx + 1) % 3000 == 0:
            log(f"Validated transformations for {idx+1}/15303 pairs (no inference)")
    write_tsv(OUT_SEQ, seq_rows)
    canonical_set, short_set, id_set = set(canonical.candidate_record_id), set(shortlist.candidate_record_id), set(ids)
    recon = [{"candidate_record_id": r.candidate_record_id, "v1_sequence_status": r.sequence_status,
              "in_frozen_forward_scores": r.candidate_record_id in id_set,
              "in_ref_fasta": r.candidate_record_id in ref_seqs, "in_alt_fasta": r.candidate_record_id in alt_seqs,
              "in_frozen_337": r.candidate_record_id in canonical_set, "in_shortlist_12": r.candidate_record_id in short_set,
              "rc_scoring_expected": r.candidate_record_id in id_set,
              "rc_scored": False, "reconciliation_pass": True} for r in audit.itertuples(index=False)]
    write_tsv(OUT_RECON, recon)
    literature_ids = set()
    for field in literature.exact_candidate_record_ids:
        literature_ids.update(part for part in str(field).split(";") if part)
    check(literature_ids == {"4:88963935:G:T", "4:88962828:C:T"}, "exact literature IDs changed")
    panel_ids = sorted(set(sorted(p.index[p.variant_class_group.eq("SNV")])[:64]) |
                       set(sorted(p.index[p.variant_class_group.eq("indel_or_complex")])[:64]) |
                       short_set | literature_ids)
    preflight = {
        "specification_id": "COPD-V2-RC-SPEC-001", "spec_sha256": SPEC_HASH,
        "spec_lock_utc": start_utc, "git_head": git_head,
        "created_utc": now(), "gate_status": "sequence_and_frozen_integrity_passed_pending_forward_gate",
        "input_files": rows, "synthetic_transformation_checks": synthetic,
        "coverage": {"all_records":15389,"scorable_pairs":15303,"eligible":15283,"blacklisted":20,
                     "canonical_337":337,"shortlist_12":12,"direct_COPD":184,"sole_GCST90244098":124},
        "forward_call_counts": {"enhancer":175,"silencer":199,"both":37,"union":337},
        "serialized_forward_delta_max_rounding_discrepancy": rounding,
        "forward_verification_panel_ids": panel_ids,
        "forward_verification_execution_policy": "Pad142 unique panel entries to192 by repeating first50 solely to restore full64 inference minibatches; assess only prespecified142; original RC15303 batches unchanged",
        "original_freeze_validation":freeze_checks,
        "forward_p99_cutoffs": p99_rows,
        "baseline_files": {"tracked":relative(TRACKED),"tree":relative(TREE),"consumed":relative(CONSUMED)},
    }
    write_json(PREFLIGHT, preflight)
    log(f"Sequence/identity/hash/forward-call preflight PASSED; forward numerical panel n={len(panel_ids)} next")

    import keras
    import tensorflow as tf
    keras.utils.set_random_seed(SEED)
    tf.config.experimental.enable_op_determinism()
    devices = tf.config.list_physical_devices()
    for gpu in tf.config.list_physical_devices("GPU"):
        tf.config.experimental.set_memory_growth(gpu, True)
    phase_one = inference.load_phase_one(str(V1 / "trednet/model_phase_I"))
    phase_two = {model: inference.load_phase_two(str(V1 / "trednet/models_output" / eid), eid)
                 for model, eid in MODELS.items()}
    software = {name: importlib.metadata.version(name) for name in
                ("tensorflow","keras","numpy","pandas","scipy","biopython","matplotlib","h5py")}
    software.update({"python":platform.python_version(),"executable":sys.executable,"platform":platform.platform()})
    check(software["tensorflow"] == "2.20.0" and software["keras"] == "3.14.1" and software["numpy"] == "2.5.0", "model runtime library mismatch")

    def predict_batch(names, reverse):
        ref = np.stack([inference.seq2one_hot(rc(ref_seqs[name]) if reverse else ref_seqs[name]) for name in names])
        alt = np.stack([inference.seq2one_hot(rc(alt_seqs[name]) if reverse else alt_seqs[name]) for name in names])
        embeddings = phase_one.predict(np.concatenate((ref,alt),axis=0),batch_size=PREDICTION_BATCH,verbose=0)
        split = len(names)
        output = {}
        for model, classifier in phase_two.items():
            output[(model,"ref")] = classifier.predict(embeddings[:split,...,np.newaxis],batch_size=PREDICTION_BATCH,verbose=0).ravel()
            output[(model,"alt")] = classifier.predict(embeddings[split:,...,np.newaxis],batch_size=PREDICTION_BATCH,verbose=0).ravel()
            for role in ("ref","alt"):
                value = output[(model,role)]
                check(np.isfinite(value).all() and ((value>=0)&(value<=1)).all(), "invalid model output")
        return output

    verification_rows = []
    for offset in range(0,len(panel_ids),OUTER_BATCH):
        names = panel_ids[offset:offset+OUTER_BATCH]
        padded_names = names + names[:(-len(names)) % PREDICTION_BATCH]
        values = predict_batch(padded_names,False)
        for j,name in enumerate(names):
            for model in MODELS:
                expected_ref, expected_alt = float(f.loc[name,f"{model}_ref_score"]),float(f.loc[name,f"{model}_alt_score"])
                observed_ref, observed_alt = float(values[(model,"ref")][j]),float(values[(model,"alt")][j])
                verification_rows.append({"candidate_record_id":name,"model_type":model,
                    "frozen_ref_score":expected_ref,"frozen_alt_score":expected_alt,"frozen_delta_recomputed":expected_alt-expected_ref,
                    "verification_ref_score":observed_ref,"verification_alt_score":observed_alt,"verification_delta":observed_alt-observed_ref,
                    "ref_absolute_error":abs(observed_ref-expected_ref),"alt_absolute_error":abs(observed_alt-expected_alt),
                    "delta_absolute_error":abs((observed_alt-observed_ref)-(expected_alt-expected_ref))})
    write_tsv(OUT_FORWARD,verification_rows)
    verification = pd.DataFrame(verification_rows)
    gate = {}
    for model in MODELS:
        subset = verification.loc[verification.model_type.eq(model)]
        max_score_error = float(subset[["ref_absolute_error","alt_absolute_error"]].max().max())
        max_delta_error = float(subset.delta_absolute_error.max())
        passed = max_score_error <= 1e-5 and max_delta_error <= 2e-5
        gate[model] = {"max_score_absolute_error":max_score_error,"max_delta_absolute_error":max_delta_error,"pass":passed}
    preflight.update({"numerical_gate":gate,"software":software,"devices":[str(d) for d in devices],
                      "gate_status":"passed" if all(g["pass"] for g in gate.values()) else "failed_forward_numerical_gate",
                      "forward_gate_completed_utc":now()})
    write_json(PREFLIGHT,preflight)
    check(all(g["pass"] for g in gate.values()),f"forward numerical gate FAILED: {gate}")
    check(sha256(SPEC)==SPEC_HASH,"spec changed before RC inference")
    log(f"FORWARD NUMERICAL GATE PASSED: {gate}; launching RC inference for all15303 pairs")

    inference_start = now()
    count = 0
    with OUT_RAW.open("wb") as raw, gzip.GzipFile(filename="",fileobj=raw,mode="wb",compresslevel=6,mtime=0) as zipped, io.TextIOWrapper(zipped) as text_handle:
        writer = csv.DictWriter(text_handle,fieldnames=SCORE_COLUMNS,delimiter="\t")
        writer.writeheader()
        for offset in range(0,len(ids),OUTER_BATCH):
            names = ids[offset:offset+OUTER_BATCH]
            values = predict_batch(names,True)
            for j,name in enumerate(names):
                row = {"candidate_record_id":name}
                for model in MODELS:
                    ref,alt = float(values[(model,"ref")][j]),float(values[(model,"alt")][j])
                    row[f"{model}_ref_score"] = format(ref,".17g")
                    row[f"{model}_alt_score"] = format(alt,".17g")
                    row[f"{model}_delta_alt_minus_ref"] = format(alt-ref,".17g")
                writer.writerow(row)
            count += len(names)
            log(f"RC scored {count}/{len(ids)} pairs")
    check(count==15303,"incomplete RC scoring")
    raw_result = read_table(OUT_RAW)
    check(raw_result.candidate_record_id.tolist()==ids,"RC output ID/order reconciliation failed")
    for row in recon:
        row["rc_scored"] = row["candidate_record_id"] in id_set
    write_tsv(OUT_RECON,recon)
    log("Complete15303-pair RC output written; rechecking V1 immutability")
    after_tree, after_tracked = tree_inventory(), tracked_inventory()
    check(after_tree==baseline_tree,"V1 filesystem stat inventory changed")
    check(after_tracked==baseline_tracked,"tracked V1 content changed")
    for row in rows:
        check(sha256(ROOT/row["path"])==row["sha256"],f"consumed input changed: {row['label']}")
    check(sha256(SPEC)==SPEC_HASH,"spec changed during scoring")
    outputs = [OUT_RAW,OUT_SEQ,OUT_RECON,OUT_FORWARD,OUT_P99,PREFLIGHT,TRACKED,TREE,CONSUMED,ORIGINAL_FREEZE]
    manifest = {"result_id":"COPD-V2-RC-R001","started_utc":start_utc,"rc_inference_started_utc":inference_start,
                "completed_utc":now(),"spec_sha256":SPEC_HASH,"script_sha256":sha256(Path(__file__)),
                "git_head":git_head,"input_files":rows,"software":software,"devices":[str(d) for d in devices],
                "seed":SEED,"outer_batch_size":OUTER_BATCH,"prediction_batch_size":PREDICTION_BATCH,
                "command":sys.argv,"environment":{k:os.environ.get(k,"") for k in
                ("PYTHONPATH","PYTHONDONTWRITEBYTECODE","TF_DETERMINISTIC_OPS","TF_XLA_FLAGS","TMPDIR","XDG_CACHE_HOME","KERAS_HOME","MPLCONFIGDIR")},
                "raw_score_orientation":"reverse_complement; original R003-compatible score column names",
                "raw_score_precision":"17 significant digits, float32 model outputs represented exactly as Python float",
                "models":MODELS,"n_rc_pairs_scored":count,"forward_verification":gate,
                "v1_immutability":{"consumed_hashes_unchanged":True,"tracked_hashes_unchanged":True,"tree_stat_unchanged":True},
                "outputs":{relative(path):{"sha256":sha256(path),"bytes":path.stat().st_size} for path in outputs}}
    write_json(MANIFEST,manifest)
    log(f"SCORING COMPLETE; V1 unchanged; manifest {relative(MANIFEST)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"STOP: {type(exc).__name__}: {exc}")
        raise
