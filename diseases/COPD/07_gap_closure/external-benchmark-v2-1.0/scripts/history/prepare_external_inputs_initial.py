#!/usr/bin/env python3
"""Prepare exact frozen benchmark identities, without labels or model selection.

Reads only the already-opened immutable benchmark snapshot, frozen prospective
specifications and the existing indexed hg38 genome. All outputs are new-only.
No original benchmark/candidate source table, model or GPU module is opened.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import platform
import resource
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

STAGE = Path(__file__).resolve().parents[1]
REPO = STAGE.parents[3]
EXPECTED_FREEZE = "d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path):
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size,
            "sha256": sha(path)}


def json_new(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def tsv_new(path, rows, columns):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns,
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def one_hot(sequence):
    encoded = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return (encoded[:, None] == np.frombuffer(b"ACGT", dtype=np.uint8)).astype(np.float32)


class IndexedFasta:
    """Read-only five-column FASTA-index reader; never creates an index."""

    def __init__(self, fasta, fai):
        self.index = {}
        with fai.open() as handle:
            for line in handle:
                name, length, offset, bases, width, *_ = line.rstrip("\n").split("\t")
                self.index[name] = tuple(map(int, (length, offset, bases, width)))
        self.handle = fasta.open("rb")

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        if start < 0 or start > length or end < start:
            raise ValueError("invalid FASTA interval")
        end = min(end, length)
        if end == start:
            return ""
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.handle.seek(first)
        sequence = self.handle.read(last - first).replace(b"\n", b"").replace(b"\r", b"")
        assert len(sequence) == end - start
        return sequence.decode("ascii").upper()

    def close(self):
        self.handle.close()


def frozen_identity(row):
    reasons = []
    if row["identity_status"] != "exact":
        reasons.append("frozen_identity_not_exact")
    if row["reference_verified"].lower() != "true":
        reasons.append("frozen_reference_not_verified")
    for key in ("exact_identity_eligible", "construct_exact_allele_identity_valid"):
        if row[key].lower() == "false":
            reasons.append("frozen_" + key + "_false")
    parts = [row[key] for key in ("grch38_chrom", "grch38_pos", "grch38_ref", "grch38_alt")]
    if not all(parts) or row["canonical_variant_id"] != ":".join(parts):
        reasons.append("missing_or_inconsistent_frozen_exact_key")
    return not reasons, ";".join(reasons)


def construct(genome, chrom, pos1, ref, alt):
    """Implement the unchanged prospective V1/RC geometry and exclusions."""
    if pos1 < 1:
        return "invalid_coordinate", "", "", ""
    if not ref or not alt or set(ref + alt) - set("ACGT"):
        return "ambiguous_or_empty_allele", "", "", ""
    if max(len(ref), len(alt)) > 1001:
        return "allele_span_exceeds_2001bp_window", "", "", ""
    fasta_chrom = "chr" + chrom if "chr" + chrom in genome.index else chrom
    if fasta_chrom not in genome.index:
        return "missing_fasta_chromosome", "", "", ""
    pos0 = pos1 - 1
    if pos0 < 1000 or pos0 + len(ref) > genome.index[fasta_chrom][0]:
        return "chromosome_boundary", "", "", ""
    observed = genome.fetch(fasta_chrom, pos0, pos0 + len(ref))
    if observed != ref:
        return "reference_mismatch", "", "", observed
    left = genome.fetch(fasta_chrom, pos0 - 1000, pos0)
    right = genome.fetch(fasta_chrom, pos0 + len(ref), pos0 + len(ref) + 2001)
    ref_sequence = (left + ref + right)[:2001]
    alt_sequence = (left + alt + right)[:2001]
    if len(ref_sequence) != 2001 or len(alt_sequence) != 2001:
        return "chromosome_boundary", "", "", observed
    if set(ref_sequence + alt_sequence) - set("ACGTN"):
        return "unsupported_sequence_symbol", "", "", observed
    return "ok", ref_sequence, alt_sequence, observed


def validate_sequence(sequence, allele):
    assert len(sequence) == 2001 and set(sequence) <= set("ACGTN")
    assert sequence[1000:1000 + len(allele)] == allele
    transformed = rc(sequence)
    assert transformed[1001 - len(allele):1001] == rc(allele)
    assert transformed[1000] == rc(sequence[1000])
    assert rc(transformed) == sequence
    assert np.array_equal(one_hot(transformed), one_hot(sequence)[::-1, ::-1])
    return {"sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "rc_sha256": hashlib.sha256(transformed.encode("ascii")).hexdigest(),
            "N_count": sequence.count("N")}


def synthetic_validation():
    results = []
    left, right = "ACGT" * 250, "TGCA" * 600
    for name, ref, alt in (("SNV", "A", "G"), ("insertion", "A", "ATGC"),
                           ("deletion", "ATGC", "A"), ("MNV", "ATGC", "GCAT"),
                           ("max_allele", "A", "T" * 1001)):
        for role, allele in (("REF", ref), ("ALT", alt)):
            sequence = (left + allele + right)[:2001]
            validate_sequence(sequence, allele)
            results.append({"case": name, "allele_role": role, "status": "PASS"})
    assert np.array_equal(one_hot("N"), np.zeros((1, 4), dtype=np.float32))
    assert np.array_equal(one_hot("ACGT"), np.eye(4, dtype=np.float32))
    return results


def main():
    start, started_utc = time.perf_counter(), datetime.now(timezone.utc).isoformat()
    outputs = [STAGE / "inputs" / x for x in
               ("identity_map.tsv", "scorable_variants.tsv", "allele_sequences.npy",
                "input_manifest.json", "sequence_qc.tsv")]
    validation_path = STAGE / "provenance/sequence_construction_validation.json"
    assert not any(p.exists() for p in outputs + [validation_path]), "Refuse any output overwrite"
    freeze_path = STAGE / "provenance/prospective_freeze.json"
    assert sha(freeze_path) == EXPECTED_FREEZE
    freeze = json.loads(freeze_path.read_text())
    for rec in freeze["files"]:
        assert record(STAGE / rec["path"]) == rec, rec["path"]
    spec = json.loads((STAGE / "specification/evaluation_specification.json").read_text())
    assert spec["sequence_length"] == 2001 and spec["allele_start_zero_based"] == 1000
    assert spec["maximum_allele_length"] == 1001 and spec["new_identity_rescue"] is False
    complete = json.loads((STAGE / "provenance/benchmark_ingestion_complete.json").read_text())
    snapshot_path = STAGE / "inputs/benchmark_snapshot.json.gz"
    assert complete["status"] == "PASS" and record(snapshot_path) == complete["snapshot"]
    snapshot = json.loads(gzip.decompress(snapshot_path.read_bytes()))
    master_name = "COPD-V2-BENCH-R003_frozen_benchmark_master.tsv"
    harmon_name = "COPD-V2-BENCH-R002_exact_allele_harmonization.tsv"
    unique_name = "COPD-V2-BENCH-R015_unique_variant_summary.tsv"
    master = snapshot["tables"][master_name]["rows"]
    harmon = snapshot["tables"][harmon_name]["rows"]
    unique = snapshot["tables"][unique_name]["rows"]
    assert len(master) == len(harmon) == 14025
    by_assay = {r["assay_id"]: r for r in harmon}
    assert len(by_assay) == len(harmon) == len({r["assay_id"] for r in master})
    for row in master:
        assert all(row[k] == v for k, v in by_assay[row["assay_id"]].items()), row["assay_id"]
    unique_keys = {r["benchmark_variant_key"] for r in unique}
    assert len(unique_keys) == len(unique) == 1731
    observations, identities = [], {}
    for i, row in enumerate(master):
        exact, exclusion = frozen_identity(row)
        key = row["canonical_variant_id"] or (
            "unresolved:" + row["rsid"] + ":" + row["tested_allele1"] + "/" + row["tested_allele2"])
        assert key in unique_keys, key
        variant_id = row["canonical_variant_id"] if exact else ""
        if exact:
            entry = (row["grch38_chrom"], row["grch38_pos"], row["grch38_ref"], row["grch38_alt"])
            assert variant_id not in identities or identities[variant_id] == entry
            identities[variant_id] = entry
        observations.append({
            "observation_index": i, "assay_id": row["assay_id"], "study_id": row["study_id"],
            "rsid": row["rsid"], "benchmark_variant_key": key,
            "canonical_variant_id": row["canonical_variant_id"], "variant_id": variant_id,
            "identity_status": row["identity_status"], "reference_verified": row["reference_verified"],
            "exact_identity_eligible": row["exact_identity_eligible"],
            "construct_exact_allele_identity_valid": row["construct_exact_allele_identity_valid"],
            "identity_exact": exact, "identity_exclusion_reason": exclusion,
            "frozen_identity_reason": row["identity_reason"],
            "sequence_status": "pending" if exact else "unresolved_identity",
            "sequence_exclusion_reason": "" if exact else exclusion,
            "sequence_row_index": ""})
    assert {r["benchmark_variant_key"] for r in observations} == unique_keys
    fasta_rec = next(r for r in spec["dependencies"] if r["path"].endswith("/hg38.fa"))
    fai_rec = next(r for r in spec["dependencies"] if r["path"].endswith("/hg38.fa.fai"))
    for rec in (fasta_rec, fai_rec):
        source = REPO / rec["path"]
        assert source.stat().st_size == rec["bytes"] and sha(source) == rec["sha256"], rec["path"]
    synthetic = synthetic_validation()
    genome = IndexedFasta(REPO / fasta_rec["path"], REPO / fai_rec["path"])
    scorable, sequence_arrays, qc = [], [], []
    observed_by_key = {}
    for key in sorted(unique_keys):
        matching = [r for r in observations if r["benchmark_variant_key"] == key]
        assert matching
        observed_by_key[key] = matching
    ordered = sorted(identities, key=lambda k: (int(identities[k][0]), int(identities[k][1]),
                                                identities[k][2], identities[k][3]))
    for variant_id in ordered:
        chrom, pos_text, ref, alt = identities[variant_id]
        pos1 = int(pos_text)
        status, ref_seq, alt_seq, observed_ref = construct(genome, chrom, pos1, ref, alt)
        qrow = {"benchmark_variant_key": variant_id, "variant_id": variant_id,
                "chrom": chrom, "pos1": pos1, "ref": ref, "alt": alt,
                "sequence_status": status, "reference_observed": observed_ref,
                "sequence_row_index": "", "ref_sha256": "", "alt_sha256": "",
                "ref_rc_sha256": "", "alt_rc_sha256": "", "ref_N_count": "", "alt_N_count": "",
                "allele_spans_valid": "", "rc_involution_valid": "", "one_hot_rc_equivalence_valid": ""}
        if status == "ok":
            row_index = len(scorable)
            checks = [validate_sequence(ref_seq, ref), validate_sequence(alt_seq, alt)]
            qrow.update(sequence_row_index=row_index, allele_spans_valid=True,
                        rc_involution_valid=True, one_hot_rc_equivalence_valid=True)
            for role, result in zip(("ref", "alt"), checks):
                qrow.update({role + "_sha256": result["sha256"],
                             role + "_rc_sha256": result["rc_sha256"],
                             role + "_N_count": result["N_count"]})
            scorable.append({"variant_id": variant_id, "chrom": chrom, "pos1": pos1,
                             "ref": ref, "alt": alt, "row_index": row_index,
                             "variant_class": "SNV" if len(ref) == len(alt) == 1 else "indel_or_complex",
                             "n_assay_observations": len(observed_by_key[variant_id]),
                             "ref_sequence_sha256": checks[0]["sha256"],
                             "alt_sequence_sha256": checks[1]["sha256"]})
            sequence_arrays.append(np.stack([np.frombuffer(s.encode("ascii"), dtype=np.uint8)
                                             for s in (ref_seq, alt_seq)]))
        for row in observed_by_key[variant_id]:
            row.update(sequence_status=status, sequence_exclusion_reason="" if status == "ok" else status,
                       sequence_row_index=qrow["sequence_row_index"])
        qc.append(qrow)
    genome.close()
    for key in sorted(unique_keys - set(identities)):
        first = observed_by_key[key][0]
        qrow = dict.fromkeys(qc[0], "")
        qrow.update(benchmark_variant_key=key, sequence_status="unresolved_identity")
        qc.append(qrow)
        assert not first["identity_exact"]
    array = np.stack(sequence_arrays)
    assert array.dtype == np.uint8 and array.shape == (len(scorable), 2, 2001)
    assert all(row["sequence_status"] != "pending" for row in observations)
    assert len(qc) == len(unique) and len(observations) == len(master)
    tsv_new(outputs[0], observations, list(observations[0]))
    tsv_new(outputs[1], scorable, list(scorable[0]))
    with outputs[2].open("xb") as handle:
        np.save(handle, array, allow_pickle=False)
    tsv_new(outputs[4], qc, list(qc[0]))
    validation = {
        "status": "PASS", "stage": spec["stage"], "started_utc": started_utc,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "prospective_freeze": record(freeze_path), "snapshot": record(snapshot_path),
        "master_harmonization_identity_equality": True, "complete_original_population_accounted": True,
        "n_observations": len(observations), "n_unique_original_identities": len(unique),
        "n_exact_identities": len(identities), "n_scorable_variants": len(scorable),
        "identity_status_counts_observations": dict(Counter(r["identity_status"] for r in observations)),
        "sequence_status_counts_unique": dict(Counter(r["sequence_status"] for r in qc)),
        "sequence_status_counts_observations": dict(Counter(r["sequence_status"] for r in observations)),
        "synthetic_checks": synthetic, "real_sequence_pair_checks": len(scorable),
        "labels_or_scope_used_for_sequence_selection": False, "new_identity_rescue": False,
        "allele_swaps_or_strand_rescue": False, "original_benchmark_reopened": False,
        "candidate_universe_accessed": False, "GPU_or_models_imported": False,
        "sources": {"hg38_fasta": fasta_rec, "hg38_fai": fai_rec},
        "script": record(Path(__file__).resolve()),
        "command_argv": [sys.executable] + sys.argv,
        "environment": {"python": platform.python_version(), "executable": sys.executable,
                        "numpy": np.__version__, "platform": platform.platform(),
                        "FASTA_reader": "stdlib read-only 5-column FAI direct byte offsets"},
        "compute": {"wall_seconds": time.perf_counter() - start,
                    "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "files": [record(path) for path in (outputs[0], outputs[1], outputs[2], outputs[4])]}
    json_new(validation_path, validation)
    manifest = {
        "status": "PASS", "stage": spec["stage"], "n_variants": len(scorable),
        "n_observations": len(observations), "n_original_variant_identities": len(unique),
        "n_unresolved_original_identities": len(unique) - len(identities),
        "dtype": "uint8", "shape": list(array.shape), "encoding": "ASCII ACGTN",
        "allele_axis_order": ["REF", "ALT"], "base_axis": "positive genomic orientation",
        "row_order": "scorable_variants.tsv row_index; numeric chromosome, position, REF, ALT",
        "identity_rule": "identity_status exact AND reference_verified True AND complete canonical exact key; explicit False eligibility/construct flags veto; blank optional flags do not veto",
        "sequence_rule": "unchanged prospective 1000-left allele substitution and independent right crop to 2001bp",
        "outcome_or_scope_selection": False, "prospective_freeze": record(freeze_path),
        "snapshot": record(snapshot_path), "sources": {"hg38_fasta": fasta_rec, "hg38_fai": fai_rec},
        "validation": record(validation_path), "files": validation["files"],
        "script": record(Path(__file__).resolve()), "created_utc": datetime.now(timezone.utc).isoformat()}
    json_new(outputs[3], manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
