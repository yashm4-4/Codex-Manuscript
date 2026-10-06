#!/usr/bin/env python3
"""Independent CPU audit of frozen external allele construction and annotations.

This validator imports neither the sequence builder nor annotation/scoring/model
code. A separate stdlib FAI reader reconstructs every allele from the immutable
hg38 file, and taxonomy checks are derived from the new-stage source snapshot.
It writes only its exclusive validation receipt, never scientific inputs.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np

STAGE = Path(__file__).resolve().parents[1]
REPO = STAGE.parents[3]
PROSPECTIVE_SHA = "d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830"
SPECIFICATION_SHA = "32432be4a7698256b8ae3cb6fe4139eac8a1ad6c83c5b8de41b880ce771e11bc"
MASTER = "COPD-V2-BENCH-R003_frozen_benchmark_master.tsv"
HARMONIZATION = "COPD-V2-BENCH-R002_exact_allele_harmonization.tsv"
CONTEXTUAL = "COPD-V2-BENCH-R004B_contextual_and_excluded_evidence.tsv"
VARIANTS = "COPD-V2-BENCH-R015_unique_variant_summary.tsv"
REPORTERS = {"MPRA_allele_effect", "conventional_reporter"}
CHECKS = Counter()


def check(condition, category, detail=""):
    CHECKS[category] += 1
    if not bool(condition):
        raise AssertionError(f"{category}: {detail}")


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def record(path):
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": digest(path)}


def load_json(path):
    with path.open() as stream:
        return json.load(stream)


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def truth(value):
    return str(value).lower() in {"true", "1"}


def explicitly_false(value):
    return str(value).lower() in {"false", "0"}


def row_digest(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


class IndependentFasta:
    """Independent whole-span read with FASTA newlines removed after byte seek."""

    def __init__(self, fasta, fai):
        self.records = {}
        with fai.open() as stream:
            for line in stream:
                name, length, offset, bases, width = line.rstrip("\n").split("\t")[:5]
                self.records[name] = tuple(map(int, (length, offset, bases, width)))
        self.stream = fasta.open("rb")

    def fetch(self, chromosome, start, end):
        length, offset, bases, width = self.records[chromosome]
        check(0 <= start < end <= length, "fasta_boundaries", (chromosome, start, end))
        byte_start = offset + (start // bases) * width + start % bases
        last = end - 1
        byte_end = offset + (last // bases) * width + last % bases + 1
        self.stream.seek(byte_start)
        raw = self.stream.read(byte_end - byte_start)
        result = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
        check(len(result) == end - start, "fasta_read_length", chromosome)
        return result.upper()


def audit():
    check(platform.python_version() == "3.13.0" and np.__version__ == "2.5.0", "frozen_CPU_environment")
    check(os.environ.get("PYTHONHASHSEED") == "104729", "frozen_CPU_environment")
    check(not (STAGE / "predictions").exists(), "pre_prediction_chronology")
    check(digest(STAGE / "provenance/prospective_freeze.json") == PROSPECTIVE_SHA, "prospective_hash")
    check(digest(STAGE / "specification/evaluation_specification.json") == SPECIFICATION_SHA, "prospective_hash")
    specification = load_json(STAGE / "specification/evaluation_specification.json")
    prospective = load_json(STAGE / "provenance/prospective_freeze.json")
    input_manifest = load_json(STAGE / "inputs/input_manifest.json")
    annotation_rules = load_json(STAGE / "provenance/annotation_rules.json")
    for group in (prospective["files"], input_manifest["files"], annotation_rules["outputs"]):
        for item in group:
            check(record(STAGE / item["path"]) == item, "input_or_annotation_file_hash", item["path"])
    check(record(STAGE / "inputs/benchmark_snapshot.json.gz") == input_manifest["snapshot"], "snapshot_hash")
    check(input_manifest["snapshot"] == annotation_rules["source"], "snapshot_hash")
    sources = input_manifest["sources"]
    dependencies = {r["path"]: r for r in specification["dependencies"]}
    for source in sources.values():
        frozen = dependencies[source["path"]]
        check(source == frozen, "frozen_reference_metadata", source["path"])
        path = REPO / source["path"]
        check(path.stat().st_size == source["bytes"] and digest(path) == source["sha256"],
              "frozen_reference_byte_hash", source["path"])
    scoring_rows = table(STAGE / "inputs/scorable_variants.tsv")
    sequence_qc = table(STAGE / "inputs/sequence_qc.tsv")
    qc_by_id = {r["variant_id"]: r for r in sequence_qc if r["variant_id"]}
    sequences = np.load(STAGE / "inputs/allele_sequences.npy", allow_pickle=False)
    check(len(scoring_rows) == 1710 == input_manifest["n_variants"], "sequence_population")
    check(sequences.shape == (1710, 2, 2001) and sequences.dtype == np.uint8, "sequence_array_schema")
    check(input_manifest["shape"] == list(sequences.shape) and input_manifest["allele_axis_order"] == ["REF", "ALT"], "sequence_array_schema")
    check(len({r["variant_id"] for r in scoring_rows}) == 1710, "sequence_population")
    complement = str.maketrans("ACGTN", "TGCAN")
    fasta = IndependentFasta(REPO / sources["hg38_fasta"]["path"], REPO / sources["hg38_fai"]["path"])
    nucleotide_counts, allele_types, chromosomes = Counter(), Counter(), Counter()
    sequence_details = []
    for i, row in enumerate(scoring_rows):
        identity = ":".join(row[k] for k in ("chrom", "pos1", "ref", "alt"))
        check(identity == row["variant_id"] and int(row["row_index"]) == i, "identity_and_row_order", identity)
        p0 = int(row["pos1"]) - 1
        chromosome = "chr" + row["chrom"]
        chromosomes[row["chrom"]] += 1
        ref, alt = row["ref"], row["alt"]
        check(all(1 <= len(a) <= 1001 and set(a) <= set("ACGT") for a in (ref, alt)), "allele_validity", identity)
        observed_ref = fasta.fetch(chromosome, p0, p0 + len(ref))
        check(observed_ref == ref, "actual_hg38_REF_match", identity)
        # A single contiguous genomic extraction provides a separately expressed
        # reconstruction of the prospective fixed-left, shared-downstream rule.
        genome = fasta.fetch(chromosome, p0 - 1000, p0 + len(ref) + 2001)
        left, right = genome[:1000], genome[1000 + len(ref):]
        expected = [(left + a + right)[:2001] for a in (ref, alt)]
        check(genome[1000:1000 + len(ref)] == ref, "actual_hg38_REF_match", identity)
        qc = qc_by_id[identity]
        check(qc["sequence_status"] == "ok" and int(qc["sequence_row_index"]) == i, "sequence_QC_consistency", identity)
        check(qc["reference_observed"] == ref, "sequence_QC_consistency", identity)
        if len(ref) == len(alt) == 1:
            allele_types["SNV"] += 1
        elif len(ref) == len(alt):
            allele_types["MNV"] += 1
        elif len(ref) < len(alt):
            allele_types["insertion"] += 1
        else:
            allele_types["deletion"] += 1
        detail = {"variant_id": identity, "row_index": i}
        for a, name in enumerate(("ref", "alt")):
            actual = sequences[i, a].tobytes().decode("ascii")
            check(actual == expected[a] and len(actual) == 2001, "all_actual_allele_genome_reconstructions", (identity, name))
            check(set(actual) <= set("ACGTN"), "ACGTN_only_context", (identity, name))
            check(actual[1000:1000 + len(row[name])] == row[name], "forward_allele_span", (identity, name))
            rc = actual.translate(complement)[::-1]
            check(rc.translate(complement)[::-1] == actual, "RC_involution", (identity, name))
            check(rc[1001 - len(row[name]):1001] == row[name].translate(complement)[::-1], "RC_allele_span", (identity, name))
            forward_onehot = np.asarray([[base == x for x in "ACGT"] for base in actual], dtype=np.uint8)
            rc_onehot = np.asarray([[base == x for x in "ACGT"] for base in rc], dtype=np.uint8)
            check(np.array_equal(rc_onehot, forward_onehot[::-1, ::-1]), "onehot_nucleotide_RC_equivalence", (identity, name))
            n_mask = sequences[i, a] == ord("N")
            check(np.all(forward_onehot[n_mask] == 0), "N_all_zero_onehot", (identity, name))
            nucleotide_counts.update(actual)
            digest_forward = hashlib.sha256(actual.encode()).hexdigest()
            digest_rc = hashlib.sha256(rc.encode()).hexdigest()
            check(digest_forward == row[name + "_sequence_sha256"] == qc[name + "_sha256"], "sequence_sha256_consistency", (identity, name))
            check(digest_rc == qc[name + "_rc_sha256"], "sequence_sha256_consistency", (identity, name))
            check(actual.count("N") == int(qc[name + "_N_count"]), "sequence_N_count_consistency", (identity, name))
            detail[name + "_sha256"] = digest_forward
        sequence_details.append(detail)
    fasta.stream.close()
    # The independent one-hot implementation must handle N even if no actual
    # frozen allele contains it. Synthetic validation changes no biological data.
    synthetic = "ACGTNN"
    synthetic_rc = synthetic.translate(complement)[::-1]
    a = np.asarray([[base == x for x in "ACGT"] for base in synthetic], dtype=np.uint8)
    b = np.asarray([[base == x for x in "ACGT"] for base in synthetic_rc], dtype=np.uint8)
    check(np.array_equal(a[::-1, ::-1], b) and not a[4:].any(), "synthetic_N_RC_encoding")
    with gzip.open(STAGE / "inputs/benchmark_snapshot.json.gz", "rt") as stream:
        snapshot = json.load(stream)
    master = snapshot["tables"][MASTER]["rows"]
    harmonized = {r["assay_id"]: r for r in snapshot["tables"][HARMONIZATION]["rows"]}
    original_variant_rows = snapshot["tables"][VARIANTS]["rows"]
    original_variants = {r["benchmark_variant_key"]: r for r in original_variant_rows}
    observations = table(STAGE / "inputs/observation_annotations.tsv")
    variants = table(STAGE / "inputs/variant_annotations.tsv")
    contextual = table(STAGE / "inputs/contextual_annotations.tsv")
    check(len(master) == len(observations) == len(harmonized) == 14025, "complete_master_population")
    check(len(variants) == len(original_variants) == 1731, "complete_variant_population")
    check(len(contextual) == len(snapshot["tables"][CONTEXTUAL]["rows"]) == 267, "complete_contextual_population")
    check({int(r["snapshot_row_index"]) for r in observations} == set(range(14025)), "master_row_bijection")
    groups, exact_ids, directional = defaultdict(list), set(), defaultdict(list)
    state_counts, scope_counts = Counter(), Counter()
    for row in observations:
        source = master[int(row["snapshot_row_index"])]
        check(row["source_table"] == MASTER and row["assay_id"] == source["assay_id"], "original_row_link", row["observation_id"])
        check(row["frozen_row_json_sha256"] == row_digest(source), "original_row_hash", row["assay_id"])
        harmony = harmonized[source["assay_id"]]
        for field in ("canonical_variant_id", "grch38_chrom", "grch38_pos", "grch38_ref", "grch38_alt", "identity_status", "reference_verified"):
            check(source.get(field, "") == harmony.get(field, ""), "master_harmonization_identity_equality", (row["assay_id"], field))
        exact = (source["identity_status"] == "exact" and truth(source["reference_verified"])
                 and not explicitly_false(source.get("exact_identity_eligible", ""))
                 and not explicitly_false(source.get("construct_exact_allele_identity_valid", "")))
        check(truth(row["identity_exact"]) == exact, "independent_identity_eligibility", row["assay_id"])
        if exact:
            key = ":".join(source[k] for k in ("grch38_chrom", "grch38_pos", "grch38_ref", "grch38_alt"))
            check(key == source["canonical_variant_id"] == row["variant_id"] == row["benchmark_variant_key"], "exact_frozen_identity_key", row["assay_id"])
            exact_ids.add(key)
        else:
            check(row["variant_id"] == "", "unresolved_not_promoted", row["assay_id"])
        positive = source["experimental_state"] == "positive"
        in_scope = positive and source["mechanism_in_model_scope"] in {"yes", "partial"}
        reporter = source["assay_class"] in REPORTERS
        endogenous = source["assay_class"] == "endogenous_allele_editing"
        for column, expected in (("experimental_positive", positive), ("positive_in_scope", in_scope),
                                 ("reporter_assay", reporter), ("endogenous_assay", endogenous),
                                 ("reporter_positive", in_scope and reporter), ("endogenous_positive", in_scope and endogenous)):
            check(truth(row[column]) == expected, "independent_observation_taxonomy", (row["assay_id"], column))
        higher = source.get("higher_activity_grch38_allele", "")
        expected_sign = 1 if higher == source["grch38_alt"] and higher else (-1 if higher == source["grch38_ref"] and higher else 0)
        eligible_direction = (exact and in_scope and (reporter or endogenous)
                              and truth(source.get("direction_identity_resolved", "")) and expected_sign != 0
                              and source.get("reported_direction_alt_minus_ref", "") == str(expected_sign))
        check(truth(row["enhancer_direction_eligible"]) == eligible_direction, "strict_direction_eligibility", row["assay_id"])
        if eligible_direction:
            check(row["enhancer_expected_sign"] == str(expected_sign), "strict_direction_allele_sign", row["assay_id"])
            directional[row["variant_id"]].append(expected_sign)
        check(not truth(row["h3k27me3_direction_eligible"]) and row["h3k27me3_expected_sign"] == "", "no_unsupported_H3_sign_inversion", row["assay_id"])
        check(row["benchmark_variant_key"] in original_variants, "frozen_reporting_group", row["assay_id"])
        groups[row["benchmark_variant_key"]].append(row)
        state_counts[source["experimental_state"]] += 1
        if in_scope:
            scope_counts[source["mechanism_in_model_scope"]] += 1
    check(exact_ids == {r["variant_id"] for r in scoring_rows}, "all_exact_identities_scored_without_outcome_selection")
    for row in variants:
        key = row["benchmark_variant_key"]
        original = original_variants[key]
        members = groups[key]
        check(int(row["observation_count"]) == int(original["assay_rows"]) == len(members), "variant_observation_accounting", key)
        check(truth(row["positive_in_scope"]) == truth(original["any_in_scope_positive_assay"]), "frozen_R015_positive_flag", key)
        for column in ("experimental_positive", "positive_in_scope", "reporter_positive", "endogenous_positive"):
            check(truth(row[column]) == any(truth(r[column]) for r in members), "independent_variant_taxonomy", (key, column))
        signs = set(directional[row["variant_id"]])
        expected = str(next(iter(signs))) if len(signs) == 1 else ""
        check(truth(row["enhancer_direction_conflict"]) == (len(signs) > 1), "variant_direction_consensus", key)
        check(truth(row["enhancer_direction_eligible"]) == (len(signs) == 1), "variant_direction_consensus", key)
        check(row["enhancer_expected_sign"] == row["enhancer_direction_consensus_sign"] == expected, "variant_direction_consensus", key)
    for row in contextual:
        original = snapshot["tables"][CONTEXTUAL]["rows"][int(row["snapshot_row_index"])]
        check(row["source_table"] == CONTEXTUAL and row["frozen_row_json_sha256"] == row_digest(original), "original_contextual_row_hash", row["observation_id"])
        check(not truth(row["identity_exact"]) and not truth(row["sequence_scoring_population"]) and row["variant_id"] == "", "contextual_not_promoted_to_variant", row["observation_id"])
        check(truth(row["experimental_positive_context"]) == (original["experimental_state"] == "positive"), "contextual_state_preserved", row["observation_id"])
    check(not (STAGE / "predictions").exists(), "pre_prediction_chronology")
    check(not any(name.startswith(("tensorflow", "keras")) for name in sys.modules), "no_model_runtime_import")
    return {
        "n_actual_variant_pairs_checked": len(scoring_rows), "n_actual_allele_sequences_checked": 2 * len(scoring_rows),
        "shape": list(sequences.shape), "dtype": str(sequences.dtype), "actual_variant_types": dict(allele_types),
        "chromosome_counts": dict(chromosomes), "actual_nucleotide_counts_both_alleles": dict(nucleotide_counts),
        "all_alleles_independently_reconstructed_from_hg38": True,
        "actual_sequence_row_digest": row_digest(sequence_details),
        "master_observations": len(observations), "contextual_observations": len(contextual),
        "original_benchmark_groups": len(variants), "exact_unique_variants": len(exact_ids),
        "experimental_state_observations": dict(state_counts), "positive_observation_scope_counts": dict(scope_counts),
        "in_scope_positive_groups": sum(truth(r["positive_in_scope"]) for r in variants),
        "exact_scorable_in_scope_positive_variants": sum(truth(r["positive_in_scope"]) and r["variant_id"] in exact_ids for r in variants),
        "strict_enhancer_direction_observations": sum(map(len, directional.values())),
        "strict_enhancer_direction_variants": sum(len(set(x)) == 1 for x in directional.values()),
        "strict_H3K27me3_direction_observations": 0,
        "files_verified": [record(STAGE / name) for name in (
            "inputs/input_manifest.json", "inputs/scorable_variants.tsv", "inputs/allele_sequences.npy",
            "inputs/sequence_qc.tsv", "inputs/observation_annotations.tsv", "inputs/variant_annotations.tsv",
            "inputs/contextual_annotations.tsv", "inputs/benchmark_snapshot.json.gz", "provenance/annotation_rules.json")],
        "reference_dependencies_independently_rehashed": sources,
    }


def main():
    receipt = STAGE / "provenance/independent_input_validation.json"
    if receipt.exists():
        raise RuntimeError("Independent validation receipt already exists; refusing overwrite")
    started = time.monotonic()
    result = {"stage": STAGE.name, "started_utc": datetime.now(timezone.utc).isoformat(),
              "script": record(Path(__file__).resolve()), "argv": sys.argv,
              "command_summary": "PYTHONHASHSEED=104729 bash scripts/runtime.sh scripts/validate_external_inputs_independently.py; CPU-only NumPy/stdlib",
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "executable": sys.executable, "platform": platform.platform(),
                              "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED")},
              "independence": "No builder, annotator, scorer, TensorFlow, Keras, or prior-stage execution imports; independent contiguous-span stdlib FAI reader and direct snapshot taxonomy derivation",
              "V2_prediction_access": False, "GPU_or_model_runtime_used": False,
              "candidate_universe_access": False, "original_benchmark_files_reopened": False,
              "scientific_inputs_modified": False}
    error = None
    try:
        result.update(audit())
        result["status"] = "PASS"
    except BaseException as exc:
        error = exc
        result.update(status="FAIL", error=repr(exc))
    result.update(completed_utc=datetime.now(timezone.utc).isoformat(),
                  check_counts_by_category=dict(CHECKS), n_checks=sum(CHECKS.values()),
                  resources={"wall_seconds": time.monotonic() - started,
                             "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    with receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "n_checks": sum(CHECKS.values()),
                      "receipt": str(receipt), "error": result.get("error")}), flush=True)
    if error:
        raise error


if __name__ == "__main__":
    main()
