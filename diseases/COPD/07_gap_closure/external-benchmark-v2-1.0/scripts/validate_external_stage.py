#!/usr/bin/env python3
"""Independent CPU-only verification of the completed frozen external stage.

Does not import scoring, evaluation, annotation, model or training modules. It
reconstructs expected arithmetic and populations from frozen snapshot fields and
stored raw network outputs. Original benchmark files are hashed as opaque bytes,
never parsed again. Original checkpoints are hashed, never loaded or copied.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
import traceback

import numpy as np

STAGE = Path(__file__).resolve().parents[1]
REPO = STAGE.parents[3]
BASELINE = "c9515770db4dc38064d4c1fd7fa323788040ec37"
FREEZE_SHA = "d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830"
CONTEXTS = ("enhancer", "h3k27me3")
SEEDS = (104729, 130363, 155921)
THRESHOLDS = {"enhancer": float.fromhex("0x1.7f39710000001p-1"),
              "h3k27me3": float.fromhex("0x1.8a0a12aaaaaacp-1")}
REPORTER = {"MPRA_allele_effect", "conventional_reporter"}
AUDIT_ROWS = []


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def record(path):
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size,
            "sha256": digest(path)}


def check(value, category, key="", detail=""):
    passed = bool(value)
    AUDIT_ROWS.append({"category": category, "key": str(key),
                       "status": "PASS" if passed else "FAIL", "detail": str(detail)})
    if not passed:
        raise AssertionError(f"{category}: {key}: {detail}")


def true(value):
    if value is True or value == "True" or value == "true" or value == "1" or value == 1:
        return True
    if value is False or value in ("False", "false", "0", "", None, 0):
        return False
    raise ValueError(f"Invalid Boolean {value!r}")


def read_json(relative):
    with (STAGE / relative).open() as handle:
        return json.load(handle)


def read_tsv(relative):
    with (STAGE / relative).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    check(all(None not in row for row in rows), "valid_TSV_schema", relative)
    return rows


def file_checks(base, records, category):
    names = [r["path"] for r in records]
    check(len(names) == len(set(names)), category + "_unique_paths")
    for r in records:
        path = base / r["path"]
        check(not Path(r["path"]).is_absolute() and ".." not in Path(r["path"]).parts,
              "scoped_manifest_path", r["path"])
        check(path.stat().st_size == r["bytes"] and digest(path) == r["sha256"], category, r["path"])


def scalar_match(text, expected):
    return text == "" if expected is None else float(text) == expected


def fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def signed(value):
    return int(value > 0) - int(value < 0)


def json_row_hash(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode()).hexdigest()


def exact_identity(row):
    fields = [row.get(k, "") for k in ("grch38_chrom", "grch38_pos", "grch38_ref", "grch38_alt")]
    return (row["identity_status"] == "exact" and true(row["reference_verified"])
            and all(fields) and row["canonical_variant_id"] == ":".join(fields)
            and row.get("exact_identity_eligible", "").lower() != "false"
            and row.get("construct_exact_allele_identity_valid", "").lower() != "false")


def direction_sign(row):
    if not exact_identity(row) or row["experimental_state"] != "positive":
        return None
    if row["mechanism_in_model_scope"] not in {"yes", "partial"}:
        return None
    if row["assay_class"] not in REPORTER | {"endogenous_allele_editing"}:
        return None
    if not true(row["direction_identity_resolved"]):
        return None
    ref, alt = row["grch38_ref"], row["grch38_alt"]
    higher, reported = row["higher_activity_grch38_allele"], row["reported_direction_alt_minus_ref"]
    if not higher or higher not in {ref, alt} or reported not in {"-1", "1"}:
        return None
    if row["allele_transform"] != "identity" or {row["tested_allele1"], row["tested_allele2"]} != {ref, alt}:
        raise AssertionError("Eligible direction lacks unchanged experimental nucleotide mapping")
    sign = 1 if higher == alt else -1
    if int(reported) != sign:
        raise AssertionError("Experimental activity allele and ALT-minus-REF direction disagree")
    return sign


def independent_arithmetic(raw):
    if raw.dtype != np.float64 or raw.ndim != 5 or raw.shape[:2] != (2, 3) or raw.shape[-2:] != (2, 4):
        raise ValueError("Invalid raw array contract")
    if not np.isfinite(raw).all() or np.any(raw < 0) or np.any(raw > 1):
        raise ValueError("Invalid probabilities")
    if not np.array_equal(raw, raw.astype(np.float32).astype(np.float64)):
        raise ValueError("Raw network outputs are not exact float32 promotions")
    seed = (raw[..., 0] + raw[..., 1]) * 0.5
    audit_seed = (raw[..., 2] + raw[..., 3]) * 0.5
    ensemble = np.add.reduce(seed, axis=1, dtype=np.float64) / np.float64(3)
    audit_ensemble = np.add.reduce(audit_seed, axis=1, dtype=np.float64) / np.float64(3)
    return seed, audit_seed, ensemble, audit_ensemble


def expected_v1_region(row, orientation):
    available = true(row["v1_model_evaluable" if orientation == "forward" else "rc_model_evaluable"])
    if not available:
        return None
    calls = {}
    for context, prefix in (("enhancer", "enhancer"), ("h3k27me3", "silencer")):
        values = [float(row[f"{prefix}_{orientation}_{role}_score"]) for role in ("ref", "alt")]
        threshold = float(row[f"{prefix}_region_threshold"])
        if not all(np.isfinite(x) for x in values + [threshold]):
            raise ValueError("Nonfinite frozen V1 probability")
        region = max(values)
        if region != float(row[f"{prefix}_{orientation}_region_score"]):
            raise ValueError("Frozen V1 region score inconsistent")
        call = region >= threshold
        if call != true(row[f"{prefix}_{orientation}_region_gate"]):
            raise ValueError("Frozen V1 region gate inconsistent")
        calls[context] = call
    calls["union"] = any(calls.values())
    return calls


def direction_totals(records):
    counts = Counter(records)
    if set(counts) - {"concordant", "discordant", "tie"}:
        raise ValueError("Unsupported direction state")
    n = sum(counts.values())
    non_tied = n - counts["tie"]
    return {"denominator": n, "concordant": counts["concordant"], "discordant": counts["discordant"],
            "tied": counts["tie"], "strict_concordance_fraction": fraction(counts["concordant"], n),
            "non_tied_denominator": non_tied, "secondary_non_tied_fraction": fraction(counts["concordant"], non_tied)}


def timestamp(value):
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        raise ValueError("Naive timestamp")
    return result


def validate_manifests_and_chronology():
    freeze = read_json("provenance/prospective_freeze.json")
    check(digest(STAGE / "provenance/prospective_freeze.json") == FREEZE_SHA, "prospective_freeze_exact_hash")
    file_checks(STAGE, freeze["files"], "prospective_payload_hash")
    annotation = read_json("provenance/annotation_freeze.json")
    input_manifest = read_json("inputs/input_manifest.json")
    inference = read_json("predictions/inference_manifest.json")
    evaluation = read_json("results/evaluation_manifest.json")
    for name, manifest in (("annotation", annotation), ("inputs", input_manifest),
                           ("inference", inference), ("evaluation", evaluation)):
        check(manifest["status"] == "PASS", "manifest_status", name)
        file_checks(STAGE, manifest["files"], name + "_payload_hash")
    completed = read_json("predictions/completed.json")
    file_checks(STAGE, [completed["inference_manifest"]], "completed_inference_receipt_hash")
    opening = read_json("provenance/benchmark_opening_event.json")
    ingestion = read_json("provenance/benchmark_ingestion_complete.json")
    file_checks(STAGE, [ingestion["snapshot"], ingestion["single_open_event"]], "opening_snapshot_hash")
    started = read_json("predictions/started.json")
    before = read_json("predictions/verified_before_inference.json")
    evaluation_started = read_json("results/evaluation_started.json")
    times = [freeze["created_utc"], opening["opened_utc"], ingestion["completed_utc"],
             annotation["created_utc"], started["started_utc"], before["verified_utc"],
             completed["resources"]["finished_utc"], evaluation_started["started_utc"], evaluation["completed_utc"]]
    parsed = list(map(timestamp, times))
    check(all(a < b for a, b in zip(parsed, parsed[1:])), "strict_prospective_open_annotation_prediction_chronology", detail=json.dumps(times))
    check(freeze["before_any_benchmark_outcome_or_label_opening"] is True
          and annotation["before_any_V2_external_prediction"] is True
          and opening["single_open_event"] is True, "prospective_single_open_flags")
    check(inference["ordered_seeds"] == list(SEEDS) and inference["contexts"] == list(CONTEXTS)
          and inference["all_six_seeds_contexts_complete"] is True, "complete_frozen_six_model_population")
    check(inference["training_or_recalibration"] is False
          and inference["candidate_universe_access"] is False
          and inference["experimental_outcome_tables_parsed"] is False
          and inference["prior_training_or_internal_test_entrypoint_execution"] is False, "inference_firewall_receipts")
    check(all(evaluation[k] is False for k in ("classification_metrics_calculated", "significance_tests_performed",
              "internal_reference_distribution_constructed", "candidate_universe_accessed", "thresholds_changed", "models_or_seeds_changed")),
          "evaluation_scope_receipts")
    for r in before["script_files"]:
        file_checks(STAGE, [r], "executed_script_unchanged")
    return inference, evaluation, before


def validate_predictions(inputs, inference, before):
    n = len(inputs)
    raw = np.load(STAGE / "predictions/raw_orientation_probabilities.npy", allow_pickle=False)
    check(raw.shape == (2, 3, n, 2, 4), "raw_prediction_shape")
    arrays = independent_arithmetic(raw)
    names = ("seed_symmetric_probabilities.npy", "seed_independent_RC_probabilities.npy",
             "ensemble_probabilities.npy", "ensemble_independent_RC_probabilities.npy")
    for name, expected in zip(names, arrays):
        actual = np.load(STAGE / "predictions" / name, allow_pickle=False)
        check(actual.dtype == np.float64 and np.array_equal(expected, actual), "independent_exact_array_arithmetic", name)
    seed, audit_seed, ensemble, audit_ensemble = arrays
    check(inference["n_variants"] == n and inference["n_context_variant_scores"] == 2 * n
          and inference["n_seed_variant_scores"] == 6 * n, "complete_score_population")
    seed_rows = read_tsv("predictions/seed_allele_scores.tsv")
    score_rows = read_tsv("predictions/variant_context_scores.tsv")
    positions = {row["variant_id"]: i for i, row in enumerate(inputs)}
    seed_keys, score_keys = set(), set()
    for rows, is_seed in ((seed_rows, True), (score_rows, False)):
        for row in rows:
            key = (row["variant_id"], row["context"], int(row["seed"])) if is_seed else (row["variant_id"], row["context"])
            seen = seed_keys if is_seed else score_keys
            check(key not in seen, "unique_prediction_key", key)
            seen.add(key)
            c, i = CONTEXTS.index(row["context"]), positions[row["variant_id"]]
            s = SEEDS.index(int(row["seed"])) if is_seed else None
            pair = seed[c, s, i] if is_seed else ensemble[c, i]
            repeated = audit_seed[c, s, i] if is_seed else audit_ensemble[c, i]
            ref, alt = map(float, pair)
            expected = {"ref_score": ref, "alt_score": alt, "delta": alt - ref,
                        "abs_delta": abs(alt - ref), "max_allele_score": max(ref, alt), "threshold": THRESHOLDS[row["context"]]}
            check(all(float(row[k]) == value and row[k + "_hex"] == value.hex()
                      and float.fromhex(row[k + "_hex"]) == float(row[k]) for k, value in expected.items()),
                  "score_decimal_hex_delta_arithmetic", key)
            threshold = THRESHOLDS[row["context"]]
            check(int(row["delta_sign"]) == signed(alt - ref)
                  and true(row["ref_region_call"]) == (ref >= threshold)
                  and true(row["alt_region_call"]) == (alt >= threshold)
                  and true(row["either_allele_region_call"]) == (max(ref, alt) >= threshold)
                  and row["region_call_operator"] == ">=", "exact_sign_and_inclusive_threshold", key)
            check(all(row[k] == inputs[i][k] for k in ("variant_id", "chrom", "pos1", "ref", "alt")), "prediction_identity_join", key)
            for a, allele in enumerate(("ref", "alt")):
                suffix = "q_independent_RC" if is_seed else "ensemble_independent_RC"
                name = allele + "_" + suffix
                check(float(row[name]) == repeated[a] and row[name + "_hex"] == float(repeated[a]).hex(), "independent_RC_table_value", (key, allele))
                if is_seed:
                    for o, orientation in enumerate(("forward", "nucleotide_RC", "independent_RC_first", "independent_forward_second")):
                        name = f"{allele}_p_{orientation}"
                        value = float(raw[c, s, i, a, o])
                        check(float(row[name]) == value and row[name + "_hex"] == value.hex(), "raw_network_table_binding", (key, allele, orientation))
    check(len(seed_keys) == 6 * n and len(score_keys) == 2 * n, "all_prediction_keys_present")
    audit_manifest = read_json("predictions/real_network_invariance.json")
    check(audit_manifest["n_seed_audits"] == 6 and audit_manifest["n_ensemble_audits"] == 2
          and audit_manifest["real_network"] is True and audit_manifest["phase_I_and_phase_II_independently_repeated"] is True,
          "six_seed_two_ensemble_real_network_audit_count")
    max_residual = 0.0
    for audit in audit_manifest["audit_rows"]:
        c = CONTEXTS.index(audit["context"])
        if audit["unit"] == "seed":
            s = SEEDS.index(int(audit["seed"]))
            a, b = seed[c, s], audit_seed[c, s]
            filename = f"{audit['context']}_seed{audit['seed']}_RC_invariance.json"
        else:
            a, b = ensemble[c], audit_ensemble[c]
            filename = f"{audit['context']}_ensemble_RC_invariance.json"
        residual = np.abs(a - b)
        tolerance = 1e-6 + 1e-6 * np.abs(b)
        max_residual = max(max_residual, float(residual.max()))
        check(audit == read_json("predictions/" + filename), "individual_combined_RC_audit_equality", filename)
        check(audit["status"] == "PASS" and audit["n_probabilities"] == a.size
              and audit["n_failed"] == int(np.sum(residual > tolerance)) == 0
              and audit["max_absolute_difference"] == float(residual.max())
              and audit["max_tolerance_normalized_difference"] == float((residual / tolerance).max())
              and audit["atol"] == audit["rtol"] == 1e-6, "independently_recomputed_RC_tolerance", filename)
        for label, array in (("original", a), ("independent_RC", b)):
            check(audit[label + "_array_sha256"] == hashlib.sha256(array.tobytes()).hexdigest(), "RC_audit_array_hash", (filename, label))
    check(audit_manifest["maximum_absolute_residual"] == max_residual, "combined_RC_maximum")
    phase_one = read_json("predictions/phase_I_provenance.json")
    check(phase_one["independent_network_calls_per_sequence"] == 4
          and phase_one["orientation_sequence_evaluations"] == n * 2 * 4
          and phase_one["audit_recomputes_phase_I_features"] is True
          and phase_one["representation_axis_reversal"] is False
          and phase_one["training"] is False and phase_one["optimizer"] is None
          and phase_one["phase_I_batch_size"] == 32 and phase_one["representation"] == "final4560sigmoid",
          "real_phase_I_call_provenance")
    feature_arrays = [np.load(STAGE / "predictions" / name, allow_pickle=False, mmap_mode="r")
                      for name in phase_one["feature_names_in_orientation_order"]]
    check(len(feature_arrays) == 4 and all(a.shape == (2 * n, 4560) and a.dtype == np.float32 for a in feature_arrays), "phase_I_four_orientation_geometry")
    for (a, b), field in (((0, 3), "repeated_forward_max_abs_difference"), ((1, 2), "repeated_RC_max_abs_difference")):
        maximum = max(float(np.max(np.abs(feature_arrays[a][i:i + 128] - feature_arrays[b][i:i + 128]))) for i in range(0, 2 * n, 128))
        check(maximum == phase_one[field], "phase_I_repeat_residual", field)
    check(before["n_variants"] == n and len(before["sequence_records"]) == 2 * n, "preinference_complete_population")
    environment = read_json("predictions/environment.json")
    check(environment["python"] == "3.13.0" and all(environment["packages"][k] == v for k, v in
          {"numpy": "2.5.0", "tensorflow": "2.20.0", "keras": "3.14.1"}.items()), "frozen_numerical_runtime")
    return ensemble, score_rows


def validate_annotations(snapshot, inputs):
    master = snapshot["tables"]["COPD-V2-BENCH-R003_frozen_benchmark_master.tsv"]["rows"]
    unique = snapshot["tables"]["COPD-V2-BENCH-R015_unique_variant_summary.tsv"]["rows"]
    contextual = snapshot["tables"]["COPD-V2-BENCH-R004B_contextual_and_excluded_evidence.tsv"]["rows"]
    observations = read_tsv("inputs/observation_annotations.tsv")
    variants = read_tsv("inputs/variant_annotations.tsv")
    contexts = read_tsv("inputs/contextual_annotations.tsv")
    identity_map = read_tsv("inputs/identity_map.tsv")
    check(len(master) == len(observations) == len(identity_map) == 14025
          and len(unique) == len(variants) == 1731 and len(contextual) == len(contexts) == 267, "complete_original_observation_and_variant_counts")
    unique_by_key = {r["benchmark_variant_key"]: r for r in unique}
    source_groups, direction_obs = defaultdict(list), {}
    positions = {r["variant_id"] for r in inputs}
    for i, (source, annotation, identity) in enumerate(zip(master, observations, identity_map)):
        key = source["canonical_variant_id"] or "unresolved:" + source["rsid"] + ":" + "/".join(sorted((source["tested_allele1"], source["tested_allele2"])))
        exact = bool(exact_identity(source))
        variant = source["canonical_variant_id"] if exact else ""
        positive = source["experimental_state"] == "positive"
        scoped = source["mechanism_in_model_scope"] in ("yes", "partial")
        assay = source["assay_class"]
        check(annotation["assay_id"] == source["assay_id"] == identity["assay_id"]
              and int(annotation["snapshot_row_index"]) == i
              and annotation["benchmark_variant_key"] == identity["benchmark_variant_key"] == key
              and annotation["variant_id"] == identity["variant_id"] == variant,
              "master_annotation_identity_alignment", source["assay_id"])
        check(annotation["frozen_row_json_sha256"] == json_row_hash(source), "exact_original_row_evidence_binding", source["assay_id"])
        original_fields = set(source) & set(annotation)
        check(all(source[k] == annotation[k] for k in original_fields), "unchanged_retained_source_fields", source["assay_id"])
        expected_flags = {"identity_exact": exact, "experimental_positive": positive,
                          "positive_in_scope": positive and scoped,
                          "reporter_positive": positive and scoped and assay in REPORTER,
                          "endogenous_positive": positive and scoped and assay == "endogenous_allele_editing",
                          "out_of_model_positive": positive and not scoped,
                          "preexisting_known_case": source["rsid"] == "rs2013701"}
        check(all(true(annotation[k]) == v for k, v in expected_flags.items()), "independent_observation_population_flags", source["assay_id"])
        direction = direction_sign(source)
        check(true(annotation["enhancer_direction_eligible"]) == (direction is not None)
              and (int(annotation["enhancer_expected_sign"]) == direction if direction is not None else annotation["enhancer_expected_sign"] == "")
              and not true(annotation["h3k27me3_direction_eligible"])
              and annotation["h3k27me3_expected_sign"] == "", "strict_observation_direction_mapping", source["assay_id"])
        if direction is not None:
            direction_obs[source["assay_id"]] = (variant, direction, source)
        source_groups[key].append(source)
        check((identity["sequence_status"] == "ok") == (variant in positions), "identity_sequence_population_alignment", source["assay_id"])
    for source, annotation in zip(contextual, contexts):
        check(annotation["frozen_row_json_sha256"] == json_row_hash(source)
              and all(annotation[k] == value for k, value in source.items())
              and not true(annotation["sequence_scoring_population"]), "all_contextual_evidence_preserved", annotation["observation_id"])
    by_key, direction_consensus, conflicts = {}, {}, set()
    for annotation in variants:
        key = annotation["benchmark_variant_key"]
        sources = source_groups[key]
        original = unique_by_key[key]
        exact_ids = {r["canonical_variant_id"] for r in sources if exact_identity(r)}
        check(len(exact_ids) <= 1, "unique_genomic_identity_per_frozen_key", key)
        exact_id = next(iter(exact_ids), "")
        scoped = [r for r in sources if r["experimental_state"] == "positive" and r["mechanism_in_model_scope"] in ("yes", "partial")]
        flags = {"positive_in_scope": bool(scoped), "experimental_positive": any(r["experimental_state"] == "positive" for r in sources),
                 "reporter_positive": any(r["assay_class"] in REPORTER for r in scoped),
                 "endogenous_positive": any(r["assay_class"] == "endogenous_allele_editing" for r in scoped),
                 "splice_only": any(r["assay_class"] == "splicing" for r in sources) and not scoped,
                 "out_of_model_positive": any(r["experimental_state"] == "positive" and r["mechanism_in_model_scope"] not in ("yes", "partial") for r in sources)}
        check(annotation["variant_id"] == exact_id == original["canonical_variant_id"]
              and true(original["any_in_scope_positive_assay"]) == flags["positive_in_scope"]
              and int(annotation["observation_count"]) == len(sources) == int(original["assay_rows"])
              and annotation["frozen_R015_row_json_sha256"] == json_row_hash(original)
              and all(true(annotation[k]) == bool(value) for k, value in flags.items()), "independent_variant_population_flags", key)
        signs = {value for row in sources if (value := direction_sign(row)) is not None}
        conflict = len(signs) > 1
        eligible = len(signs) == 1
        check(true(annotation["enhancer_direction_conflict"]) == conflict
              and true(annotation["enhancer_direction_eligible"]) == eligible
              and not true(annotation["h3k27me3_direction_eligible"]), "independent_variant_direction_consensus", key)
        if conflict:
            conflicts.add(exact_id)
        if eligible:
            consensus = next(iter(signs))
            check(int(annotation["enhancer_expected_sign"]) == consensus, "unique_direction_sign", key)
            direction_consensus[exact_id] = consensus
        by_key[key] = dict(annotation)
    check(set(positions) == {v["variant_id"] for v in variants if v["variant_id"]}
          and len(positions) == 1710, "all_exact_identities_scored_regardless_outcome_scope")
    check(len(direction_obs) == 11 and len(direction_consensus) == 8 and not conflicts, "strict_direction_denominators_11_observations_8_variants_no_conflicts")
    return observations, by_key, unique_by_key, direction_obs, direction_consensus


def validate_evaluation(snapshot, inputs, ensemble, observations, annotations, unique, direction_obs, consensus):
    positions = {r["variant_id"]: i for i, r in enumerate(inputs)}
    methods = {"V1_forward": {}, "V1_RC": {}, "V2_C_symmetric": {}}
    for key, row in unique.items():
        for name, orientation in (("V1_forward", "forward"), ("V1_RC", "rc")):
            result = expected_v1_region(row, orientation)
            if result is not None:
                methods[name][key] = result
    for key, index in positions.items():
        calls = {c: bool(np.max(ensemble[j, index]) >= THRESHOLDS[c]) for j, c in enumerate(CONTEXTS)}
        calls["union"] = any(calls.values())
        methods["V2_C_symmetric"][key] = calls
    populations = {name: {k for k, row in annotations.items() if true(row[field])} for name, field in
                   (("all_in_scope_positive", "positive_in_scope"), ("reporter_positive", "reporter_positive"),
                    ("endogenous_positive", "endogenous_positive"))}
    positive = populations["all_in_scope_positive"]
    primary = positive & set(positions)
    check(len(positive) == 40 and len(primary) == 39 and len(positive - primary) == 1, "primary_recovery_denominator_40_to39")
    recovery = read_tsv("results/region_recovery_summary.tsv")
    check(len(recovery) == 27 and len({(r["subset"], r["method"], r["context"]) for r in recovery}) == 27, "complete_recovery_grid")
    for row in recovery:
        population, mapping = populations[row["subset"]], methods[row["method"]]
        usable = population.intersection(mapping)
        recovered = sum(mapping[k][row["context"]] for k in usable)
        check(int(row["positive_keys_before_method_coverage"]) == len(population)
              and int(row["evaluable_positive_denominator"]) == len(usable)
              and int(row["recovered"]) == recovered
              and int(row["unevaluable_count"]) == len(population) - len(usable)
              and scalar_match(row["recovery_fraction"], fraction(recovered, len(usable))), "independent_region_recovery", (row["subset"], row["method"], row["context"]))
    common = positive.intersection(*[set(v) for v in methods.values()])
    common_rows = read_tsv("results/three_method_common_identity_recovery.tsv")
    check(len(common_rows) == 9 and len({(r["method"], r["context"]) for r in common_rows}) == 9, "complete_three_method_grid")
    for row in common_rows:
        recovered = sum(methods[row["method"]][k][row["context"]] for k in common)
        check(int(row["denominator"]) == len(common) and int(row["recovered"]) == recovered
              and scalar_match(row["fraction"], fraction(recovered, len(common))), "three_method_common_identity_region_recovery", (row["method"], row["context"]))
    expected_transitions = {}
    expected_pairs = {}
    for baseline in ("V1_forward", "V1_RC"):
        shared = positive.intersection(methods[baseline], methods["V2_C_symmetric"])
        for context in (*CONTEXTS, "union"):
            counts = Counter()
            for key in shared:
                old, new = methods[baseline][key][context], methods["V2_C_symmetric"][key][context]
                transition = {(True, True): "retained_recovery", (False, True): "newly_recovered", (True, False): "lost", (False, False): "still_missed"}[old, new]
                expected_transitions[baseline, context, key] = (old, new, transition)
                counts[transition] += 1
            expected_pairs[baseline, context] = (len(shared), counts)
    transitions = read_tsv("results/paired_V1_V2_variant_transitions.tsv")
    check(len(transitions) == len(expected_transitions)
          and {(r["baseline"], r["context"], r["benchmark_variant_key"]) for r in transitions} == set(expected_transitions), "all_paired_transition_rows")
    for row in transitions:
        key = (row["baseline"], row["context"], row["benchmark_variant_key"])
        old, new, transition = expected_transitions[key]
        check(true(row["V1_region_call"]) == old and true(row["V2_region_call"]) == new and row["transition"] == transition,
              "V1_region_only_paired_transition_not_combined_delta_call", key)
    pairs = read_tsv("results/paired_V1_V2_region_comparison.tsv")
    check(len(pairs) == 6, "paired_summary_count")
    for row in pairs:
        n, counts = expected_pairs[row["baseline"], row["context"]]
        check(int(row["shared_positive_denominator"]) == n and all(int(row[name]) == counts[name] for name in
              ("retained_recovery", "newly_recovered", "lost", "still_missed")), "paired_summary_counts", (row["baseline"], row["context"]))
    expansion = read_tsv("results/newly_V2_scorable_positive_variants.tsv")
    check({r["benchmark_variant_key"] for r in expansion} == primary - set(methods["V1_forward"]), "coverage_expansion_not_paired_gain")
    for row in expansion:
        check(all(true(row[c + "_recovered"]) == methods["V2_C_symmetric"][row["benchmark_variant_key"]][c]
                  for c in (*CONTEXTS, "union")), "new_coverage_region_calls", row["benchmark_variant_key"])
    all_results = read_tsv("results/all_variant_evaluation.tsv")
    check(len(all_results) == 1731 and {r["benchmark_variant_key"] for r in all_results} == set(annotations), "all_variants_retained_in_results")
    for row in all_results:
        key, variant = row["benchmark_variant_key"], row["variant_id"]
        for method in ("V1_forward", "V1_RC"):
            check(true(row[method + "_available"]) == (key in methods[method]), "V1_exact_coverage_flag", (key, method))
            for context in (*CONTEXTS, "union"):
                field = method + "_" + context + "_region_call"
                check(true(row[field]) == methods[method][key][context] if key in methods[method] else row[field] == "", "V1_region_only_variant_calls", (key, method, context))
        check(true(row["v2_scorable"]) == (variant in positions), "V2_variant_coverage", key)
        if variant in positions:
            i = positions[variant]
            for c, context in enumerate(CONTEXTS):
                ref, alt = map(float, ensemble[c, i])
                values = {"ref_score": ref, "alt_score": alt, "delta": alt - ref, "abs_delta": abs(alt - ref),
                          "delta_sign": signed(alt - ref), "max_allele_score": max(ref, alt)}
                check(all(float(row["V2_" + context + "_" + k]) == v for k, v in values.items()), "variant_table_numerics", (variant, context))
            check(true(row["V2_union_region_call"]) == methods["V2_C_symmetric"][variant]["union"], "variant_union_call", variant)
    missing = read_tsv("results/unevaluable_variants.tsv")
    check({r["benchmark_variant_key"] for r in missing} == set(annotations) - set(positions) and len(missing) == 21, "explicit_21_unresolved_variants")
    check({r["benchmark_variant_key"] for r in read_tsv("results/region_recovery_per_variant.tsv")} == positive,
          "all_positive_variants_listed_including_unresolved")
    check(read_tsv("results/contextual_and_excluded_evidence.tsv") == read_tsv("inputs/contextual_annotations.tsv"),
          "contextual_and_excluded_outputs_exactly_retained")
    check(read_tsv("results/frozen_denominator_rules.tsv") == snapshot["tables"]["COPD-V2-BENCH-R005_assayed_denominator_audit.tsv"]["rows"], "source_denominator_audit_unchanged")
    direction_rows = read_tsv("results/direction_resolved_allelic_effects.tsv")
    consensus_rows = read_tsv("results/direction_unique_variant_consensus.tsv")
    check(len(direction_rows) == 11 and {r["assay_id"] for r in direction_rows} == set(direction_obs), "exact_11_eligible_direction_observations")
    check(len(consensus_rows) == 8 and {r["variant_id"] for r in consensus_rows} == set(consensus), "exact_8_nonconflicting_consensus_variants")
    for row in direction_rows + consensus_rows:
        variant = row["variant_id"]
        expected = direction_obs[row["assay_id"]][1] if "assay_id" in row else consensus[variant]
        delta = float(ensemble[0, positions[variant], 1] - ensemble[0, positions[variant], 0])
        sign = signed(delta)
        outcome = "tie" if sign == 0 else "concordant" if sign == expected else "discordant"
        check(row["context"] == "enhancer" and int(row["expected_ALT_minus_REF_sign"]) == expected
              and int(row["predicted_sign"]) == sign and float(row["delta"]) == delta and row["concordance"] == outcome,
              "unrounded_direction_with_strict_ties", row.get("assay_id", variant))
    summaries = read_tsv("results/direction_concordance_summary.tsv")
    check(len(summaries) == 8, "direction_summary_grid")
    for row in summaries:
        source = consensus_rows if row["unit"] == "unique_variant" else direction_rows
        selected = [r for r in source if r["context"] == row["context"]]
        if row["subset"] == "reporter":
            selected = [r for r in selected if r["assay_class"] in REPORTER]
        elif row["subset"] == "endogenous_editing":
            selected = [r for r in selected if r["assay_class"] == "endogenous_allele_editing"]
        totals = direction_totals([r["concordance"] for r in selected])
        check(all(scalar_match(row[k], value) for k, value in totals.items()), "direction_counts_and_strict_nontied_fractions", (row["unit"], row["subset"], row["context"]))
    locus_rows = read_tsv("results/direction_locus_summary.tsv")
    locus_groups = defaultdict(list)
    for row in direction_rows:
        locus_groups[row["context"], row["locus"]].append(row)
    check(len(locus_rows) == len(locus_groups)
          and {(r["context"], r["locus"]) for r in locus_rows} == set(locus_groups), "all_direction_locus_groups")
    for row in locus_rows:
        group = locus_groups[row["context"], row["locus"]]
        totals = direction_totals([r["concordance"] for r in group])
        check(all(scalar_match(row[k], value) for k, value in totals.items())
              and int(row["unique_variants"]) == len({r["variant_id"] for r in group}),
              "direction_locus_counts_not_independent_observations", (row["context"], row["locus"]))
    check(not read_tsv("results/direction_conflicts.tsv"), "no_direction_conflict_silently_majority_voted")
    direction_missing = read_tsv("results/direction_unevaluable_observations.tsv")
    expected_missing = {(r["assay_id"], context) for r in observations
                        if true(r["positive_in_scope"]) or r["assay_class"] == "splicing"
                        for context in CONTEXTS
                        if not true(r[context + "_direction_eligible"]) or r["variant_id"] not in positions}
    check(len(direction_missing) == len(expected_missing)
          and {(r["assay_id"], r["context"]) for r in direction_missing} == expected_missing,
          "all_direction_exclusions_and_H3_unevaluable_rows_accounted")
    distribution_sets = {"all_scorable_exact": set(positions), "in_scope_positive": primary,
                         "reporter_positive": populations["reporter_positive"] & primary,
                         "endogenous_positive": populations["endogenous_positive"] & primary}
    distributions = read_tsv("results/continuous_score_distributions.tsv")
    check(len(distributions) == 40, "continuous_distribution_grid")
    for row in distributions:
        c = CONTEXTS.index(row["context"])
        values = []
        for key in sorted(distribution_sets[row["subset"]]):
            ref, alt = map(float, ensemble[c, positions[key]])
            values.append({"ref_score": ref, "alt_score": alt, "max_allele_score": max(ref, alt),
                           "delta": alt - ref, "abs_delta": abs(alt - ref)}[row["score"]])
        quantiles = np.quantile(np.array(values, dtype=np.float64), [0, .25, .5, .75, 1], method="linear") if values else [None] * 5
        check(int(row["n"]) == len(values) and all(scalar_match(row[k], v) for k, v in zip(("minimum", "q25", "median", "q75", "maximum"), quantiles)), "independent_float64_linear_quantiles", (row["subset"], row["context"], row["score"]))
    strata = defaultdict(set)
    for row in observations:
        variant = row["variant_id"]
        if not true(row["positive_in_scope"]) or variant not in primary:
            continue
        for axis in ("study_id", "cell_context", "assay_class", "mechanism_in_model_scope", "locus"):
            strata[axis, row.get(axis, "") or "unspecified"].add(variant)
        evidence = "reporter" if row["assay_class"] in REPORTER else "endogenous_editing" if row["assay_class"] == "endogenous_allele_editing" else "other_in_scope_evidence"
        mechanism = {"reporter": "transcriptional_reporter_not_H3_specific", "endogenous_editing": "endogenous_expression_partial_scope", "other_in_scope_evidence": "TF_binding_partial_scope"}[evidence]
        strata["evidence_unit", evidence].add(variant)
        strata["regulatory_mechanism", mechanism].add(variant)
    for variant in primary:
        ref, alt = variant.split(":")[-2:]
        variant_type = "SNV" if len(ref) == len(alt) == 1 else "indel_or_multibase"
        v1_presence = "exact_forward_and_RC" if variant in methods["V1_forward"] and variant in methods["V1_RC"] else "forward_only" if variant in methods["V1_forward"] else "RC_only" if variant in methods["V1_RC"] else "absent_unavailable"
        strata["variant_type", variant_type].add(variant)
        strata["V1_presence", v1_presence].add(variant)
    stratified = read_tsv("results/source_mechanism_stratified_recovery.tsv")
    check(len(stratified) == len(strata) * 3, "complete_stratifier_grid")
    for row in stratified:
        keys = strata[row["axis"], row["stratum"]]
        recovered = sum(methods["V2_C_symmetric"][k][row["context"]] for k in keys)
        check(int(row["unique_positive_scorable_variants"]) == len(keys) and int(row["recovered"]) == recovered
              and scalar_match(row["fraction"], fraction(recovered, len(keys))), "independent_stratified_recovery", (row["axis"], row["stratum"], row["context"]))
    state_sets = defaultdict(set)
    for row in observations:
        state_sets[row["experimental_state"], row["mechanism_in_model_scope"]].add(row["benchmark_variant_key"])
    states = read_tsv("results/context_state_populations.tsv")
    check(len(states) == len(state_sets), "all_context_state_populations")
    for row in states:
        keys = state_sets[row["experimental_state"], row["mechanism_in_model_scope"]]
        check(int(row["unique_benchmark_keys"]) == len(keys) and int(row["exact_V2_scorable_keys"]) == len(keys & set(positions)),
              "context_null_positive_populations_not_universal_negatives", (row["experimental_state"], row["mechanism_in_model_scope"]))
    population_expected = {
        "original_master_observations": len(observations),
        "contextual_excluded_source_rows": len(snapshot["tables"]["COPD-V2-BENCH-R004B_contextual_and_excluded_evidence.tsv"]["rows"]),
        "all_benchmark_identity_keys": len(annotations),
        "exact_benchmark_identity_keys": sum(bool(r["variant_id"]) for r in annotations.values()),
        "exact_V2_scorable_identity_keys_all_states_mechanisms": len(positions),
        "experimentally_positive_keys_all_mechanisms_including_unresolved": sum(true(r["experimental_positive"]) for r in annotations.values()),
        "exact_experimentally_positive_keys_all_mechanisms": sum(true(r["experimental_positive"]) and bool(r["variant_id"]) for r in annotations.values()),
        "splice_only_positive_keys_out_of_model_not_false_negatives": sum(true(r["splice_only"]) for r in annotations.values()),
        "positive_keys_with_out_of_model_observations": sum(true(r["out_of_model_positive"]) for r in annotations.values()),
        "in_scope_positive_identity_keys_including_unresolved": len(positive),
        "exact_in_scope_positive_identity_keys": len(primary),
        "exact_V2_scorable_in_scope_positive_primary": len(primary),
        "exact_V2_scorable_reporter_positive": len(populations["reporter_positive"] & primary),
        "exact_V2_scorable_endogenous_positive": len(populations["endogenous_positive"] & primary),
        "V1_forward_RC_V2_common_positive": len(common),
        "strict_enhancer_direction_observations": len(direction_rows),
        "strict_enhancer_direction_unique_variants_nonconflicting": len(consensus_rows),
        "strict_H3K27me3_direction_observations": 0}
    accounting = read_tsv("results/population_accounting.tsv")
    check(len(accounting) == len(population_expected)
          and {r["population"]: int(r["count"]) for r in accounting} == population_expected,
          "complete_independent_population_accounting_table")
    evaluation = read_json("results/evaluation_manifest.json")
    check({r["population"]: int(r["count"]) for r in evaluation["population_accounting"]} == population_expected,
          "evaluation_manifest_population_summary_matches")
    for field, table in (("region_recovery", recovery), ("common_identity_recovery", common_rows),
                         ("paired_comparisons", pairs), ("direction_summaries", summaries)):
        serialized = [{k: "" if value is None else str(value) for k, value in row.items()}
                      for row in evaluation[field]]
        check(serialized == table, "evaluation_manifest_summary_matches_independently_verified_TSV", field)
    return {"n_original_observations": 14025, "n_original_identity_keys": 1731,
            "n_exact_scorable_variants": 1710, "n_unresolved_identity_keys": 21,
            "n_in_scope_positive_keys": len(positive), "n_primary_positive_variants": len(primary),
            "n_three_method_common_positive_variants": len(common),
            "enhancer_direction_observations": len(direction_rows), "enhancer_direction_consensus_variants": len(consensus_rows),
            "h3k27me3_direction_observations": 0}


def git(*args):
    return subprocess.check_output(["git", "-C", str(REPO), *args])


def validate_final_reports():
    prefinal_path = STAGE / "provenance/independent_validation_pre-final.json"
    prefinal = read_json("provenance/independent_validation_pre-final.json")
    historical_script = STAGE / "scripts/history/validate_external_stage_pre_final.py"
    check(prefinal["status"] == "PASS" and prefinal["failed_checks"] == 0,
          "pre_final_independent_validation_receipt_PASS")
    check(digest(historical_script) == prefinal["script"]["sha256"]
          and historical_script.stat().st_size == prefinal["script"]["bytes"],
          "pre_final_executed_validator_source_preserved_exactly")
    file_checks(STAGE, [prefinal["check_table"]], "pre_final_validation_check_table_unchanged")
    report_command = ("provenance/commands/write_external_report_final.completed.json"
                      if (STAGE / "provenance/commands/write_external_report_final.completed.json").exists()
                      else "provenance/commands/write_external_report.completed.json")
    completed = read_json(report_command)
    check(completed["returncode"] == 0, "report_command_success")
    file_checks(STAGE, [completed["stdout"], completed["stderr"]], "report_command_log_hash")
    with (STAGE / completed["stdout"]["path"]).open() as handle:
        report_receipt = json.load(handle)
    check(report_receipt["status"] == "PASS" and report_receipt["scientific_results_recomputed"] is False
          and report_receipt["original_benchmark_reopened"] is False, "report_generation_nonanalytic_status")
    check({r["path"] for r in report_receipt["reports"]} == {"EXTERNAL_BENCHMARK_V2_REPORT.md", "LIMITATIONS.md"},
          "complete_expected_reports")
    file_checks(STAGE, report_receipt["reports"], "final_reports_exact_hashes")
    evaluation = read_json("results/evaluation_manifest.json")
    headline = report_receipt["report_headline_metrics"]
    primary = {r["context"]: r for r in evaluation["region_recovery"]
               if r["method"] == "V2_C_symmetric" and r["subset"] == "all_in_scope_positive"}
    common_union = {r["method"]: r for r in evaluation["common_identity_recovery"] if r["context"] == "union"}
    strict = next(r for r in evaluation["direction_summaries"] if r["context"] == "enhancer" and r["subset"] == "all_strict_observations")
    consensus = next(r for r in evaluation["direction_summaries"] if r["context"] == "enhancer" and r["unit"] == "unique_variant")
    check(headline["primary_region_recovery"] == primary and headline["common_identity_union"] == common_union
          and headline["strict_enhancer_direction"] == strict and headline["unique_variant_enhancer_direction"] == consensus,
          "report_headline_metrics_match_independently_verified_results")
    check(headline["exact_scorable_variants"] == 1710 and headline["exact_scorable_in_scope_positive_variants"] == 39
          and headline["unresolved_reporting_groups"] == 21
          and headline["coverage_expansion_positive_count"] == evaluation["new_V2_coverage_positive_count"],
          "report_population_and_coverage_headlines")
    report_text = (STAGE / "EXTERNAL_BENCHMARK_V2_REPORT.md").read_text()
    limitations_text = (STAGE / "LIMITATIONS.md").read_text()
    check("provenance/independent_final_validation.json" in report_text
          and "independent_validation_pre-final.json" in report_text,
          "report_references_exact_independent_validation_receipts")
    check("rs2013701" in report_text and "NPNT" in report_text
          and "not sensitivity" in report_text and "region" in limitations_text.lower()
          and "allelic" in limitations_text.lower(), "report_contains_required_case_and_interpretation_boundaries")
    check(timestamp(prefinal["completed_utc"]) < timestamp(report_receipt["generated_utc"]),
          "finished_report_postdates_pre_final_validation")
    return {"reports": report_receipt["reports"], "report_headline_metrics": headline,
            "final_report_command_receipt": record(STAGE / report_command),
            "pre_final_validation": record(prefinal_path),
            "pre_final_executed_source_preserved": record(historical_script),
            "final_validator_changes": "Added final report/headline/source-lineage and direction-locus/manifest consistency checks only; no scientific method or stored result changed.",
            "no_report_edits_after_final_validation_required": True}


def validate_preservation(phase):
    sources = read_json("provenance/preopen_source_manifest.json")
    check(sources["baseline_commit"] == BASELINE and len(sources["benchmark_records"]) == 218, "preservation_baseline_and_benchmark_record_count")
    file_checks(REPO, sources["benchmark_records"], "original_benchmark_opaque_hash_unchanged")
    file_checks(REPO, sources["dependencies"], "frozen_dependencies_unchanged")
    check(git("rev-parse", "HEAD").decode().strip() == BASELINE, "no_commit_or_push_baseline_HEAD_unchanged")
    registers = {r["path"]: r for r in sources["register_baselines"]}
    changed = set(git("diff", "--name-only", BASELINE).decode().splitlines())
    check(changed <= set(registers), "all_historical_tracked_artifacts_unchanged_except_authorized_registers", detail=json.dumps(sorted(changed)))
    status = git("status", "--porcelain=v1", "-z", "--untracked-files=all").decode().split("\0")
    stage_prefix = str(STAGE.relative_to(REPO)) + "/"
    scoped = []
    for entry in status:
        if not entry:
            continue
        path = entry[3:]
        scoped.append(path)
        check(path in registers or path.startswith(stage_prefix), "working_tree_change_scope", path)
    updates = []
    for path, baseline in registers.items():
        original = git("show", BASELINE + ":" + path)
        current = (REPO / path).read_bytes()
        check(len(original) == baseline["bytes"] and hashlib.sha256(original).hexdigest() == baseline["sha256"], "register_exact_baseline", path)
        check(current.startswith(original), "register_append_only_prefix", path)
        suffix = current[len(original):]
        if phase == "final":
            check(bool(suffix), "final_stage_register_updated", path)
        check(not suffix or b"external-benchmark-v2-1.0" in suffix, "register_append_stage_identifier", path)
        updates.append({"path": path, "before_bytes": len(original), "appended_bytes": len(suffix),
                        "current_sha256": hashlib.sha256(current).hexdigest()})
    if phase == "final":
        receipt = read_json("provenance/shared_register_updates.json")
        check(receipt["status"] == "PASS" and receipt["prior_bytes_changed"] is False,
              "final_shared_register_receipt_status")
        check({r["path"] for r in receipt["registers"]} == set(registers), "all_three_registers_receipted")
        file_checks(REPO, receipt["registers"], "final_shared_register_receipt_hash")
        check(all(r["append_only"] is True for r in receipt["registers"]), "all_register_appends_receipted")
    checkpoint_prefix = "diseases/COPD/07_gap_closure/internal-training-1.0/runs"
    paths = [x for x in git("ls-tree", "-r", "--name-only", BASELINE, "--", checkpoint_prefix).decode().splitlines()
             if x.endswith("/selected_checkpoint.keras")]
    check(len(paths) == 18, "all_18_original_checkpoint_paths_preserved")
    checkpoints = []
    for path in paths:
        pointer = git("show", BASELINE + ":" + path).decode()
        lines = pointer.splitlines()
        check(lines[0] == "version https://git-lfs.github.com/spec/v1", "baseline_checkpoint_LFS_pointer", path)
        oid = next(x.removeprefix("oid sha256:") for x in lines if x.startswith("oid sha256:"))
        size = int(next(x.removeprefix("size ") for x in lines if x.startswith("size ")))
        original = REPO / path
        check(original.stat().st_size == size and digest(original) == oid, "original_checkpoint_archive_SHA256_unchanged", path)
        checkpoints.append({"path": path, "bytes": size, "sha256": oid})
    check(not list(STAGE.rglob("*.keras")), "no_checkpoint_duplicate_in_new_stage")
    return {"original_benchmark_records_opaque_verified": 218, "original_checkpoints_verified": checkpoints,
            "register_updates": updates, "historical_Git_preservation": "PASS",
            "candidate_universe_opened": False, "original_benchmark_parsed_again": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True, choices=("pre-final", "final"))
    parser.add_argument("--label", default="", help="Unique receipt suffix for a separately authorized correction/revalidation")
    args = parser.parse_args()
    stem = "independent_final_validation" if args.phase == "final" else "independent_validation_pre-final"
    if args.label:
        if not args.label.replace("-", "").replace("_", "").isalnum():
            raise ValueError("Unsafe receipt label")
        stem += "_" + args.label
    json_path = STAGE / "provenance" / (stem + ".json")
    tsv_path = json_path.with_suffix(".tsv")
    if json_path.exists() or tsv_path.exists():
        raise RuntimeError("Refuse validation receipt overwrite")
    started = time.perf_counter()
    result = {"stage": STAGE.name, "phase": args.phase, "started_utc": datetime.now(timezone.utc).isoformat(),
              "script": record(Path(__file__).resolve()), "command": [sys.executable] + sys.argv,
              "independent_implementation": True, "scorer_evaluator_or_model_imported": False,
              "scientific_inference_rerun": False}
    exit_code = 0
    try:
        inference, evaluation, before = validate_manifests_and_chronology()
        inputs = read_tsv("inputs/scorable_variants.tsv")
        check(len(inputs) == 1710 and [int(r["row_index"]) for r in inputs] == list(range(1710)), "frozen_sequence_row_order")
        ensemble, _ = validate_predictions(inputs, inference, before)
        with gzip.open(STAGE / "inputs/benchmark_snapshot.json.gz", "rt") as handle:
            snapshot = json.load(handle)
        annotations = validate_annotations(snapshot, inputs)
        result["population_summary"] = validate_evaluation(snapshot, inputs, ensemble, *annotations)
        result["preservation"] = validate_preservation(args.phase)
        if args.phase == "final":
            result["final_reports_and_validation_lineage"] = validate_final_reports()
        result["status"] = "PASS"
    except BaseException as exc:
        result.update(status="FAIL", exception_type=type(exc).__name__, error=str(exc), traceback=traceback.format_exc())
        exit_code = 1
    with tsv_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("category", "key", "status", "detail"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(AUDIT_ROWS)
    result.update(completed_utc=datetime.now(timezone.utc).isoformat(), checks=len(AUDIT_ROWS),
                  failed_checks=sum(row["status"] != "PASS" for row in AUDIT_ROWS), check_table=record(tsv_path),
                  environment={"python": platform.python_version(), "numpy": np.__version__, "CPU_only": True},
                  wall_seconds=time.perf_counter() - started, peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with json_path.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
