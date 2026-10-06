#!/usr/bin/env python3
"""Independent read-only scientific QC; writes only new validation evidence.

No model runtime, training helper, inference helper or evaluator is imported.
All predictions are audited from already-generated frozen test outputs only.
External benchmark/candidate payloads are never opened, even for hashing.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
SEEDS = (104729, 130363, 155921)
CONFIGS = ("V2-A", "V2-B", "V2-C")
MODELS = ("enhancer", "h3k27me3")
EXPECTED = {"enhancer": (25390, 12695, 8754), "h3k27me3": (3508, 1754, 1809)}
THRESHOLDS = {"enhancer": "0x1.7f39710000001p-1", "h3k27me3": "0x1.8a0a12aaaaaacp-1"}
METRICS = ("AP", "AUROC", "Brier", "BrierSkill", "AP_gain", "prevalence")
RATES = ("sensitivity", "specificity", "precision", "NPV", "accuracy", "negative_row_FPR")
PHASE_I_SHA = "483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6"
BASELINE = "07d72dad0f57d605f5f50bbc7fa6d2c7843a94ba"
GAP = "diseases/COPD/07_gap_closure"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def js(path):
    with Path(path).open() as stream:
        return json.load(stream)


def table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def array_hash(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


class Audit:
    def __init__(self, stage, repo):
        self.stage, self.repo = stage, repo
        self.rows, self.files = [], {}

    def check(self, name, passed, detail=""):
        self.rows.append({"check": name, "status": "PASS" if bool(passed) else "FAIL", "detail": str(detail)})
        if not passed:
            raise RuntimeError(name + ": " + str(detail))

    def close(self, name, actual, expected, atol=2e-12):
        self.check(name, np.isfinite(actual) and np.isfinite(expected) and abs(float(actual)-float(expected)) <= atol,
                   f"observed={actual!r}; independent={expected!r}; atol={atol}")

    def record(self, base, item):
        path = Path(item["path"])
        path = path if path.is_absolute() else base / path
        path = path.resolve()
        # The caller selects authorized stage/test/reference records. Fail closed
        # if a manifest tries to redirect a record into benchmark/candidate data.
        roots = [self.stage, self.stage.parent / "internal-training-1.0", self.stage.parent / "pretraining-1.1"]
        refs = [(self.repo / "diseases/COPD/04_modeling/trednet" / r).resolve()
                for r in ("fasta/hg38.fa", "fasta/hg38.fa.fai", "model_phase_I/phase_one_weights.h5")]
        self.check("authorized_hash_path", any(path.is_relative_to(root) for root in roots) or path in refs, path)
        key = str(path)
        if key not in self.files:
            self.files[key] = {"path": key, "bytes": path.stat().st_size, "sha256": sha(path)}
        observed = self.files[key]
        self.check("file_hash_unchanged", observed["sha256"] == item["sha256"], path)
        if "bytes" in item:
            self.check("file_size_unchanged", observed["bytes"] == int(item["bytes"]), path)
        return path

    def prospective(self):
        freeze = js(self.stage / "provenance/prospective_specification_freeze.json")
        self.check("prospective_status_and_order_attestation", freeze["status"] == "PASS" and freeze["before_any_test_model_inference"] is True)
        for key in ("specification", "markdown", "authorization"):
            self.record(self.stage, freeze[key])
        integrity = js(self.stage / "provenance/input_integrity_verified_before.json")
        self.check("independent_initial_source_integrity_PASS", integrity["status"] == "PASS")
        self.record(self.stage, integrity["prospective_specification_freeze"])
        return freeze

    def inputs(self):
        self.prospective()
        manifest = js(self.stage / "inputs/input_manifest.json")
        self.check("input_manifest_complete_PASS", manifest["status"] == "PASS" and manifest["completed"] is True)
        for item in manifest["files"]:
            self.record(self.stage, item)
        for item in manifest["frozen_sources"]:
            self.record(self.repo, item)
        cache_rows = table(self.stage / "inputs/cache_index.tsv.gz")
        by_id = {r["interval_id"]: r for r in cache_rows}
        self.check("unique_union_rows", len(cache_rows) == len(by_id) == 26225)
        self.check("union_components", len({r["component_id"] for r in cache_rows}) == 8926)
        self.check("canonical_cache_one_to_one", {int(r["cache_row"]) for r in cache_rows} == set(range(26225)))
        expected_order = sorted(cache_rows, key=lambda r: (int(r["chrom"][3:]), int(r["core_start"]), r["interval_id"]))
        self.check("frozen_genomic_union_order", cache_rows == expected_order)
        union = set()
        for model in MODELS:
            frozen = self.stage.parent / f"pretraining-1.1/data/evaluation/{model}_common_challenge_panel.tsv.gz"
            source = [r for r in table(frozen) if r["validation_role"] == "test"]
            prepared = table(self.stage / f"inputs/common/{model}_test.tsv.gz")
            configuration = {r["interval_id"]: r for r in table(self.stage.parent / f"pretraining-1.1/data/configurations/V2-C_{model}_interval_manifest.tsv.gz") if r["validation_role"] == "test"}
            n, positive, components = EXPECTED[model]
            self.check(model + ":exact_population", len(source) == len(prepared) == len(configuration) == n)
            self.check(model + ":class_balance", sum(int(r["label"]) for r in prepared) == positive)
            self.check(model + ":exact_components", len({r["component_id"] for r in prepared}) == components)
            for old, new in zip(source, prepared):
                self.check(model + ":all_frozen_cells_order_unchanged", all(new[k] == value for k, value in old.items()) and old == configuration[old["interval_id"]], old["interval_id"])
                self.check(model + ":role_label_eligibility", new["chrom"] in ("chr8", "chr9") and new["partition"] == new["validation_role"] == "test" and new["sequence_available"] == "1" and new["sequence_length"] == "2001" and new["label"] in ("0", "1"), new["interval_id"])
                mapped = by_id[new["interval_id"]]
                self.check(model + ":cache_mapping", all(new[k] == mapped[k] for k in mapped), new["interval_id"])
                union.add(new["interval_id"])
        self.check("test_panel_union_exact", union == set(by_id))
        sequence = np.load(self.stage / "inputs/canonical_sequences.npy", mmap_mode="r", allow_pickle=False)
        self.check("sequence_array_geometry", sequence.dtype == np.uint8 and sequence.shape == (26225, 2001))
        fasta = self.repo / "diseases/COPD/04_modeling/trednet/fasta/hg38.fa"
        fai = {}
        for line in Path(str(fasta)+".fai").read_text().splitlines():
            name, length, offset, bases, width = line.split("\t")[:5]
            fai[name] = tuple(map(int, (length, offset, bases, width)))
        normalize = bytes(c if c in b"ACGT" else 78 for c in range(256))
        complement = bytes.maketrans(b"ACGTN", b"TGCAN")
        with fasta.open("rb") as stream:
            for row in cache_rows:
                chrom = row["chrom"]
                start, end, core_start, core_end = [int(row[k]) for k in ("input_start", "input_end", "core_start", "core_end")]
                self.check("test_only_exact_input_geometry", chrom in ("chr8", "chr9") and end-start == 2001 and core_end-core_start == 1000 and start == core_start-501 and end == core_end+500, row["interval_id"])
                length, offset, bases, width = fai[chrom if chrom in fai else chrom[3:]]
                self.check("reference_bounds", 0 <= start < end <= length)
                first = offset + (start // bases)*width + start % bases
                last = offset + ((end-1)//bases)*width + (end-1) % bases
                stream.seek(first)
                raw = stream.read(last-first+1).replace(b"\n", b"").replace(b"\r", b"").upper()
                normalized = raw.translate(normalize)
                reverse = normalized.translate(complement)[::-1]
                canonical = min(normalized, reverse)
                stored = sequence[int(row["cache_row"])].tobytes()
                self.check("raw_and_encoded_sequence_hashes", len(raw) == 2001 and hashlib.sha256(raw).hexdigest() == row["sequence_sha256"] and stored == canonical and hashlib.sha256(stored).hexdigest() == row["canonical_rc_sequence_sha256"], row["interval_id"])
                self.check("genomic_forward_cache_orientation", int(row["forward_is_canonical"]) == int(normalized == canonical), row["interval_id"])
                self.close("sequence_GC", float(row["gc_fraction"]), (raw.count(b"C")+raw.count(b"G"))/2001, 1e-9)
                self.close("sequence_ambiguity", float(row["non_acgt_fraction"]), normalized.count(b"N")/2001, 1e-9)
        self.check("exact_sequence_count", manifest["sequence_intervals"] == manifest["encoded_identity_count"] == 26225)

    def cache(self):
        self.prospective()
        manifest = js(self.stage / "cache/cache_manifest.json")
        self.check("cache_complete_PASS", manifest["completed"] is True and manifest["status"] == "PASS")
        self.check("exact_phaseI_weight", manifest["phase_I_weights_sha256"] == PHASE_I_SHA)
        self.check("phaseI_geometry_and_orientation_contract", manifest["shape"] == [26225,4560] and manifest["dtype"] == "float32" and manifest["representation"] == "final4560sigmoid" and manifest["nucleotide_rc_independently_evaluated"] is True and manifest["orientation_network_evaluations"] == 52450)
        for item in manifest["files"]:
            self.record(self.stage, item)
        self.record(self.stage, manifest["prospective_freeze"])
        arrays = []
        for name in ("features_canonical.npy", "features_rc.npy"):
            arr = np.load(self.stage / "cache" / name, mmap_mode="r", allow_pickle=False)
            self.check("cache_array_geometry", arr.dtype == np.float32 and arr.shape == (26225,4560), name)
            for start in range(0, len(arr), 512):
                block = arr[start:start+512]
                self.check("cache_finite_sigmoid_values", np.isfinite(block).all() and np.all((block >= 0) & (block <= 1)), f"{name}:{start}")
            arrays.append(arr)
        self.check("orientations_not_identical_or_column_reversal", not np.array_equal(arrays[0], arrays[1]) and not np.array_equal(arrays[0][:,::-1], arrays[1]))
        env = js(self.stage / "cache/environment.json")
        self.check("exact_runtime", env["python"] == "3.13.0" and all(env["packages"][k] == v for k,v in (("tensorflow","2.20.0"),("keras","3.14.1"),("numpy","2.5.0"))))
        self.check("phaseI_architecture_exact", js(self.stage / "cache/architecture.json") == js(self.stage.parent / "internal-training-1.0/cache/architecture.json"))
        started = js(self.stage / "cache/started.json")
        self.check("nucleotide_RC_independent_calls_recorded", started["batch_size"] == 32 and started["two_independent_nucleotide_orientation_passes"] is True and started["representation_axis_reversal"] is False and started["optimizer"] is None)
        self.record(self.stage, started["input_manifest"])
        self.record(self.stage, started["phase_I_gate"])

    def prediction_arrays(self, model, independent=True):
        metadata = table(self.stage / f"inputs/common/{model}_test.tsv.gz")
        rows = table(self.stage / f"predictions/{model}_chr8_9.tsv.gz")
        self.check(model+":exact_prediction_order_and_membership", [(r["interval_id"],r["chrom"]) for r in rows] == [(r["interval_id"],r["chrom"]) for r in metadata])
        repeats = table(self.stage / f"predictions/{model}_independent_rc_pass.tsv.gz") if independent else None
        audit_units = {(r["configuration"],r["unit"],str(r["seed"])):r
                       for r in js(self.stage/"predictions/real_network_invariance.json")["audit_rows"]
                       if r["model"]==model} if independent else {}
        if repeats is not None:
            self.check(model+":repeated_prediction_order", [(r["interval_id"],r["chrom"]) for r in repeats] == [(r["interval_id"],r["chrom"]) for r in metadata])
        def values(data, col):
            v = np.array([float(r[col]) for r in data], dtype=np.float64)
            self.check("finite_range:"+col, np.isfinite(v).all() and np.all((v>=0)&(v<=1)))
            return v
        output = {}
        for config in CONFIGS:
            seeds, reverse_seeds = [], []
            for seed in SEEDS:
                fwd = values(rows, f"{config}_p_forward_{seed}")
                rc = values(rows, f"{config}_p_rc_{seed}")
                q = (fwd+rc)/2
                self.check("exact_symmetric_float64", np.array_equal(q, values(rows, f"{config}_q_{seed}")), f"{model}:{config}:{seed}")
                output[(config,f"seed:{seed}")] = q
                seeds.append(q)
                if repeats is not None:
                    first = values(repeats, f"{config}_rc_input_first_{seed}")
                    second = values(repeats, f"{config}_rc_input_second_{seed}")
                    qrc = (first+second)/2
                    self.check("exact_repeated_symmetric_float64", np.array_equal(qrc, values(repeats, f"{config}_q_rc_input_{seed}")))
                    self.check("seed_numeric_RC_tolerance", np.all(np.abs(q-qrc) <= 1e-6+1e-6*np.abs(qrc)), f"{model}:{config}:{seed}")
                    unit = audit_units[config,"seed",str(seed)]
                    for name,value in (("original_forward",fwd),("original_rc",rc),
                                       ("reverse_pass_first",first),("reverse_pass_second",second),
                                       ("original_wrapper",q),("reverse_wrapper",qrc)):
                        self.check("RC_audit_array_digest",unit[name+"_array_sha256"]==array_hash(value),f"{model}:{config}:{seed}:{name}")
                    self.close("RC_recorded_max_residual",unit["max_absolute_difference"],np.max(np.abs(q-qrc)),0)
                    reverse_seeds.append(qrc)
            ensemble = np.mean(np.stack(seeds), axis=0,dtype=np.float64)
            self.check("equal_fixed_three_seed_float64_ensemble", np.array_equal(ensemble, values(rows, f"{config}_ensemble_q")), f"{model}:{config}")
            output[(config,"ensemble")] = ensemble
            if repeats is not None:
                reverse = np.mean(np.stack(reverse_seeds),axis=0,dtype=np.float64)
                self.check("exact_repeated_ensemble", np.array_equal(reverse,values(repeats,f"{config}_ensemble_q_rc_input")))
                self.check("ensemble_numeric_RC_tolerance", np.all(np.abs(ensemble-reverse)<=1e-6+1e-6*np.abs(reverse)))
                unit = audit_units[config,"ensemble","all3"]
                self.check("ensemble_RC_array_digests",unit["original_wrapper_array_sha256"]==array_hash(ensemble) and unit["reverse_wrapper_array_sha256"]==array_hash(reverse))
                self.close("ensemble_RC_recorded_max_residual",unit["max_absolute_difference"],np.max(np.abs(ensemble-reverse)),0)
        return metadata, output

    def predictions(self):
        self.prospective()
        record = js(self.stage / "predictions/real_network_invariance.json")
        self.check("actual_network_RC_gate_PASS", record["status"] == "PASS" and record["real_network"] is True and record["all_configurations_all_seeds_and_ensembles"] is True)
        for item in record["files"]:
            self.record(self.stage, item)
        self.check("exact_24_RC_units", record["n_seed_audits"] == 18 and record["n_ensemble_audits"] == 6 and len(record["audit_rows"]) == 24)
        units = set()
        for row in record["audit_rows"]:
            key = row["configuration"],row["model"],row["unit"],str(row["seed"])
            self.check("RC_audit_unique_unit", key not in units)
            units.add(key)
            self.check("RC_full_sequence_count_and_tolerance", row["status"] == "PASS" and row["n_sequences"] == EXPECTED[row["model"]][0] and row["n_failed"] == 0 and row["atol"] == row["rtol"] == 1e-6 and row["max_tolerance_normalized_difference"] <= 1)
        expected = {(c,m,"seed",str(s)) for c in CONFIGS for m in MODELS for s in SEEDS} | {(c,m,"ensemble","all3") for c in CONFIGS for m in MODELS}
        self.check("all_prescribed_RC_units", units == expected)
        for model in MODELS:
            self.prediction_arrays(model)
        checkpoint_gate = js(self.stage / "provenance/checkpoint_verification.json")
        for checkpoint in checkpoint_gate["checkpoints"]:
            self.record(self.repo,checkpoint)

    def evaluation(self):
        self.prospective()
        # Independent implementation, plus sklearn cross-check of all point
        # estimates; no imports from the executed evaluator or frozen helper.
        from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss
        out = self.stage / "results/chr8_9_evaluation"
        manifest = js(out / "evaluation_manifest.json")
        for rec in manifest["artifacts"]:
            self.record(self.repo,rec)
        metrics = table(out / "metrics.tsv")
        boot = table(out / "bootstrap_replicates.tsv.gz")
        summaries = table(out / "bootstrap_summary.tsv")
        operations = table(out / "fixed_threshold_operating.tsv")
        operation_replicates = table(out / "fixed_threshold_bootstrap.tsv.gz")
        operation_intervals = table(out / "fixed_threshold_intervals.tsv")
        reliability = table(out / "reliability_bins.tsv")
        strata_output = table(out / "stratified_metrics.tsv")
        paired = table(out / "paired_ablation.tsv")
        stability = table(out / "seed_stability.tsv")
        gates = table(out / "C_holdout_generalization.tsv")
        audit = js(out / "bootstrap_audit.json")
        decisions = js(out / "holdout_decision.json")
        cuts = js(self.stage.parent / "internal-training-1.0/results/chr7_evaluation/stratification_cutpoints.json")
        self.check("complete_metric_replicate_output_counts", len(metrics)==24 and len(boot)==48000 and len(reliability)==240 and len(operation_replicates)==4000)
        independent_decisions = {}
        for model in MODELS:
            metadata,predictors = self.prediction_arrays(model,independent=False)
            y = np.array([int(r["label"]) for r in metadata])
            components = np.array([r["component_id"] for r in metadata])
            names,inverse = np.unique(components,return_inverse=True)
            prepared = {key: IndependentMetrics(y,p) for key,p in predictors.items()}
            points = {key:calc.calculate(np.ones(len(y))) for key,calc in prepared.items()}
            for row in [r for r in metrics if r["model"]==model]:
                key = row["configuration"],row["predictor"]
                p = predictors[key]
                for metric,value in points[key].items():
                    self.close(model+":point:"+metric,float(row[metric]),value)
                self.close("sklearn_AP",points[key]["AP"],average_precision_score(y,p))
                self.close("sklearn_AUROC",points[key]["AUROC"],roc_auc_score(y,p))
                self.close("sklearn_Brier",points[key]["Brier"],brier_score_loss(y,p))
            br = {(int(r["valid_replicate"]),r["configuration"],r["predictor"]):r for r in boot if r["model"]==model}
            op_rows = {int(r["valid_replicate"]):r for r in operation_replicates if r["model"]==model}
            independent_draws = {key:{metric:[] for metric in METRICS} for key in predictors}
            independent_rates = {metric:[] for metric in RATES}
            rng = np.random.default_rng(314159)
            attempt,valid,invalid = 0,0,0
            while valid < 2000 and attempt < 20000:
                attempt += 1
                mult = np.bincount(rng.integers(0,len(names),size=len(names)),minlength=len(names))
                weights = mult[inverse]
                try:
                    current = {key:calc.calculate(weights) for key,calc in prepared.items()}
                except ValueError:
                    invalid += 1
                    continue
                valid += 1
                for key,values in current.items():
                    stored = br[valid,*key]
                    self.check("exact_bootstrap_attempt_number",int(stored["attempt"])==attempt)
                    for metric,value in values.items():
                        self.close("independent_bootstrap:"+metric,float(stored[metric]),value)
                        independent_draws[key][metric].append(value)
                operating = independent_operating(y,predictors[("V2-C","ensemble")],float.fromhex(THRESHOLDS[model]),weights)
                for metric,value in operating.items():
                    if value is None:
                        self.check("undefined_ratio_left_empty",op_rows[valid][metric]=="")
                    else:
                        self.close("independent_threshold_replicate:"+metric,float(op_rows[valid][metric]),value)
                    if metric in RATES and value is not None:
                        independent_rates[metric].append(value)
            self.check("bootstrap_validity_exact",valid==2000 and audit[model]["valid"]==valid and audit[model]["attempted"]==attempt and audit[model]["invalid"]==invalid)
            self.close("invalid_fraction",float(audit[model]["invalid_fraction"]),invalid/attempt)
            for row in [r for r in summaries if r["model"]==model]:
                values = independent_draws[row["configuration"],row["predictor"]][row["metric"]]
                for label,q in (("lower95",2.5),("upper95",97.5),("lower90",5),("upper90",95)):
                    self.close("independent_CI:"+label,float(row[label]),np.percentile(values,q,method="linear"))
            for row in [r for r in paired if r["model"]==model]:
                lhs,rhs = row["comparison"].split(" minus ")
                metric = row["metric"]
                values = np.array(independent_draws[lhs,"ensemble"][metric])-np.array(independent_draws[rhs,"ensemble"][metric])
                self.close("paired_point",float(row["point"]),points[lhs,"ensemble"][metric]-points[rhs,"ensemble"][metric])
                for label,q in (("lower95",2.5),("upper95",97.5),("lower90",5),("upper90",95)):
                    self.close("paired_CI",float(row[label]),np.percentile(values,q,method="linear"))
            for row in [r for r in stability if r["model"]==model]:
                ap = [points[row["configuration"],f"seed:{s}"]["AP"] for s in SEEDS]
                self.close("seed_sample_SD",float(row["AP_sample_SD"]),np.std(ap,ddof=1))
                self.close("seed_range",float(row["AP_range"]),max(ap)-min(ap))
            expected_operation = independent_operating(y,predictors["V2-C","ensemble"],float.fromhex(THRESHOLDS[model]),np.ones(len(y)))
            row = next(r for r in operations if r["model"]==model)
            self.check("unchanged_threshold_exact_encoding",row["threshold_hex"]==THRESHOLDS[model] and float(row["threshold_decimal_17g"])==float.fromhex(THRESHOLDS[model]) and row["call_rule"]=="score>=threshold")
            for metric,value in expected_operation.items():
                self.close("threshold_operating_point",float(row[metric]),value)
            for row in [r for r in operation_intervals if r["model"]==model]:
                values = independent_rates[row["metric"]]
                self.check("rate_valid_undefined_counts",int(row["valid_ratio_replicates"])==len(values) and int(row["undefined_ratio_replicates"])==2000-len(values))
                for label,q in (("lower95",2.5),("upper95",97.5),("lower90",5),("upper90",95)):
                    self.close("rate_CI",float(row[label]),np.percentile(values,q,method="linear"))
            for row in [r for r in reliability if r["model"]==model]:
                p = predictors[row["configuration"],row["predictor"]]
                i = int(row["bin"])
                mask = (p>=i/10)&((p<(i+1)/10) if i<9 else (p<=1))
                self.check("fixed_reliability_bin_count",int(row["n"])==int(mask.sum()))
                if mask.any():
                    self.close("reliability_mean_prediction",float(row["mean_probability"]),np.mean(p[mask]))
                    self.close("reliability_observed_fraction",float(row["observed_positive_fraction"]),np.mean(y[mask]))
                else:
                    self.check("empty_reliability_bins_retained",row["mean_probability"]==row["observed_positive_fraction"]=="")
            covariates = ("gc_fraction","atac_signal_percentile_max_train_only","repeat_fraction_2001")
            definitions = {v:np.array([f"quintile_{i+1}" for i in np.searchsorted(cuts[model]["cutpoints"][v],[float(r[v]) for r in metadata],side="right")]) for v in covariates}
            definitions.update({"atac_lobe_signature":np.array([r["atac_lobes"] or "none" for r in metadata]),"same_lobe_mark_support":np.array([r[model+"_same_lobe_support_lobes"] or "none" for r in metadata]),"chromosome":np.array([r["chrom"] for r in metadata]),"ambiguous_base_status":np.array(["present" if float(r["non_acgt_fraction"])>0 else "absent" for r in metadata])})
            for row in [r for r in strata_output if r["model"]==model]:
                self.check("only_frozen_stratum_definitions",row["stratifier"] in definitions)
                mask = definitions[row["stratifier"]]==row["stratum"]
                self.check("stratum_exact_population",int(row["n"])==int(mask.sum()) and int(row["n_positive"])==int(y[mask].sum()) and int(row["n_components"])==len(set(components[mask])))
                self.check("small_stratum_flag",row["cells_lt100_descriptive_only"]==str(mask.sum()<100))
                if mask.any() and len(set(y[mask]))==2:
                    values = IndependentMetrics(y[mask],predictors[row["configuration"],row["predictor"]][mask]).calculate(np.ones(mask.sum()))
                    for metric,value in values.items():
                        self.close("stratum_metric",float(row[metric]),value)
                else:
                    self.check("undefined_one_class_metrics_not_imputed",all(row[m]=="" for m in METRICS if m!="prevalence"))
            c_values = independent_draws["V2-C","ensemble"]
            aps = [points["V2-C",f"seed:{s}"]["AP"] for s in SEEDS]
            conditions = {"all_three_C_seeds_evaluable":len(aps)==3,"real_network_invariance":True,"finite_valid_probabilities":True,"exact_test_positive":int(y.sum())==EXPECTED[model][1],"exact_test_control":int((1-y).sum())==EXPECTED[model][1],"exact_test_components":len(names)==EXPECTED[model][2],"bootstrap_valid_replicates":valid==2000,"bootstrap_invalid_fraction":invalid/attempt<.1,"lower95_AUROC":np.percentile(c_values["AUROC"],2.5)>0.5,"lower95_AP_gain":np.percentile(c_values["AP_gain"],2.5)>0,"lower95_BrierSkill":np.percentile(c_values["BrierSkill"],2.5)>0,"seed_AP_sample_SD":np.std(aps,ddof=1)<=.03,"seed_AP_range":max(aps)-min(aps)<=.1}
            for row in [r for r in gates if r["model"]==model]:
                self.check("independent_release_gate_decision",row["status"]==("PASS" if conditions[row["gate"]] else "FAIL"),row["gate"])
            independent_decisions[model] = "PASS" if all(conditions.values()) else "FAIL"
            self.check("context_release_decision",decisions["model_contexts"][model]["status"]==independent_decisions[model])
        self.check("both_context_release_decision",decisions["both_C_contexts_pass"]==all(x=="PASS" for x in independent_decisions.values()))
        self.check("no_downstream_authorization_inferred",decisions["external_access_authorized"] is False and decisions["candidate_scoring_authorized"] is False and decisions["thresholds_changed"] is False)

    def preservation(self):
        self.prospective()
        source = js(self.stage / "provenance/input_integrity_verified_before.json")
        for rec in source["files"]:
            if not rec["path"].startswith(GAP+"/internal-test-1.0/"):
                self.record(self.repo,rec)
        historical = js(self.stage / "provenance/historical_preservation_baseline.json")
        head = subprocess.check_output(["git","-C",str(self.repo),"rev-parse","HEAD"],text=True).strip()
        self.check("published_git_commit_unchanged",head==historical["baseline_commit"]==BASELINE)
        tree = []
        raw = subprocess.check_output(["git","-C",str(self.repo),"ls-tree","-r","-l","-z","HEAD"])
        for entry in raw.split(b"\0"):
            if entry:
                meta,path = entry.split(b"\t",1)
                mode,kind,oid,size = meta.decode().split()
                tree.append({"path":path.decode(),"mode":mode,"type":kind,"git_object":oid,"git_blob_bytes":size})
        self.check("historical_git_tree_identity_exact",tree==table(self.stage/"provenance/historical_git_tree_before.tsv"))
        registers = {r["path"] for r in historical["registers"]}
        for rec in historical["registers"]:
            data = (self.repo/rec["path"]).read_bytes()
            prefix = data[:rec["bytes"]]
            self.check("register_prior_bytes_unchanged",len(data)>=rec["bytes"] and hashlib.sha256(prefix).hexdigest()==rec["sha256"],rec["path"])
        status = subprocess.check_output(["git","-C",str(self.repo),"status","--porcelain=v1","-z","--untracked-files=all","--ignore-submodules=all"])
        for entry in status.split(b"\0"):
            if entry:
                path = entry.decode()[3:]
                self.check("only_new_stage_or_append_only_register_changes",path.startswith(GAP+"/internal-test-1.0/") or path in registers,path)


class IndependentMetrics:
    """Separate weighted empirical ROC/PR implementation, no evaluator reuse."""
    def __init__(self,y,p):
        self.y,self.p = np.asarray(y,dtype=np.int64),np.asarray(p,dtype=np.float64)
        self.order = np.argsort(self.p,kind="stable")
        scores = self.p[self.order]
        self.starts = np.r_[0,np.flatnonzero(scores[1:]!=scores[:-1])+1]
        self.error = (self.y-self.p)**2

    def calculate(self,weights):
        weights = np.asarray(weights,dtype=np.float64)
        positives = float(np.sum(weights*self.y)); negatives = float(np.sum(weights*(1-self.y)))
        if positives<=0 or negatives<=0:
            raise ValueError("One-class component draw")
        w,y = weights[self.order],self.y[self.order]
        positive_groups = np.add.reduceat(w*y,self.starts)
        negative_groups = np.add.reduceat(w*(1-y),self.starts)
        less_negative = np.cumsum(negative_groups)-negative_groups
        auroc = float(np.sum(positive_groups*(less_negative+.5*negative_groups))/(positives*negatives))
        reverse_tp = np.cumsum(positive_groups[::-1]); reverse_n = np.cumsum((positive_groups+negative_groups)[::-1])
        precision = np.divide(reverse_tp,reverse_n,out=np.zeros_like(reverse_tp),where=reverse_n>0)
        ap = float(np.sum(positive_groups[::-1]*precision)/positives)
        prevalence = positives/(positives+negatives)
        brier = float(np.sum(weights*self.error)/(positives+negatives))
        return {"AP":ap,"AUROC":auroc,"Brier":brier,"BrierSkill":1-brier/(prevalence*(1-prevalence)),"AP_gain":ap-prevalence,"prevalence":prevalence}


def independent_operating(y,p,threshold,w):
    called = p>=threshold
    tp,fp,tn,fn = [float(np.sum(w[mask])) for mask in ((y==1)&called,(y==0)&called,(y==0)&~called,(y==1)&~called)]
    ratio = lambda a,b: a/b if b else None
    return {"TP":tp,"FP":fp,"TN":tn,"FN":fn,"sensitivity":ratio(tp,tp+fn),"specificity":ratio(tn,tn+fp),"precision":ratio(tp,tp+fp),"NPV":ratio(tn,tn+fn),"accuracy":ratio(tp+tn,tp+fp+tn+fn),"negative_row_FPR":ratio(fp,fp+tn)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--stage",type=Path,required=True)
    parser.add_argument("--mode",choices=("inputs","cache","predictions","evaluation","preservation"),required=True)
    args = parser.parse_args()
    repo,stage = args.repo.resolve(),args.stage.resolve()
    if stage != repo/GAP/"internal-test-1.0":
        raise RuntimeError("Output stage must be the authorized new internal-test-1.0")
    prefix = stage/"provenance"/(args.mode+"_independent_validation")
    if prefix.with_suffix(".json").exists() or prefix.with_suffix(".tsv").exists():
        raise RuntimeError("Independent evidence already exists; refusing replacement")
    audit = Audit(stage,repo)
    started,started_utc = time.monotonic(),utc()
    error = None
    try:
        getattr(audit,args.mode)()
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
        if not audit.rows or audit.rows[-1]["status"] != "FAIL":
            audit.rows.append({"check":"execution_exception","status":"FAIL","detail":error})
    with prefix.with_suffix(".tsv").open("x",newline="") as stream:
        writer = csv.DictWriter(stream,["check","status","detail"],delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(audit.rows)
    result = {"status":"FAIL" if error else "PASS","mode":args.mode,"started_utc":started_utc,"completed_utc":utc(),"wall_seconds":time.monotonic()-started,"checks":len(audit.rows),"failures":sum(r["status"]!="PASS" for r in audit.rows),"error":error,"script":{"path":str(Path(__file__).relative_to(stage)),"sha256":sha(__file__)},"checked_files":list(audit.files.values()),"model_execution_performed":False,"external_benchmark_opened":False,"candidate_payloads_opened":False,"historical_sources_written":False}
    with prefix.with_suffix(".json").open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
    print(json.dumps({k:result[k] for k in ("status","mode","checks","failures","error","wall_seconds")}))
    return 1 if error else 0


if __name__ == "__main__":
    raise SystemExit(main())
