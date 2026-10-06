#!/usr/bin/env python3
"""Render saved, independently validated chr8-9 summaries; never evaluate.

No model, checkpoint, sequence array, row-level score table, historical stage,
external benchmark, or candidate input is opened. Only this new stage's frozen
specification, summary tables, manifests, and execution/audit records are read.
The report is exclusively created and cannot overwrite an earlier report.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

MODELS = ("enhancer", "h3k27me3")
CONFIGURATIONS = ("V2-A", "V2-B", "V2-C")
SEEDS = (104729, 130363, 155921)
METRICS = ("AP", "AUROC", "Brier", "BrierSkill")
QC_PHASES = ("inputs", "cache", "predictions", "evaluation", "preservation")


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def read_table(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stage_path(stage, relative):
    path = (stage / relative).resolve()
    if not path.is_relative_to(stage):
        raise RuntimeError("Report input is outside the new test stage")
    return path


def checked(stage, record):
    path = stage_path(stage, record["path"])
    if ("bytes" in record and path.stat().st_size != int(record["bytes"])) or sha256(path) != record["sha256"]:
        raise RuntimeError("Report summary hash/size mismatch: " + str(path))
    return path


def escape(value):
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value).replace("|", "\\|").replace("\n", " ")


def number(value, digits=7):
    if value is None or value == "":
        return "—"
    try:
        return format(float(value), f".{digits}g")
    except (ValueError, TypeError):
        return escape(value)


def table(headers, rows):
    rows = list(rows)
    if not rows:
        return "No records.\n"
    return "\n".join(["| " + " | ".join(map(escape, headers)) + " |",
                       "| " + " | ".join("---" for _ in headers) + " |",
                       *("| " + " | ".join(map(escape, row)) + " |" for row in rows)]) + "\n"


def confidence(row):
    return f"{number(row['point'])} [{number(row.get('lower95'))}, {number(row.get('upper95'))}]"


def context(model):
    return "Enhancer-associated" if model == "enhancer" else "H3K27me3-associated"


def disposition(decision, qc):
    if not decision["both_C_contexts_pass"]:
        return ("STOP — at least one C context is FAIL or INCONCLUSIVE under the prospective internal-holdout rule. "
                "Investigator review of the recorded conditions is required; no redesign, retraining, external benchmark access, or candidate scoring is authorized.")
    if any(item.get("status") != "PASS" for item in qc.values()):
        return ("Both C scientific holdout gates passed, but independent validation is non-PASS. "
                "The stage is not ready for downstream release; stop for investigator review of the validation findings.")
    return ("READY FOR INVESTIGATOR REVIEW BEFORE ANY EXTERNAL FUNCTIONAL EVALUATION. "
            "Both C contexts pass the prospective holdout rule and all five recorded independent validation phases pass. "
            "This is a review disposition only: separate authorization is still required for external evaluation or candidate/allele scoring.")


def paired_description(rows, model):
    parts = []
    for comparison in ("V2-B minus V2-A", "V2-C minus V2-B", "V2-C minus V2-A"):
        found = [row for row in rows if row["model"] == model and row["comparison"] == comparison]
        summaries = []
        for metric in ("AP", "AUROC", "Brier"):
            matches = [row for row in found if row["metric"] == metric]
            if not matches:
                summaries.append(metric + ": interval unavailable")
                continue
            row = matches[0]
            if row.get("lower95", "") == "" or row.get("upper95", "") == "":
                direction = "interval unavailable"
            elif float(row["lower95"]) > 0:
                direction = "entire interval above zero"
            elif float(row["upper95"]) < 0:
                direction = "entire interval below zero"
            else:
                direction = "interval includes zero"
            summaries.append(f"{metric} {number(row['point'])} ({direction})")
        parts.append(comparison + ": " + "; ".join(summaries) + ".")
    return " ".join(parts)


def report(stage, evaluation_relative="results/chr8_9_evaluation", predictions_relative="predictions"):
    stage = Path(stage).resolve()
    if stage.name != "internal-test-1.0":
        raise RuntimeError("Only the new internal-test-1.0 stage is permitted")
    evaluation = stage_path(stage, evaluation_relative)
    predictions = stage_path(stage, predictions_relative)
    freeze = read_json(stage / "provenance/prospective_specification_freeze.json")
    if freeze.get("status") != "PASS" or freeze.get("before_any_test_model_inference") is not True:
        raise RuntimeError("Prospective test specification was not frozen before inference")
    specification = read_json(checked(stage, freeze["specification"]))
    checked(stage, freeze["markdown"])
    if specification.get("version") != "internal-test-1.0":
        raise RuntimeError("Unexpected test-stage specification")
    checkpoints = read_json(stage / "provenance/checkpoint_verification.json")
    expected_checkpoints = {(c, m, s) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}
    checkpoint_rows = checkpoints.get("checkpoints", [])
    observed_checkpoints = {(r["configuration"], r["model"], int(r["seed"])) for r in checkpoint_rows}
    if checkpoints.get("status") != "PASS" or len(checkpoint_rows) != 18 or observed_checkpoints != expected_checkpoints:
        raise RuntimeError("Complete all-18 frozen-checkpoint verification is required")
    manifest = read_json(evaluation / "evaluation_manifest.json")
    artifact_map = {stage_path(stage, record["path"]): record for record in manifest["artifacts"]}

    def evaluation_file(name, as_json=False):
        path = evaluation / name
        if path not in artifact_map:
            raise RuntimeError("Summary is not bound by the evaluation manifest: " + name)
        verified = checked(stage, artifact_map[path])
        return read_json(verified) if as_json else read_table(verified)

    decision = evaluation_file("holdout_decision.json", True)
    if set(decision["model_contexts"]) != set(MODELS):
        raise RuntimeError("Both saved C-context decisions are required")
    for model in MODELS:
        if decision["model_contexts"][model]["status"] not in ("PASS", "FAIL", "INCONCLUSIVE"):
            raise RuntimeError("Unrecognized C-context decision")
    both_pass = all(decision["model_contexts"][model]["status"] == "PASS" for model in MODELS)
    if decision["both_C_contexts_pass"] is not both_pass:
        raise RuntimeError("Contradictory saved C-context decisions")
    if decision.get("external_access_authorized") is not False or decision.get("candidate_scoring_authorized") is not False:
        raise RuntimeError("Report cannot broaden the test-stage authorization")
    qc = {phase: read_json(stage / f"provenance/{phase}_independent_validation.json") for phase in QC_PHASES}
    inputs = read_json(stage / "inputs/input_manifest.json")
    input_audit = read_json(stage / "inputs/test_panel_integrity_audit.json")
    cache = read_json(stage / "cache/cache_manifest.json")
    invariance = read_json(predictions / "real_network_invariance.json")
    metrics = evaluation_file("metrics.tsv")
    summaries = evaluation_file("bootstrap_summary.tsv")
    paired = evaluation_file("paired_ablation.tsv")
    stability = evaluation_file("seed_stability.tsv")
    gates = evaluation_file("C_holdout_generalization.tsv")
    bootstrap = evaluation_file("bootstrap_audit.json", True)
    operating = evaluation_file("fixed_threshold_operating.tsv")
    operating_intervals = evaluation_file("fixed_threshold_intervals.tsv")
    reliability = evaluation_file("reliability_bins.tsv")
    stratified = evaluation_file("stratified_metrics.tsv")
    expected_predictors = {(m, c, p) for m in MODELS for c in CONFIGURATIONS
                           for p in ("ensemble", *(f"seed:{s}" for s in SEEDS))}
    if (len(metrics) != 24 or {(r["model"], r["configuration"], r["predictor"]) for r in metrics} != expected_predictors
            or len(stability) != 6 or len(operating) != 2 or len(reliability) != 240):
        raise RuntimeError("Incomplete saved test summary population")
    intervals = {(r["model"], r["configuration"], r["predictor"], r["metric"]): r for r in summaries}

    def metric_cells(row):
        result = []
        for metric in METRICS:
            entry = intervals.get((row["model"], row["configuration"], row["predictor"], metric))
            result.append(confidence(entry) if entry else f"{number(row[metric])} [interval unavailable]")
        return result

    output = ["# COPD V2 frozen chr8–9 internal holdout evaluation", "",
              "Version: `internal-test-1.0`. Models: unchanged published `internal-training-1.0`; test construction: unchanged `pretraining-1.1`.", "",
              "## Decision and stop boundary", "", disposition(decision, qc), "",
              table(["C context", "Prospective holdout decision", "Failed or inconclusive conditions"],
                    ((context(m), decision["model_contexts"][m]["status"],
                      "; ".join(f"{r['gate']}: {r['status']} (observed {r.get('observed')}; {r.get('criterion')})"
                                for r in decision["model_contexts"][m]["failed_or_inconclusive_gates"]) or "none") for m in MODELS)),
              "Chr8–9 was held out from V2 training, V2 checkpoint selection and chr7 calibration. It is **not historically pristine**: "
              "V1 used these chromosomes and V2 construction/QC inspected outcome-blind covariates. "
              "This is a **fixed internal V2 region-label holdout, not independent external validation**. "
              "No result changes the model configuration, matching, seed list, checkpoint, threshold or prospective contract.", "",
              "The external functional benchmark remains unopened and unparsed. No COPD variant universe, 337 V1 candidates, "
              "12-shortlist variants, rs2013701 or other special functional variants, or arbitrary REF/ALT sequences were scored. "
              "No retraining, fine-mapping, target-gene analysis, manuscript revision or other V2 module was performed. "
              "No Git commit or push is authorized by this stage request.", "",
              "## Prospective contract and provenance", "",
              table(["Record", "Value"], [("Prospective freeze UTC", freeze["created_utc"]),
                    ("Before any test inference", freeze["before_any_test_model_inference"]),
                    ("Baseline published commit", specification["baseline_commit"]),
                    ("Machine specification SHA-256", freeze["specification"]["sha256"]),
                    ("Readable specification SHA-256", freeze["markdown"]["sha256"]),
                    ("Prior training freeze SHA-256", next(r["sha256"] for r in freeze["sources"] if r["path"].endswith("internal-training-1.0/provenance/freeze.json")))]),
              "The [prospective specification](TEST_SPECIFICATION.md), [machine contract](specification/test_specification.json), "
              "[pre-inference freeze](provenance/prospective_specification_freeze.json), and [preserved authorization](provenance/user_authorization.txt) "
              "define this stage. Hash-bound gates precede input preparation, phase-I extraction, checkpoint inference and one-time evaluation. "
              "Model-performance outcomes were not used to revise the contract.", "",
              "## Frozen population, membership and sequence eligibility", "",
              table(["Context", "Positive", "Control", "Rows", "Global components", "chr8", "chr9", "Input audit"],
                    ((context(m), specification["test_population"][m]["n_positive"], specification["test_population"][m]["n_control"],
                      specification["test_population"][m]["n_rows"], specification["test_population"][m]["n_components"],
                      specification["test_population"][m]["chromosome_counts"]["chr8"], specification["test_population"][m]["chromosome_counts"]["chr9"],
                      next(a["status"] for a in input_audit["panels"] if a["model"] == m)) for m in MODELS)),
              f"The union contains {inputs['sequence_intervals']:,} unique genomic intervals and {inputs['encoded_identity_count']:,} canonical encoded identities "
              f"in {input_audit['union_components']:,} existing global overlap/encoded-sequence components. "
              f"The two contexts share {input_audit['cross_model_shared_intervals']:,} interval identities. "
              "All selected source rows, cell strings, order, labels, roles, components and train-derived covariates are preserved; no matching was rebuilt and no row was substituted.", "",
              table(["Frozen consumed input", "Bytes", "SHA-256"], ((r["path"], r["bytes"], r["sha256"]) for r in inputs["frozen_sources"])),
              "[Input integrity details](inputs/test_panel_integrity_audit.json) and the [input manifest](inputs/input_manifest.json) record "
              "1,000-bp cores, `[core_start−501, core_end+500)` 2,001-bp inputs, exact hg38/FAI hashes, raw sequence hashes, "
              "canonical nucleotide/RC hashes, GC and ambiguous-base checks. FASTA access is limited to the exact frozen test intervals on chr8 and chr9.", "",
              "## Frozen phase-I extraction and unchanged checkpoint inference", "",
              table(["Record", "Observed"], [("Sequence geometry / encoding", f"{inputs['sequence_shape']} / {inputs['sequence_dtype']}"),
                    ("Per-orientation feature geometry / dtype", f"{cache['shape']} / {cache['dtype']}"),
                    ("Representation", cache["representation"]),
                    ("Independent nucleotide-orientation phase-I evaluations", cache["orientation_network_evaluations"]),
                    ("Missing sequences / extraction failures", f"{cache['missing_sequences']} / {cache['failures']}"),
                    ("Phase-I weight SHA-256", cache["phase_I_weights_sha256"])]),
              "Non-ACGT bases are normalized to N and all-zero one-hot encoding. Each canonical nucleotide sequence and its nucleotide reverse complement "
              "are independently evaluated through the original frozen phase-I network; representation columns are never reversed. "
              "The orientation map restores genomic-forward and nucleotide-RC inputs for each interval.", "",
              table(["Test cache", "Bytes", "SHA-256"], ((r["path"], r["bytes"], r["sha256"]) for r in cache["files"] if Path(r["path"]).suffix == ".npy")),
              "All 18 original selected `.keras` archives were loaded unchanged with `compile=False` and verified before and after inference. "
              "No retraining, checkpoint replacement, re-saving, optimizer-state editing or favorable-seed selection occurred. "
              "Each seed uses `q_s(x) = [float64(p_s(x)) + float64(p_s(RC(x)))] / 2`; the final score is the equal-three-seed float64 mean "
              "in order 104729, 130363, 155921. Actual network operations remain float32.", "",
              table(["Configuration", "Context", "Seed", "Original checkpoint SHA-256"],
                    ((r["configuration"], r["model"], r["seed"], r["sha256"]) for r in checkpoint_rows)),
              "## Actual-network reverse-complement invariance", "",
              f"Saved invariance status: **{invariance['status']}**; {invariance['n_seed_audits']} seed and {invariance['n_ensemble_audits']} ensemble audits. "
              "Every frozen test-panel sequence is covered. Additional real-network calls with swapped nucleotide-orientation inputs use the same row order and batching; "
              "the second wrapper is not obtained by algebraically recycling saved probabilities. "
              "The frozen tolerance is `abs(Q(x)−Q(RC(x))) <= 1e−6 + 1e−6*abs(Q(RC(x)))`. "
              "All raw, repeated, seed-symmetric and ensemble predictions must remain finite within [0,1].", "",
              table(["Configuration", "Context", "Unit", "Seed", "Sequences", "Status", "Failed", "Maximum absolute residual"],
                    ((r["configuration"], r["model"], r["unit"], r["seed"], r["n_sequences"], r["status"], r["n_failed"],
                      number(r["max_absolute_difference"], 17)) for r in invariance["audit_rows"])),
              "## Primary common-panel A/B/C performance", "",
              "All A/B/C models use the identical frozen C-task test panel within each context. AP is stepwise tie-grouped average precision, "
              "not trapezoidal PR-AUC; AUROC is tie-aware. Brier is mean squared probability error and Brier skill is "
              "`1 − Brier/[prevalence × (1−prevalence)]`. Display rounding does not replace full-precision TSV values. "
              "Entries below are point estimates [95% component-bootstrap confidence intervals].", "",
              table(["Context", "Configuration", "AP [95% CI]", "AUROC [95% CI]", "Brier [95% CI]", "Brier skill [95% CI]"],
                    ((r["model"], r["configuration"], *metric_cells(r)) for r in metrics if r["predictor"] == "ensemble")),
              "### All frozen seed-specific symmetric predictors", "",
              table(["Context", "Configuration", "Predictor", "AP [95% CI]", "AUROC [95% CI]", "Brier [95% CI]", "Brier skill [95% CI]"],
                    ((r["model"], r["configuration"], r["predictor"], *metric_cells(r)) for r in metrics if r["predictor"] != "ensemble")),
              "### Component-bootstrap validity and scope", "",
              table(["Context", "Components", "Attempted", "Valid", "Invalid", "Invalid fraction", "Status"],
                    ((m, bootstrap[m]["represented_components"], bootstrap[m]["attempted"], bootstrap[m]["valid"],
                      bootstrap[m]["invalid"], number(bootstrap[m]["invalid_fraction"]), bootstrap[m]["status"]) for m in MODELS)),
              "The bootstrap uses the existing global full-input overlap/encoded-sequence component IDs, not matching pairs. "
              "Within each context all seed and ensemble A/B/C predictors and fixed-threshold C operating rates share the same draws: lexical component ordering, "
              "fresh NumPy default_rng(314159), 2,000 valid replicates, at most 20,000 attempts, and invalid/attempted strictly below 0.10. "
              "CIs use linear 2.5/97.5 percentiles. AP-minus-prevalence and Brier skill use each replicate's own weighted prevalence. "
              "Uncertainty is conditional on this internal genomic-component/fixed-matching design, not independent-donor or population uncertainty.", "",
              "### Paired diagnostic ablations", "",
              table(["Context", "Paired difference", "Metric", "Point [95% paired CI]"],
                    ((r["model"], r["comparison"], r["metric"], confidence(r)) for r in paired)),
              *[f"{context(m)}: {paired_description(paired, m)}\n" for m in MODELS],
              "For AP/AUROC, positive differences mean higher discrimination; for Brier, negative differences mean lower squared error. "
              "These comparisons are diagnostic ablations only. B→C changes both positive-label definition and independently rematched controls, "
              "so it is not a label-only contrast. C is not required to outperform B; no ordering changes the already-frozen model choice.", "",
              "## C seed stability and prospective holdout-generalization safeguards", "",
              table(["Context", "Configuration", "AP seed 104729", "AP seed 130363", "AP seed 155921", "AP sample SD", "AP range"],
                    ((r["model"], r["configuration"], *(number(r[f"AP_seed_{s}"]) for s in SEEDS),
                      number(r["AP_sample_SD"]), number(r["AP_range"])) for r in stability)),
              "For C, sample SD (`ddof=1`) ≤0.03 and range ≤0.10 are conjunctive release safeguards. "
              "These are internal generalization safeguards, not clinically meaningful effect-size thresholds. A/B stability is descriptive.", "",
              table(["Context", "Prospective condition", "Status", "Observed", "Frozen criterion"],
                    ((r["model"], r["gate"], r["status"], r["observed"], r["criterion"]) for r in gates)),
              "Any FAIL or INCONCLUSIVE context requires stopping for investigator review without redesign, retraining or external access. "
              "Even two PASS decisions do not themselves authorize the next module.", "",
              "## Frozen chr7 region-score thresholds on the test panels", "",
              "Apply the original C ensemble thresholds using `score >= threshold`; no chr8–9 score is used for calibration or adjustment.", "",
              table(["Context", "Original decimal threshold", "Exact float64 hex", "TP", "FP", "TN", "FN"],
                    ((r["model"], r["threshold_decimal_17g"], r["threshold_hex"], r["TP"], r["FP"], r["TN"], r["FN"]) for r in operating)),
              table(["Context", "Operating metric", "Point [95% component CI]", "Defined draws", "Undefined draws"],
                    ((r["model"], r["metric"], confidence(r), r["valid_ratio_replicates"], r["undefined_ratio_replicates"]) for r in operating_intervals)),
              table(["Context", "Observed chr7 negative-row FPR", "Observed test negative-row FPR", "Test minus chr7", "Test FPR >5%", "Recalibrated"],
                    ((r["model"], number(r["observed_chr7_calibration_FPR"]), number(r["negative_row_FPR"]),
                      number(r["test_minus_chr7_calibration_FPR"]), r["test_FPR_above_0_05"], r["threshold_recalibrated"]) for r in operating)),
              "Any increase is reported as calibration drift. A test FPR above 5% is not automatically a protocol failure: "
              "the chr7 empirical 5% rule never guaranteed a test/population FPR. No recalibration is performed. "
              "These are region-label thresholds only—not variant-level FPR control, allele-effect thresholds or causal probabilities. "
              "Precision/NPV intervals omit undefined denominators with explicit counts; those omissions do not invalidate otherwise valid metric-bootstrap draws. "
              "Confusion counts are descriptive; rate CIs retain component dependence rather than assuming independent rows.", "",
              "## Fixed reliability bins", "",
              "All 24 predictors retain ten fixed bins `[0,.1), …, [.9,1]`, including empty bins. "
              "The ensemble bins are shown below; the full per-seed and ensemble table is linked in the artifact index. "
              "Observed fractions on this balanced, matched region panel are not population disease-risk calibration.", "",
              table(["Context", "Configuration", "Bin", "Lower", "Upper", "Upper inclusive", "Rows", "Mean score", "Observed positive fraction"],
                    ((r["model"], r["configuration"], r["bin"], r["lower"], r["upper"], r["upper_inclusive"], r["n"],
                      number(r["mean_probability"]), number(r["observed_positive_fraction"])) for r in reliability if r["predictor"] == "ensemble")),
              "## Frozen stratified diagnostics", "",
              "GC, train-derived ATAC rank and repeat quintiles reuse the exact frozen training-positive cutpoints with `searchsorted(side='right')`; "
              "no test-derived cutpoints or new subgroup definitions are introduced. Exact ATAC-lobe signatures, same-lobe histone-support signatures, "
              "chromosome and ambiguous-base status retain their frozen definitions. The table shows C ensembles; all seeds and A/B/C ensembles are retained in the full artifact. "
              "Every subgroup is descriptive. Cells below 100 rows are explicitly flagged, and one-class/undefined metrics remain empty (shown as —).", "",
              table(["Context", "Stratifier", "Stratum", "Rows", "Positive/control", "Components", "n<100", "AP", "AUROC", "Brier", "Brier skill"],
                    ((r["model"], r["stratifier"], r["stratum"], r["n"], f"{r['n_positive']}/{r['n_control']}", r["n_components"],
                      r["cells_lt100_descriptive_only"], *(number(r[k]) for k in METRICS)) for r in stratified
                     if r["configuration"] == "V2-C" and r["predictor"] == "ensemble")),
              "## Independent validation and historical preservation", "",
              table(["Validation phase", "Status", "Checks", "Failures", "Evidence"],
                    ((phase, qc[phase]["status"], qc[phase].get("checks", qc[phase].get("check_count", "recorded in evidence")),
                      qc[phase].get("failures", "recorded in evidence"), f"provenance/{phase}_independent_validation.json") for phase in QC_PHASES)),
              ("Independent preservation validation is PASS: historical tracked Git object identities and worktree scope are unchanged. "
               "Authorized pretraining-1.1 and internal-training-1.0 payloads and reference inputs are also directly rehashed. "
               "Benchmark payloads were not opened, parsed or newly hashed; their preservation evidence is Git object identity and worktree status only."
               if qc["preservation"].get("status") == "PASS" else
               "Historical preservation is not independently PASS; see the preservation findings. No clean-preservation claim is made."), "",
              "## Environment, compute and execution ledger", ""]
    for label, path in (("Phase-I extraction", stage / "cache/environment.json"), ("Checkpoint inference", predictions / "environment.json")):
        environment = read_json(path)
        output.extend([f"### {label}", "",
                       table(["Runtime field", "Recorded value"], [("Python", environment.get("python")), ("Executable", environment.get("executable")),
                             *( (name, version) for name, version in environment.get("packages", {}).items()),
                             ("GPU/driver record", environment.get("nvidia_smi", {}).get("stdout")),
                             ("Logical devices", "; ".join(environment.get("logical_devices", []))),
                             *( (key, value) for key, value in environment.get("environment", {}).items())])])
    cache_resources = cache.get("resources", {})
    output.extend([table(["Stage operation", "Wall seconds", "User CPU seconds", "System CPU seconds", "Peak RSS KiB"],
                        [("Sequence preparation", number(inputs.get("elapsed_seconds")), "not separately recorded", "not separately recorded", inputs.get("peak_rss_KiB")),
                         ("Paired phase-I extraction", number(cache_resources.get("wall_seconds")), number(cache_resources.get("user_cpu_seconds")),
                          number(cache_resources.get("system_cpu_seconds")), cache_resources.get("peak_RSS_KiB")),
                         ("All frozen checkpoint/RC inference", number(invariance.get("wall_seconds")), number(invariance.get("user_cpu_seconds")),
                          number(invariance.get("system_cpu_seconds")), invariance.get("peak_RSS_KiB"))]),
                   "The frozen runtime is CPython 3.13.0 / TensorFlow 2.20.0 / Keras 3.14.1 / NumPy 2.5.0, "
                   "with deterministic settings and suppressed bytecode writes to historical stages. Phase-I batch size is 32; phase-II inference batch size is 256. "
                   "Compute records describe the actual allocated execution, not a projected retraining cost.", ""])
    commands = [read_json(path) for path in sorted((stage / "provenance/commands").glob("*.completed.json"))]
    output.extend([table(["Command", "Return code", "Wall seconds", "Started UTC", "Completed UTC"],
                        ((r["label"], r["returncode"], number(r.get("elapsed_seconds")), r.get("started_utc"), r.get("completed_utc")) for r in commands)),
                   "Exact commands, environment settings, stdout/stderr and return codes are preserved under [provenance/commands](provenance/commands/). "
                   "Nonzero command returns are retained below rather than suppressed. Evaluator return code 3 denotes the saved FAIL/INCONCLUSIVE scientific stop decision, "
                   "not permission for an outcome-driven rerun.", "",
                   table(["Nonzero command", "Return code", "Log"], ((r["label"], r["returncode"], f"provenance/commands/{r['label']}.log") for r in commands if r["returncode"] != 0)),
                   "## Artifact index and final stop", ""])
    names = ("metrics.tsv", "bootstrap_summary.tsv", "paired_ablation.tsv", "seed_stability.tsv", "C_holdout_generalization.tsv",
             "holdout_decision.json", "fixed_threshold_operating.tsv", "fixed_threshold_intervals.tsv", "reliability_bins.tsv",
             "stratified_metrics.tsv", "bootstrap_audit.json", "evaluation_manifest.json")
    output.extend([table(["Artifact", "Stage-relative path"], ((name, f"[{name}]({evaluation_relative}/{name})") for name in names)),
                   "Complete per-seed/ensemble scores, repeated independent RC outputs, 2,000-replicate bootstrap records and score-linked components "
                   "are preserved in the new-stage prediction/result tables. This report reads summary artifacts only and does not recompute model scores, "
                   "bootstrap statistics or threshold calibration. The final checksum/freeze manifest records this report and all completed payloads; "
                   "shared V2 activity, decision and result registers receive append-only stage-completion entries.", "",
                   disposition(decision, qc), "",
                   "**Stop. Do not open the external functional benchmark or begin candidate/allele scoring without another explicit authorization.**", ""])
    return "\n".join(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=Path)
    parser.add_argument("--evaluation-relative", default="results/chr8_9_evaluation")
    parser.add_argument("--predictions-relative", default="predictions")
    args = parser.parse_args()
    stage = args.stage.resolve()
    path = stage / "INTERNAL_TEST_REPORT.md"
    text = report(stage, args.evaluation_relative, args.predictions_relative)
    with path.open("x") as stream:
        stream.write(text)
    print(json.dumps({"status": "REPORT_CREATED", "path": str(path), "bytes": path.stat().st_size,
                      "sha256": sha256(path), "model_or_evaluation_execution": False}))


if __name__ == "__main__":
    main()
