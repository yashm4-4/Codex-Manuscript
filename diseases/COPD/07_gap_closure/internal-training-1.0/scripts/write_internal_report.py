#!/usr/bin/env python3
"""Render completed internal-stage summaries; never run or recompute science.

The all-18 checkpoint freeze is checked before reading any outcome summary.
Only existing manifests, histories' summary records, evaluation summary tables,
audit records and optional calibration summaries are read. No network weights,
sequence arrays, row-level scores, benchmark or chr8/9 records are opened.
Output is exclusively created; the script cannot replace an existing report.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


MODELS = ("enhancer", "h3k27me3")
SEEDS = (104729, 130363, 155921)
CONFIGURATIONS = ("V2-A", "V2-B", "V2-C")


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def read_table(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_local_reference(stage, record):
    path = (stage / record["path"]).resolve()
    if not path.is_relative_to(stage):
        raise RuntimeError("Summary reference points outside the internal execution bundle")
    if sha256(path) != record["sha256"]:
        raise RuntimeError(f"Summary metadata hash mismatch: {path}")
    return path


def escape(value):
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value).replace("|", "\\|").replace("\n", " ")


def number(value, precision=7):
    if value is None or value == "":
        return "—"
    try:
        return format(float(value), f".{precision}g")
    except (ValueError, TypeError):
        return str(value)


def table(headers, rows):
    rows = list(rows)
    if not rows:
        return "No records.\n"
    return "\n".join([
        "| " + " | ".join(map(escape, headers)) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *("| " + " | ".join(map(escape, row)) + " |" for row in rows),
    ]) + "\n"


def context_label(model):
    return "Enhancer-associated" if model == "enhancer" else "H3K27me3-associated"


def confidence(row):
    return f"{number(row['point'])} [{number(row['lower95'])}, {number(row['upper95'])}]"


def derive_disposition(decision, calibration, qc_complete):
    """Describe saved decisions only; no performance threshold is recalculated."""
    if decision["both_C_contexts_pass"] is not True:
        return "STOP: at least one V2-C context failed or was inconclusive. No calibration and no A/B fallback."
    if calibration is None or calibration.get("status") != "PASS":
        return "C internal adequacy passed in both contexts; calibration is not complete/PASS, so model release remains blocked."
    if not qc_complete:
        return "C adequacy and calibration passed; final independent QC remains incomplete or non-PASS. No holdout release."
    return "READY FOR INVESTIGATOR REVIEW BEFORE OPENING CHR8–9. Both C contexts, conditional calibration and the recorded independent QC passed. Separate test/external authorization remains required."


def qc_phase_statuses(both_c_pass, latest):
    phases = ("inputs", "prepared", "cache", "runs", "packaging", "predictions", "evaluation", "calibration", "preservation")
    statuses = {
        phase: ("NOT_APPLICABLE (scientifically gated off)" if phase == "calibration" and not both_c_pass
                else latest.get(phase, {}).get("status", "NOT PRESENT"))
        for phase in phases
    }
    complete = all(status == "PASS" or status.startswith("NOT_APPLICABLE") for status in statuses.values())
    return statuses, complete


def report(stage):
    stage = Path(stage).resolve()
    freeze_path = stage / "provenance/checkpoint_freeze.json"
    freeze = read_json(freeze_path)
    expected = {(c, m, s) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}
    keys = {(r["configuration"], r["model"], int(r["seed"])) for r in freeze.get("checkpoints", [])}
    if (freeze.get("status") != "PASS" or freeze.get("all_18_selected_checkpoints_frozen") is not True
            or len(freeze.get("checkpoints", [])) != 18 or keys != expected
            or any(r.get("completed") is not True for r in freeze["checkpoints"])):
        raise RuntimeError("Refusing outcome summary access before the complete all-18 checkpoint freeze")
    # Outcome summaries are first opened ONLY after the gate above.
    evaluation = stage / "results/chr7_evaluation"
    decision_path = evaluation / "adequacy_decision.json"
    decision = read_json(decision_path)
    if set(decision["model_contexts"]) != set(MODELS):
        raise RuntimeError("Incomplete saved adequacy decision")
    if decision["both_C_contexts_pass"] != all(decision["model_contexts"][m]["status"] == "PASS" for m in MODELS):
        raise RuntimeError("Contradictory saved adequacy summary")
    for context in decision["model_contexts"].values():
        if context["status"] not in ("PASS", "FAIL", "INCONCLUSIVE"):
            raise RuntimeError("Unknown saved C adequacy status")
    calibration_path = stage / "results/chr7_calibration/C_region_thresholds.json"
    calibration = read_json(calibration_path) if calibration_path.exists() else None
    if calibration is not None and decision["both_C_contexts_pass"] is not True:
        raise RuntimeError("Calibration present despite C non-PASS; cannot issue a compliant report")
    if calibration is not None and set(calibration.get("models", {})) != set(MODELS):
        raise RuntimeError("Incomplete saved calibration summary")
    input_manifest_path = checked_local_reference(stage, freeze["input_manifest"])
    inputs = read_json(input_manifest_path)
    cache_path = stage / "cache/cache_manifest.json"
    cache = read_json(cache_path)
    invariance_path = stage / "predictions/real_network_invariance.json"
    invariance = read_json(invariance_path)
    outcomes_path = checked_local_reference(stage, freeze["training_outcomes"])
    outcomes = read_table(outcomes_path)
    metrics = read_table(evaluation / "metrics.tsv")
    stability = read_table(evaluation / "seed_stability.tsv")
    summaries = read_table(evaluation / "bootstrap_summary.tsv")
    paired = read_table(evaluation / "paired_ablation.tsv")
    gates = read_table(evaluation / "C_absolute_adequacy.tsv")
    bootstrap = read_json(evaluation / "bootstrap_audit.json")
    if len(outcomes) != 18 or len(metrics) != 24 or len(stability) != 6:
        raise RuntimeError("Completed report tables do not contain the expected prescribed fits/predictors")
    audits = []
    for path in sorted((stage / "provenance").glob("*.json")):
        record = read_json(path)
        if isinstance(record, dict) and record.get("validator") == "independent-execution-audit":
            audits.append((path, record))
    latest = {}
    for path, audit in sorted(audits, key=lambda item: item[1].get("created_utc", "")):
        latest[audit.get("phase", "unspecified")] = audit
    qc_statuses, qc_complete = qc_phase_statuses(decision["both_C_contexts_pass"], latest)
    packaging_path = stage / "provenance/checkpoint_packaging_final_validation.json"
    packaging_section = "### Selected-checkpoint container packaging\n\n"
    if packaging_path.exists():
        packaging = read_json(packaging_path)
        packaging_section += (
            f"Independent final archive audit: **{escape(packaging['status'])}**; "
            f"{len(packaging.get('checkpoints', []))} checkpoint containers recorded. "
            "[checkpoint_packaging_final_validation.json](provenance/checkpoint_packaging_final_validation.json) "
            "and its linked tensor-level ledger document the original, unchanged containers.\n\n"
        )
        packaging_section += table(
            ["Run", "Packaging checks", "Model tensors", "Optimizer tensors", "Optimizer slots", "Iteration", "All slots zero"],
            ((r["run_id"], "PASS" if r["packaging_checks_pass"] else "FAIL", r["model_tensor_count"],
              r["optimizer_tensor_count"], r["optimizer_slot_count"], r["optimizer_iterations"],
              r["checks"]["all_optimizer_slots_zero"]) for r in packaging.get("checkpoints", [])),
        )
    else:
        packaging_section += "Final all-18 archive-packaging validation is not yet present; packaging QC is not complete.\n"
    lines = ["# COPD V2 internally controlled training and chr7 evaluation", "",
             "Version: `internal-training-1.0`; immutable design: `pretraining-1.1`.", "",
             "## Decision and stop boundary", "", derive_disposition(decision, calibration, qc_complete), "",
             table(["V2-C context", "Saved absolute-adequacy decision", "Failed/inconclusive conditions"],
                   ((context_label(m), decision["model_contexts"][m]["status"],
                     "; ".join(row["gate"] for row in decision["model_contexts"][m]["failed_or_inconclusive_gates"]) or "none") for m in MODELS)),
             "Only region-label modeling was evaluated. A and B are diagnostic ablations, never eligible fallbacks or competitors for final-model selection. "
             "Neither chr8–9 performance nor the external COPD-V2-BENCH benchmark was opened in this stage. "
             "The COPD variant universe and 337 historical candidates were not scored. No new candidate list, fine-mapping, target-gene work or manuscript revision was performed.", "",
             "## Frozen design and execution", "",
             "A preserves historical labels and controls while correcting orientation handling. B replaces historical controls with frozen donor-accessible matched controls. "
             "C uses same-lobe positives and independently rematched controls. Consequently B→C changes both label membership and control matching, not labels alone. "
             "All comparisons below use the identical frozen C-task chr7 selection panel within each context; native-task metrics across different class definitions are not used to choose a model.", "",
             "All 18 fits use seeds `104729`, `130363`, `155921`, the unchanged V1 phase-II architecture, unweighted binary cross-entropy, "
             "Adadelta (learning rate 0.001, rho 0.95, epsilon 1e-7), batch size 256, maximum 50 epochs and patience 15. "
             "One nucleotide orientation per biological interval per epoch is sampled with probability 0.5. Frozen seed substreams assign orientations in genomic order before the independent epoch permutation. "
             "The earliest epoch with minimum native chr7 checkpoint-role symmetric BCE is retained. Selection/calibration rows cannot choose a checkpoint. "
             "The selected checkpoint contains all selected-epoch trainable and BatchNorm weights and is inference-only, not an optimizer-resumption checkpoint.", "",
             "Checkpoint-packaging qualification: Keras preserved `compile_config` during JSON reconstruction, so the `.keras` containers are not literally uncompiled or optimizer-free. "
             "They include fresh ancillary Adadelta state, not the selected epoch's training-optimizer state. Actual prediction loads use `compile=False`, which ignores that ancillary optimizer state. "
             "Selected network weights were restored and checked by bitwise serialization round-trip. These artifacts are frozen inference checkpoints only and must not be used to resume training. "
             "Any original implementation comment describing the saved container as uncompiled is superseded by the append-only packaging clarification; the checkpoints and scientific procedure were not rewritten or rerun.", "",
             "The final symmetric score is the float64 equal-weight mean, in the frozen seed order, of each seed's `(p_forward+p_nucleotide_RC)/2`. "
             "No representation-axis reversal, seed selection, seed weighting, mixed precision or hyperparameter search is permitted.", "",
             "## Phase-I cache and sequence integrity", "",
             table(["Record", "Value"], [
                 ("Unique genomic intervals used", cache["sequence_intervals"]),
                 ("Encoded sequence/RC identities", cache["encoded_identity_count"]),
                 ("Independent nucleotide-orientation network evaluations", cache["orientation_network_evaluations"]),
                 ("Per-orientation cache geometry", " × ".join(map(str, cache["shape"]))),
                 ("Cache dtype / representation", f"{cache['dtype']} / {cache['representation']}"),
                 ("Missing sequences / extraction failures", f"{cache['missing_sequences']} / {cache['failures']}"),
                 ("chr8–9 sequence extractions", cache["test_sequence_extractions"]),
                 ("Phase-I weight SHA-256", cache["phase_I_weights_sha256"]),
             ]),
             "Both canonical nucleotide sequence and its nucleotide reverse complement independently passed through the same frozen, inference-mode phase-I network. "
             "Each interval's forward-orientation flag maps the two caches back to genomic forward/RC. Ambiguous bases retain the frozen all-zero encoding after normalization; no new interval exclusion was introduced.", "",
             table(["Cache artifact", "Bytes", "SHA-256"],
                   ((r["path"], r["bytes"], r["sha256"]) for r in cache["files"] if Path(r["path"]).name in ("features_canonical.npy", "features_rc.npy"))),
             "Chromosome interval counts (training chromosomes plus chr7 only):", "",
             table(["Chromosome", "Intervals"], inputs["chromosome_counts"].items()),
             "## All 18 fit outcomes", "",
             table(["Configuration", "Context", "Seed", "Outcome", "Epochs", "Selected epoch", "Native checkpoint BCE", "Attempt", "Retries"],
                   ((r["configuration"], r["model"], r["seed"], r["status"], r["epochs_completed"], r["selected_epoch_one_based"],
                     number(r["checkpoint_symmetric_BCE"]), r["successful_attempt"], r["infrastructure_retries"]) for r in outcomes)),
             "Complete histories, epoch RNG hashes, checkpoint-improvement ledgers, selected-checkpoint records and serialized-weight round-trip checks remain under [runs](runs/). "
             "The all-18 checkpoint release and exact hashes are in [checkpoint_freeze.json](provenance/checkpoint_freeze.json).", "",
             packaging_section, "",
             "## Actual-network reverse-complement invariance", "",
             f"Saved gate: **{escape(invariance['status'])}**. Actual-network evidence: {escape(invariance['real_network'])}; "
             f"{invariance['n_seed_audits']} seed-level audits and {invariance['n_ensemble_audits']} ensemble audits. "
             "Every frozen internal chr7 validation sequence in the context-specific native/common union is covered. "
             "The second wrapper pass calls the trained network again on swapped nucleotide-orientation inputs, using the same frozen row order and batch size 256; original probabilities are not reused as the purported second network pass. "
             "The frozen comparison is `abs(Q(x)-Q(RC(x))) <= 1e-6 + 1e-6*abs(Q(RC(x)))`, with finite [0,1] checks on all probabilities.", "",
             table(["Configuration", "Context", "Unit", "Seed", "Sequences", "Status", "Failed", "Maximum absolute difference"],
                   ((r["configuration"], r["model"], r["unit"], r["seed"], r["n_sequences"], r["status"], r["n_failed"],
                     number(r["max_absolute_difference"], 17)) for r in invariance["audit_rows"])),
             "Detailed actual per-orientation and independent reverse-wrapper outputs are in [predictions](predictions/); "
             "[real_network_invariance.json](predictions/real_network_invariance.json) binds the checkpoints, inputs and implementation.", "",
             "## Common C-task chr7 selection performance", "",
             "The following are descriptive comparisons on the same frozen panel per context. AP is stepwise average precision; Brier skill uses the constant-prevalence baseline. "
             "Only C's predeclared absolute gates determine adequacy. Rounded display values below do not replace the full-precision TSV/JSON records.", "",
             table(["Context", "Configuration", "Predictor", "Positive / control", "Components", "AP", "AUROC", "Brier", "Brier skill"],
                   ((r["model"], r["configuration"], r["predictor"], f"{r['n_positive']} / {r['n_control']}", r["n_components"],
                     *[number(r[k]) for k in ("AP", "AUROC", "Brier", "BrierSkill")]) for r in metrics)),
             "### Fixed-seed stability", "",
             table(["Context", "Configuration", "AP 104729", "AP 130363", "AP 155921", "AP sample SD", "AP range"],
                   ((r["model"], r["configuration"], *[number(r[f"AP_seed_{s}"]) for s in SEEDS],
                     number(r["AP_sample_SD"]), number(r["AP_range"])) for r in stability)),
             "C requires sample SD ≤0.03 (`ddof=1`) and range ≤0.10 for each context. These limits were not revised after outcomes. A/B seed behavior is descriptive only.", "",
             "### Paired component-bootstrap uncertainty", "",
             table(["Context", "Components", "Attempts", "Valid", "Invalid", "Invalid fraction", "Status"],
                   ((m, bootstrap[m]["represented_components"], bootstrap[m]["attempted"], bootstrap[m]["valid"], bootstrap[m]["invalid"],
                     number(bootstrap[m]["invalid_fraction"]), bootstrap[m]["status"]) for m in MODELS)),
             "The frozen bootstrap uses seed 314159 and 2,000 valid genomic overlap/encoded-identity component draws, at most 20,000 attempts and invalid fraction <0.10. "
             "A/B/C share each draw and component multiplicities. The uncertainty is conditional on this donor and fixed matching design; it is not independent-donor or population generalization. "
             "AP minus prevalence and Brier skill are computed within each replicate before confidence limits, not by substituting the original prevalence afterward.", "",
             table(["Context", "Configuration", "Metric", "Point [95% component-bootstrap CI]"],
                   ((r["model"], r["configuration"], r["metric"], confidence(r)) for r in summaries if r["metric"] in ("AP", "AUROC", "Brier", "BrierSkill", "AP_gain"))),
             "### Diagnostic A→B, B→C and A→C ablations", "",
             table(["Context", "Difference", "Metric", "Point [95% paired CI]"],
                   ((r["model"], r["comparison"], r["metric"], confidence(r)) for r in paired)),
             "Negative Brier differences mean lower error for the first configuration. Comparative superiority, equivalence or noninferiority is not a retention gate. "
             "A favorable A/B result cannot rescue a failed or inconclusive C context.", "",
             "### Reliability and fixed strata", "",
             "[reliability_bins.tsv](results/chr7_evaluation/reliability_bins.tsv) reports the ten prespecified bins `[0,.1), …, [.9,1]` for every seed and ensemble. "
             "[stratified_metrics.tsv](results/chr7_evaluation/stratified_metrics.tsv) reports GC, train-derived ATAC rank and repeat quintiles, ATAC lobe signature, same-lobe mark support, chromosome and ambiguous-base status. "
             "Quintile cutpoints come only from retained C training positives and are fixed across configurations. All strata are descriptive; cells with fewer than 100 rows are explicitly flagged. "
             "No fitted recalibration, subgroup selection or threshold tuning was performed from these displays.", "",
             "## V2-C absolute-adequacy gates", "",
             "All conditions are conjunctive within each C context, and both contexts must pass. Boundary equality cannot pass the strict lower-confidence-bound gates. "
             "An inconclusive support/bootstrap result blocks release just as a failed gate does.", "",
             table(["Context", "Gate", "Status", "Saved observation", "Frozen criterion"],
                   ((r["model"], r["gate"], r["status"], r["observed"], r["criterion"]) for r in gates)),
             "Saved next action: `" + str(decision["next_action"]) + "`.", "",
             "## Conditional chr7 negative-control calibration", ""]
    if calibration is None:
        lines += ["No region thresholds were released. " + (
            "Calibration is prohibited because at least one C context did not pass; no A/B fallback is permitted."
            if not decision["both_C_contexts_pass"] else "Both C contexts passed, but a completed PASS calibration record is not present."), ""]
    else:
        lines += [f"Saved calibration status: **{escape(calibration['status'])}**. Calibration used only dedicated common chr7 calibration negative controls, after both C contexts and the six C checkpoints/ensemble were frozen. "
                  "No positive score chose a threshold. With `k=floor(0.05*n)` and descending negative scores `d`, the rule is `t=nextafter(d[k], +infinity)` and calls use `score>=t`.", "",
                  table(["Context", "Controls", "Control components", "Boundary (17 digits)", "Threshold (17 digits)", "Threshold (hex)", "Called / controls", "Observed row FPR", "Boundary ties"],
                        ((m, calibration["models"][m]["n_controls"], calibration["models"][m]["n_control_components"],
                          calibration["models"][m]["boundary_decimal_17g"], calibration["models"][m]["threshold_decimal_17g"],
                          calibration["models"][m]["threshold_hex"],
                          f"{calibration['models'][m]['n_called_controls']} / {calibration['models'][m]['n_controls']}",
                          number(calibration["models"][m]["observed_calibration_row_FPR"], 17), calibration["models"][m]["n_at_boundary"]) for m in MODELS)),
                  "All observations tied at the boundary are excluded; ties are never split. Exact decimal/hex serialization, empirical-row FPR checks and any threshold-above-one edge case are recorded in "
                  "[C_region_thresholds.json](results/chr7_calibration/C_region_thresholds.json). "
                  "This bounds the observed calibration-row fraction only; it does not establish a population-FPR guarantee, variant FPR, allele-effect FDR or causal probability. No allele-delta threshold was defined.", ""]
    environments = [("phase_I", read_json(stage / "cache/environment.json"))]
    resources = [("Phase-I cache extraction", cache["resources"]["wall_seconds"], cache["resources"].get("peak_RSS_KiB"))]
    for checkpoint in freeze["checkpoints"]:
        completed_path = checked_local_reference(stage, checkpoint["completed_record"])
        completed = read_json(completed_path)
        environments.append((checkpoint["run_id"], read_json(completed_path.parent / "environment.json")))
        resources.append((checkpoint["run_id"], completed["wall_seconds"], completed.get("peak_RSS_KiB")))
    environments.append(("chr7_inference", read_json(stage / "predictions/environment.json")))
    resources.append(("All-context chr7 inference and RC audit", invariance["wall_seconds"], invariance.get("peak_RSS_KiB")))
    software_rows = []
    for name, environment in environments:
        packages = environment["packages"]
        software_rows.append((name, environment["python"], packages.get("tensorflow"), packages.get("keras"), packages.get("numpy"),
                              environment.get("environment", {}).get("SLURM_JOB_ID")))
    lines += ["## Software, command provenance and compute", "",
              "[provenance/commands](provenance/commands/) retains exact commands, launcher hashes, environment fields, streamed logs, return codes and elapsed times. "
              "Per-fit/cache/inference environment records preserve Python/package versions, deterministic settings, GPU identity and TensorFlow build metadata. "
              "The following are recorded process elapsed times, including I/O, hashing and serialization—not measured active GPU-kernel time or full scheduler billing.", "",
              table(["Process", "Recorded wall seconds", "Peak RSS (KiB)"],
                    ((name, number(seconds, 10), rss) for name, seconds, rss in resources)),
              f"Sum of recorded successful one-GPU process walltimes: {number(sum(float(r[1]) for r in resources) / 3600, 10)} GPU-process-hours. "
              "This accounting label does not imply 100% GPU utilization; failed-attempt allocations, queue wait and scheduler overhead are not included.", "",
              table(["Process", "Python", "TensorFlow", "Keras", "NumPy", "Slurm job"], software_rows),
              ]
    accounting_path = stage / "provenance/slurm_accounting.tsv"
    if accounting_path.exists():
        with accounting_path.open(newline="") as stream:
            header = stream.readline()
            stream.seek(0)
            accounting = list(csv.DictReader(stream, delimiter="\t" if "\t" in header else "|"))
        wanted = ("JobIDRaw", "JobID", "JobName", "State", "ExitCode", "ElapsedRaw", "Elapsed", "AllocTRES", "TotalCPU", "MaxRSS", "TRESUsageInAve", "TRESUsageInTot")
        columns = [key for key in wanted if accounting and key in accounting[0]]
        lines += ["### Scheduler accounting", "",
                  "[slurm_accounting.tsv](provenance/slurm_accounting.tsv) preserves scheduler allocation and job-step accounting. "
                  "Allocation and `.batch`/other step rows overlap and must not be summed together. Allocated GPU time is reserved resource time, not measured GPU-kernel utilization.", "",
                  table(columns, ([row.get(key) for key in columns] for row in accounting)) if columns else "The raw accounting record is linked above; its fields are not summarized by this formatter.", ""]
    else:
        lines += ["Scheduler accounting was not yet present when this report was formatted; the recorded process-walltime table is not a substitute for final allocation accounting.", ""]
    lines += ["### Infrastructure failures and retries", ""]
    retries = freeze.get("infrastructure_failures_and_retries", [])
    if not retries:
        lines += ["The checkpoint freeze records no phase-II infrastructure retries; each prescribed seed/configuration/context has exactly one successful fit.", ""]
    else:
        lines += [table(["Run", "Failed attempt", "Infrastructure reason", "Identical-seed retry authorized"],
                        ((r["run_id"], r["attempt"], r["retry_authorization"]["reason"], r["retry_authorization"]["same_seed_identical_procedure"]) for r in retries)),
                  "Failures and retries are preserved; no seed was replaced based on performance. The checkpoint freeze binds each retry authorization and available failure record.", ""]
    other_failures = []
    for path in sorted(stage.rglob("failure.json")):
        if "runs" in path.relative_to(stage).parts:
            continue
        failure = read_json(path)
        other_failures.append((str(path.relative_to(stage)), failure.get("status"), failure.get("error")))
    if other_failures:
        lines += ["Other recorded execution failures:", "", table(["Record", "Status", "Error"], other_failures)]
    nonzero_commands = []
    for path in sorted((stage / "provenance/commands").glob("*.completed.json")):
        command = read_json(path)
        if command.get("returncode", 0) != 0:
            nonzero_commands.append((command.get("label"), command.get("returncode"), str(path.relative_to(stage))))
    if nonzero_commands:
        lines += ["Nonzero command exits (including administrative failures or scientific STOP codes, distinct from failed model fits):", "",
                  table(["Command", "Exit code", "Preserved completion record"], nonzero_commands)]
    if (stage / "provenance/calibration_gate_path_resolution.json").exists():
        lines += ["The original calibration-gate command stopped before creating a release or calibration output because the evaluation manifest used repository-relative artifact paths while the gate helper resolved them from the stage directory. "
                  "An independently reviewed, append-only metadata-path adapter verified the original bytes and hashes and resolved those references without modifying the frozen helper, evaluation manifest, metrics, selected checkpoints or scientific design. "
                  "Training, inference and evaluation were not rerun. The correction and exact original-to-resolved mappings are preserved in "
                  "[calibration_gate_path_resolution.json](provenance/calibration_gate_path_resolution.json).", ""]
    lines += ["## Independent validation and preservation", "",
              table(["Audit record", "Phase", "Status", "Checks", "Failures"],
                    ((str(path.relative_to(stage)), audit.get("phase"), audit["status"], audit["checks"], audit["failures"]) for path, audit in audits)),
              "Phase-level audit applicability and status used for the report's review-readiness wording. Calibration QC is not required when either C context fails or is inconclusive because calibration is scientifically prohibited:", "",
              table(["Phase", "Latest applicable independent audit"], qc_statuses.items()),
              f"Applicable execution QC complete/PASS: {escape(qc_complete)}. A valid scientific STOP bundle can have complete execution QC without being ready for holdout evaluation.", "",
              "Historical pretraining-1.0 and pretraining-1.1 remain immutable inputs. Shared-register activity, decision and result updates are recorded separately from the frozen scientific artifacts. "
              "The independent validator checks exact interval/role membership, source hashes, all prescribed seeds, RNG histories, earliest-minimum checkpoints, serialized outputs, probability reductions, "
              "component-bootstrap calculations, frozen adequacy rules and conditional calibration. A PASS artifact audit does not convert failed scientific adequacy into PASS.", "",
              "## Interpretation limits and required stop", "",
              "This is an internally controlled, single-donor region-label experiment, not independent population validation. Bulk same-lobe provenance does not establish matched cells or aliquots. "
              "H3K27me3-associated labels are not experimentally proven silencer labels. Matching reduces specified measured imbalances but cannot remove unmeasured confounding, and the repaired design does not imply local-caliper support for every matched pair. "
              "Improved region AP/AUROC/Brier does not establish biologically valid REF–ALT effects, causal variants, fine-mapping probabilities or variant-level false-positive control.", "",
              "Stop after this internal stage. Do not score chr8–9, open the external functional benchmark, score the candidate universe, generate a candidate list or redesign/rerun this version in response to model outcomes. "
              "Any future test/external stage and any scientific redesign require separate investigator authorization. The frozen checkpoints, ensemble implementation and any permissible region thresholds must remain fixed.", "",
              "## Summary-source integrity", "",
              "This report formats existing stored summaries only. It does not rerun models, recompute metrics, bootstrap, recalibrate thresholds or read row-level scores. "
              "Display rounding is descriptive; all decisions and exact threshold values originate from the linked frozen records.", ""]
    sources = [freeze_path, input_manifest_path, cache_path, invariance_path, outcomes_path, decision_path,
               *[evaluation / name for name in ("metrics.tsv", "seed_stability.tsv", "bootstrap_summary.tsv", "paired_ablation.tsv", "C_absolute_adequacy.tsv", "bootstrap_audit.json")]]
    if calibration is not None:
        sources.append(calibration_path)
    if packaging_path.exists():
        sources.append(packaging_path)
    lines += [table(["Summary source", "SHA-256"], ((str(path.relative_to(stage)), sha256(path)) for path in sources))]
    return "\n".join(lines).rstrip() + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    args = parser.parse_args()
    stage = args.stage.resolve()
    output = stage / "INTERNAL_TRAINING_REPORT.md"
    if output.exists():
        raise FileExistsError("Report already exists; refusing automatic overwrite")
    rendered = report(stage)
    with output.open("x") as stream:
        stream.write(rendered)
    print(f"Wrote summary-only report: {output}")


if __name__ == "__main__":
    main()
