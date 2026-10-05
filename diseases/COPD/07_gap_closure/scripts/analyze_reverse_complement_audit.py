#!/usr/bin/env python3
"""Analyze the locked COPD frozen-model orientation diagnostic without editing V1.

The authoritative forward values are serialized V1 R003 REF/ALT scores. All
decisions preserve V1 blacklist eligibility, fixed thresholds, and R010 order.
No metric estimates a causal probability or nominates replacement candidates.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import pearsonr, spearmanr


SCRIPT = Path(__file__).resolve()
V2 = SCRIPT.parents[1]
COPD = V2.parent
ROOT = COPD.parents[1]
RESULTS = V2 / "results"
PROVENANCE = V2 / "provenance"
SPEC = PROVENANCE / "COPD-V2-RC_analysis_specification.md"
SPEC_SHA = "09d3656f3e06c24fa5402c28b976d6be06fff91604ad0b7b29c3150120d76ac1"
MODELS = ("enhancer", "silencer")
CONTEXTS = ("neither", "enhancer_only", "H3K27me3_only", "both")
CLASSES = ("SNV", "indel_or_complex")
EXPECTED_THRESHOLDS = {
    ("enhancer", "SNV"): (0.643623, 0.05706318769999998),
    ("enhancer", "indel_or_complex"): (0.643623, 0.04936093850000001),
    ("silencer", "SNV"): (0.58505, 0.028802613899999996),
    ("silencer", "indel_or_complex"): (0.58505, 0.026183359500000003),
}
INPUT_PATHS = {
    "forward": COPD / "04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz",
    "prioritized": COPD / "04_modeling/results/COPD-S4-R004_prioritized_candidates.tsv.gz",
    "cutoffs": COPD / "04_modeling/results/COPD-S4-R004_delta_thresholds.tsv",
    "frozen337": COPD / "04_modeling/results/COPD-S4-R010_THE_LIST.tsv",
    "shortlist": COPD / "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv",
    "phenotype": RESULTS / "COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv",
    "phenotype_shortlist": RESULTS / "COPD-V2-PHENO-R007_shortlist_phenotype_support.tsv",
    "literature": COPD / "04_modeling/results/COPD-S4-R009_literature_functional_variant_recovery.tsv",
    "sequence_qc": RESULTS / "COPD-V2-RC-R002_sequence_transformation_qc.tsv.gz",
    "rc": RESULTS / "COPD-V2-RC-R001_reverse_complement_scores.tsv.gz",
    "p99": PROVENANCE / "COPD-V2-RC_forward_p99_cutoffs.tsv",
    "preflight": PROVENANCE / "COPD-V2-RC_preflight_manifest.json",
    "specification": SPEC,
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, low_memory=False, float_precision="round_trip")


def bools(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin(("true", "1", "yes"))


def tokens(value):
    return [x for x in str(value).split(";") if x and x.lower() != "nan"]


def fraction(a, b):
    return float(a / b) if b else None


def correlations(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2:
        return {"pearson_r": None, "spearman_rho": None, "correlation_status": "not_evaluable_n_lt_2"}
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return {"pearson_r": None, "spearman_rho": None, "correlation_status": "not_evaluable_constant_vector"}
    return {"pearson_r": float(pearsonr(x, y).statistic),
            "spearman_rho": float(spearmanr(x, y).statistic),
            "correlation_status": "evaluable"}


def distribution(values, prefix=""):
    values = np.asarray(values, float)
    names = ("mean", "median", "p90", "p95", "p99", "maximum")
    if not len(values):
        return {prefix + name: None for name in names}
    vals = [float(values.mean())] + [float(x) for x in np.quantile(values, [.5, .9, .95, .99], method="linear")] + [float(values.max())]
    return dict(zip([prefix + name for name in names], vals))


def score_metrics(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    difference = y - x
    absolute = np.abs(difference)
    return {"n": len(x), **correlations(x, y),
            "signed_mean_difference": float(difference.mean()) if len(x) else None,
            **distribution(absolute, "abs_difference_"),
            **{f"fraction_abs_difference_gt_{label}": fraction(int((absolute > value).sum()), len(x))
               for label, value in (("1e_5", 1e-5), ("0_02", .02), ("0_10", .10))}}


def signs(delta, epsilon=0):
    delta = np.asarray(delta, float)
    return np.where(delta > epsilon, 1, np.where(delta < -np.asarray(epsilon), -1, 0))


def sign_metrics(fwd, rev, cutoff):
    fwd, rev, cutoff = np.asarray(fwd, float), np.asarray(rev, float), np.asarray(cutoff, float)
    raw_f, raw_r = signs(fwd), signs(rev)
    nonzero = (raw_f != 0) & (raw_r != 0)
    epsilon = np.maximum(2e-5, .1 * cutoff)
    meaningful_f, meaningful_r = signs(fwd, epsilon), signs(rev, epsilon)
    selected = meaningful_f != 0
    same = selected & (meaningful_f == meaningful_r)
    collapse = selected & (meaningful_r == 0)
    opposite = selected & (meaningful_f == -meaningful_r)
    return {
        "raw_both_nonzero_n": int(nonzero.sum()),
        "raw_sign_agree_n": int(((raw_f == raw_r) & nonzero).sum()),
        "raw_nonzero_sign_agreement": fraction(int(((raw_f == raw_r) & nonzero).sum()), int(nonzero.sum())),
        "raw_gain_to_loss_n": int(((raw_f == 1) & (raw_r == -1)).sum()),
        "raw_loss_to_gain_n": int(((raw_f == -1) & (raw_r == 1)).sum()),
        "forward_meaningful_n": int(selected.sum()),
        "meaningful_direction_retained_n": int(same.sum()),
        "meaningful_direction_retention": fraction(int(same.sum()), int(selected.sum())),
        "meaningful_near_zero_collapse_n": int(collapse.sum()),
        "meaningful_near_zero_collapse_fraction": fraction(int(collapse.sum()), int(selected.sum())),
        "meaningful_opposite_direction_n": int(opposite.sum()),
        "meaningful_opposite_direction_fraction": fraction(int(opposite.sum()), int(selected.sum())),
    }


def delta_metrics(fwd, rev, cutoff):
    fwd, rev, cutoff = np.asarray(fwd, float), np.asarray(rev, float), np.asarray(cutoff, float)
    difference = rev - fwd
    return {"n": len(fwd), **correlations(fwd, rev),
            "signed_mean_difference": float(difference.mean()) if len(fwd) else None,
            **distribution(np.abs(difference), "delta_abs_difference_"),
            **distribution(np.abs(np.abs(rev) - np.abs(fwd)), "absolute_effect_magnitude_difference_"),
            **distribution(np.abs(difference) / cutoff, "normalized_delta_disagreement_"),
            **sign_metrics(fwd, rev, cutoff)}


def call_summary(fwd, rev, well_separated):
    fwd, rev, well_separated = np.asarray(fwd, bool), np.asarray(rev, bool), np.asarray(well_separated, bool)
    both = int((fwd & rev).sum())
    forward_only, reverse_only = int((fwd & ~rev).sum()), int((~fwd & rev).sum())
    neither = int((~fwd & ~rev).sum())
    well_n, well_lost = int(well_separated.sum()), int((well_separated & ~rev).sum())
    return {"n": len(fwd), "forward_positive_n": int(fwd.sum()), "rc_positive_n": int(rev.sum()),
            "both_positive_n": both, "forward_only_n": forward_only, "reverse_only_n": reverse_only,
            "neither_n": neither, "transition_1_to_1": both, "transition_1_to_0": forward_only,
            "transition_0_to_1": reverse_only, "transition_0_to_0": neither,
            "jaccard": fraction(both, both + forward_only + reverse_only),
            "forward_retention": fraction(both, int(fwd.sum())),
            "well_separated_forward_positive_n": well_n,
            "well_separated_lost_n": well_lost,
            "well_separated_retention": fraction(well_n - well_lost, well_n)}


def contexts(enhancer, silencer):
    return np.select([enhancer & silencer, enhancer, silencer], ["both", "enhancer_only", "H3K27me3_only"], default="neither")


def positive_call(eligible, ref, alt, region_threshold, delta_threshold):
    return np.asarray(eligible, bool) & (np.maximum(ref, alt) >= region_threshold) & (np.abs(np.asarray(alt) - np.asarray(ref)) >= delta_threshold)


def severity_label(material, negligible):
    return "scientifically_material" if material else "negligible" if negligible else "modest"


def input_frames(paths):
    tables = {name: read_table(path) for name, path in paths.items() if str(path).endswith((".tsv", ".tsv.gz"))}
    fwd, rc, prio = tables["forward"], tables["rc"], tables["prioritized"]
    for name in ("forward", "rc", "prioritized", "frozen337", "shortlist", "phenotype", "phenotype_shortlist", "sequence_qc"):
        table = tables[name]
        assert "candidate_record_id" in table, (name, "missing ID")
        assert not table.candidate_record_id.duplicated().any(), (name, "duplicate ID")
    assert len(fwd) == len(rc) == 15303
    assert fwd.candidate_record_id.tolist() == rc.candidate_record_id.tolist(), "forward/RC record order mismatch"
    assert len(prio) == 15389
    assert set(fwd.candidate_record_id) == set(prio.loc[bools(prio.scored), "candidate_record_id"])
    assert len(tables["sequence_qc"]) == 15303
    assert set(tables["sequence_qc"].candidate_record_id) == set(fwd.candidate_record_id)
    assert len(tables["frozen337"]) == len(tables["phenotype"]) == 337
    assert tables["frozen337"].candidate_record_id.tolist() == tables["phenotype"].candidate_record_id.tolist()
    assert tables["frozen337"].predicted_causal_priority_rank.astype(int).tolist() == list(range(1, 338))
    assert len(tables["shortlist"]) == len(tables["phenotype_shortlist"]) == 12
    assert tables["shortlist"].candidate_record_id.tolist() == tables["phenotype_shortlist"].candidate_record_id.tolist()
    assert tables["shortlist"].experimental_shortlist_rank.astype(int).tolist() == list(range(1, 13))
    assert len(tables["cutoffs"]) == len(tables["p99"]) == 4
    for row in tables["cutoffs"].itertuples():
        expected = EXPECTED_THRESHOLDS[(row.model_type, row.variant_class_group)]
        assert abs(float(row.region_score_cutoff) - expected[0]) < 1e-14
        assert abs(float(row.abs_delta_cutoff) - expected[1]) < 1e-14
    for name in ("forward", "rc"):
        numeric = tables[name].drop(columns="candidate_record_id").astype(float)
        assert np.isfinite(numeric.to_numpy()).all(), name
        assert numeric[[c for c in numeric if c.endswith("_score")]].ge(0).all().all(), name
        assert numeric[[c for c in numeric if c.endswith("_score")]].le(1).all().all(), name
    return tables


def build_comparisons(tables):
    forward = tables["forward"].set_index("candidate_record_id", drop=False)
    ids = forward.index
    prio = tables["prioritized"].set_index("candidate_record_id").loc[ids]
    rc = tables["rc"].set_index("candidate_record_id").loc[ids]
    sequence = tables["sequence_qc"].set_index("candidate_record_id").loc[ids]
    keep = ["chromosome_grch38", "position_grch38", "ref", "alt", "variant_class_group", "gws_tag_ids", "source_focal_tags", "ld_panels"]
    out = prio[keep].copy()
    out.insert(0, "candidate_record_id", ids)
    out.insert(0, "v2_record_id", [f"COPD-V2-RC-C{i:05d}" for i in range(1, len(out) + 1)])
    out["blacklisted"] = bools(prio.blacklisted)
    out["causal_call_eligible"] = bools(prio.causal_call_eligible)
    assert int(out.blacklisted.sum()) == 20
    assert int(out.causal_call_eligible.sum()) == 15283
    assert (out.causal_call_eligible == ~out.blacklisted).all()
    assert dict(out.loc[out.causal_call_eligible].variant_class_group.value_counts()) == {"SNV": 13747, "indel_or_complex": 1536}
    out["ref_N_count"] = sequence.ref_N_count.astype(int)
    out["alt_N_count"] = sequence.alt_N_count.astype(int)
    out["has_ambiguous_N"] = (out.ref_N_count + out.alt_N_count) > 0
    for column in ("ref_allele_length", "alt_allele_length"):
        out[column] = sequence[column].astype(int)
    frozen = tables["frozen337"].set_index("candidate_record_id")
    out["in_frozen_337"] = out.index.isin(frozen.index)
    out["phenotype_stratum_applicable"] = out["in_frozen_337"]
    out["v1_priority_rank"] = out.index.to_series().map(frozen.predicted_causal_priority_rank).astype("Int64")
    pheno = tables["phenotype"].set_index("candidate_record_id")
    for column in ("primary_retained", "GCST90244098_only_study"):
        out[column] = out.index.to_series().map(bools(pheno[column])).fillna(False).astype(bool)
    for column in ("support_state", "gws_supporting_studies", "primary_supporting_studies", "v1_component_id"):
        out["phenotype_" + column] = out.index.to_series().map(pheno[column]).fillna("not_applicable_outside_frozen_337")
    assert int(out.primary_retained.sum()) == 184
    assert int(out.GCST90244098_only_study.sum()) == 124
    shortlist = tables["shortlist"].set_index("candidate_record_id")
    out["in_shortlist_12"] = out.index.isin(shortlist.index)
    out["v1_shortlist_rank"] = out.index.to_series().map(shortlist.experimental_shortlist_rank).astype("Int64")
    out["v1_selection_anchor"] = out.index.to_series().map(shortlist.selection_anchor_criterion).fillna("not_applicable")
    literal_ids = set()
    for value in tables["literature"].exact_candidate_record_ids:
        literal_ids.update(tokens(value))
    assert literal_ids == {"4:88963935:G:T", "4:88962828:C:T"}
    assert literal_ids.issubset(set(ids))
    out["in_exact_v1_literature"] = out.index.isin(literal_ids)
    p99 = tables["p99"].set_index(["model_type", "variant_class_group"])
    cutoffs = tables["cutoffs"].set_index(["model_type", "variant_class_group"])
    for model in MODELS:
        T = float(cutoffs.loc[(model, "SNV"), "region_score_cutoff"])
        D = out.variant_class_group.map({c: float(cutoffs.loc[(model, c), "abs_delta_cutoff"]) for c in CLASSES})
        out[f"{model}_region_threshold"] = T
        out[f"{model}_abs_delta_threshold"] = D
        out[f"{model}_meaningful_epsilon"] = np.maximum(2e-5, .1 * D)
        for orientation, scores in (("forward", forward), ("rc", rc)):
            for allele in ("ref", "alt"):
                out[f"{model}_{orientation}_{allele}_score"] = scores[f"{model}_{allele}_score"].astype(float)
            ref, alt = out[f"{model}_{orientation}_ref_score"], out[f"{model}_{orientation}_alt_score"]
            delta = alt - ref
            region = np.maximum(ref, alt)
            out[f"{model}_{orientation}_serialized_delta"] = scores[f"{model}_delta_alt_minus_ref"].astype(float)
            out[f"{model}_{orientation}_delta"] = delta
            out[f"{model}_{orientation}_serialization_delta_rounding"] = delta - out[f"{model}_{orientation}_serialized_delta"]
            out[f"{model}_{orientation}_abs_delta"] = np.abs(delta)
            out[f"{model}_{orientation}_region_score"] = region
            out[f"{model}_{orientation}_region_margin"] = region - T
            out[f"{model}_{orientation}_delta_margin"] = np.abs(delta) - D
            out[f"{model}_{orientation}_normalized_delta_margin"] = np.abs(delta) / D - 1
            out[f"{model}_{orientation}_region_gate"] = region >= T
            out[f"{model}_{orientation}_delta_gate"] = np.abs(delta) >= D
            out[f"{model}_{orientation}_call"] = positive_call(out.causal_call_eligible, ref, alt, T, D)
            out[f"{model}_{orientation}_raw_sign"] = signs(delta)
            out[f"{model}_{orientation}_meaningful_sign"] = signs(delta, out[f"{model}_meaningful_epsilon"])
        for quantity in ("ref_score", "alt_score", "region_score", "delta", "abs_delta"):
            difference = out[f"{model}_rc_{quantity}"] - out[f"{model}_forward_{quantity}"]
            out[f"{model}_{quantity}_signed_disagreement"] = difference
            out[f"{model}_{quantity}_absolute_disagreement"] = difference.abs()
        out[f"{model}_abs_delta_disagreement"] = out[f"{model}_delta_absolute_disagreement"]
        out[f"{model}_normalized_delta_disagreement"] = out[f"{model}_delta_absolute_disagreement"] / D
        out[f"{model}_raw_sign_changed"] = out[f"{model}_forward_raw_sign"] != out[f"{model}_rc_raw_sign"]
        out[f"{model}_meaningful_sign_changed"] = out[f"{model}_forward_meaningful_sign"] != out[f"{model}_rc_meaningful_sign"]
        out[f"{model}_forward_score_borderline"] = out[f"{model}_forward_region_margin"].abs() <= .02
        out[f"{model}_forward_delta_borderline"] = out[f"{model}_forward_normalized_delta_margin"].abs() <= .10
        out[f"{model}_forward_any_borderline"] = out[f"{model}_forward_score_borderline"] | out[f"{model}_forward_delta_borderline"]
        out[f"{model}_forward_well_separated"] = out.causal_call_eligible & (out[f"{model}_forward_region_score"] >= T + .05) & (out[f"{model}_forward_abs_delta"] >= 1.5 * D)
        out[f"{model}_forward_tail95"] = out.causal_call_eligible & out[f"{model}_forward_delta_gate"]
        for c in CLASSES:
            eligible = out.causal_call_eligible & out.variant_class_group.eq(c)
            q = float(np.quantile(out.loc[eligible, f"{model}_forward_abs_delta"], .99, method="linear"))
            assert abs(q - float(p99.loc[(model, c), "forward_abs_delta_p99"])) < 1e-13, (model, c, q)
            assert int(p99.loc[(model, c), "n_eligible_background"]) == int(eligible.sum())
        out[f"{model}_forward_abs_delta_p99_cutoff"] = out.variant_class_group.map({c: float(p99.loc[(model, c), "forward_abs_delta_p99"]) for c in CLASSES})
        out[f"{model}_forward_tail99"] = out.causal_call_eligible & (out[f"{model}_forward_abs_delta"] >= out[f"{model}_forward_abs_delta_p99_cutoff"])
        F, R = out[f"{model}_forward_call"], out[f"{model}_rc_call"]
        out[f"{model}_call_changed"] = F != R
        out[f"{model}_lost_forward_call"] = F & ~R
        out[f"{model}_gained_rc_only_call"] = ~F & R
        out[f"{model}_lost_well_separated_call"] = out[f"{model}_forward_well_separated"] & ~R
        assert (F == bools(prio[f"predicted_causal_{model}"])).all(), model
    for orientation in ("forward", "rc"):
        E, S = out[f"enhancer_{orientation}_call"], out[f"silencer_{orientation}_call"]
        out[f"{orientation}_union_call"] = E | S
        out[f"{orientation}_model_context"] = contexts(E, S)
    out["forward_union_well_separated"] = out.enhancer_forward_well_separated | out.silencer_forward_well_separated
    out["union_call_changed"] = out.forward_union_call != out.rc_union_call
    out["lost_union_call"] = out.forward_union_call & ~out.rc_union_call
    out["retained_union_call"] = out.forward_union_call & out.rc_union_call
    out["lost_any_original_model_call"] = out.enhancer_lost_forward_call | out.silencer_lost_forward_call
    out["gained_any_additional_rc_only_model_call"] = out.enhancer_gained_rc_only_call | out.silencer_gained_rc_only_call
    out["lost_any_well_separated_model_call"] = out.enhancer_lost_well_separated_call | out.silencer_lost_well_separated_call
    out["exact_model_context_stable"] = out.forward_model_context == out.rc_model_context
    out["support_switching"] = out.lost_any_original_model_call & out.gained_any_additional_rc_only_model_call
    out["pure_enhancer_to_H3K27me3_switch"] = out.forward_model_context.eq("enhancer_only") & out.rc_model_context.eq("H3K27me3_only")
    out["pure_H3K27me3_to_enhancer_switch"] = out.forward_model_context.eq("H3K27me3_only") & out.rc_model_context.eq("enhancer_only")
    assert int(out.enhancer_forward_call.sum()) == 175
    assert int(out.silencer_forward_call.sum()) == 199
    assert int((out.enhancer_forward_call & out.silencer_forward_call).sum()) == 37
    assert int(out.forward_union_call.sum()) == 337
    assert set(out.index[out.forward_union_call]) == set(frozen.index)
    return out.copy()


def strata(out):
    masks = {"all_15303": np.ones(len(out), bool), "SNV": out.variant_class_group.eq("SNV"),
             "indel_or_complex": out.variant_class_group.eq("indel_or_complex"),
             "call_eligible": out.causal_call_eligible, "blacklisted": out.blacklisted,
             "frozen_337": out.in_frozen_337, "direct_COPD_184": out.primary_retained,
             "sole_GCST90244098_124": out.GCST90244098_only_study,
             "shortlist_12": out.in_shortlist_12, "exact_v1_literature": out.in_exact_v1_literature}
    if out.has_ambiguous_N.any():
        masks["N_containing"] = out.has_ambiguous_N
        masks["N_free"] = ~out.has_ambiguous_N
    return masks


def support_summary(frame):
    n = len(frame)
    fields = ("exact_model_context_stable", "lost_union_call", "retained_union_call",
              "lost_any_original_model_call", "gained_any_additional_rc_only_model_call",
              "lost_any_well_separated_model_call", "support_switching",
              "pure_enhancer_to_H3K27me3_switch", "pure_H3K27me3_to_enhancer_switch")
    result = {"n": n, "forward_union_positive_n": int(frame.forward_union_call.sum()),
              "rc_union_positive_n": int(frame.rc_union_call.sum()),
              "either_model_call_changed_n": int((~frame.exact_model_context_stable).sum()),
              "both_model_calls_changed_n": int((frame.enhancer_call_changed & frame.silencer_call_changed).sum())}
    for field in fields:
        result[field + "_n"] = int(frame[field].sum())
        result[field + "_fraction"] = fraction(int(frame[field].sum()), n)
    for model in MODELS:
        changed = frame[f"{model}_call_changed"]
        result[f"{model}_changed_call_n"] = int(changed.sum())
        result[f"{model}_changed_call_forward_borderline_n"] = int((changed & frame[f"{model}_forward_any_borderline"]).sum())
        result[f"{model}_changed_call_forward_borderline_fraction"] = fraction(result[f"{model}_changed_call_forward_borderline_n"], result[f"{model}_changed_call_n"])
    changed_n = result["enhancer_changed_call_n"] + result["silencer_changed_call_n"]
    borderline_n = result["enhancer_changed_call_forward_borderline_n"] + result["silencer_changed_call_forward_borderline_n"]
    result["changed_model_decisions_n"] = changed_n
    result["changed_model_decisions_forward_borderline_n"] = borderline_n
    result["changed_model_decisions_forward_borderline_fraction"] = fraction(borderline_n, changed_n)
    return result


def metric_tables(out):
    rows = {name: [] for name in ("score", "delta", "sign", "tail", "call", "context", "proximity", "bins", "support", "gate")}
    labels = {-1: "loss", 0: "zero", 1: "gain"}
    for stratum, mask in strata(out).items():
        frame = out.loc[mask]
        rows["support"].append({"stratum": stratum, **support_summary(frame)})
        for model in MODELS:
            D = frame[f"{model}_abs_delta_threshold"].to_numpy()
            fdelta, rdelta = frame[f"{model}_forward_delta"].to_numpy(), frame[f"{model}_rc_delta"].to_numpy()
            for quantity in ("ref", "alt", "region"):
                rows["score"].append({"stratum": stratum, "model_type": model, "score_type": quantity,
                                      **score_metrics(frame[f"{model}_forward_{quantity}_score"], frame[f"{model}_rc_{quantity}_score"])})
            rows["delta"].append({"stratum": stratum, "model_type": model, **delta_metrics(fdelta, rdelta, D)})
            for sign_mode in ("raw", "meaningful"):
                fs = frame[f"{model}_forward_{sign_mode}_sign"]
                rs = frame[f"{model}_rc_{sign_mode}_sign"]
                for source in (-1, 0, 1):
                    for target in (-1, 0, 1):
                        count = int(((fs == source) & (rs == target)).sum())
                        rows["sign"].append({"stratum": stratum, "model_type": model, "sign_mode": sign_mode,
                                             "forward_sign": labels[source] if sign_mode == "raw" or source else "near_zero",
                                             "rc_sign": labels[target] if sign_mode == "raw" or target else "near_zero",
                                             "n": count, "stratum_n": len(frame),
                                             "forward_sign_n": int((fs == source).sum()),
                                             "fraction_within_forward_sign": fraction(count, int((fs == source).sum()))})
            for percentile in (95, 99):
                tail = frame.loc[frame[f"{model}_forward_tail{percentile}"]]
                F, R = tail[f"{model}_forward_call"], tail[f"{model}_rc_call"]
                rows["tail"].append({"stratum": stratum, "model_type": model,
                                     "forward_tail": f"matched_class_p{percentile}",
                                     **delta_metrics(tail[f"{model}_forward_delta"], tail[f"{model}_rc_delta"], tail[f"{model}_abs_delta_threshold"]),
                                     "forward_call_n": int(F.sum()), "rc_call_n": int(R.sum()),
                                     "retained_forward_call_n": int((F & R).sum()),
                                     "forward_call_retention": fraction(int((F & R).sum()), int(F.sum())),
                                     "rc_positive_fraction_of_whole_forward_tail": fraction(int(R.sum()), len(tail)),
                                     "retained_call_fraction_of_whole_forward_tail": fraction(int((F & R).sum()), len(tail)),
                                     "tail_selection": "eligible; frozen forward-only model/class threshold; includes RC collapse"})
            rows["call"].append({"stratum": stratum, "model_type": model,
                                 **call_summary(frame[f"{model}_forward_call"], frame[f"{model}_rc_call"], frame[f"{model}_forward_well_separated"])})
            sb, db = frame[f"{model}_forward_score_borderline"], frame[f"{model}_forward_delta_borderline"]
            for score_borderline in (False, True):
                for delta_borderline in (False, True):
                    block = frame.loc[(sb == score_borderline) & (db == delta_borderline)]
                    score_cross = block[f"{model}_forward_region_gate"] != block[f"{model}_rc_region_gate"]
                    delta_cross = block[f"{model}_forward_delta_gate"] != block[f"{model}_rc_delta_gate"]
                    row = {"stratum": stratum, "model_type": model, "forward_score_borderline": score_borderline,
                           "forward_delta_borderline": delta_borderline,
                           **call_summary(block[f"{model}_forward_call"], block[f"{model}_rc_call"], block[f"{model}_forward_well_separated"]),
                           "eligible_n": int(block.causal_call_eligible.sum()),
                           "region_gate_crossing_n": int(score_cross.sum()), "region_gate_crossing_fraction": fraction(int(score_cross.sum()), len(block)),
                           "delta_gate_crossing_n": int(delta_cross.sum()), "delta_gate_crossing_fraction": fraction(int(delta_cross.sum()), len(block)),
                           "call_change_n": int(block[f"{model}_call_changed"].sum()),
                           "call_change_fraction": fraction(int(block[f"{model}_call_changed"].sum()), len(block)),
                           **distribution(block[f"{model}_region_score_absolute_disagreement"], "region_abs_difference_")}
                    rows["proximity"].append(row)
            T = EXPECTED_THRESHOLDS[(model, "SNV")][0]
            boundaries = [0., T - .10, T - .02, T + .02, T + .10, 1.]
            for i, (lower, upper) in enumerate(zip(boundaries[:-1], boundaries[1:])):
                selected = (frame[f"{model}_forward_region_score"] >= lower) & ((frame[f"{model}_forward_region_score"] <= upper) if i == 4 else (frame[f"{model}_forward_region_score"] < upper))
                block = frame.loc[selected]
                rows["bins"].append({"stratum": stratum, "model_type": model, "forward_region_bin": i + 1,
                                     "lower_inclusive": lower, "upper": upper, "upper_inclusive": i == 4,
                                     **score_metrics(block[f"{model}_forward_region_score"], block[f"{model}_rc_region_score"]),
                                     "changed_call_n": int(block[f"{model}_call_changed"].sum()),
                                     "changed_call_fraction": fraction(int(block[f"{model}_call_changed"].sum()), len(block))})
            for direction, selected, opposite_orientation in (
                    ("lost_forward_call", frame[f"{model}_lost_forward_call"], "rc"),
                    ("gained_rc_only_call", frame[f"{model}_gained_rc_only_call"], "forward")):
                for failure in ("region_only", "delta_only", "both"):
                    region_failed, delta_failed = ~frame[f"{model}_{opposite_orientation}_region_gate"], ~frame[f"{model}_{opposite_orientation}_delta_gate"]
                    failure_mask = region_failed & ~delta_failed if failure == "region_only" else ~region_failed & delta_failed if failure == "delta_only" else region_failed & delta_failed
                    rows["gate"].append({"stratum": stratum, "model_type": model, "transition": direction,
                                         "failed_orientation": opposite_orientation, "failed_gate": failure,
                                         "n": int((selected & failure_mask).sum()), "transition_n": int(selected.sum()),
                                         "fraction_of_transition": fraction(int((selected & failure_mask).sum()), int(selected.sum()))})
        rows["call"].append({"stratum": stratum, "model_type": "union", **call_summary(frame.forward_union_call, frame.rc_union_call, frame.forward_union_well_separated)})
        for fcontext in CONTEXTS:
            for rcontext in CONTEXTS:
                n = int((frame.forward_model_context.eq(fcontext) & frame.rc_model_context.eq(rcontext)).sum())
                forward_n = int(frame.forward_model_context.eq(fcontext).sum())
                rows["context"].append({"stratum": stratum, "forward_context": fcontext, "rc_context": rcontext,
                                        "n": n, "forward_context_n": forward_n, "fraction_within_forward_context": fraction(n, forward_n), "stratum_n": len(frame)})
    return {key: pd.DataFrame(value) for key, value in rows.items()}


def severity_table(tables):
    result = []
    scores = tables["score"].loc[tables["score"].stratum.eq("all_15303")]
    deltas = tables["delta"].loc[tables["delta"].stratum.eq("all_15303")].set_index("model_type")
    tails = tables["tail"].loc[tables["tail"].stratum.eq("all_15303") & tables["tail"].forward_tail.eq("matched_class_p95")].set_index("model_type")
    calls = tables["call"].loc[tables["call"].stratum.eq("all_15303")].set_index("model_type")
    supports = tables["support"].set_index("stratum")
    for model in MODELS:
        selected = scores.loc[scores.model_type.eq(model)]
        material = (selected.abs_difference_mean >= .05).any() or (selected.abs_difference_p95 >= .10).any()
        negligible = (selected.abs_difference_mean <= .005).all() and (selected.abs_difference_p95 <= .02).all()
        result.append({"domain": "scores", "model_type": model, "severity": severity_label(material, negligible),
                       "material_condition_met": bool(material), "all_negligible_conditions_met": bool(negligible),
                       "metrics": json.dumps({r.score_type: {"MAE": r.abs_difference_mean, "p95": r.abs_difference_p95} for r in selected.itertuples()}, sort_keys=True),
                       "rule": "material: any REF/ALT/region MAE>=0.05 or p95>=0.10; negligible: all MAE<=0.005 and p95<=0.02"})
        d, t = deltas.loc[model], tails.loc[model]
        retention = t.meaningful_direction_retention
        material = d.normalized_delta_disagreement_median >= .50 or d.normalized_delta_disagreement_p95 >= 1 or retention < .90
        negligible = d.normalized_delta_disagreement_median <= .05 and d.normalized_delta_disagreement_p95 <= .25 and retention >= .99
        result.append({"domain": "allele_deltas", "model_type": model, "severity": severity_label(material, negligible),
                       "material_condition_met": bool(material), "all_negligible_conditions_met": bool(negligible),
                       "metrics": json.dumps({"median_q": d.normalized_delta_disagreement_median, "p95_q": d.normalized_delta_disagreement_p95,
                                               "forward_p95_tail_meaningful_direction_retention": retention,
                                               "tail_denominator": int(t.forward_meaningful_n), "tail_direction_retained_n": int(t.meaningful_direction_retained_n)}, sort_keys=True),
                       "rule": "material: median q>=0.50 OR p95 q>=1 OR tail meaningful retention<0.90; negligible: median q<=0.05 AND p95 q<=0.25 AND tail retention>=0.99"})
    for model in (*MODELS, "union"):
        c = calls.loc[model]
        well_loss_material = c.well_separated_lost_n >= 3 and c.well_separated_retention < .95
        material = c.jaccard < .80 or c.forward_retention < .90 or well_loss_material
        negligible = c.jaccard >= .98 and c.forward_retention >= .99 and c.well_separated_lost_n == 0
        result.append({"domain": "candidate_calls", "model_type": model, "severity": severity_label(material, negligible),
                       "material_condition_met": bool(material), "all_negligible_conditions_met": bool(negligible),
                       "metrics": json.dumps({key: float(c[key]) for key in ("jaccard", "forward_retention", "well_separated_retention", "well_separated_forward_positive_n", "well_separated_lost_n")}, sort_keys=True),
                       "rule": "material: Jaccard<0.80 OR retention<0.90 OR well-separated retention<0.95 with >=3 losses; negligible: Jaccard>=0.98 AND retention>=0.99 AND no well-separated loss"})
    s = supports.loc["frozen_337"]
    material = s.lost_union_call_fraction > .05 or s.lost_any_original_model_call_fraction > .10
    negligible = s.exact_model_context_stable_fraction >= .99 and s.lost_any_well_separated_model_call_n == 0
    result.append({"domain": "frozen_337_support_vectors", "model_type": "both_models", "severity": severity_label(material, negligible),
                   "material_condition_met": bool(material), "all_negligible_conditions_met": bool(negligible),
                   "metrics": json.dumps({key: float(s[key]) for key in ("exact_model_context_stable_fraction", "lost_union_call_n", "lost_union_call_fraction", "lost_any_original_model_call_n", "lost_any_original_model_call_fraction", "lost_any_well_separated_model_call_n")}, sort_keys=True),
                   "rule": "material: >5% lose union OR >10% lose any original model; negligible: >=99% exact support vector retained AND no well-separated model loss"})
    severity_order = {"negligible": 0, "modest": 1, "scientifically_material": 2}
    worst = max((r["severity"] for r in result), key=severity_order.get)
    result.append({"domain": "overall", "model_type": "both_models_and_union", "severity": worst,
                   "material_condition_met": worst == "scientifically_material", "all_negligible_conditions_met": worst == "negligible",
                   "metrics": json.dumps({"changed_decisions_n": int(supports.loc["all_15303", "changed_model_decisions_n"]),
                                           "forward_borderline_changed_decisions_n": int(supports.loc["all_15303", "changed_model_decisions_forward_borderline_n"]),
                                           "forward_borderline_fraction": float(supports.loc["all_15303", "changed_model_decisions_forward_borderline_fraction"])}, sort_keys=True),
                   "rule": "worst prespecified domain; modest predominantly borderline only when >=90% of changed model decisions fall in forward proximity band"})
    return pd.DataFrame(result)


def literature_tables(out, literature):
    coverage, exact = [], []
    for row in literature.to_dict("records"):
        ids = tokens(row["exact_candidate_record_ids"])
        evaluated = [candidate for candidate in ids if candidate in out.index]
        row.update({"v2_result_id": "COPD-V2-RC-R007", "rc_evaluable_exact_n": len(evaluated),
                    "rc_evaluable_exact_ids": ";".join(evaluated),
                    "rc_coverage_status": "evaluable_exact_frozen_identity" if evaluated else "not_evaluable_no_exact_frozen_identity",
                    "scope": "orientation behavior of already cataloged V1 identities; no new functional benchmark or LD substitution"})
        coverage.append(row)
        for candidate in evaluated:
            values = out.loc[candidate].to_dict()
            values.update({"variant_evidence_id": row["variant_evidence_id"], "source_id": row["source_id"],
                           "tested_rsids": row["tested_rsids"], "conclusion_tier": row["conclusion_tier"],
                           "endogenous_allele_tested": row["endogenous_allele_tested"], "literature_caveat": row["literature_caveat"]})
            exact.append(values)
    assert len(coverage) == 8 and len(exact) == 2
    return pd.DataFrame(coverage), pd.DataFrame(exact)


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean_json(v) for v in value]
    if value is pd.NA or value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_table(table, filename):
    path = RESULTS / filename
    kwargs = {"sep": "\t", "index": False, "na_rep": "not_evaluable", "float_format": "%.17g"}
    if filename.endswith(".gz"):
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                import io
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                    table.to_csv(text, **kwargs)
    else:
        table.to_csv(path, **kwargs)
    return path


def self_test():
    # Inclusive threshold equality, independent blacklist eligibility, class cutoff.
    assert positive_call([True, False, True], np.array([.25, .25, .25]), np.array([.5, .5, .5]), .5, np.array([.25, .25, .3])).tolist() == [True, False, False]
    assert signs([.1, -.1, .05, -.05, 0], .05).tolist() == [1, -1, 0, 0, 0]
    m = sign_metrics([.2, -.2, .2], [.3, .2, 0], [.1, .1, .1])
    assert m["forward_meaningful_n"] == 3 and m["meaningful_direction_retained_n"] == 1
    assert m["meaningful_near_zero_collapse_n"] == 1 and m["meaningful_opposite_direction_n"] == 1
    assert m["meaningful_direction_retention"] == 1/3 and m["raw_nonzero_sign_agreement"] == .5
    d = delta_metrics([.2, -.2], [-.2, .2], [.1, .2])
    assert d["normalized_delta_disagreement_median"] == 3.
    assert d["absolute_effect_magnitude_difference_maximum"] == 0.
    assert correlations([1, 1], [2, 3])["correlation_status"] == "not_evaluable_constant_vector"
    assert correlations([1], [1])["pearson_r"] is None
    assert abs(correlations([1, 2, 2, 3], [4, 2, 2, 1])["spearman_rho"] + 1) < 1e-12
    c = call_summary([1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 0, 0])
    assert [c[k] for k in ("both_positive_n", "forward_only_n", "reverse_only_n", "neither_n")] == [1, 1, 1, 1]
    assert c["jaccard"] == 1/3 and c["forward_retention"] == .5 and c["well_separated_lost_n"] == 1
    assert call_summary([], [], [])["jaccard"] is None
    assert distribution([0, 1, 2, 3], "x_")["x_p95"] == 2.8499999999999996
    assert severity_label(True, True) == "scientifically_material"
    print(json.dumps({"self_tests": "passed", "tests": "inclusive_fixed_calls; blacklist; class_cutoff; meaningful_zero_boundaries; collapse_denominator; normalized_delta_vs_amplitude; constant_correlation; tied_Spearman; 2x2; empty_denominator; linear_quantile; severity_precedence"}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="Only test synthetic metric/rule behavior; do not read or produce RC results")
    parser.add_argument("--sequence-qc", type=Path, default=INPUT_PATHS["sequence_qc"])
    args = parser.parse_args()
    assert sha256(SPEC) == SPEC_SHA, "locked specification hash mismatch"
    if args.self_test:
        self_test()
        return
    start = datetime.now(timezone.utc).isoformat()
    paths = dict(INPUT_PATHS)
    paths["sequence_qc"] = args.sequence_qc
    for path in paths.values():
        assert path.is_file(), f"Required input missing: {path}"
    inputs_before = {name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "size_bytes": path.stat().st_size} for name, path in paths.items()}
    preflight = json.loads(paths["preflight"].read_text())
    assert preflight["spec_sha256"] == SPEC_SHA
    assert str(preflight.get("gate_status", "")).lower() in ("pass", "passed", "complete", "passed_before_rc"), preflight.get("gate_status")
    tables = input_frames(paths)
    frame = build_comparisons(tables)
    metrics = metric_tables(frame)
    severity = severity_table(metrics)
    coverage, exact = literature_tables(frame, tables["literature"])
    frozen = frame.loc[tables["frozen337"].candidate_record_id].copy()
    shortlist = frame.loc[tables["shortlist"].candidate_record_id].copy()
    assert frozen.v1_priority_rank.tolist() == list(range(1, 338))
    assert shortlist.v1_shortlist_rank.tolist() == list(range(1, 13))
    output_specs = [
        (frame, "R004_all_scorable_comparison.tsv.gz"),
        (frozen, "R005_frozen_337_comparison.tsv"),
        (shortlist, "R006_shortlist_12_comparison.tsv"),
        (coverage, "R007_literature_coverage.tsv"),
        (exact, "R008_exact_literature_comparison.tsv"),
        (metrics["score"], "R009_score_metrics.tsv"),
        (metrics["delta"], "R010_delta_metrics.tsv"),
        (metrics["sign"], "R011_sign_transitions.tsv"),
        (metrics["tail"], "R012_forward_tail_metrics.tsv"),
        (metrics["call"], "R013_call_transitions.tsv"),
        (metrics["context"], "R014_model_context_transitions.tsv"),
        (metrics["proximity"], "R015_threshold_proximity.tsv"),
        (metrics["bins"], "R016_forward_region_bins.tsv"),
        (metrics["support"], "R017_subset_support_summary.tsv"),
        (severity, "R018_severity.tsv"),
        (metrics["gate"], "R019_gate_failure_transitions.tsv"),
    ]
    checks = []
    def check(name, condition, detail):
        assert condition, (name, detail)
        checks.append({"check_id": f"COPD-V2-RC-AQC-{len(checks)+1:03d}", "check": name, "status": "PASS", "detail": detail})
    check("locked_specification", sha256(SPEC) == SPEC_SHA, SPEC_SHA)
    check("complete_scorable_record_identity", len(frame) == 15303 and frame.index.is_unique, "15303 exact paired IDs; no duplicate or excluded score")
    check("frozen_call_eligibility", int(frame.causal_call_eligible.sum()) == 15283 and int(frame.blacklisted.sum()) == 20, "15283 eligible and 20 scorable blacklisted")
    check("forward_175_199_37_337", (int(frame.enhancer_forward_call.sum()),int(frame.silencer_forward_call.sum()),int((frame.enhancer_forward_call&frame.silencer_forward_call).sum()),int(frame.forward_union_call.sum())) == (175,199,37,337), "inclusive original thresholds and eligibility; authoritative REF/ALT subtraction")
    check("frozen_337_order", frozen.candidate_record_id.tolist() == tables["frozen337"].candidate_record_id.tolist(), "original R010 identity and 1..337 rank")
    check("shortlist_12_order", shortlist.candidate_record_id.tolist() == tables["shortlist"].candidate_record_id.tolist(), "original S6 identity and 1..12 shortlist rank")
    check("phenotype_subsets", (int(frame.primary_retained.sum()), int(frame.GCST90244098_only_study.sum())) == (184,124), "frozen phenotype flags used descriptively")
    check("literature_coverage", len(coverage) == 8 and len(exact) == 2, "all8 prior register rows; only2 exact extant IDs; no LD substitution")
    check("blacklisted_calls_stay_negative", not frame.loc[frame.blacklisted,["enhancer_forward_call","silencer_forward_call","enhancer_rc_call","silencer_rc_call"]].any().any(), "all20 included in score metrics and excluded from calls")
    for stratum, mask in strata(frame).items():
        n = int(np.asarray(mask).sum())
        c = metrics["call"].loc[metrics["call"].stratum.eq(stratum)]
        check("call_transition_reconciliation_" + stratum, (c[["transition_0_to_0","transition_0_to_1","transition_1_to_0","transition_1_to_1"]].sum(axis=1) == n).all(), f"each model and union 2x2 sums to {n}")
        ctx = metrics["context"].loc[metrics["context"].stratum.eq(stratum)]
        check("context_reconciliation_" + stratum, int(ctx.n.sum()) == n, f"4x4 sums to {n}")
        sign = metrics["sign"].loc[metrics["sign"].stratum.eq(stratum)]
        check("sign_reconciliation_" + stratum, (sign.groupby(["model_type","sign_mode"]).n.sum() == n).all(), f"each raw/meaningful model 3x3 sums to {n}")
        gates = metrics["gate"].loc[metrics["gate"].stratum.eq(stratum)]
        check("gate_loss_reconciliation_" + stratum, all(int(g.n.sum()) == int(g.transition_n.iloc[0]) for _,g in gates.groupby(["model_type","transition"])), "region-only/delta-only/both failures exactly exhaust each lost/gained call")
    check("inputs_unchanged_during_analysis", all(sha256(path) == inputs_before[name]["sha256"] for name,path in paths.items()), "all consumed input hashes rechecked; broader V1 freeze checked by final validator")
    output_specs.append((pd.DataFrame(checks), "R020_analysis_validation.tsv"))
    written = []
    for table, name in output_specs:
        path = write_table(table, "COPD-V2-RC-" + name)
        written.append({"result_id": "COPD-V2-RC-" + name.split("_")[0], "path": str(path.relative_to(V2)),
                        "rows": len(table), "columns": len(table.columns), "sha256": sha256(path), "size_bytes": path.stat().st_size})
    summary = {
        "analysis_id": "COPD-V2-RC-ANALYSIS-001", "started_utc": start, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "spec_sha256": SPEC_SHA, "overall_severity": severity.loc[severity.domain.eq("overall"),"severity"].iloc[0],
        "complete_universe_score_metrics": metrics["score"].loc[metrics["score"].stratum.eq("all_15303")].to_dict("records"),
        "complete_universe_delta_metrics": metrics["delta"].loc[metrics["delta"].stratum.eq("all_15303")].to_dict("records"),
        "complete_universe_calls": metrics["call"].loc[metrics["call"].stratum.eq("all_15303")].to_dict("records"),
        "all_subset_support": metrics["support"].to_dict("records"), "all_subset_calls": metrics["call"].to_dict("records"),
        "complete_universe_tail_metrics": metrics["tail"].loc[metrics["tail"].stratum.eq("all_15303")].to_dict("records"),
        "class_delta_metrics": metrics["delta"].loc[metrics["delta"].stratum.isin(CLASSES)].to_dict("records"),
        "class_score_metrics": metrics["score"].loc[metrics["score"].stratum.isin(CLASSES)].to_dict("records"),
        "severity": severity.to_dict("records"), "exact_literature": exact.to_dict("records"),
        "shortlist": shortlist.to_dict("records"), "forward_serialized_delta_max_abs_rounding": {
            model: float(frame[f"{model}_forward_serialization_delta_rounding"].abs().max()) for model in MODELS},
        "ambiguous_N_pairs": int(frame.has_ambiguous_N.sum()), "analysis_qc_passed": len(checks),
        "interpretation": "diagnostic only; authoritative forward scores and all V1 ranks unchanged; no replacement candidate set"}
    summary_path = V2 / "logs/COPD-V2-RC_analysis_summary.json"
    summary_path.write_text(json.dumps(clean_json(summary), indent=2, sort_keys=True, allow_nan=False) + "\n")
    manifest = {"analysis_id": "COPD-V2-RC-ANALYSIS-001", "started_utc": start, "completed_utc": datetime.now(timezone.utc).isoformat(),
                "spec_sha256": SPEC_SHA, "script": {"path": str(SCRIPT.relative_to(ROOT)), "sha256": sha256(SCRIPT)},
                "software": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
                "command": sys.argv, "inputs": inputs_before, "outputs": written,
                "summary": {"path": str(summary_path.relative_to(V2)), "sha256": sha256(summary_path)},
                "analysis_qc_passed": len(checks), "input_hashes_unchanged": True}
    (PROVENANCE / "COPD-V2-RC_analysis_manifest.json").write_text(json.dumps(clean_json(manifest), indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"analysis_qc_passed": len(checks), "overall_severity": summary["overall_severity"], "outputs_written": len(written), "summary": str(summary_path)}))


if __name__ == "__main__":
    main()

