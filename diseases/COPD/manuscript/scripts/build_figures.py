#!/usr/bin/env python3
"""Build COPD manuscript figures and their plotted source-data tables.

All values are read from frozen primary result tables. The script does not
recompute upstream analyses or alter their ranking.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


SCRIPT = Path(__file__).resolve()
REPO = SCRIPT.parents[4]
COPD = REPO / "diseases" / "COPD"
OUT = SCRIPT.parents[1] / "figures"
DATA = SCRIPT.parents[1] / "figure_data"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

COLORS = {
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "green": "#009E73",
    "orange": "#E69F00",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "gray": "#6B7280",
    "lightgray": "#D1D5DB",
    "dark": "#1F2937",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
    }
)


def read(rel: str) -> pd.DataFrame:
    return pd.read_csv(REPO / rel, sep="\t", low_memory=False)


def bool_col(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.07,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)


def figure_1() -> None:
    ld_path = "diseases/COPD/02_gwas/results/COPD-S2-R006F_summary.tsv"
    consequence_path = (
        "diseases/COPD/02_gwas/results/COPD-S2-R004_coding_noncoding_summary.tsv"
    )
    nonexclusive_path = (
        "diseases/COPD/03_regulatory_landscape/results/"
        "COPD-S3-R003_nonexclusive_summary.tsv"
    )
    exclusive_path = (
        "diseases/COPD/03_regulatory_landscape/results/"
        "COPD-S3-R003_exclusive_summary.tsv"
    )
    ld = read(ld_path).set_index("metric")
    cons = read(consequence_path).query("subset == 'gws_unique_tag_variants'").iloc[0]
    nonex_all = read(nonexclusive_path)
    nonex = nonex_all.query(
        "stratum_type == 'overall' and stratum == 'all_candidate_records'"
    ).set_index("feature")
    excl = read(exclusive_path)
    excl = excl.query(
        "stratum_type == 'overall' and stratum == 'all_candidate_records' "
        "and hierarchy == 'comprehensive'"
    ).sort_values("precedence_rank")

    gws_tag_records = int(
        nonex_all.query(
            "stratum_type == 'role' and stratum == 'gws_tag_record' "
            "and feature == 'coding_CDS'"
        ).iloc[0].n_total_records
    )
    role_overlap = int(
        nonex_all.query(
            "stratum_type == 'candidate_origin' and stratum == 'gws_tag_and_ld_proxy' "
            "and feature == 'coding_CDS'"
        ).iloc[0].n_total_records
    )
    union_metrics = [
        ("Genome-wide significant tag identifiers", int(ld.loc["gws_tag_variants", "value"]), ld_path),
        ("Reference-matched tag identifiers", int(ld.loc["strictly_matched_tag_variants", "value"]), ld_path),
        ("GWS tag records retained", gws_tag_records, nonexclusive_path),
        ("Unique LD proxy records", int(ld.loc["unique_ld_proxy_reference_records", "value"]), ld_path),
        ("Records with both roles", role_overlap, nonexclusive_path),
        ("Frozen record union", int(ld.loc["expanded_candidate_records", "value"]), ld_path),
    ]
    consequence = [
        ("Coding", int(cons.n_coding)),
        ("Splice", int(cons.n_splice)),
        ("Noncoding", int(cons.n_noncoding)),
        ("Unknown", int(cons.n_unknown)),
    ]
    annotation_features = [
        ("Coding CDS", "coding_CDS"),
        ("Refined donor regulatory", "donor_refined_any"),
        ("Reference regulatory", "known_regulatory_any"),
        ("RepeatMasker", "repeatmasker"),
        ("No tested annotation", "unannotated_across_tested_biological_features"),
    ]
    annotations = [
        (label, int(nonex.loc[key, "n_overlapping"]), float(nonex.loc[key, "percent_of_all_records"]))
        for label, key in annotation_features
    ]

    records = []
    for label, n, source_path in union_metrics:
        records.append(
            {
                "panel": "A",
                "metric": label,
                "n": n,
                "denominator": "",
                "percent": "",
                "source_file": source_path,
            }
        )
    for label, n in consequence:
        records.append(
            {
                "panel": "B",
                "metric": label,
                "n": n,
                "denominator": 660,
                "percent": n / 660 * 100,
                "source_file": consequence_path,
            }
        )
    for label, n, pct in annotations:
        records.append(
            {
                "panel": "C",
                "metric": label,
                "n": n,
                "denominator": 15389,
                "percent": pct,
                "source_file": nonexclusive_path,
            }
        )
    for row in excl.itertuples():
        records.append(
            {
                "panel": "D",
                "metric": row.exclusive_class,
                "n": int(row.n_records),
                "denominator": int(row.n_total_records),
                "percent": float(row.percent_of_all_records),
                "source_file": exclusive_path,
            }
        )
    pd.DataFrame(records).to_csv(DATA / "figure_1_source_data.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.4))
    ax = axes[0, 0]
    ax.axis("off")
    boxes = [
        (0.01, 0.50, 0.35, 0.31, COLORS["blue"], "660 GWS tags", "575 matched\nfor LD calculation"),
        (0.01, 0.10, 0.35, 0.31, COLORS["orange"], "14,913 proxies", "r² ≥ 0.8\nwithin ±500 kb"),
        (0.57, 0.30, 0.41, 0.35, COLORS["green"], "15,389 union records", "652 tag records + 14,913 proxies\n- 176 shared identities"),
    ]
    for bx, by, bw, bh, color, main, sub in boxes:
        rect = mpl.patches.FancyBboxPatch((bx, by), bw, bh, boxstyle="round,pad=0.015,rounding_size=0.015", transform=ax.transAxes, facecolor=color, edgecolor="white")
        ax.add_patch(rect)
        ax.text(bx + bw / 2, by + bh * 0.65, main, color="white", ha="center", va="center", fontsize=8.4, fontweight="bold", transform=ax.transAxes)
        ax.text(bx + bw / 2, by + bh * 0.28, sub, color="white", ha="center", va="center", fontsize=6.4, linespacing=1.05, transform=ax.transAxes)
    for y0 in [0.655, 0.255]:
        ax.annotate("", xy=(0.56, 0.475), xytext=(0.37, y0), xycoords=ax.transAxes, arrowprops={"arrowstyle": "->", "color": COLORS["gray"], "lw": 1.6})
    ax.set_title("Tag records and ancestry-matched LD proxies form a union", loc="left")
    ax.text(0.01, -0.005, "All tags retained; tag and proxy identities can overlap.\nLD expansion is not statistical fine-mapping.", transform=ax.transAxes, fontsize=6.8, color=COLORS["dark"], va="bottom")
    panel_label(ax, "A")

    ax = axes[0, 1]
    colors = [COLORS["vermillion"], COLORS["purple"], COLORS["blue"], COLORS["lightgray"]]
    left = 0
    for (label, n), color in zip(consequence, colors):
        ax.barh([0], [n], left=left, color=color, height=0.42, label=f"{label}: {n}")
        if n > 20:
            ax.text(left + n / 2, 0, str(n), ha="center", va="center", color="white" if label != "Unknown" else COLORS["dark"], fontsize=8, fontweight="bold")
        left += n
    ax.set_xlim(0, 660)
    ax.set_yticks([])
    ax.set_xlabel("Genome-wide significant tag variants")
    ax.set_title("Catalog consequence among 660 significant tags", loc="left")
    ax.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.47), fontsize=8)
    panel_label(ax, "B")

    ax = axes[1, 0]
    labels = [x[0] for x in annotations][::-1]
    pcts = [x[2] for x in annotations][::-1]
    ns = [x[1] for x in annotations][::-1]
    bars = ax.barh(labels, pcts, color=[COLORS["gray"], COLORS["orange"], COLORS["blue"], COLORS["green"], COLORS["vermillion"]])
    for bar, pct, n in zip(bars, pcts, ns):
        ax.text(pct + 1, bar.get_y() + bar.get_height() / 2, f"{pct:.1f}% ({n:,})", va="center", fontsize=8)
    ax.set_xlim(0, 62)
    ax.set_xlabel("Candidate records (%)")
    ax.set_title("Nonexclusive annotation of 15,389 candidates", loc="left")
    panel_label(ax, "C")

    ax = axes[1, 1]
    class_label = {
        "bed_ineligible": "BED ineligible",
        "coding_CDS": "Coding CDS",
        "refined_enhancer": "Refined enhancer",
        "refined_silencer": "Refined silencer",
        "preliminary_enhancer": "Preliminary enhancer",
        "preliminary_silencer": "Preliminary silencer",
        "other_known_regulatory": "Other reference regulatory",
        "repetitive_only_in_hierarchy": "Repeat only",
        "other": "Other",
    }
    labels = [class_label.get(x, x) for x in excl.exclusive_class]
    vals = excl.n_records.astype(int).to_numpy()
    palette = [COLORS["lightgray"], COLORS["vermillion"], COLORS["green"], COLORS["purple"], COLORS["sky"], COLORS["yellow"], COLORS["blue"], COLORS["orange"], COLORS["gray"]]
    y = np.arange(len(labels))
    ax.barh(y, vals, color=palette)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Candidate records")
    ax.set_title("Mutually exclusive comprehensive classification", loc="left")
    for yi, value in zip(y, vals):
        ax.text(value + 100, yi, f"{value:,}", va="center", fontsize=7)
    ax.set_xlim(0, max(vals) * 1.17)
    panel_label(ax, "D")

    fig.suptitle("Figure 1. From COPD GWAS signals to a frozen regulatory candidate set", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout(h_pad=2.1, w_pad=2.6)
    save(fig, "figure_1_gwas_regulatory_landscape")


def figure_2() -> None:
    performance_path = "diseases/COPD/04_modeling/results/COPD-S4-R001_model_performance.tsv"
    variant_path = "diseases/COPD/04_modeling/results/COPD-S4-R004_variant_summary.tsv"
    locus_path = "diseases/COPD/04_modeling/results/COPD-S4-R004_locus_summary.tsv"
    validation_path = "diseases/COPD/05_computational_validation/results/COPD-S5-R005_integrated_summary.tsv"
    perf = read(performance_path)
    var = read(variant_path).set_index("metric")
    loci = read(locus_path).query("analysis_set == 'replicated_gwas_genes'").set_index("metric")
    validation = read(validation_path).set_index("metric")

    model_calls = [
        ("Enhancer only", int(var.loc["predicted_causal_enhancer", "n"] - var.loc["predicted_causal_both_models", "n"])),
        ("Both models", int(var.loc["predicted_causal_both_models", "n"])),
        ("Silencer only", int(var.loc["predicted_causal_silencer", "n"] - var.loc["predicted_causal_both_models", "n"])),
    ]
    locus_classes = [
        ("Coding + regulatory", int(loci.loc["coding_plus_predicted_causal_regulatory", "n"])),
        ("Coding only", int(loci.loc["coding_only", "n"])),
        ("Regulatory only", int(loci.loc["predicted_causal_regulatory_only", "n"])),
        ("Other", int(loci.loc["other", "n"])),
    ]
    validation_metrics = [
        ("Exact GTEx Lung eQTL", int(validation.loc["exact_significant_GTEx_v10_Lung_cis_eqtl", "n"]), 337, "variant-level association"),
        ("MPRAbase element", int(validation.loc["any_MPRAbase_element_evidence", "n"]), 337, "regional coverage"),
        ("No primary-query evidence", int(validation.loc["no_variant_or_regional_evidence_in_primary_queries", "n"]), 337, "not a negative causal call"),
        ("Open Targets gene context", int(validation.loc["candidate_with_OpenTargets_COPD_gene_context", "n"]), 337, "gene-level context only"),
    ]

    records = []
    for row in perf.itertuples():
        records.extend(
            [
                {"panel": "A", "metric": "ROC AUC", "group": row.model_type, "n": row.roc_auc, "denominator": row.n_test, "percent": "", "note": f"95% CI {row.roc_auc_ci95_low:.6f}-{row.roc_auc_ci95_high:.6f}", "source_file": performance_path},
                {"panel": "A", "metric": "PR AUC", "group": row.model_type, "n": row.pr_auc_trapezoid, "denominator": row.n_test, "percent": "", "note": f"95% CI {row.pr_auc_ci95_low:.6f}-{row.pr_auc_ci95_high:.6f}", "source_file": performance_path},
            ]
        )
    for label, n in model_calls:
        records.append({"panel": "B", "metric": label, "group": "model-nominated candidate", "n": n, "denominator": 337, "percent": n / 337 * 100, "note": "", "source_file": variant_path})
    for label, n in locus_classes:
        records.append({"panel": "C", "metric": label, "group": "replicated GWAS genes", "n": n, "denominator": 140, "percent": n / 140 * 100, "note": "", "source_file": locus_path})
    for label, n, denominator, note in validation_metrics:
        records.append({"panel": "D", "metric": label, "group": "public validation", "n": n, "denominator": denominator, "percent": n / denominator * 100, "note": note, "source_file": validation_path})
    pd.DataFrame(records).to_csv(DATA / "figure_2_source_data.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.5))
    ax = axes[0, 0]
    point_x = []
    point_labels = []
    for metric_i, metric_name in enumerate(["ROC AUC", "PR AUC"]):
        for model_i, row in enumerate(perf.itertuples()):
            xval = metric_i * 1.25 + model_i * 0.34
            if metric_name == "ROC AUC":
                value, low, high = row.roc_auc, row.roc_auc_ci95_low, row.roc_auc_ci95_high
            else:
                value, low, high = row.pr_auc_trapezoid, row.pr_auc_ci95_low, row.pr_auc_ci95_high
            color = [COLORS["blue"], COLORS["purple"]][model_i]
            ax.errorbar(xval, value, yerr=[[value - low], [high - value]], fmt="o", color=color, capsize=4, markersize=7, lw=1.5)
            ax.text(xval, value + 0.009, f"{value:.3f}", ha="center", fontsize=7.5)
            point_x.append(xval)
            model_short = "Enh." if row.model_type == "enhancer" else "Sil."
            metric_short = "ROC" if metric_name == "ROC AUC" else "PR"
            point_labels.append(f"{model_short}\n{metric_short}")
    ax.set_xticks(point_x, point_labels)
    ax.set_xlim(-0.3, 1.9)
    ax.set_ylim(0.86, 0.985)
    ax.set_ylabel("Held-out performance")
    ax.set_title("Chromosome-held-out discrimination (95% CI)", loc="left")
    panel_label(ax, "A")

    ax = axes[0, 1]
    labels, vals = zip(*model_calls)
    ax.bar(labels, vals, color=[COLORS["green"], COLORS["orange"], COLORS["purple"]])
    for i, value in enumerate(vals):
        ax.text(i, value + 4, f"{value}\n({value/337*100:.1f}%)", ha="center", fontsize=8)
    ax.set_ylim(0, max(vals) * 1.24)
    ax.set_ylabel("Candidate records")
    ax.set_title("337 predicted regulatory candidates", loc="left")
    panel_label(ax, "B")

    ax = axes[1, 0]
    labels, vals = zip(*locus_classes)
    wedges, _ = ax.pie(vals, startangle=90, colors=[COLORS["vermillion"], COLORS["gray"], COLORS["green"], COLORS["lightgray"]], wedgeprops={"linewidth": 1, "edgecolor": "white"})
    ax.legend(wedges, [f"{l}: {v}" for l, v in locus_classes], loc="center left", bbox_to_anchor=(0.94, 0.5), fontsize=8)
    ax.set_title("Replicated GWAS-gene locus classes (n = 140)", loc="left")
    panel_label(ax, "C")

    ax = axes[1, 1]
    labels = [x[0] for x in validation_metrics][::-1]
    vals = [x[1] for x in validation_metrics][::-1]
    colors = [COLORS["sky"], COLORS["lightgray"], COLORS["orange"], COLORS["blue"]]
    bars = ax.barh(labels, vals, color=colors)
    for bar, value in zip(bars, vals):
        ax.text(value + 5, bar.get_y() + bar.get_height() / 2, f"{value} ({value/337*100:.1f}%)", va="center", fontsize=8)
    ax.set_xlim(0, 337)
    ax.set_xlabel("Predicted candidates")
    ax.set_title("Nonexclusive evidence summaries of different scope", loc="left")
    panel_label(ax, "D")

    fig.suptitle("Figure 2. Predictive performance, candidate nomination, and public evidence", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout(h_pad=2.0, w_pad=2.5)
    save(fig, "figure_2_modeling_validation")


def figure_3() -> None:
    shortlist_path = "diseases/COPD/06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv"
    constructs_path = "diseases/COPD/06_experimental_validation/results/COPD-S6-R003_MPRA_constructs.tsv"
    plans_path = "diseases/COPD/06_experimental_validation/results/COPD-S6-R004_candidate_validation_plan.tsv"
    tf_path = "diseases/COPD/06_experimental_validation/results/COPD-S6-R004_TF_first_plan.tsv"
    shortlist = read(shortlist_path).sort_values("experimental_shortlist_rank")
    constructs = read(constructs_path)
    plans = read(plans_path)
    tfplans = read(tf_path)

    enh = shortlist.enhancer_abs_delta_percentile_within_class.astype(float)
    sil = shortlist.silencer_abs_delta_percentile_within_class.astype(float)
    evidence = pd.DataFrame(
        {
            "Exact GTEx Lung eQTL": bool_col(shortlist.gtex_lung_exact_significant_eqtl),
            "Refined donor interval": bool_col(shortlist.donor_refined_any),
            "MPRAbase element": bool_col(shortlist.mprabase_any_element_evidence),
            "Motif created/disrupted": (
                bool_col(shortlist.enhancer_tf_disrupts_any_motif_compatible_site)
                | bool_col(shortlist.enhancer_tf_creates_any_motif_compatible_site)
                | bool_col(shortlist.silencer_tf_disrupts_any_motif_compatible_site)
                | bool_col(shortlist.silencer_tf_creates_any_motif_compatible_site)
            ),
            "Indel/complex": shortlist.variant_class.eq("indel_or_complex"),
        }
    )
    evidence.index = shortlist.candidate_record_id.to_numpy()

    source = shortlist[
        [
            "experimental_shortlist_rank",
            "predicted_causal_priority_rank",
            "candidate_record_id",
            "model_context",
            "variant_class",
            "enhancer_delta_alt_minus_ref",
            "silencer_delta_alt_minus_ref",
            "enhancer_abs_delta_percentile_within_class",
            "silencer_abs_delta_percentile_within_class",
            "donor_refined_any",
            "gtex_lung_exact_significant_eqtl",
            "mprabase_any_element_evidence",
            "enhancer_tf_disrupts_any_motif_compatible_site",
            "enhancer_tf_creates_any_motif_compatible_site",
            "silencer_tf_disrupts_any_motif_compatible_site",
            "silencer_tf_creates_any_motif_compatible_site",
            "section6_exact_GTEx_R006_target_genes",
            "nearest_protein_coding_gene",
        ]
    ].copy()
    source = source.rename(columns={"predicted_causal_priority_rank": "computational_priority_rank"})
    source["plotted_enhancer_matched_class_percentile"] = enh
    source["plotted_silencer_matched_class_percentile"] = sil
    source["source_file"] = shortlist_path
    source.to_csv(DATA / "figure_3_shortlist_source_data.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"design_item": "shortlisted variants", "n": len(shortlist), "source_file": shortlist_path},
            {"design_item": "MPRA inserts", "n": len(constructs), "source_file": constructs_path},
            {"design_item": "candidate validation plans", "n": len(plans), "source_file": plans_path},
            {"design_item": "TF-first plans", "n": len(tfplans), "source_file": tf_path},
        ]
    ).to_csv(DATA / "figure_3_design_counts.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(12.3, 9.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 0.8], width_ratios=[1.32, 1], hspace=0.56, wspace=0.32)
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(shortlist))
    height = 0.36
    ax.barh(y - height / 2, enh, height, color=COLORS["green"], label="Enhancer matched-class percentile")
    ax.barh(y + height / 2, sil, height, color=COLORS["purple"], label="Silencer matched-class percentile")
    ax.set_yticks(y, shortlist.candidate_record_id)
    ax.invert_yaxis()
    ax.set_xlim(0, 102)
    ax.set_xlabel("Within-model |ALT - REF| percentile in matched variant class")
    ax.set_title("Calibrated allelic effects in shortlist rank order", loc="left")
    ax.legend(loc="lower right", fontsize=8)
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    matrix = evidence.astype(int).to_numpy()
    ax.imshow(matrix, cmap=ListedColormap(["#F3F4F6", COLORS["blue"]]), vmin=0, vmax=1, aspect="auto")
    matrix_labels = ["GTEx\neQTL", "Donor\ninterval", "MPRAbase\nelement", "Motif\nchange", "Indel or\ncomplex"]
    ax.set_xticks(np.arange(evidence.shape[1]), matrix_labels, rotation=0, ha="center", fontsize=7.5)
    ax.set_yticks(np.arange(len(shortlist)), shortlist.experimental_shortlist_rank.astype(str))
    ax.set_ylabel("Experimental shortlist rank")
    ax.set_title("Evidence and design features", loc="left")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, "●" if matrix[i, j] else "", ha="center", va="center", color="white", fontsize=8)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    composition = shortlist.model_context.value_counts()
    order = ["both", "enhancer_only", "silencer_only"]
    labels = ["Both models", "Enhancer only", "Silencer only"]
    vals = [int(composition.get(k, 0)) for k in order]
    bars = ax.bar(labels, vals, color=[COLORS["orange"], COLORS["green"], COLORS["purple"]])
    for bar, value in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.18, str(value), ha="center", fontsize=9)
    ax.set_ylim(0, 10.5)
    ax.set_ylabel("Shortlisted candidates")
    ax.set_title("Model-context diversity", loc="left")
    panel_label(ax, "C")

    ax = fig.add_subplot(gs[1, 1])
    items = ["Shortlisted\nvariants", "MPRA\ninserts", "Candidate\nplans", "TF-first\nplans"]
    vals = [len(shortlist), len(constructs), len(plans), len(tfplans)]
    bars = ax.bar(items, vals, color=[COLORS["blue"], COLORS["sky"], COLORS["green"], COLORS["orange"]])
    for bar, value in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1, str(value), ha="center", fontsize=9)
    ax.set_ylim(0, 54)
    ax.set_ylabel("Records")
    ax.set_title("Proposed experimental package", loc="left")
    panel_label(ax, "D")

    fig.suptitle("Figure 3. Nonredundant candidate panel for prospective experimental testing", fontsize=14, fontweight="bold", y=0.99)
    fig.text(0.52, 0.018, "All plans are proposed and not performed; inserts are not ordering-ready oligonucleotides.", ha="center", fontsize=8)
    fig.subplots_adjust(bottom=0.10, top=0.91)
    save(fig, "figure_3_experimental_shortlist")


def main() -> None:
    figure_1()
    figure_2()
    figure_3()
    print(f"Wrote figures to {OUT}")
    print(f"Wrote plotted source data to {DATA}")


if __name__ == "__main__":
    main()
