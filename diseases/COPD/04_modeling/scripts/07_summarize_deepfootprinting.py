#!/usr/bin/env python3
"""Summarize COPD TREDNet DeepExplainer and TF-MoDISco runs.

The result is deliberately model-aware: enhancer and silencer feature SHAP
values, discovered patterns, and JASPAR matches are kept in one tidy schema
with an explicit ``model_type`` column.  TF names are resolved only from the
local JASPAR 2024 MEME file used by ``modisco report``.

This script interprets sequence patterns as computational hypotheses.  A
TF-MoDISco-to-JASPAR motif match is sequence similarity, not evidence that a
named TF binds the COPD locus in vivo.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases" / "COPD" / "04_modeling"
RESULTS = SECTION / "results"
DEEP_ROOT = SECTION / "deepexplainer" / "output" / "COPD_Section4"
FEATURES = SECTION / "trednet" / "model_phase_I" / "features_list.txt"
MEME = (
    SECTION
    / "deepexplainer"
    / "backup"
    / "JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme"
)

MODEL_IDS = {
    "enhancer": "COPD_SevereEmphysema_Lung_Enhancer_DHS_x2",
    "silencer": "COPD_SevereEmphysema_Lung_Silencer_DHS_x2",
}

OUT_FEATURES = RESULTS / "COPD-S4-R007_phase1_feature_shap.tsv.gz"
OUT_PATTERNS = RESULTS / "COPD-S4-R007_tfmodisco_patterns.tsv"
OUT_MOTIFS = RESULTS / "COPD-S4-R007_jaspar_motif_matches.tsv"
OUT_RUNS = RESULTS / "COPD-S4-R007_run_metadata.tsv"
OUT_FILES = RESULTS / "COPD-S4-R007_run_file_checksums.tsv"
OUT_SUMMARY = RESULTS / "COPD-S4-R007_summary.tsv"
OUT_REPORT = RESULTS / "COPD-S4-R007_interpretation.md"
OUT_MANIFEST = RESULTS / "COPD-S4-R007_analysis_manifest.json"

REQUIRED_RUN_FILES = (
    "selection_manifest.json",
    "submission_manifest.tsv",
    "shap_values_features_TREDNet.npy",
    "actual_shap.npy",
    "raw_shap_explanations.npy",
    "ohe_seqs_to_explain.npy",
    "tfmodisco_output/modisco_results.h5",
    "tfmodisco_output/motifs.html",
)


def sha256(path: Path, cache: dict[Path, str]) -> str:
    """Return a cached SHA-256 for a real file (symlinks are dereferenced)."""
    resolved = path.resolve()
    if resolved not in cache:
        digest = hashlib.sha256()
        with resolved.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        cache[resolved] = digest.hexdigest()
    return cache[resolved]


def file_record(
    model_type: str,
    run_dir: Path,
    role: str,
    path: Path,
    cache: dict[Path, str],
) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"missing {role}: {path}")
    return {
        "model_type": model_type,
        "run_directory": str(run_dir.resolve()),
        "file_role": role,
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path, cache),
    }


def read_key_value_tsv(path: Path) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(table.columns) != ["field", "value"]:
        raise ValueError(f"expected field/value table: {path}")
    if table["field"].duplicated().any():
        raise ValueError(f"duplicate fields in {path}")
    return dict(zip(table["field"], table["value"]))


def read_meme_catalog(path: Path) -> tuple[dict[str, dict[str, Any]], dict[str, float]]:
    """Parse motif names/widths and background frequencies from MEME text."""
    lines = path.read_text().splitlines()
    background: dict[str, float] = {}
    motifs: dict[str, dict[str, Any]] = {}
    for index, line in enumerate(lines):
        if line.startswith("Background letter frequencies") and index + 1 < len(lines):
            tokens = lines[index + 1].split()
            background = {
                tokens[i]: float(tokens[i + 1]) for i in range(0, len(tokens), 2)
            }
        if not line.startswith("MOTIF "):
            continue
        parts = line.split(maxsplit=2)
        motif_id = parts[1]
        tf_name = parts[2] if len(parts) == 3 else motif_id
        cursor = index + 1
        while cursor < len(lines) and "letter-probability matrix" not in lines[cursor]:
            if lines[cursor].startswith("MOTIF "):
                break
            cursor += 1
        if cursor >= len(lines) or "letter-probability matrix" not in lines[cursor]:
            raise ValueError(f"missing probability matrix for {motif_id}")
        width_match = re.search(r"\bw\s*=\s*(\d+)", lines[cursor])
        if not width_match:
            raise ValueError(f"missing width for {motif_id}")
        motifs[motif_id] = {"tf_name": tf_name, "motif_width": int(width_match.group(1))}
    if not motifs:
        raise ValueError(f"no motifs parsed from {path}")
    if set(background) != {"A", "C", "G", "T"}:
        raise ValueError(f"unexpected MEME background: {background}")
    return motifs, background


def complete_run(path: Path) -> bool:
    return path.is_dir() and all((path / item).is_file() for item in REQUIRED_RUN_FILES)


def select_run(model_type: str, explicit: str | None) -> Path:
    model_id = MODEL_IDS[model_type]
    if explicit:
        run = Path(explicit).resolve()
        if not complete_run(run):
            missing = [item for item in REQUIRED_RUN_FILES if not (run / item).is_file()]
            raise FileNotFoundError(f"incomplete {model_type} run {run}; missing {missing}")
        return run
    candidates = sorted(
        (path for path in DEEP_ROOT.glob(f"{model_id}_*") if complete_run(path)),
        key=lambda path: path.name,
    )
    if not candidates:
        raise FileNotFoundError(f"no complete DeepFootprinting run for {model_type} in {DEEP_ROOT}")
    return candidates[-1]


def deterministic_ranks(values: np.ndarray, feature_indices: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rank absolute, positive signed, and negative signed effects with IND tie-breaks."""
    if not np.isfinite(values).all():
        raise ValueError("cannot rank non-finite SHAP summaries")
    ranks = []
    for key in (-np.abs(values), -values, values):
        order = np.lexsort((feature_indices, key))
        rank = np.empty(len(values), dtype=np.int64)
        rank[order] = np.arange(1, len(values) + 1)
        ranks.append(rank)
    return ranks[0], ranks[1], ranks[2]


def summarize_features(
    model_type: str, run: Path, features: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = run / "shap_values_features_TREDNet.npy"
    values = np.load(path, mmap_mode="r")
    values = np.squeeze(values)
    if values.ndim != 2:
        raise ValueError(f"feature SHAP array must be 2D after squeeze, got {values.shape}")
    if values.shape[1] != len(features):
        raise ValueError(
            f"feature SHAP/features mismatch for {model_type}: {values.shape[1]} != {len(features)}"
        )
    if not np.isfinite(values).all():
        raise ValueError(f"non-finite feature SHAP values in {path}")

    mean_abs = np.abs(values).mean(axis=0)
    mean_signed = values.mean(axis=0)
    abs_rank, signed_desc_rank, signed_asc_rank = deterministic_ranks(
        mean_signed, features["IND"].to_numpy(dtype=np.int64)
    )
    # The absolute rank must use mean absolute SHAP, rather than abs(mean SHAP).
    absolute_order = np.lexsort((features["IND"].to_numpy(dtype=np.int64), -mean_abs))
    abs_rank = np.empty(len(features), dtype=np.int64)
    abs_rank[absolute_order] = np.arange(1, len(features) + 1)

    output = features.copy()
    output.insert(0, "model_type", model_type)
    output.insert(1, "model_id", MODEL_IDS[model_type])
    output["rank_mean_abs_shap"] = abs_rank
    output["rank_mean_signed_shap_desc"] = signed_desc_rank
    output["rank_mean_signed_shap_asc"] = signed_asc_rank
    output["mean_abs_shap"] = mean_abs
    output["mean_signed_shap"] = mean_signed
    output["median_signed_shap"] = np.median(values, axis=0)
    output["sd_signed_shap"] = values.std(axis=0, ddof=0)
    output["fraction_positive_shap"] = (values > 0).mean(axis=0)
    output["fraction_negative_shap"] = (values < 0).mean(axis=0)
    output["n_explained_sequences"] = values.shape[0]
    output = output.sort_values(
        ["rank_mean_abs_shap", "IND"], kind="mergesort"
    ).reset_index(drop=True)
    metadata = {
        "feature_shap_shape": "x".join(map(str, values.shape)),
        "feature_shap_dtype": str(values.dtype),
        "n_feature_attributions": int(values.shape[1]),
        "n_explained_sequences": int(values.shape[0]),
    }
    return output, metadata


def parse_html_matches(path: Path) -> pd.DataFrame:
    tables = pd.read_html(io.StringIO(path.read_text()))
    if len(tables) != 1:
        raise ValueError(f"expected one motif table in {path}, found {len(tables)}")
    table = tables[0]
    required = {"pattern", "num_seqlets"}
    if not required.issubset(table.columns):
        raise ValueError(f"motif report lacks {required - set(table.columns)}")
    return table


def h5_pattern_metadata(path: Path) -> tuple[pd.DataFrame, set[int]]:
    rows: list[dict[str, Any]] = []
    all_examples: set[int] = set()
    with h5py.File(path, "r") as handle:
        for group_name, polarity in (("pos_patterns", "positive"), ("neg_patterns", "negative")):
            if group_name not in handle:
                continue
            group = handle[group_name]
            for pattern_key in sorted(group.keys()):
                pattern = group[pattern_key]
                if "seqlets" not in pattern:
                    continue
                seqlets = pattern["seqlets"]
                examples = np.asarray(seqlets["example_idx"][:], dtype=np.int64)
                starts = np.asarray(seqlets["start"][:], dtype=np.int64)
                ends = np.asarray(seqlets["end"][:], dtype=np.int64)
                reverse = np.asarray(seqlets["is_revcomp"][:], dtype=bool)
                lengths = ends - starts
                n_seqlets = len(examples)
                if not (len(starts) == len(ends) == len(reverse) == n_seqlets):
                    raise ValueError(f"inconsistent seqlet arrays in {group_name}/{pattern_key}")
                declared = (
                    int(np.asarray(seqlets["n_seqlets"][:]).reshape(-1)[0])
                    if "n_seqlets" in seqlets
                    else n_seqlets
                )
                if declared != n_seqlets:
                    raise ValueError(
                        f"declared/observed seqlet mismatch for {group_name}/{pattern_key}: "
                        f"{declared} != {n_seqlets}"
                    )
                all_examples.update(int(value) for value in examples)
                rows.append(
                    {
                        "pattern": f"{group_name}.{pattern_key}",
                        "pattern_polarity": polarity,
                        "n_seqlets": n_seqlets,
                        "n_unique_sequence_examples": len(set(examples.tolist())),
                        "n_reverse_complement_seqlets": int(reverse.sum()),
                        "fraction_reverse_complement_seqlets": (
                            float(reverse.mean()) if n_seqlets else np.nan
                        ),
                        "median_seqlet_start": (
                            float(np.median(starts)) if n_seqlets else np.nan
                        ),
                        "median_seqlet_end": float(np.median(ends)) if n_seqlets else np.nan,
                        "median_seqlet_width": (
                            float(np.median(lengths)) if n_seqlets else np.nan
                        ),
                        "pattern_matrix_width": int(pattern["sequence"].shape[0]),
                    }
                )
    return pd.DataFrame(rows), all_examples


def summarize_patterns(
    model_type: str,
    run: Path,
    meme_catalog: dict[str, dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    html = parse_html_matches(run / "tfmodisco_output" / "motifs.html")
    h5_patterns, unique_examples = h5_pattern_metadata(
        run / "tfmodisco_output" / "modisco_results.h5"
    )
    if h5_patterns.empty:
        raise ValueError(f"no patterns in {run}")
    html_patterns = set(html["pattern"].astype(str))
    h5_names = set(h5_patterns["pattern"].astype(str))
    if not html_patterns.issubset(h5_names):
        raise ValueError(f"HTML patterns absent from HDF5: {sorted(html_patterns - h5_names)}")

    motif_rows: list[dict[str, Any]] = []
    match_columns = sorted(
        (
            int(match.group(1))
            for column in html.columns
            if (match := re.fullmatch(r"match(\d+)", str(column)))
        )
    )
    html_seqlets: dict[str, int] = {}
    for record in html.to_dict("records"):
        pattern = str(record["pattern"])
        html_seqlets[pattern] = int(record["num_seqlets"])
        for rank in match_columns:
            motif_value = record.get(f"match{rank}")
            q_value = record.get(f"qval{rank}")
            if pd.isna(motif_value) or str(motif_value).strip() in {"", "nan"}:
                continue
            motif_id = str(motif_value).strip()
            if motif_id not in meme_catalog:
                raise ValueError(f"{motif_id} from {pattern} is absent from local JASPAR MEME")
            if pd.isna(q_value):
                raise ValueError(f"missing q-value for {pattern} match {rank}")
            q_float = float(q_value)
            motif_rows.append(
                {
                    "model_type": model_type,
                    "model_id": MODEL_IDS[model_type],
                    "pattern": pattern,
                    "match_rank_zero_based": rank,
                    "motif_id": motif_id,
                    "tf_name": meme_catalog[motif_id]["tf_name"],
                    "motif_width": meme_catalog[motif_id]["motif_width"],
                    "q_value": q_float,
                    "significant_q_le_0_05": q_float <= 0.05,
                    "match_interpretation": (
                        "JASPAR sequence-similarity hypothesis; not evidence of TF binding"
                    ),
                }
            )
    motifs = pd.DataFrame(motif_rows)
    if motifs.empty:
        motifs = pd.DataFrame(
            columns=[
                "model_type",
                "model_id",
                "pattern",
                "match_rank_zero_based",
                "motif_id",
                "tf_name",
                "motif_width",
                "q_value",
                "significant_q_le_0_05",
                "match_interpretation",
            ]
        )

    patterns = h5_patterns.copy()
    patterns.insert(0, "model_type", model_type)
    patterns.insert(1, "model_id", MODEL_IDS[model_type])
    patterns["reported_num_seqlets"] = patterns["pattern"].map(html_seqlets)
    mismatch = patterns[
        patterns["reported_num_seqlets"].notna()
        & patterns["reported_num_seqlets"].ne(patterns["n_seqlets"])
    ]
    if not mismatch.empty:
        raise ValueError(
            "motifs.html and HDF5 seqlet counts differ: "
            + ", ".join(mismatch["pattern"].astype(str))
        )
    match_counts = motifs.groupby("pattern").size() if not motifs.empty else pd.Series(dtype=int)
    sig = motifs[motifs["significant_q_le_0_05"]].copy()
    sig_counts = sig.groupby("pattern").size() if not sig.empty else pd.Series(dtype=int)
    patterns["n_reported_jaspar_matches"] = (
        patterns["pattern"].map(match_counts).fillna(0).astype(int)
    )
    patterns["n_significant_jaspar_matches_q_le_0_05"] = (
        patterns["pattern"].map(sig_counts).fillna(0).astype(int)
    )
    best = (
        motifs.sort_values(["q_value", "match_rank_zero_based"], kind="mergesort")
        .drop_duplicates("pattern")
        .set_index("pattern")
        if not motifs.empty
        else pd.DataFrame()
    )
    for output_col, source_col in (
        ("best_match_motif_id", "motif_id"),
        ("best_match_tf_name", "tf_name"),
        ("best_match_q_value", "q_value"),
    ):
        patterns[output_col] = (
            patterns["pattern"].map(best[source_col]) if not best.empty else pd.NA
        )
    patterns["interpretation"] = (
        "de novo sequence pattern with JASPAR similarity hypotheses; not binding evidence"
    )
    patterns = patterns.sort_values(
        ["pattern_polarity", "n_seqlets", "pattern"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    motifs = motifs.sort_values(
        ["pattern", "q_value", "match_rank_zero_based"], kind="mergesort"
    ).reset_index(drop=True)
    metadata = {
        "n_tfmodisco_patterns": len(patterns),
        "n_positive_patterns": int(patterns["pattern_polarity"].eq("positive").sum()),
        "n_negative_patterns": int(patterns["pattern_polarity"].eq("negative").sum()),
        "n_pattern_seqlets": int(patterns["n_seqlets"].sum()),
        "n_unique_sequences_with_pattern_seqlets": len(unique_examples),
        "n_jaspar_matches_reported": len(motifs),
        "n_significant_jaspar_matches_q_le_0_05": int(
            motifs["significant_q_le_0_05"].sum()
        ),
        "n_distinct_significant_motifs": int(
            motifs.loc[motifs["significant_q_le_0_05"], "motif_id"].nunique()
        ),
        "n_distinct_significant_tf_names": int(
            motifs.loc[motifs["significant_q_le_0_05"], "tf_name"].nunique()
        ),
    }
    return patterns, motifs, metadata


def log_seqlet_count(run: Path, model_id: str) -> int | None:
    log = run / f"{model_id}_log.out"
    if not log.is_file():
        return None
    matches = re.findall(r"Using\s+(\d+)\s+positive seqlets", log.read_text(errors="replace"))
    return int(matches[-1]) if matches else None


def input_and_artifact_paths(run: Path, selection: dict[str, Any]) -> dict[str, Path]:
    phase_one_weights = Path(selection["phase_one_weights"])
    phase_two_weights = Path(selection["phase_two_weights"])
    roles = {
        "input_feature_list": FEATURES,
        "input_phase_one_model": phase_one_weights.parent / "phase_one_model.h5",
        "input_phase_one_weights": phase_one_weights,
        "input_phase_two_model": phase_two_weights.parent / "phase_two_model.keras",
        "input_phase_two_weights": phase_two_weights,
        "input_positive_bed": Path(selection["positive_bed"]),
        "input_control_bed": Path(selection["control_bed"]),
        "input_reference_fasta": Path(selection["fasta"]),
        "input_jaspar_meme": Path(selection["meme_file"]),
        "run_selection_manifest": run / "selection_manifest.json",
        "run_submission_manifest": run / "submission_manifest.tsv",
        "run_configured_script": run / "DeepExplainer_TREDNet.py",
        "output_feature_shap": run / "shap_values_features_TREDNet.npy",
        "output_actual_shap": run / "actual_shap.npy",
        "output_raw_shap": run / "raw_shap_explanations.npy",
        "output_explained_ohe": run / "ohe_seqs_to_explain.npy",
        "output_modisco_ohe": run / "ohe_tfmodisco.npz",
        "output_modisco_hypothetical_scores": run / "hyp_score_tfmodisco.npz",
        "output_modisco_h5": run / "tfmodisco_output" / "modisco_results.h5",
        "output_modisco_report": run / "tfmodisco_output" / "motifs.html",
    }
    model_id = str(selection.get("model_id", ""))
    if not model_id:
        # Selection manifests predate the explicit model_id field; the run basename is authoritative.
        model_id = next((value for value in MODEL_IDS.values() if run.name.startswith(value + "_")), "")
    for suffix, role in (("_log.out", "run_stdout"), ("_log.err", "run_stderr")):
        candidate = run / f"{model_id}{suffix}"
        if candidate.is_file():
            roles[role] = candidate
    return roles


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--models",
        nargs="+",
        choices=sorted(MODEL_IDS),
        default=sorted(MODEL_IDS),
        help="Models to summarize; default requires both enhancer and silencer.",
    )
    parser.add_argument("--enhancer-run", help="Explicit complete enhancer run directory")
    parser.add_argument("--silencer-run", help="Explicit complete silencer run directory")
    args = parser.parse_args()

    requested_models = list(dict.fromkeys(args.models))
    explicit = {"enhancer": args.enhancer_run, "silencer": args.silencer_run}
    RESULTS.mkdir(parents=True, exist_ok=True)

    features = pd.read_csv(FEATURES, sep="\t", dtype=str, keep_default_na=False)
    expected_columns = {"IND", "SRC_TYPE", "TARGET", "CL_ID", "FILE_ACCESSION", "CL_NAME"}
    if not expected_columns.issubset(features.columns):
        raise ValueError(f"feature list lacks {expected_columns - set(features.columns)}")
    features["IND"] = pd.to_numeric(features["IND"], errors="raise").astype(int)
    if len(features) != 4560 or features["IND"].duplicated().any():
        raise ValueError(f"expected 4,560 unique phase-I features, observed {len(features)}")

    meme_catalog, meme_background = read_meme_catalog(MEME)
    feature_tables: list[pd.DataFrame] = []
    pattern_tables: list[pd.DataFrame] = []
    motif_tables: list[pd.DataFrame] = []
    run_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    checksum_cache: dict[Path, str] = {}

    for model_type in requested_models:
        run = select_run(model_type, explicit[model_type])
        model_id = MODEL_IDS[model_type]
        submission = read_key_value_tsv(run / "submission_manifest.tsv")
        selection = json.loads((run / "selection_manifest.json").read_text())
        if submission.get("eid") != model_id:
            raise ValueError(
                f"run/model mismatch for {model_type}: {submission.get('eid')} != {model_id}"
            )
        if Path(selection["meme_file"]).resolve() != MEME.resolve():
            raise ValueError(f"run did not use the expected local JASPAR file: {run}")

        feature_table, feature_meta = summarize_features(model_type, run, features)
        patterns, motifs, pattern_meta = summarize_patterns(model_type, run, meme_catalog)
        feature_tables.append(feature_table)
        pattern_tables.append(patterns)
        motif_tables.append(motifs)

        run_row = {
            "model_type": model_type,
            "model_id": model_id,
            "run_directory": str(run.resolve()),
            "submitted_utc": submission.get("submitted_utc", ""),
            "slurm_job_id": submission.get("job_id", ""),
            "seed": int(selection["seed"]),
            "positive_sequences_available": int(selection["positive_sequences_available"]),
            "positive_sequences_explained": int(selection["positive_sequences_explained"]),
            "background_controls": int(selection["background_controls"]),
            "control_sequences_available": int(selection["control_sequences_available"]),
            "tfmodisco_input_positive_seqlets": log_seqlet_count(run, model_id),
            **feature_meta,
            **pattern_meta,
        }
        if run_row["positive_sequences_explained"] != run_row["n_explained_sequences"]:
            raise ValueError(f"selection/SHAP sequence count mismatch in {run}")
        run_rows.append(run_row)

        for role, path in input_and_artifact_paths(run, selection).items():
            file_rows.append(file_record(model_type, run, role, path, checksum_cache))

    feature_output = pd.concat(feature_tables, ignore_index=True)
    pattern_output = pd.concat(pattern_tables, ignore_index=True)
    motif_output = pd.concat(motif_tables, ignore_index=True)
    run_output = pd.DataFrame(run_rows).sort_values("model_type").reset_index(drop=True)
    file_output = pd.DataFrame(file_rows).sort_values(
        ["model_type", "file_role"], kind="mergesort"
    ).reset_index(drop=True)

    if len(feature_output) != 4560 * len(requested_models):
        raise ValueError("combined feature output does not have 4,560 rows per model")
    feature_output.to_csv(
        OUT_FEATURES,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    pattern_output.to_csv(OUT_PATTERNS, sep="\t", index=False, float_format="%.12g")
    motif_output.to_csv(OUT_MOTIFS, sep="\t", index=False, float_format="%.12g")
    run_output.to_csv(OUT_RUNS, sep="\t", index=False)
    file_output.to_csv(OUT_FILES, sep="\t", index=False)

    summary_rows: list[dict[str, Any]] = []
    for row in run_rows:
        for metric in (
            "n_feature_attributions",
            "n_explained_sequences",
            "tfmodisco_input_positive_seqlets",
            "n_tfmodisco_patterns",
            "n_pattern_seqlets",
            "n_unique_sequences_with_pattern_seqlets",
            "n_jaspar_matches_reported",
            "n_significant_jaspar_matches_q_le_0_05",
            "n_distinct_significant_motifs",
            "n_distinct_significant_tf_names",
        ):
            summary_rows.append(
                {
                    "model_type": row["model_type"],
                    "metric": metric,
                    "value": row.get(metric),
                    "denominator_or_note": (
                        "q <= 0.05" if "significant" in metric else "run-level count"
                    ),
                }
            )
        top = feature_output[
            (feature_output["model_type"] == row["model_type"])
            & (feature_output["rank_mean_abs_shap"] == 1)
        ].iloc[0]
        summary_rows.append(
            {
                "model_type": row["model_type"],
                "metric": "top_mean_abs_shap_feature",
                "value": f"IND={top['IND']};{top['SRC_TYPE']}:{top['TARGET']}:{top['CL_NAME']}",
                "denominator_or_note": f"mean_abs_shap={top['mean_abs_shap']:.12g}",
            }
        )
    summary_output = pd.DataFrame(summary_rows)
    summary_output.to_csv(OUT_SUMMARY, sep="\t", index=False)

    report_lines = [
        "# COPD DeepFootprinting interpretation (COPD-S4-R007)",
        "",
        (
            f"This result summarizes {len(requested_models)} model run(s) over "
            f"{len(features):,} TREDNet phase-I features per model. Feature SHAP values are "
            "reported as mean absolute and mean signed contributions across the explained "
            "positive sequences, with deterministic feature-index tie breaks."
        ),
        "",
        "## Run summary",
        "",
        "| Model | Explained sequences | Input positive seqlets | Retained patterns | Pattern seqlets | Significant JASPAR matches | Distinct TF names |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in run_rows:
        report_lines.append(
            f"| {row['model_type']} | {row['n_explained_sequences']:,} | "
            f"{row['tfmodisco_input_positive_seqlets'] or 0:,} | "
            f"{row['n_tfmodisco_patterns']:,} | {row['n_pattern_seqlets']:,} | "
            f"{row['n_significant_jaspar_matches_q_le_0_05']:,} | "
            f"{row['n_distinct_significant_tf_names']:,} |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            (
                "TF-MoDISco discovers recurring sequence patterns in model attributions. "
                "The reported JASPAR q-values quantify sequence similarity between those "
                "patterns and reference motifs. A match nominates a TF-family hypothesis; "
                "it does not establish occupancy, direction of regulation, or COPD-specific binding."
            ),
            "",
            (
                "Feature SHAP summaries explain these trained model predictions relative to the "
                "selected background controls. They are not causal effects and should not be "
                "interpreted as differential activity between people with and without COPD."
            ),
            "",
        ]
    )
    OUT_REPORT.write_text("\n".join(report_lines))

    output_paths = {
        "feature_shap": OUT_FEATURES,
        "tfmodisco_patterns": OUT_PATTERNS,
        "jaspar_motif_matches": OUT_MOTIFS,
        "run_metadata": OUT_RUNS,
        "run_file_checksums": OUT_FILES,
        "summary": OUT_SUMMARY,
        "interpretation": OUT_REPORT,
    }
    manifest = {
        "result_id": "COPD-S4-R007",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "models": requested_models,
        "status": (
            "complete_two_model_summary"
            if set(requested_models) == set(MODEL_IDS)
            else "partial_model_summary"
        ),
        "feature_count_per_model": len(features),
        "feature_ranking": {
            "primary": "descending mean absolute SHAP across explained positive sequences",
            "signed": "both descending and ascending mean signed SHAP ranks",
            "tie_break": "ascending phase-I feature IND",
        },
        "motif_match_significance_q_value": 0.05,
        "meme_background_frequencies": meme_background,
        "interpretation_limit": (
            "TF-MoDISco/JASPAR matches are sequence-similarity hypotheses and do not prove TF binding"
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "h5py": h5py.__version__,
        },
        "runs": run_rows,
        "outputs": {
            label: {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path, checksum_cache),
            }
            for label, path in output_paths.items()
        },
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__), checksum_cache),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"COPD-S4-R007: summarized {len(requested_models)} model(s), "
        f"{len(feature_output):,} feature rows, {len(pattern_output):,} patterns, "
        f"and {len(motif_output):,} JASPAR matches."
    )


if __name__ == "__main__":
    main()
