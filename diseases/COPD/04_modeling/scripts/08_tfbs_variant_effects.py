#!/usr/bin/env python3
"""Score allele-specific JASPAR motif effects for COPD candidate variants.

For each enhancer or silencer model, this script uses only significant
(q <= 0.05 by default) JASPAR matches from COPD-S4-R007.  Every sequence
window that overlaps the complete REF or ALT allele span is scored on both
strands.  A relative PWM score of 0.8 is the prespecified motif-compatible
site threshold.  Crossing that threshold from REF to ALT is called a
computational disruption or creation.

These calls are sequence-model hypotheses.  They do not demonstrate TF
occupancy, TF identity at a locus, regulatory direction, or molecular binding.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, TextIO

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases" / "COPD" / "04_modeling"
RESULTS = SECTION / "results"

VARIANTS = RESULTS / "COPD-S4-R004_prioritized_candidates.tsv.gz"
MOTIF_MATCHES = RESULTS / "COPD-S4-R007_jaspar_motif_matches.tsv"
PATTERNS = RESULTS / "COPD-S4-R007_tfmodisco_patterns.tsv"
REF_FASTA = SECTION / "data" / "COPD_candidate_variants_ref_2001bp.fa"
ALT_FASTA = SECTION / "data" / "COPD_candidate_variants_alt_2001bp.fa"
MEME = (
    SECTION
    / "deepexplainer"
    / "backup"
    / "JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme"
)

OUT_SCORES = RESULTS / "COPD-S4-R008_allele_specific_motif_scores.tsv.gz"
OUT_CANDIDATES = RESULTS / "COPD-S4-R008_candidate_tfbs_summary.tsv.gz"
OUT_FRACTIONS = RESULTS / "COPD-S4-R008_variant_fractions.tsv"
OUT_TF_RANKS = RESULTS / "COPD-S4-R008_tf_rankings.tsv"
OUT_MOTIF_RANKS = RESULTS / "COPD-S4-R008_motif_rankings.tsv"
OUT_METHODS = RESULTS / "COPD-S4-R008_methods.md"
OUT_MANIFEST = RESULTS / "COPD-S4-R008_analysis_manifest.json"

MODEL_TYPES = ("enhancer", "silencer")
SEQUENCE_CENTER = 1000
PWM_PSEUDOCOUNT = 1e-6
DEFAULT_RELATIVE_THRESHOLD = 0.8
DEFAULT_MOTIF_Q_THRESHOLD = 0.05

BASE_TO_CODE = np.full(256, 255, dtype=np.uint8)
for base, code in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    BASE_TO_CODE[ord(base)] = code


@dataclass(frozen=True)
class Motif:
    motif_id: str
    tf_name: str
    probabilities: np.ndarray
    forward_log_odds: np.ndarray
    reverse_log_odds: np.ndarray
    forward_min: float
    forward_max: float
    reverse_min: float
    reverse_max: float

    @property
    def width(self) -> int:
        return int(self.probabilities.shape[0])


@dataclass(frozen=True)
class BestSite:
    relative_score: float
    start_zero_based: int | None
    strand: str
    n_valid_overlapping_windows: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna(False).astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes"}
    )


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    sequence: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier in {path}")
                sequence = []
            else:
                if name is None:
                    raise ValueError(f"sequence precedes FASTA header in {path}")
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence).upper()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for name, sequence in fasta_records(path):
        if name in records:
            raise ValueError(f"duplicate FASTA identifier {name} in {path}")
        records[name] = sequence
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def parse_meme(path: Path) -> tuple[dict[str, tuple[str, np.ndarray]], np.ndarray]:
    lines = path.read_text().splitlines()
    background: dict[str, float] = {}
    for index, line in enumerate(lines):
        if line.startswith("Background letter frequencies") and index + 1 < len(lines):
            fields = lines[index + 1].split()
            background = {
                fields[i]: float(fields[i + 1]) for i in range(0, len(fields), 2)
            }
            break
    if set(background) != {"A", "C", "G", "T"} or any(
        value <= 0 for value in background.values()
    ):
        raise ValueError(f"invalid MEME background in {path}: {background}")
    background_array = np.asarray(
        [background[base] for base in ("A", "C", "G", "T")], dtype=float
    )

    parsed: dict[str, tuple[str, np.ndarray]] = {}
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].startswith("MOTIF "):
            cursor += 1
            continue
        parts = lines[cursor].split(maxsplit=2)
        motif_id = parts[1]
        tf_name = parts[2] if len(parts) == 3 else motif_id
        cursor += 1
        while cursor < len(lines) and "letter-probability matrix" not in lines[cursor]:
            if lines[cursor].startswith("MOTIF "):
                raise ValueError(f"missing matrix for {motif_id}")
            cursor += 1
        if cursor >= len(lines):
            raise ValueError(f"missing matrix for {motif_id}")
        width_match = re.search(r"\bw\s*=\s*(\d+)", lines[cursor])
        if not width_match:
            raise ValueError(f"missing matrix width for {motif_id}")
        width = int(width_match.group(1))
        cursor += 1
        matrix: list[list[float]] = []
        while cursor < len(lines) and len(matrix) < width:
            fields = lines[cursor].split()
            if fields:
                try:
                    row = [float(value) for value in fields[:4]]
                except ValueError as exc:
                    raise ValueError(f"invalid matrix row for {motif_id}: {lines[cursor]}") from exc
                if len(row) != 4:
                    raise ValueError(f"short matrix row for {motif_id}")
                matrix.append(row)
            cursor += 1
        if len(matrix) != width:
            raise ValueError(f"matrix width mismatch for {motif_id}")
        if motif_id in parsed:
            raise ValueError(f"duplicate motif ID in MEME: {motif_id}")
        probabilities = np.asarray(matrix, dtype=float)
        if np.any(probabilities < 0) or np.any(probabilities.sum(axis=1) <= 0):
            raise ValueError(f"invalid probabilities for {motif_id}")
        parsed[motif_id] = (tf_name, probabilities)
    return parsed, background_array


def build_motif(
    motif_id: str,
    tf_name: str,
    probabilities: np.ndarray,
    background: np.ndarray,
) -> Motif:
    smoothed = probabilities + PWM_PSEUDOCOUNT
    smoothed = smoothed / smoothed.sum(axis=1, keepdims=True)
    reverse_probabilities = smoothed[::-1, :][:, [3, 2, 1, 0]]
    forward = np.log2(smoothed / background[None, :])
    reverse = np.log2(reverse_probabilities / background[None, :])
    return Motif(
        motif_id=motif_id,
        tf_name=tf_name,
        probabilities=probabilities,
        forward_log_odds=forward,
        reverse_log_odds=reverse,
        forward_min=float(forward.min(axis=1).sum()),
        forward_max=float(forward.max(axis=1).sum()),
        reverse_min=float(reverse.min(axis=1).sum()),
        reverse_max=float(reverse.max(axis=1).sum()),
    )


def encode_sequence(sequence: str) -> np.ndarray:
    raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return BASE_TO_CODE[raw]


def best_overlapping_site(
    encoded: np.ndarray,
    motif: Motif,
    allele_start: int,
    allele_length: int,
) -> BestSite:
    """Return the best two-strand site among windows intersecting an allele span."""
    if allele_length < 1:
        raise ValueError("VCF-normalized allele length must be positive")
    width = motif.width
    first_start = max(0, allele_start - width + 1)
    last_start = min(len(encoded) - width, allele_start + allele_length - 1)
    if first_start > last_start:
        return BestSite(np.nan, None, "", 0)
    starts = np.arange(first_start, last_start + 1, dtype=np.int64)
    positions = starts[:, None] + np.arange(width, dtype=np.int64)[None, :]
    windows = encoded[positions]
    valid = np.all(windows < 4, axis=1)
    if not valid.any():
        return BestSite(np.nan, None, "", 0)
    starts = starts[valid]
    windows = windows[valid]
    columns = np.arange(width, dtype=np.int64)[None, :]
    forward_raw = motif.forward_log_odds[columns, windows].sum(axis=1)
    reverse_raw = motif.reverse_log_odds[columns, windows].sum(axis=1)
    forward_denominator = motif.forward_max - motif.forward_min
    reverse_denominator = motif.reverse_max - motif.reverse_min
    if forward_denominator <= 0 or reverse_denominator <= 0:
        raise ValueError(f"degenerate PWM for {motif.motif_id}")
    forward_relative = (forward_raw - motif.forward_min) / forward_denominator
    reverse_relative = (reverse_raw - motif.reverse_min) / reverse_denominator
    combined = np.column_stack((forward_relative, reverse_relative))
    flat_index = int(np.argmax(combined))
    window_index, strand_index = np.unravel_index(flat_index, combined.shape)
    score = float(np.clip(combined[window_index, strand_index], 0.0, 1.0))
    return BestSite(
        relative_score=score,
        start_zero_based=int(starts[window_index]),
        strand="+" if strand_index == 0 else "-",
        n_valid_overlapping_windows=len(starts),
    )


def consolidate_motif_matches(
    motif_matches: pd.DataFrame,
    patterns: pd.DataFrame,
    q_threshold: float,
) -> pd.DataFrame:
    required = {"model_type", "pattern", "motif_id", "tf_name", "q_value"}
    if not required.issubset(motif_matches.columns):
        raise ValueError(f"R007 motif table lacks {required - set(motif_matches.columns)}")
    selected = motif_matches[pd.to_numeric(motif_matches["q_value"]) <= q_threshold].copy()
    selected = selected[selected["model_type"].isin(MODEL_TYPES)]
    if selected.empty:
        raise ValueError(f"no JASPAR motif matches at q <= {q_threshold}")
    selected["q_value"] = pd.to_numeric(selected["q_value"], errors="raise")

    if {"model_type", "pattern", "n_seqlets"}.issubset(patterns.columns):
        pattern_counts = patterns[["model_type", "pattern", "n_seqlets"]].drop_duplicates()
        selected = selected.merge(
            pattern_counts, on=["model_type", "pattern"], how="left", validate="many_to_one"
        )
    else:
        selected["n_seqlets"] = np.nan

    rows: list[dict[str, Any]] = []
    for (model_type, motif_id), group in selected.groupby(
        ["model_type", "motif_id"], sort=True
    ):
        names = sorted(set(group["tf_name"].astype(str)))
        if len(names) != 1:
            raise ValueError(f"inconsistent TF names for {motif_id}: {names}")
        rows.append(
            {
                "model_type": model_type,
                "motif_id": motif_id,
                "tf_name": names[0],
                "minimum_match_q_value": float(group["q_value"].min()),
                "matched_patterns": ";".join(sorted(set(group["pattern"].astype(str)))),
                "n_matched_patterns": group["pattern"].nunique(),
                "sum_matched_pattern_seqlets": int(group["n_seqlets"].dropna().sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["model_type", "minimum_match_q_value", "motif_id"], kind="mergesort"
    ).reset_index(drop=True)


def validate_variants_and_fastas(
    variants: pd.DataFrame,
    ref_sequences: dict[str, str],
    alt_sequences: dict[str, str],
) -> pd.DataFrame:
    required = {"candidate_record_id", "ref", "alt"}
    for model_type in MODEL_TYPES:
        required.update(
            {f"predicted_{model_type}_fpr5", f"predicted_causal_{model_type}"}
        )
    if not required.issubset(variants.columns):
        raise ValueError(f"R004 candidate table lacks {required - set(variants.columns)}")
    if variants["candidate_record_id"].duplicated().any():
        raise ValueError("R004 contains duplicate candidate_record_id values")
    if set(ref_sequences) != set(alt_sequences):
        raise ValueError("REF and ALT FASTA identifier sets differ")

    candidate_ids = set(variants["candidate_record_id"].astype(str))
    extra = set(ref_sequences) - candidate_ids
    if extra:
        raise ValueError(f"FASTA contains {len(extra)} IDs absent from R004")
    scorable = variants[
        variants["candidate_record_id"].astype(str).isin(ref_sequences)
    ].copy()
    scorable["candidate_record_id"] = scorable["candidate_record_id"].astype(str)
    for model_type in MODEL_TYPES:
        scorable[f"predicted_{model_type}_fpr5"] = as_bool(
            scorable[f"predicted_{model_type}_fpr5"]
        )
        scorable[f"predicted_causal_{model_type}"] = as_bool(
            scorable[f"predicted_causal_{model_type}"]
        )

    for row in scorable.itertuples(index=False):
        candidate_id = row.candidate_record_id
        ref = str(row.ref).upper()
        alt = str(row.alt).upper()
        if not ref or not alt or set(ref + alt) - set("ACGT"):
            raise ValueError(f"non-ACGT allele in scorable candidate {candidate_id}")
        ref_sequence = ref_sequences[candidate_id]
        alt_sequence = alt_sequences[candidate_id]
        if len(ref_sequence) != 2001 or len(alt_sequence) != 2001:
            raise ValueError(f"non-2001-bp FASTA sequence for {candidate_id}")
        if ref_sequence[SEQUENCE_CENTER : SEQUENCE_CENTER + len(ref)] != ref:
            raise ValueError(f"REF allele is not centered for {candidate_id}")
        if alt_sequence[SEQUENCE_CENTER : SEQUENCE_CENTER + len(alt)] != alt:
            raise ValueError(f"ALT allele is not centered for {candidate_id}")
    return scorable


def score_all(
    variants: pd.DataFrame,
    ref_sequences: dict[str, str],
    alt_sequences: dict[str, str],
    motif_definitions: pd.DataFrame,
    motifs: dict[str, Motif],
    relative_threshold: float,
) -> pd.DataFrame:
    memberships: dict[str, list[dict[str, Any]]] = {}
    for record in motif_definitions.to_dict("records"):
        memberships.setdefault(str(record["motif_id"]), []).append(record)

    rows: list[dict[str, Any]] = []
    for index, variant in enumerate(variants.itertuples(index=False), start=1):
        candidate_id = str(variant.candidate_record_id)
        ref_allele = str(variant.ref).upper()
        alt_allele = str(variant.alt).upper()
        ref_encoded = encode_sequence(ref_sequences[candidate_id])
        alt_encoded = encode_sequence(alt_sequences[candidate_id])
        for motif_id in sorted(memberships):
            motif = motifs[motif_id]
            ref_site = best_overlapping_site(
                ref_encoded, motif, SEQUENCE_CENTER, len(ref_allele)
            )
            alt_site = best_overlapping_site(
                alt_encoded, motif, SEQUENCE_CENTER, len(alt_allele)
            )
            if np.isnan(ref_site.relative_score) or np.isnan(alt_site.relative_score):
                ref_hit = False
                alt_hit = False
                score_delta = np.nan
                score_status = "no_valid_overlapping_window"
            else:
                ref_hit = ref_site.relative_score >= relative_threshold
                alt_hit = alt_site.relative_score >= relative_threshold
                score_delta = alt_site.relative_score - ref_site.relative_score
                score_status = "ok"
            for membership in memberships[motif_id]:
                model_type = str(membership["model_type"])
                predicted = bool(getattr(variant, f"predicted_{model_type}_fpr5"))
                causal = bool(getattr(variant, f"predicted_causal_{model_type}"))
                rows.append(
                    {
                        "model_type": model_type,
                        "candidate_record_id": candidate_id,
                        "chromosome_grch38": getattr(variant, "chromosome_grch38", pd.NA),
                        "position_grch38": getattr(variant, "position_grch38", pd.NA),
                        "ref": ref_allele,
                        "alt": alt_allele,
                        "variant_class": getattr(variant, "variant_class", pd.NA),
                        "model_predicted_fpr5": predicted,
                        "model_predicted_causal": causal,
                        "motif_id": motif_id,
                        "tf_name": motif.tf_name,
                        "minimum_match_q_value": membership["minimum_match_q_value"],
                        "matched_patterns": membership["matched_patterns"],
                        "n_matched_patterns": membership["n_matched_patterns"],
                        "sum_matched_pattern_seqlets": membership[
                            "sum_matched_pattern_seqlets"
                        ],
                        "motif_width": motif.width,
                        "relative_score_threshold": relative_threshold,
                        "ref_relative_score": ref_site.relative_score,
                        "alt_relative_score": alt_site.relative_score,
                        "delta_alt_minus_ref_relative_score": score_delta,
                        "absolute_delta_relative_score": abs(score_delta),
                        "ref_best_window_start_zero_based": ref_site.start_zero_based,
                        "ref_best_window_strand": ref_site.strand,
                        "alt_best_window_start_zero_based": alt_site.start_zero_based,
                        "alt_best_window_strand": alt_site.strand,
                        "ref_valid_overlapping_windows": ref_site.n_valid_overlapping_windows,
                        "alt_valid_overlapping_windows": alt_site.n_valid_overlapping_windows,
                        "ref_motif_compatible_site": ref_hit,
                        "alt_motif_compatible_site": alt_hit,
                        "either_allele_motif_compatible_site": ref_hit or alt_hit,
                        "disrupts_motif_compatible_site": ref_hit and not alt_hit,
                        "creates_motif_compatible_site": alt_hit and not ref_hit,
                        "retains_motif_compatible_site": ref_hit and alt_hit,
                        "score_status": score_status,
                    }
                )
        if index % 1000 == 0:
            print(f"Scored {index:,}/{len(variants):,} candidates", flush=True)
    output = pd.DataFrame(rows)
    return output.sort_values(
        ["model_type", "candidate_record_id", "minimum_match_q_value", "motif_id"],
        kind="mergesort",
    ).reset_index(drop=True)


def aggregate_candidates(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (model_type, candidate_id), group in scores.groupby(
        ["model_type", "candidate_record_id"], sort=True
    ):
        successful = group["score_status"].eq("ok")
        if successful.any():
            largest_index = group.loc[successful, "absolute_delta_relative_score"].idxmax()
            largest = group.loc[largest_index]
            largest_motif = largest["motif_id"]
            largest_tf = largest["tf_name"]
            largest_delta = largest["delta_alt_minus_ref_relative_score"]
            largest_abs_delta = largest["absolute_delta_relative_score"]
        else:
            largest_motif = pd.NA
            largest_tf = pd.NA
            largest_delta = np.nan
            largest_abs_delta = np.nan
        first = group.iloc[0]
        rows.append(
            {
                "model_type": model_type,
                "candidate_record_id": candidate_id,
                "model_predicted_fpr5": bool(first["model_predicted_fpr5"]),
                "model_predicted_causal": bool(first["model_predicted_causal"]),
                "n_significant_matched_motifs_tested": len(group),
                "n_motifs_successfully_scored": int(successful.sum()),
                "n_ref_motif_compatible_sites": int(
                    group["ref_motif_compatible_site"].sum()
                ),
                "n_alt_motif_compatible_sites": int(
                    group["alt_motif_compatible_site"].sum()
                ),
                "n_disrupted_motif_compatible_sites": int(
                    group["disrupts_motif_compatible_site"].sum()
                ),
                "n_created_motif_compatible_sites": int(
                    group["creates_motif_compatible_site"].sum()
                ),
                "n_retained_motif_compatible_sites": int(
                    group["retains_motif_compatible_site"].sum()
                ),
                "any_motif_compatible_site": bool(
                    group["either_allele_motif_compatible_site"].any()
                ),
                "disrupts_any_motif_compatible_site": bool(
                    group["disrupts_motif_compatible_site"].any()
                ),
                "creates_any_motif_compatible_site": bool(
                    group["creates_motif_compatible_site"].any()
                ),
                "retains_any_motif_compatible_site": bool(
                    group["retains_motif_compatible_site"].any()
                ),
                "largest_effect_motif_id": largest_motif,
                "largest_effect_tf_name": largest_tf,
                "largest_signed_delta_relative_score": largest_delta,
                "largest_absolute_delta_relative_score": largest_abs_delta,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["model_type", "candidate_record_id"], kind="mergesort"
    ).reset_index(drop=True)


def scope_masks(frame: pd.DataFrame) -> dict[str, pd.Series]:
    return {
        "all_scorable_candidates": pd.Series(True, index=frame.index),
        "model_predicted_fpr5": as_bool(frame["model_predicted_fpr5"]),
        "model_predicted_causal": as_bool(frame["model_predicted_causal"]),
    }


def variant_fractions(candidate_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for model_type, frame in candidate_summary.groupby("model_type", sort=True):
        for scope, mask in scope_masks(frame).items():
            subset = frame[mask]
            denominator = len(subset)
            any_site = int(subset["any_motif_compatible_site"].sum())
            disrupting = int(subset["disrupts_any_motif_compatible_site"].sum())
            creating = int(subset["creates_any_motif_compatible_site"].sum())
            both = int(
                (
                    subset["disrupts_any_motif_compatible_site"]
                    & subset["creates_any_motif_compatible_site"]
                ).sum()
            )
            altered = int(
                (
                    subset["disrupts_any_motif_compatible_site"]
                    | subset["creates_any_motif_compatible_site"]
                ).sum()
            )
            retained = int(subset["retains_any_motif_compatible_site"].sum())
            fully_scored = int(
                subset["n_motifs_successfully_scored"].eq(
                    subset["n_significant_matched_motifs_tested"]
                ).sum()
            )
            rows.append(
                {
                    "model_type": model_type,
                    "analysis_scope": scope,
                    "denominator_candidate_records": denominator,
                    "n_fully_scored_for_all_significant_motifs": fully_scored,
                    "n_with_ref_or_alt_motif_compatible_site": any_site,
                    "fraction_with_site_of_scope": any_site / denominator if denominator else np.nan,
                    "n_disrupting_at_least_one_site": disrupting,
                    "fraction_disrupting_of_scope": (
                        disrupting / denominator if denominator else np.nan
                    ),
                    "fraction_disrupting_among_variants_with_site": (
                        disrupting / any_site if any_site else np.nan
                    ),
                    "n_creating_at_least_one_site": creating,
                    "fraction_creating_of_scope": creating / denominator if denominator else np.nan,
                    "fraction_creating_among_variants_with_site": (
                        creating / any_site if any_site else np.nan
                    ),
                    "n_disrupting_or_creating_site": altered,
                    "fraction_disrupting_or_creating_of_scope": (
                        altered / denominator if denominator else np.nan
                    ),
                    "n_both_disrupting_and_creating_across_motifs": both,
                    "n_retaining_at_least_one_site": retained,
                    "n_without_any_ref_or_alt_site": denominator - any_site,
                    "site_definition": (
                        "best overlapping two-strand relative PWM score >= 0.8; "
                        "creation/disruption requires threshold crossing"
                    ),
                }
            )
    return pd.DataFrame(rows)


def ranking_table(
    scores: pd.DataFrame,
    candidate_summary: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Rank TFs or motifs after collapsing to one row per candidate and group."""
    identity_columns = ["model_type", group_column]
    if group_column == "tf_name":
        per_candidate = (
            scores.groupby(["model_type", "tf_name", "candidate_record_id"], sort=True)
            .agg(
                model_predicted_fpr5=("model_predicted_fpr5", "first"),
                model_predicted_causal=("model_predicted_causal", "first"),
                any_site=("either_allele_motif_compatible_site", "any"),
                disrupts=("disrupts_motif_compatible_site", "any"),
                creates=("creates_motif_compatible_site", "any"),
                retains=("retains_motif_compatible_site", "any"),
                max_abs_delta=("absolute_delta_relative_score", "max"),
                n_motif_models=("motif_id", "nunique"),
            )
            .reset_index()
        )
        metadata = (
            scores.groupby(["model_type", "tf_name"], sort=True)
            .agg(
                motif_ids=("motif_id", lambda values: ";".join(sorted(set(values)))),
                n_distinct_motif_ids=("motif_id", "nunique"),
                minimum_match_q_value=("minimum_match_q_value", "min"),
            )
            .reset_index()
        )
    else:
        per_candidate = scores.rename(
            columns={
                "either_allele_motif_compatible_site": "any_site",
                "disrupts_motif_compatible_site": "disrupts",
                "creates_motif_compatible_site": "creates",
                "retains_motif_compatible_site": "retains",
                "absolute_delta_relative_score": "max_abs_delta",
            }
        )[
            [
                "model_type",
                "motif_id",
                "candidate_record_id",
                "model_predicted_fpr5",
                "model_predicted_causal",
                "any_site",
                "disrupts",
                "creates",
                "retains",
                "max_abs_delta",
            ]
        ].copy()
        metadata = (
            scores.groupby(["model_type", "motif_id"], sort=True)
            .agg(
                tf_name=("tf_name", "first"),
                minimum_match_q_value=("minimum_match_q_value", "min"),
                matched_patterns=(
                    "matched_patterns",
                    lambda values: ";".join(
                        sorted(
                            {
                                token
                                for value in values
                                for token in str(value).split(";")
                                if token
                            }
                        )
                    ),
                ),
            )
            .reset_index()
        )

    rows: list[dict[str, Any]] = []
    for model_type, model_candidates in candidate_summary.groupby("model_type", sort=True):
        model_groups = per_candidate[per_candidate["model_type"].eq(model_type)]
        group_values = sorted(model_groups[group_column].unique())
        masks = scope_masks(model_candidates)
        for scope, candidate_mask in masks.items():
            included_ids = set(
                model_candidates.loc[candidate_mask, "candidate_record_id"].astype(str)
            )
            denominator = len(included_ids)
            for value in group_values:
                group = model_groups[
                    model_groups[group_column].eq(value)
                    & model_groups["candidate_record_id"].isin(included_ids)
                ]
                rows.append(
                    {
                        "model_type": model_type,
                        "analysis_scope": scope,
                        group_column: value,
                        "denominator_candidate_records": denominator,
                        "n_variants_with_site": int(group["any_site"].sum()),
                        "fraction_with_site_of_scope": (
                            float(group["any_site"].sum()) / denominator
                            if denominator
                            else np.nan
                        ),
                        "n_variants_disrupting_site": int(group["disrupts"].sum()),
                        "fraction_disrupting_of_scope": (
                            float(group["disrupts"].sum()) / denominator
                            if denominator
                            else np.nan
                        ),
                        "n_variants_creating_site": int(group["creates"].sum()),
                        "fraction_creating_of_scope": (
                            float(group["creates"].sum()) / denominator
                            if denominator
                            else np.nan
                        ),
                        "n_variants_retaining_site": int(group["retains"].sum()),
                        "n_variants_disrupting_or_creating_site": int(
                            (group["disrupts"] | group["creates"]).sum()
                        ),
                        "maximum_absolute_delta_relative_score": (
                            float(group["max_abs_delta"].max()) if len(group) else np.nan
                        ),
                        "mean_absolute_delta_relative_score": (
                            float(group["max_abs_delta"].mean()) if len(group) else np.nan
                        ),
                    }
                )
    output = pd.DataFrame(rows).merge(metadata, on=identity_columns, how="left")
    output = output.sort_values(
        [
            "model_type",
            "analysis_scope",
            "n_variants_disrupting_or_creating_site",
            "n_variants_with_site",
            "maximum_absolute_delta_relative_score",
            group_column,
        ],
        ascending=[True, True, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    output.insert(
        2,
        "priority_rank_within_model_scope",
        output.groupby(["model_type", "analysis_scope"]).cumcount() + 1,
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", type=Path, default=VARIANTS)
    parser.add_argument("--motif-matches", type=Path, default=MOTIF_MATCHES)
    parser.add_argument("--patterns", type=Path, default=PATTERNS)
    parser.add_argument("--ref-fasta", type=Path, default=REF_FASTA)
    parser.add_argument("--alt-fasta", type=Path, default=ALT_FASTA)
    parser.add_argument("--meme", type=Path, default=MEME)
    parser.add_argument(
        "--relative-score-threshold", type=float, default=DEFAULT_RELATIVE_THRESHOLD
    )
    parser.add_argument("--motif-q-threshold", type=float, default=DEFAULT_MOTIF_Q_THRESHOLD)
    args = parser.parse_args()
    if not 0 < args.relative_score_threshold <= 1:
        raise ValueError("relative score threshold must be in (0, 1]")
    if not 0 < args.motif_q_threshold <= 1:
        raise ValueError("motif q-value threshold must be in (0, 1]")
    for path in (
        args.variants,
        args.motif_matches,
        args.patterns,
        args.ref_fasta,
        args.alt_fasta,
        args.meme,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    RESULTS.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(args.variants, sep="\t", low_memory=False)
    motif_matches = pd.read_csv(args.motif_matches, sep="\t")
    patterns = pd.read_csv(args.patterns, sep="\t")
    motif_definitions = consolidate_motif_matches(
        motif_matches, patterns, args.motif_q_threshold
    )

    parsed_motifs, background = parse_meme(args.meme)
    selected_ids = set(motif_definitions["motif_id"].astype(str))
    missing = selected_ids - set(parsed_motifs)
    if missing:
        raise ValueError(f"selected motifs absent from MEME: {sorted(missing)}")
    motifs = {
        motif_id: build_motif(
            motif_id,
            parsed_motifs[motif_id][0],
            parsed_motifs[motif_id][1],
            background,
        )
        for motif_id in sorted(selected_ids)
    }
    for row in motif_definitions.itertuples(index=False):
        if row.tf_name != motifs[row.motif_id].tf_name:
            raise ValueError(
                f"R007/local MEME TF-name mismatch for {row.motif_id}: "
                f"{row.tf_name} != {motifs[row.motif_id].tf_name}"
            )

    print("Loading and validating paired candidate FASTA records", flush=True)
    ref_sequences = read_fasta(args.ref_fasta)
    alt_sequences = read_fasta(args.alt_fasta)
    scorable = validate_variants_and_fastas(variants, ref_sequences, alt_sequences)
    print(
        f"Scoring {len(scorable):,} candidates against {len(selected_ids):,} distinct "
        f"significant JASPAR motifs across {motif_definitions['model_type'].nunique()} model(s)",
        flush=True,
    )
    scores = score_all(
        scorable,
        ref_sequences,
        alt_sequences,
        motif_definitions,
        motifs,
        args.relative_score_threshold,
    )
    candidate_summary = aggregate_candidates(scores)
    fractions = variant_fractions(candidate_summary)
    tf_ranks = ranking_table(scores, candidate_summary, "tf_name")
    motif_ranks = ranking_table(scores, candidate_summary, "motif_id")

    scores.to_csv(
        OUT_SCORES,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    candidate_summary.to_csv(
        OUT_CANDIDATES,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    fractions.to_csv(OUT_FRACTIONS, sep="\t", index=False, float_format="%.12g")
    tf_ranks.to_csv(OUT_TF_RANKS, sep="\t", index=False, float_format="%.12g")
    motif_ranks.to_csv(OUT_MOTIF_RANKS, sep="\t", index=False, float_format="%.12g")

    methods = f"""# COPD allele-specific motif analysis (COPD-S4-R008)

## Tested records and denominators

The analysis tested all {len(scorable):,} candidate records with validated paired 2,001-bp REF and ALT sequences. Results remain separated by enhancer and silencer model. Each model has three explicitly reported denominators: all sequence-scorable candidates, candidates passing that model's held-out 5% FPR prediction threshold, and candidates satisfying that model's R004 predicted-causal definition. Candidate identity is the stable `candidate_record_id`.

## Motif selection and scoring

Only model-specific TF-MoDISco matches to JASPAR 2024 motifs with q <= {args.motif_q_threshold:g} were tested. JASPAR probability matrices were read from the same local MEME file used to generate the R007 motif report. A pseudocount of {PWM_PSEUDOCOUNT:g} was added to every base probability and each PWM row was renormalized. Log2 odds were calculated against the MEME background frequencies (A={background[0]:g}, C={background[1]:g}, G={background[2]:g}, T={background[3]:g}).

For each allele and motif, every sequence window that overlaps any base of the complete allele span was evaluated on both strands. This full-span rule covers SNVs, multinucleotide variants, insertions, and deletions. The best log-odds score was min-max normalized to the PWM's theoretical range. A relative score >= {args.relative_score_threshold:g} defines a motif-compatible sequence site. A disruption is REF >= {args.relative_score_threshold:g} and ALT < {args.relative_score_threshold:g}; a creation is REF < {args.relative_score_threshold:g} and ALT >= {args.relative_score_threshold:g}. Score changes that do not cross the threshold are retained as quantitative deltas but are not labeled creation or disruption.

## Interpretation boundary

The {args.relative_score_threshold:g} threshold is a prespecified computational heuristic, not an empirical binding probability. TF-MoDISco-to-JASPAR similarity nominates a motif or TF-family hypothesis. Neither a motif match nor an allele-specific PWM threshold crossing proves that the named TF binds in lung tissue, that binding changes in COPD, or that the variant is causal. Experimental occupancy and perturbation assays are required.
"""
    OUT_METHODS.write_text(methods)

    input_paths = {
        "R004_prioritized_candidates": args.variants,
        "R007_motif_matches": args.motif_matches,
        "R007_patterns": args.patterns,
        "reference_allele_fasta": args.ref_fasta,
        "alternate_allele_fasta": args.alt_fasta,
        "JASPAR2024_MEME": args.meme,
    }
    output_paths = {
        "allele_specific_motif_scores": OUT_SCORES,
        "candidate_tfbs_summary": OUT_CANDIDATES,
        "variant_fractions": OUT_FRACTIONS,
        "TF_rankings": OUT_TF_RANKS,
        "motif_rankings": OUT_MOTIF_RANKS,
        "methods": OUT_METHODS,
    }
    manifest = {
        "result_id": "COPD-S4-R008",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "models_analyzed": sorted(motif_definitions["model_type"].unique()),
        "candidate_records_in_R004": len(variants),
        "sequence_scorable_candidate_records": len(scorable),
        "significant_model_motif_pairs": len(motif_definitions),
        "distinct_jaspar_motifs": len(selected_ids),
        "motif_match_q_value_threshold": args.motif_q_threshold,
        "relative_pwm_score_threshold": args.relative_score_threshold,
        "pwm_pseudocount": PWM_PSEUDOCOUNT,
        "pwm_background": {
            base: float(value)
            for base, value in zip(("A", "C", "G", "T"), background)
        },
        "window_policy": (
            "score both strands of every motif window intersecting the complete allele span; "
            "select the maximum relative PWM score separately for REF and ALT"
        ),
        "site_definition": f"relative PWM score >= {args.relative_score_threshold:g}",
        "disruption_definition": "REF site and no ALT site at the threshold",
        "creation_definition": "ALT site and no REF site at the threshold",
        "reported_denominators": [
            "all_scorable_candidates",
            "model_predicted_fpr5",
            "model_predicted_causal",
        ],
        "interpretation_limit": (
            "motif-compatible sequence predictions do not prove TF binding or causality"
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {
            label: {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for label, path in input_paths.items()
        },
        "outputs": {
            label: {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for label, path in output_paths.items()
        },
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__)),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"COPD-S4-R008: wrote {len(scores):,} model-candidate-motif scores, "
        f"{len(candidate_summary):,} model-candidate summaries, and explicit "
        "creation/disruption denominators."
    )


if __name__ == "__main__":
    main()
