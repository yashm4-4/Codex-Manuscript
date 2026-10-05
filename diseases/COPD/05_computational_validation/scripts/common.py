#!/usr/bin/env python3
"""Shared, dependency-light helpers for COPD Section 5 validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[4]
COPD = ROOT / "diseases" / "COPD"
SECTION = COPD / "05_computational_validation"
DATA = SECTION / "data"
RESULTS = SECTION / "results"
LOGS = SECTION / "logs"

SECTION4_INPUTS = (
    COPD
    / "04_modeling/results/COPD-S4-R005_predicted_causal_population_genetics.tsv",
    COPD / "04_modeling/results/COPD-S4-R004_predicted_causal_regulatory_variants.tsv",
    COPD / "04_modeling/results/COPD-S4-R004_prioritized_candidates.tsv.gz",
)
SECTION4_ANCESTRAL_AUDIT = (
    COPD / "04_modeling/results/COPD-S4-R005_ancestral_allele_audit.tsv"
)

RSID_RE = re.compile(r"(?<![A-Za-z0-9])rs\d+(?![A-Za-z0-9])", re.IGNORECASE)


def ensure_directories() -> None:
    for path in (DATA, RESULTS, LOGS):
        path.mkdir(parents=True, exist_ok=True)


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(
        "r", newline=""
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def clean_chromosome(value: object) -> str:
    text = str(value).strip().removeprefix("chr")
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    if text == "23":
        text = "X"
    if text == "24":
        text = "Y"
    return text


def integer(value: object) -> int | None:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "."}:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def numeric(value: object) -> float | None:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "."}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def tokens(value: object) -> list[str]:
    return sorted(
        {
            item.strip()
            for item in str(value or "").split(";")
            if item.strip() and item.strip().lower() not in {"nan", "none", "."}
        }
    )


def extract_rsids(value: object) -> list[str]:
    return sorted(
        {match.lower() for match in RSID_RE.findall(str(value or ""))},
        key=lambda item: int(item[2:]),
    )


def candidate_rsids(row: Mapping[str, object]) -> list[str]:
    """Return rsIDs that identify this candidate, never its linked focal tags.

    ``source_focal_tags`` names the GWAS tag from which an LD proxy was found and
    therefore must not be treated as an identifier for the proxy itself.
    """

    result: set[str] = set()
    for column in (
        "ensembl_variation_ids",
        "selected_ensembl_variation_id",
        "candidate_rsid",
        "rsid",
    ):
        result.update(extract_rsids(row.get(column, "")))
    if is_true(row.get("is_gws_tag", False)):
        result.update(extract_rsids(row.get("gws_tag_ids", "")))
    result.update(extract_rsids(row.get("candidate_record_id", "")))
    return sorted(result, key=lambda item: int(item[2:]))


def _ancestral_rsid_map() -> dict[str, str]:
    if not SECTION4_ANCESTRAL_AUDIT.exists():
        return {}
    return {
        row["candidate_record_id"]: row.get("ensembl_variation_ids", "")
        for row in read_tsv(SECTION4_ANCESTRAL_AUDIT)
    }


def resolve_candidate_input() -> Path:
    override = os.environ.get("COPD_S5_CANDIDATE_INPUT", "").strip()
    if override:
        path = Path(override).resolve()
        if not path.exists():
            raise FileNotFoundError(
                f"COPD_S5_CANDIDATE_INPUT does not exist: {path}"
            )
        return path
    for path in SECTION4_INPUTS:
        if path.exists():
            return path
    expected = "\n".join(f"  - {path}" for path in SECTION4_INPUTS)
    raise FileNotFoundError(
        "Section 4 has not yet produced a predicted-causal/prioritized table. "
        f"Expected one of:\n{expected}"
    )


def load_candidates(path: Path | None = None) -> tuple[Path, list[dict[str, str]]]:
    """Load the final Section 4 predicted-causal candidate set.

    The population-annotated table is preferred because it contains exact
    Ensembl/dbSNP identifiers. The prioritized all-candidate table is accepted
    only after filtering the explicit predicted-causal flag.
    """

    input_path = path or resolve_candidate_input()
    rows = read_tsv(input_path)
    if not rows:
        raise ValueError(f"candidate input is empty: {input_path}")
    if "candidate_record_id" not in rows[0]:
        raise ValueError(f"candidate_record_id is absent from {input_path}")
    if "predicted_causal_regulatory" in rows[0]:
        rows = [row for row in rows if is_true(row["predicted_causal_regulatory"])]
    if not rows:
        raise ValueError(f"no predicted-causal candidates in {input_path}")

    supplemental_rsids = _ancestral_rsid_map()
    observed: set[str] = set()
    for row in rows:
        candidate_id = row["candidate_record_id"]
        if candidate_id in observed:
            raise ValueError(f"duplicate candidate_record_id: {candidate_id}")
        observed.add(candidate_id)
        if not row.get("ensembl_variation_ids") and candidate_id in supplemental_rsids:
            row["ensembl_variation_ids"] = supplemental_rsids[candidate_id]
        row["candidate_rsids"] = ";".join(candidate_rsids(row))
    return input_path, rows


def candidate_sort_key(row: Mapping[str, object]) -> tuple[int, str]:
    rank = integer(row.get("predicted_causal_priority_rank"))
    if rank is None:
        rank = integer(row.get("priority_rank_all_candidates"))
    return (rank if rank is not None else 10**12, str(row.get("candidate_record_id", "")))


def candidate_core(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "candidate_record_id": row.get("candidate_record_id", ""),
        "predicted_causal_priority_rank": row.get(
            "predicted_causal_priority_rank",
            row.get("priority_rank_all_candidates", ""),
        ),
        "chromosome_grch38": clean_chromosome(row.get("chromosome_grch38", "")),
        "position_grch38": integer(row.get("position_grch38", "")) or "",
        "ref": str(row.get("ref", "")).upper(),
        "alt": str(row.get("alt", "")).upper(),
        "variant_class": row.get("variant_class", ""),
        "candidate_rsids": row.get("candidate_rsids", ""),
        "linked_gwas_tag_ids": row.get(
            "linked_gwas_tag_ids",
            row.get("gws_tag_ids", ""),
        ),
        "linked_gwas_genes": row.get("linked_gwas_genes", ""),
        "assigned_selected_loci": row.get("assigned_selected_loci", ""),
        "priority_tier": row.get("priority_tier", ""),
        "predicted_causal_enhancer": row.get("predicted_causal_enhancer", ""),
        "predicted_causal_silencer": row.get("predicted_causal_silencer", ""),
    }


def _format_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return format(value, ".12g")
    return value


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    binary = path.suffix == ".gz"
    raw = path.open("wb") if binary else path.open("w", newline="")
    try:
        if binary:
            compressed = gzip.GzipFile(
                filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0
            )
            handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
        else:
            handle = raw
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: _format_value(row.get(name, "")) for name in fieldnames})
        if binary:
            handle.flush()
            handle.detach().close()
    finally:
        if not raw.closed:
            raw.close()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_manifest(path: Path, manifest: Mapping[str, object]) -> None:
    payload = {
        "created_utc": utc_now(),
        "software": {"python": platform.python_version()},
        **manifest,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
