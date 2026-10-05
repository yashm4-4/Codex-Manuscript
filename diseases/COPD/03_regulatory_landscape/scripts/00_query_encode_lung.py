#!/usr/bin/env python3
"""Inventory released human lung accessibility and H3K27 peak experiments.

This is a metadata-only ENCODE REST query. It does not download signal or peak
files. The output supports an explicit, auditable biosample choice for COPD
regulatory mapping and TREDNet model training.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
import requests


OUT = Path(__file__).resolve().parents[1] / "data"
OUT.mkdir(parents=True, exist_ok=True)
BASE = "https://www.encodeproject.org"


def search(assay_title: str, target: str = "") -> list[dict]:
    params: list[tuple[str, str]] = [
        ("type", "Experiment"),
        ("status", "released"),
        ("replicates.library.biosample.donor.organism.scientific_name", "Homo sapiens"),
        ("biosample_ontology.organ_slims", "lung"),
        ("assay_title", assay_title),
        ("limit", "all"),
        ("format", "json"),
    ]
    if target:
        params.append(("target.label", target))
    response = requests.get(
        f"{BASE}/search/",
        params=params,
        headers={"Accept": "application/json"},
        timeout=120,
    )
    response.raise_for_status()
    return response.json().get("@graph", [])


def audit_categories(exp: dict, level: str) -> str:
    return ";".join(
        sorted({item.get("category", "") for item in exp.get("audit", {}).get(level, []) if item.get("category")})
    )


def biosample_accessions(exp: dict) -> str:
    values = {
        rep.get("library", {}).get("biosample", {}).get("accession", "")
        for rep in exp.get("replicates", [])
    }
    return ";".join(sorted(v for v in values if v))


def main() -> None:
    queries = [
        ("ATAC-seq", ""),
        ("DNase-seq", ""),
        ("Histone ChIP-seq", "H3K27ac"),
        ("Histone ChIP-seq", "H3K27me3"),
    ]
    rows: list[dict[str, object]] = []
    for assay, target in queries:
        for exp in search(assay, target):
            ontology = exp.get("biosample_ontology", {})
            rows.append(
                {
                    "experiment": exp.get("accession", ""),
                    "assay": assay,
                    "target": target,
                    "biosample_term": ontology.get("term_name", ""),
                    "biosample_classification": ontology.get("classification", ""),
                    "biosample_summary": exp.get("biosample_summary", ""),
                    "biosample_accessions": biosample_accessions(exp),
                    "description": exp.get("description", ""),
                    "lab": exp.get("lab", {}).get("title", ""),
                    "dbxrefs": ";".join(exp.get("dbxrefs", [])),
                    "audit_error": audit_categories(exp, "ERROR"),
                    "audit_not_compliant": audit_categories(exp, "NOT_COMPLIANT"),
                    "audit_warning": audit_categories(exp, "WARNING"),
                    "url": BASE + exp.get("@id", ""),
                }
            )

    frame = pd.DataFrame(rows).sort_values(
        ["biosample_term", "assay", "target", "experiment"], kind="stable"
    )
    frame.to_csv(OUT / "encode_lung_experiments.tsv", sep="\t", index=False)

    summary = (
        frame.groupby(
            ["biosample_term", "biosample_classification", "assay", "target"],
            dropna=False,
        )
        .agg(
            n_experiments=("experiment", "nunique"),
            n_biosamples=("biosample_accessions", "nunique"),
            n_with_error=("audit_error", lambda x: sum(bool(v) for v in x)),
            n_not_compliant=("audit_not_compliant", lambda x: sum(bool(v) for v in x)),
        )
        .reset_index()
    )
    summary.to_csv(OUT / "encode_lung_experiment_summary.tsv", sep="\t", index=False)

    print(f"experiments\t{len(frame)}")
    print(f"biosample terms\t{frame['biosample_term'].nunique()}")
    print("assay counts")
    for key, count in Counter(zip(frame["assay"], frame["target"])).items():
        print(f"  {key[0]} {key[1] or 'accessibility'}\t{count}")
    print("terms with accessibility and both histone marks")
    present = frame.assign(mark=frame["target"].where(frame["target"].ne(""), frame["assay"]))
    for term, group in present.groupby("biosample_term"):
        marks = set(group["mark"])
        if (marks & {"ATAC-seq", "DNase-seq"}) and {"H3K27ac", "H3K27me3"}.issubset(marks):
            print(f"  {term}\t{','.join(sorted(marks))}")


if __name__ == "__main__":
    main()
