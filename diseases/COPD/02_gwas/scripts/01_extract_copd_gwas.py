#!/usr/bin/env python3
"""Extract auditable COPD GWAS Catalog subsets from the local release.

Core estimates use only records whose primary mapped trait is chronic
obstructive pulmonary disease (MONDO:0005002 or legacy EFO:0000341).
COPD subtypes and records that merely mention COPD are exported separately and
must not be merged into core counts.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
CATALOG = ROOT / "data" / "gwas_catalog"
OUT = Path(__file__).resolve().parents[1] / "data"
OUT.mkdir(parents=True, exist_ok=True)

MONDO = "http://purl.obolibrary.org/obo/MONDO_"
EFO = "http://www.ebi.ac.uk/efo/EFO_"

CORE_URIS = {
    MONDO + "0005002",  # chronic obstructive pulmonary disease
    EFO + "0000341",  # obsolete/legacy COPD URI retained by older records
}

# These clinically interpretable descendants are reported separately. We do
# not automatically treat every ontology descendant as COPD because the EFO
# release also places bronchiectasis below COPD, despite its distinct clinical
# and genetic definition.
SUBTYPE_URIS = {
    MONDO + "0004849",  # pulmonary emphysema
    EFO + "0000464",
    MONDO + "0005607",  # chronic bronchitis
    EFO + "0006505",
    MONDO + "0011751",  # COPD, severe early onset
}

COPD_TEXT = re.compile(
    r"\bCOPD\b|chronic\s+obstructive\s+(?:pulmonary|lung|airways?)\s+disease",
    flags=re.IGNORECASE,
)


def uri_set(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {part.strip() for part in str(value).split(",") if part.strip()}


def primary_category(uri_value: object) -> str:
    uris = uri_set(uri_value)
    if uris and uris.issubset(CORE_URIS):
        return "core_copd"
    if uris and uris.issubset(SUBTYPE_URIS):
        return "copd_subtype"
    if uris & (CORE_URIS | SUBTYPE_URIS):
        return "copd_related"
    return ""


def mentions_copd(frame: pd.DataFrame) -> pd.Series:
    text_columns = [
        col
        for col in ("DISEASE/TRAIT", "STUDY", "MAPPED BACKGROUND TRAIT")
        if col in frame.columns
    ]
    if not text_columns:
        return pd.Series(False, index=frame.index)
    combined = frame[text_columns].fillna("").astype(str).agg(" ".join, axis=1)
    return combined.str.contains(COPD_TEXT, na=False)


def extract_associations() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = CATALOG / "gwas-catalog-download-associations-alt-full.tsv"
    selected: list[pd.DataFrame] = []
    related: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=100_000):
        chunk["copd_category"] = chunk["MAPPED_TRAIT_URI"].map(primary_category)
        direct_mask = chunk["copd_category"].isin(["core_copd", "copd_subtype"])
        if direct_mask.any():
            selected.append(chunk.loc[direct_mask].copy())
        related_mask = chunk["copd_category"].eq("copd_related") | (
            chunk["copd_category"].eq("") & mentions_copd(chunk)
        )
        if related_mask.any():
            rel = chunk.loc[related_mask].copy()
            rel["copd_category"] = "copd_related"
            related.append(rel)

    direct = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()
    rel = pd.concat(related, ignore_index=True) if related else pd.DataFrame()
    core = direct.loc[direct["copd_category"].eq("core_copd")].copy()
    subtype = direct.loc[direct["copd_category"].eq("copd_subtype")].copy()
    core.to_csv(OUT / "copd_associations.tsv", sep="\t", index=False)
    subtype.to_csv(OUT / "copd_subtype_associations.tsv", sep="\t", index=False)
    rel.to_csv(OUT / "copd_related_associations.tsv", sep="\t", index=False)
    return core, subtype, rel


def extract_studies() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = CATALOG / "gwas-catalog-download-studies-v1.0.3.1.txt"
    studies = pd.read_csv(path, sep="\t", low_memory=False)
    studies["copd_category"] = studies["MAPPED_TRAIT_URI"].map(primary_category)
    related_mask = studies["copd_category"].eq("copd_related") | (
        studies["copd_category"].eq("") & mentions_copd(studies)
    )
    studies.loc[related_mask, "copd_category"] = "copd_related"
    selected = studies.loc[studies["copd_category"].ne("")].copy()
    selected.to_csv(OUT / "copd_studies.tsv", sep="\t", index=False)

    ancestry_path = CATALOG / "gwas-catalog-download-ancestries-v1.0.3.1.txt"
    ancestry = pd.read_csv(ancestry_path, sep="\t", low_memory=False)
    keep = selected[["STUDY ACCESSION", "copd_category"]].drop_duplicates()
    selected_ancestry = ancestry.loc[
        ancestry["STUDY ACCESSION"].isin(set(keep["STUDY ACCESSION"]))
    ].merge(keep, on="STUDY ACCESSION", how="left")
    selected_ancestry.to_csv(OUT / "copd_ancestries.tsv", sep="\t", index=False)
    return selected, selected_ancestry


def write_ontology_audit() -> None:
    rows = [
        (MONDO + "0005002", "core_copd", "current COPD term"),
        (EFO + "0000341", "core_copd", "legacy/obsolete COPD URI"),
        (MONDO + "0004849", "copd_subtype", "pulmonary emphysema"),
        (EFO + "0000464", "copd_subtype", "legacy pulmonary emphysema URI"),
        (MONDO + "0005607", "copd_subtype", "chronic bronchitis"),
        (EFO + "0006505", "copd_subtype", "legacy chronic bronchitis URI"),
        (MONDO + "0011751", "copd_subtype", "COPD, severe early onset"),
    ]
    pd.DataFrame(rows, columns=["trait_uri", "category", "rationale"]).to_csv(
        OUT / "copd_ontology_audit.tsv", sep="\t", index=False
    )


def main() -> None:
    write_ontology_audit()
    core, subtype, related = extract_associations()
    studies, ancestry = extract_studies()
    core_studies = studies.loc[studies["copd_category"].eq("core_copd")]

    print(f"core association rows\t{len(core)}")
    print(f"core unique study accessions\t{core_studies['STUDY ACCESSION'].nunique()}")
    print(f"core unique publications\t{core_studies['PUBMED ID'].nunique()}")
    print(f"core unique SNPS values\t{core['SNPS'].dropna().nunique()}")
    print(f"subtype association rows\t{len(subtype)}")
    print(f"related association rows\t{len(related)}")
    print(f"selected ancestry rows\t{len(ancestry)}")
    print("study categories")
    print(studies["copd_category"].value_counts().to_string())
    print("core mapped traits")
    print(core["MAPPED_TRAIT"].value_counts().to_string())


if __name__ == "__main__":
    main()
