#!/usr/bin/env python3
"""Identify non-COPD disease links for replicated COPD genes in Open Targets.

The analysis uses the local Open Targets Platform direct target-disease
association tables downloaded from the official ``latest`` endpoint on
2026-09-21.  Only MONDO disease entities are retained, excluding COPD itself,
all COPD ancestors, all COPD descendants, and disease labels containing an
explicit COPD/emphysema/chronic-bronchitis/bronchiectasis term.  The lexical
safeguard catches related compound diseases that are not represented inside
the Open Targets COPD ontology closure.  This deliberately removes COPD
subtypes and circular umbrella terms and excludes quantitative traits and
generic phenotypes from the primary disease-link table.

Association scores aggregate heterogeneous evidence and are ranking scores,
not p-values or proof of causal pleiotropy.  All direct links are exported;
the compact top table requires at least two evidence records.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd
import pyarrow
import pyarrow.dataset as ds
import pyarrow.parquet as pq


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
REPO_ROOT = SECTION_ROOT.parents[2]
SHARED_DATA = REPO_ROOT / "data"
RESULTS = SECTION_ROOT / "results"
LOGS = SECTION_ROOT / "logs"

COPD_ID = "MONDO_0005002"
OPEN_TARGETS_ENDPOINT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/"
OPEN_TARGETS_LOCAL_DATE = "2026-09-21"
COPD_RELATED_LABEL_PATTERN = (
    r"\bCOPD\b|chronic obstructive|\bemphysema\b|chronic bronchitis|bronchiectasis"
)

HUMAN_GENETIC_OR_CURATED_SOURCES = {
    "clingen",
    "eva",
    "eva_somatic",
    "gene2phenotype",
    "gene_burden",
    "genomics_england",
    "gwas_credible_sets",
    "orphanet",
    "uniprot_literature",
    "uniprot_variants",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes",
        type=Path,
        default=RESULTS / "COPD-S2-R007A_locus_conservation_per_gene.tsv",
    )
    parser.add_argument(
        "--open-targets-dir",
        type=Path,
        default=SHARED_DATA / "open_targets",
    )
    parser.add_argument("--top-per-gene", type=int, default=5)
    parser.add_argument("--minimum-top-evidence", type=int, default=2)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument(
        "--log", type=Path, default=LOGS / "COPD-S2-other_disease_links.log"
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2-other-diseases")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ\t%(levelname)s\t%(message)s")
    formatter.converter = lambda *args: datetime.now(timezone.utc).timetuple()
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def list_values(value: object) -> list[str]:
    if value is None:
        return []
    try:
        return [str(item) for item in value if item is not None]
    except TypeError:
        return []


def datasource_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (target, disease), group in frame.groupby(
        ["targetId", "diseaseId"], sort=False
    ):
        sources = sorted(set(group["aggregationValue"].dropna().astype(str)))
        details: list[str] = []
        for source in sources:
            subset = group.loc[group["aggregationValue"].eq(source)]
            details.append(
                f"{source}:score={subset['associationScore'].max():.6g},"
                f"evidence={int(subset['evidenceCount'].sum())}"
            )
        genetic_sources = sorted(set(sources) & HUMAN_GENETIC_OR_CURATED_SOURCES)
        rows.append(
            {
                "targetId": target,
                "diseaseId": disease,
                "datasources": ";".join(sources),
                "n_datasources": len(sources),
                "datasource_details": ";".join(details),
                "human_genetic_or_curated_datasources": ";".join(genetic_sources),
                "has_human_genetic_or_curated_source": bool(genetic_sources),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    if args.top_per_gene < 1 or args.minimum_top_evidence < 1:
        raise ValueError("top-per-gene and minimum evidence must be positive")
    logger = setup_logging(args.log)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(args.genes, sep="\t", low_memory=False)
    require_columns(
        genes, ["gene", "gencode_gene_id", "gene_type"], "replicated genes"
    )
    if len(genes) != 140 or genes["gencode_gene_id"].duplicated().any():
        raise AssertionError("expected 140 unique unambiguous replicated genes")
    target_ids = sorted(set(genes["gencode_gene_id"].astype(str)))

    disease_path = args.open_targets_dir / "disease" / "disease.parquet"
    disease_columns = [
        "id",
        "name",
        "description",
        "parents",
        "ancestors",
        "descendants",
        "therapeuticAreas",
    ]
    diseases = pq.read_table(disease_path, columns=disease_columns).to_pandas()
    if diseases["id"].duplicated().any():
        raise AssertionError("Open Targets disease table has duplicate IDs")
    copd_rows = diseases.loc[diseases["id"].eq(COPD_ID)]
    if len(copd_rows) != 1:
        raise AssertionError(f"expected exactly one Open Targets COPD row, found {len(copd_rows)}")
    copd = copd_rows.iloc[0]
    ancestors = set(list_values(copd["ancestors"]))
    descendants = set(list_values(copd["descendants"]))
    circular_ids = {COPD_ID} | ancestors | descendants
    lexical_safeguard_ids = set(
        diseases.loc[
            diseases["name"].str.contains(
                COPD_RELATED_LABEL_PATTERN, case=False, regex=True, na=False
            ),
            "id",
        ].astype(str)
    ) - circular_ids
    excluded_ids = circular_ids | lexical_safeguard_ids

    overall_dir = args.open_targets_dir / "association_overall_direct"
    overall_dataset = ds.dataset(str(overall_dir), format="parquet")
    overall = overall_dataset.to_table(
        columns=[
            "diseaseId",
            "targetId",
            "aggregationType",
            "aggregationValue",
            "associationScore",
            "evidenceCount",
        ],
        filter=ds.field("targetId").isin(target_ids),
    ).to_pandas()
    if not overall["aggregationType"].eq("overall").all():
        raise AssertionError("unexpected aggregation type in overall direct table")
    if overall.duplicated(["targetId", "diseaseId"]).any():
        raise AssertionError("overall direct table is not unique by target-disease pair")

    circular_associations = overall.loc[overall["diseaseId"].isin(excluded_ids)].copy()
    mondo = overall.loc[
        overall["diseaseId"].str.startswith("MONDO_", na=False)
        & ~overall["diseaseId"].isin(excluded_ids)
        & overall["associationScore"].gt(0)
        & overall["evidenceCount"].gt(0)
    ].copy()
    disease_lookup = diseases[
        ["id", "name", "description", "parents", "therapeuticAreas"]
    ].rename(columns={"id": "diseaseId", "name": "disease_name"})
    links = mondo.merge(disease_lookup, on="diseaseId", how="left", validate="many_to_one")
    if links["disease_name"].isna().any():
        raise AssertionError("MONDO association lacks Open Targets disease metadata")

    datasource_dir = args.open_targets_dir / "association_by_datasource_direct"
    datasource_dataset = ds.dataset(str(datasource_dir), format="parquet")
    by_source = datasource_dataset.to_table(
        columns=[
            "diseaseId",
            "targetId",
            "aggregationType",
            "aggregationValue",
            "associationScore",
            "evidenceCount",
        ],
        filter=ds.field("targetId").isin(target_ids),
    ).to_pandas()
    if not by_source["aggregationType"].eq("datasourceId").all():
        raise AssertionError("unexpected aggregation type in datasource direct table")
    selected_pairs = links[["targetId", "diseaseId"]].drop_duplicates()
    by_source = by_source.merge(
        selected_pairs, on=["targetId", "diseaseId"], how="inner", validate="many_to_one"
    )
    source_summary = datasource_summary(by_source)
    links = links.merge(
        source_summary,
        on=["targetId", "diseaseId"],
        how="left",
        validate="one_to_one",
    )
    links = links.merge(
        genes[["gencode_gene_id", "gene", "gene_type"]],
        left_on="targetId",
        right_on="gencode_gene_id",
        how="left",
        validate="many_to_one",
    )
    links["disease_parents"] = links["parents"].map(
        lambda value: ";".join(sorted(list_values(value)))
    )
    links["therapeutic_areas"] = links["therapeuticAreas"].map(
        lambda value: ";".join(sorted(list_values(value)))
    )
    links = links.drop(columns=["parents", "therapeuticAreas", "gencode_gene_id"])
    links.insert(0, "phenotype_scope", "core_copd")
    links["open_targets_table"] = "association_overall_direct"
    links["open_targets_local_access_date"] = OPEN_TARGETS_LOCAL_DATE
    links["link_interpretation"] = (
        "direct Open Targets target-disease evidence; association score is not a p-value"
    )
    links = links.sort_values(
        ["gene", "associationScore", "evidenceCount", "disease_name"],
        ascending=[True, False, False, True],
        kind="stable",
    ).reset_index(drop=True)

    eligible_top = links.loc[
        links["evidenceCount"].ge(args.minimum_top_evidence)
    ].copy()
    eligible_top["rank_within_gene"] = (
        eligible_top.groupby("gene", sort=False).cumcount() + 1
    )
    top = eligible_top.loc[
        eligible_top["rank_within_gene"].le(args.top_per_gene)
    ].copy()

    summary_rows: list[dict[str, object]] = []
    grouped = {key: group for key, group in links.groupby("targetId", sort=False)}
    for gene_row in genes[["gene", "gencode_gene_id", "gene_type"]].itertuples(index=False):
        group = grouped.get(gene_row.gencode_gene_id)
        if group is None or group.empty:
            summary_rows.append(
                {
                    "phenotype_scope": "core_copd",
                    "gene": gene_row.gene,
                    "gencode_gene_id": gene_row.gencode_gene_id,
                    "gene_type": gene_row.gene_type,
                    "has_other_mondo_disease_link": False,
                    "n_other_mondo_diseases_any_direct_evidence": 0,
                    "n_other_mondo_diseases_at_least_2_evidence": 0,
                    "n_other_mondo_diseases_human_genetic_or_curated": 0,
                    "top_disease_id": "",
                    "top_disease_name": "",
                    "top_association_score": pd.NA,
                    "top_evidence_count": pd.NA,
                    "top_datasources": "",
                }
            )
            continue
        strongest = group.iloc[0]
        summary_rows.append(
            {
                "phenotype_scope": "core_copd",
                "gene": gene_row.gene,
                "gencode_gene_id": gene_row.gencode_gene_id,
                "gene_type": gene_row.gene_type,
                "has_other_mondo_disease_link": True,
                "n_other_mondo_diseases_any_direct_evidence": group["diseaseId"].nunique(),
                "n_other_mondo_diseases_at_least_2_evidence": group.loc[
                    group["evidenceCount"].ge(2), "diseaseId"
                ].nunique(),
                "n_other_mondo_diseases_human_genetic_or_curated": group.loc[
                    group["has_human_genetic_or_curated_source"].fillna(False), "diseaseId"
                ].nunique(),
                "top_disease_id": strongest["diseaseId"],
                "top_disease_name": strongest["disease_name"],
                "top_association_score": strongest["associationScore"],
                "top_evidence_count": strongest["evidenceCount"],
                "top_datasources": strongest["datasources"],
            }
        )
    gene_summary = pd.DataFrame(summary_rows).sort_values(
        [
            "n_other_mondo_diseases_human_genetic_or_curated",
            "n_other_mondo_diseases_at_least_2_evidence",
            "top_association_score",
            "gene",
        ],
        ascending=[False, False, False, True],
        na_position="last",
        kind="stable",
    )

    closure_rows: list[dict[str, object]] = []
    closure_lookup = diseases.loc[diseases["id"].isin(excluded_ids), ["id", "name"]]
    for row in closure_lookup.itertuples(index=False):
        if row.id == COPD_ID:
            reason = "COPD index term"
        elif row.id in descendants:
            reason = "COPD ontology descendant"
        elif row.id in lexical_safeguard_ids:
            reason = "COPD-related lexical safeguard outside ontology closure"
        else:
            reason = "COPD ontology ancestor"
        closure_rows.append(
            {
                "excluded_disease_id": row.id,
                "excluded_disease_name": row.name,
                "exclusion_reason": reason,
                "n_target_association_rows_removed": int(
                    circular_associations["diseaseId"].eq(row.id).sum()
                ),
            }
        )
    closure = pd.DataFrame(closure_rows).sort_values(
        ["exclusion_reason", "excluded_disease_id"], kind="stable"
    )

    genes_any = int(gene_summary["has_other_mondo_disease_link"].sum())
    genes_multiple = int(
        gene_summary["n_other_mondo_diseases_at_least_2_evidence"].gt(0).sum()
    )
    genes_genetic = int(
        gene_summary["n_other_mondo_diseases_human_genetic_or_curated"].gt(0).sum()
    )
    overall_summary = pd.DataFrame(
        [
            ("replicated_unambiguous_genes", len(genes), "R007 analysis set"),
            (
                "genes_with_any_other_mondo_direct_link",
                genes_any,
                "At least one positive direct Open Targets association after circular exclusion",
            ),
            (
                "genes_with_other_mondo_link_at_least_2_evidence",
                genes_multiple,
                "At least one retained disease with evidenceCount >= 2",
            ),
            (
                "genes_with_human_genetic_or_curated_other_disease_source",
                genes_genetic,
                "At least one retained disease supported by a predefined human genetic/curated source",
            ),
            ("unique_other_mondo_diseases", links["diseaseId"].nunique(), "Across all 140 genes"),
            ("direct_gene_disease_links", len(links), "Unique targetId-diseaseId pairs"),
            (
                "direct_links_at_least_2_evidence",
                int(links["evidenceCount"].ge(2).sum()),
                "Compact top-table eligibility threshold",
            ),
            (
                "circular_target_disease_rows_removed",
                len(circular_associations),
                "COPD ontology closure plus explicit COPD-related lexical safeguards",
            ),
            (
                "excluded_copd_related_disease_terms",
                len(excluded_ids),
                "Ontology closure IDs plus compound disease labels caught lexically",
            ),
            (
                "lexical_safeguard_disease_terms_outside_ontology_closure",
                len(lexical_safeguard_ids),
                "Labels matching COPD, chronic obstructive, emphysema, chronic bronchitis, or bronchiectasis",
            ),
        ],
        columns=["metric", "value", "definition"],
    )
    overall_summary.insert(0, "phenotype_scope", "core_copd")
    overall_summary["open_targets_local_access_date"] = OPEN_TARGETS_LOCAL_DATE

    gene_summary_path = args.results_dir / "COPD-S2-R008A_other_disease_gene_summary.tsv"
    top_path = args.results_dir / "COPD-S2-R008B_top_other_disease_links.tsv"
    all_path = args.results_dir / "COPD-S2-R008C_all_other_disease_links.tsv.gz"
    closure_path = args.results_dir / "COPD-S2-R008D_circular_disease_exclusions.tsv"
    summary_path = args.results_dir / "COPD-S2-R008E_summary.tsv"
    gene_summary.to_csv(gene_summary_path, sep="\t", index=False)
    top.to_csv(top_path, sep="\t", index=False)
    links.to_csv(
        all_path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    closure.to_csv(closure_path, sep="\t", index=False)
    overall_summary.to_csv(summary_path, sep="\t", index=False)

    output_rows = {
        gene_summary_path: len(gene_summary),
        top_path: len(top),
        all_path: len(links),
        closure_path: len(closure),
        summary_path: len(overall_summary),
    }
    overall_files = sorted(overall_dir.glob("*.parquet"))
    datasource_files = sorted(datasource_dir.glob("*.parquet"))
    manifest_rows: list[tuple[str, object]] = [
        ("run_utc", datetime.now(timezone.utc).isoformat()),
        ("python", platform.python_version()),
        ("pandas", pd.__version__),
        ("pyarrow", pyarrow.__version__),
        ("script", str(SCRIPT)),
        ("script_sha256", sha256_file(SCRIPT)),
        ("gene_input", str(args.genes.resolve())),
        ("gene_input_sha256", sha256_file(args.genes)),
        ("open_targets_endpoint", OPEN_TARGETS_ENDPOINT),
        ("open_targets_local_access_date", OPEN_TARGETS_LOCAL_DATE),
        ("disease_table", str(disease_path.resolve())),
        ("disease_table_sha256", sha256_file(disease_path)),
        ("overall_direct_parquet_files", len(overall_files)),
        ("overall_direct_total_size_bytes", sum(path.stat().st_size for path in overall_files)),
        ("datasource_direct_parquet_files", len(datasource_files)),
        (
            "datasource_direct_total_size_bytes",
            sum(path.stat().st_size for path in datasource_files),
        ),
        ("copd_ontology_id", COPD_ID),
        ("n_circular_ontology_ids", len(circular_ids)),
        ("copd_related_label_pattern", COPD_RELATED_LABEL_PATTERN),
        ("n_lexical_safeguard_ids_outside_ontology_closure", len(lexical_safeguard_ids)),
        ("n_total_excluded_copd_related_ids", len(excluded_ids)),
        ("primary_disease_identifier_filter", "MONDO_"),
        ("top_per_gene", args.top_per_gene),
        ("minimum_top_evidence", args.minimum_top_evidence),
        (
            "human_genetic_or_curated_sources",
            ";".join(sorted(HUMAN_GENETIC_OR_CURATED_SOURCES)),
        ),
    ]
    for path, row_count in output_rows.items():
        manifest_rows.extend(
            [
                (f"output_{path.name}_rows", row_count),
                (f"output_{path.name}_sha256", sha256_file(path)),
            ]
        )
    pd.DataFrame(manifest_rows, columns=["field", "value"]).to_csv(
        LOGS / "COPD-S2-other_disease_links_manifest.tsv", sep="\t", index=False
    )

    logger.info(
        "Other-disease links: %d pairs, %d MONDO diseases, %d/140 genes",
        len(links),
        links["diseaseId"].nunique(),
        genes_any,
    )
    logger.info(
        "At least 2 evidence records: %d links across %d genes; human genetic/curated: %d genes",
        int(links["evidenceCount"].ge(2).sum()),
        genes_multiple,
        genes_genetic,
    )
    logger.info(
        "Removed %d circular association rows spanning the COPD ontology closure and lexical safeguard",
        len(circular_associations),
    )
    logger.info("Wrote COPD-S2-R008 outputs to %s", args.results_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
