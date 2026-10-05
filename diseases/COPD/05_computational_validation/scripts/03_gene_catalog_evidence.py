#!/usr/bin/env python3
"""Resolve candidate-linked genes in HGNC and query exact COPD Open Targets evidence.

Open Targets evidence is restricted to the exact MONDO_0005002 disease node.
It is contextual gene-disease evidence and is never labeled variant validation.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import pyarrow
import pyarrow.parquet as pq

from common import (
    COPD,
    LOGS,
    RESULTS,
    ROOT,
    candidate_core,
    candidate_sort_key,
    ensure_directories,
    file_record,
    integer,
    load_candidates,
    numeric,
    read_tsv,
    tokens,
    write_manifest,
    write_tsv,
)


RESULT_ID = "COPD-S5-R003"
DISEASE_ID = "MONDO_0005002"
DISEASE_NAME = "chronic obstructive pulmonary disease"
HGNC = ROOT / "data/hgnc/hgnc_complete_set.txt"
OPEN_TARGETS = ROOT / "data/open_targets"
TARGET_SUMMARY = (
    COPD / "04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv"
)
TARGET_EVIDENCE = (
    COPD / "04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv"
)

OUT_RESOLUTION = RESULTS / f"{RESULT_ID}_HGNC_gene_resolution.tsv"
OUT_ASSOCIATIONS = RESULTS / f"{RESULT_ID}_OpenTargets_COPD_gene_associations.tsv"
OUT_DATASOURCES = RESULTS / f"{RESULT_ID}_OpenTargets_COPD_datasources.tsv"
OUT_CONTEXT = RESULTS / f"{RESULT_ID}_candidate_gene_context.tsv.gz"
OUT_SUMMARY = RESULTS / f"{RESULT_ID}_gene_catalog_summary.tsv"
OUT_MANIFEST = RESULTS / f"{RESULT_ID}_manifest.json"
LOG = LOGS / "03_gene_catalog_evidence.log"

GENE_SPLIT = re.compile(r"\s+-\s+|,\s*")
ENSEMBL_RE = re.compile(r"^ENSG\d+(?:\.\d+)?$")

TARGET_COLUMNS = (
    "nearest_gene",
    "nearest_protein_coding_gene",
    "genes_with_tss_within_100kb",
    "coding_CDS_gene_names",
    "target_gene",
    "target_genes",
    "mapped_target_gene",
    "mapped_target_genes",
    "putative_target_gene",
    "putative_target_genes",
    "all_target_genes",
    "proximity_genes",
    "gtex_lung_eqtl_genes",
    "assigned_genes",
)

RESOLUTION_FIELDS = [
    "input_gene_token",
    "resolution_status",
    "resolution_basis",
    "approved_symbol",
    "hgnc_id",
    "approved_name",
    "locus_group",
    "locus_type",
    "chromosomal_location",
    "ensembl_gene_id",
    "entrez_id",
    "omim_id",
]

ASSOCIATION_FIELDS = [
    "approved_symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "open_targets_disease_id",
    "open_targets_disease_name",
    "exact_disease_node",
    "overall_association_present",
    "overall_association_score",
    "overall_evidence_count",
    "datasource_count",
    "datasources",
    "datasource_evidence_count_sum",
    "interpretation",
]

DATASOURCE_FIELDS = [
    "approved_symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "open_targets_disease_id",
    "datasource_id",
    "association_score",
    "evidence_count",
]

CONTEXT_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "variant_class",
    "candidate_rsids",
    "input_gene_token",
    "mapping_origins",
    "gene_mapping_scope",
    "resolution_status",
    "approved_symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "open_targets_exact_COPD_association_present",
    "open_targets_overall_score",
    "open_targets_evidence_count",
    "open_targets_datasources",
]


def split_gene_field(value: object) -> list[str]:
    result: set[str] = set()
    for item in tokens(value):
        for part in GENE_SPLIT.split(item):
            cleaned = part.strip().strip('"')
            if cleaned and cleaned.lower() not in {
                "intergenic",
                "not reported",
                "unknown",
                "na",
            }:
                result.add(cleaned)
    return sorted(result)


def candidate_gene_links(
    candidates: list[dict[str, str]],
) -> tuple[dict[str, dict[str, set[str]]], dict[str, object]]:
    links: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for candidate in candidates:
        candidate_id = candidate["candidate_record_id"]
        for gene in split_gene_field(candidate.get("linked_gwas_genes", "")):
            links[candidate_id][gene].add("GWAS_Catalog_mapped_gene")
        for gene in split_gene_field(candidate.get("assigned_selected_loci", "")):
            links[candidate_id][gene].add("selected_GWAS_gene_locus")

    target_summary_status: dict[str, object] = {
        "path": str(TARGET_SUMMARY),
        "status": "not_available",
        "recognized_target_columns": [],
        "candidate_gene_links_added": 0,
    }
    if TARGET_SUMMARY.exists():
        target_rows = read_tsv(TARGET_SUMMARY)
        recognized = [
            column
            for column in TARGET_COLUMNS
            if target_rows and column in target_rows[0]
        ]
        added = 0
        for row in target_rows:
            candidate_id = row.get("candidate_record_id", "")
            if candidate_id not in links and not any(
                candidate_id == candidate["candidate_record_id"]
                for candidate in candidates
            ):
                continue
            for column in recognized:
                for gene in split_gene_field(row.get(column, "")):
                    before = len(links[candidate_id][gene])
                    links[candidate_id][gene].add(
                        f"provisional_Section4_target:{column}"
                    )
                    if len(links[candidate_id][gene]) > before:
                        added += 1
        target_summary_status = {
            "path": str(TARGET_SUMMARY),
            "status": "consumed" if recognized else "schema_not_recognized",
            "recognized_target_columns": recognized,
            "candidate_gene_links_added": added,
            "input": file_record(TARGET_SUMMARY),
        }
        if TARGET_EVIDENCE.exists():
            evidence_rows = read_tsv(TARGET_EVIDENCE)
            target_summary_status["evidence_path"] = str(TARGET_EVIDENCE)
            target_summary_status["evidence_rows"] = len(evidence_rows)
            target_summary_status["mapping_methods"] = sorted(
                {row.get("mapping_method", "") for row in evidence_rows}
                - {""}
            )
            target_summary_status["evidence_input"] = file_record(TARGET_EVIDENCE)
    return links, target_summary_status


def load_hgnc():
    rows: list[dict[str, str]] = []
    by_symbol: dict[str, dict[str, str]] = {}
    by_ensembl: dict[str, dict[str, str]] = {}
    aliases: dict[str, list[dict[str, str]]] = defaultdict(list)
    with HGNC.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
            by_symbol[row["symbol"].upper()] = row
            ensembl = row.get("ensembl_gene_id", "").split(".")[0]
            if ensembl:
                by_ensembl[ensembl] = row
            for column in ("alias_symbol", "prev_symbol"):
                for alias in row.get(column, "").split("|"):
                    if alias.strip():
                        aliases[alias.strip().upper()].append(row)
    return rows, by_symbol, by_ensembl, aliases


def resolve_gene(
    token: str,
    by_symbol: dict[str, dict[str, str]],
    by_ensembl: dict[str, dict[str, str]],
    aliases: dict[str, list[dict[str, str]]],
) -> dict[str, object]:
    key = token.upper()
    record = None
    status = "unresolved"
    basis = "none"
    if ENSEMBL_RE.fullmatch(token):
        record = by_ensembl.get(token.split(".")[0])
        if record:
            status = "resolved"
            basis = "exact_Ensembl_gene_id"
    if record is None and key in by_symbol:
        record = by_symbol[key]
        status = "resolved"
        basis = "exact_HGNC_approved_symbol"
    if record is None and key in aliases:
        possible = {row["hgnc_id"]: row for row in aliases[key]}
        if len(possible) == 1:
            record = next(iter(possible.values()))
            status = "resolved"
            basis = "unambiguous_HGNC_alias_or_previous_symbol"
        else:
            status = "ambiguous_alias"
            basis = ";".join(sorted(possible))
    return {
        "input_gene_token": token,
        "resolution_status": status,
        "resolution_basis": basis,
        "approved_symbol": record.get("symbol", "") if record else "",
        "hgnc_id": record.get("hgnc_id", "") if record else "",
        "approved_name": record.get("name", "") if record else "",
        "locus_group": record.get("locus_group", "") if record else "",
        "locus_type": record.get("locus_type", "") if record else "",
        "chromosomal_location": record.get("location", "") if record else "",
        "ensembl_gene_id": (
            record.get("ensembl_gene_id", "").split(".")[0] if record else ""
        ),
        "entrez_id": record.get("entrez_id", "") if record else "",
        "omim_id": record.get("omim_id", "") if record else "",
    }


def open_targets_rows(
    folder: Path,
    columns: list[str],
    target_ids: set[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(folder.glob("*.parquet")):
        table = pq.read_table(
            path,
            columns=columns,
            filters=[("diseaseId", "=", DISEASE_ID)],
        )
        values = table.to_pydict()
        for index, target_id in enumerate(values["targetId"]):
            if target_id not in target_ids:
                continue
            rows.append({column: values[column][index] for column in columns})
    return rows


def verify_disease_node() -> dict[str, object]:
    table = pq.read_table(
        OPEN_TARGETS / "disease/disease.parquet",
        columns=["id", "name", "code", "exactSynonyms", "descendants"],
        filters=[("id", "=", DISEASE_ID)],
    )
    values = table.to_pydict()
    if table.num_rows != 1:
        raise ValueError(
            f"expected one Open Targets disease row for {DISEASE_ID}, got {table.num_rows}"
        )
    return {column: values[column][0] for column in values}


def main() -> None:
    ensure_directories()
    candidate_input, candidates = load_candidates()
    candidates.sort(key=candidate_sort_key)
    candidates_by_id = {row["candidate_record_id"]: row for row in candidates}
    links, target_summary_status = candidate_gene_links(candidates)
    unique_tokens = sorted({gene for genes in links.values() for gene in genes})

    _, by_symbol, by_ensembl, aliases = load_hgnc()
    resolution_rows = [
        resolve_gene(token, by_symbol, by_ensembl, aliases)
        for token in unique_tokens
    ]
    resolution_by_token = {row["input_gene_token"]: row for row in resolution_rows}
    target_ids = {
        str(row["ensembl_gene_id"])
        for row in resolution_rows
        if row["resolution_status"] == "resolved" and row["ensembl_gene_id"]
    }

    disease = verify_disease_node()
    overall_rows = open_targets_rows(
        OPEN_TARGETS / "association_overall_direct",
        [
            "diseaseId",
            "targetId",
            "aggregationType",
            "aggregationValue",
            "associationScore",
            "evidenceCount",
        ],
        target_ids,
    )
    datasource_rows_raw = open_targets_rows(
        OPEN_TARGETS / "association_by_datasource_direct",
        [
            "diseaseId",
            "targetId",
            "aggregationType",
            "aggregationValue",
            "associationScore",
            "evidenceCount",
        ],
        target_ids,
    )
    overall_by_target = {str(row["targetId"]): row for row in overall_rows}
    datasources_by_target: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in datasource_rows_raw:
        datasources_by_target[str(row["targetId"])].append(row)

    resolved_by_target: dict[str, dict[str, object]] = {}
    for row in resolution_rows:
        if row["ensembl_gene_id"]:
            resolved_by_target[str(row["ensembl_gene_id"])] = row

    association_rows = []
    for target_id, resolution in sorted(
        resolved_by_target.items(), key=lambda item: str(item[1]["approved_symbol"])
    ):
        overall = overall_by_target.get(target_id)
        sources = sorted(
            datasources_by_target.get(target_id, []),
            key=lambda row: str(row["aggregationValue"]),
        )
        association_rows.append(
            {
                "approved_symbol": resolution["approved_symbol"],
                "hgnc_id": resolution["hgnc_id"],
                "ensembl_gene_id": target_id,
                "open_targets_disease_id": DISEASE_ID,
                "open_targets_disease_name": disease["name"],
                "exact_disease_node": True,
                "overall_association_present": overall is not None,
                "overall_association_score": (
                    overall["associationScore"] if overall else None
                ),
                "overall_evidence_count": overall["evidenceCount"] if overall else 0,
                "datasource_count": len(sources),
                "datasources": ";".join(
                    str(row["aggregationValue"]) for row in sources
                ),
                "datasource_evidence_count_sum": sum(
                    int(row["evidenceCount"] or 0) for row in sources
                ),
                "interpretation": (
                    "exact COPD gene-disease association in Open Targets; contextual "
                    "gene evidence, not variant-level validation"
                    if overall
                    else "no exact COPD association row for this resolved candidate-linked gene"
                ),
            }
        )
    association_by_target = {
        row["ensembl_gene_id"]: row for row in association_rows
    }

    datasource_rows = []
    for target_id, values in sorted(datasources_by_target.items()):
        resolution = resolved_by_target[target_id]
        for row in sorted(values, key=lambda item: str(item["aggregationValue"])):
            datasource_rows.append(
                {
                    "approved_symbol": resolution["approved_symbol"],
                    "hgnc_id": resolution["hgnc_id"],
                    "ensembl_gene_id": target_id,
                    "open_targets_disease_id": DISEASE_ID,
                    "datasource_id": row["aggregationValue"],
                    "association_score": row["associationScore"],
                    "evidence_count": row["evidenceCount"],
                }
            )

    context_rows = []
    for candidate_id, genes in links.items():
        candidate = candidates_by_id[candidate_id]
        for token, origins in genes.items():
            resolution = resolution_by_token[token]
            association = association_by_target.get(
                str(resolution["ensembl_gene_id"]), {}
            )
            context_rows.append(
                {
                    **candidate_core(candidate),
                    "input_gene_token": token,
                    "mapping_origins": ";".join(sorted(origins)),
                    "gene_mapping_scope": (
                        "provisional regulatory target hypothesis"
                        if any(origin.startswith("provisional_Section4_target") for origin in origins)
                        else "GWAS mapped gene or physically assigned selected locus; not a confirmed regulatory target"
                    ),
                    "resolution_status": resolution["resolution_status"],
                    "approved_symbol": resolution["approved_symbol"],
                    "hgnc_id": resolution["hgnc_id"],
                    "ensembl_gene_id": resolution["ensembl_gene_id"],
                    "open_targets_exact_COPD_association_present": association.get(
                        "overall_association_present", False
                    ),
                    "open_targets_overall_score": association.get(
                        "overall_association_score", ""
                    ),
                    "open_targets_evidence_count": association.get(
                        "overall_evidence_count", 0
                    ),
                    "open_targets_datasources": association.get("datasources", ""),
                }
            )
    context_rows.sort(
        key=lambda row: (
            integer(row["predicted_causal_priority_rank"]) or 10**12,
            str(row["candidate_record_id"]),
            str(row["approved_symbol"] or row["input_gene_token"]),
        )
    )

    resolved_count = sum(row["resolution_status"] == "resolved" for row in resolution_rows)
    associated_genes = sum(
        bool(row["overall_association_present"]) for row in association_rows
    )
    candidates_with_gene = {row["candidate_record_id"] for row in context_rows}
    candidates_with_ot = {
        row["candidate_record_id"]
        for row in context_rows
        if row["open_targets_exact_COPD_association_present"]
    }
    summary_rows = [
        {
            "metric": "predicted_causal_candidates",
            "n": len(candidates),
            "denominator": len(candidates),
            "fraction": 1.0,
            "note": "Section 4 predicted-causal union",
        },
        {
            "metric": "candidates_with_at_least_one_gene_mapping",
            "n": len(candidates_with_gene),
            "denominator": len(candidates),
            "fraction": len(candidates_with_gene) / len(candidates),
            "note": "union of GWAS mapped genes, selected loci, and available provisional targets",
        },
        {
            "metric": "unique_gene_tokens",
            "n": len(unique_tokens),
            "denominator": "",
            "fraction": "",
            "note": "compound GWAS strings split only on comma or space-hyphen-space",
        },
        {
            "metric": "HGNC_resolved_gene_tokens",
            "n": resolved_count,
            "denominator": len(unique_tokens),
            "fraction": resolved_count / len(unique_tokens) if unique_tokens else None,
            "note": "approved symbol, Ensembl ID, or unambiguous alias/previous symbol",
        },
        {
            "metric": "resolved_genes_with_exact_OpenTargets_COPD_association",
            "n": associated_genes,
            "denominator": len(association_rows),
            "fraction": associated_genes / len(association_rows) if association_rows else None,
            "note": f"exact disease node {DISEASE_ID}, direct associations",
        },
        {
            "metric": "candidates_with_at_least_one_OpenTargets_COPD_associated_gene",
            "n": len(candidates_with_ot),
            "denominator": len(candidates),
            "fraction": len(candidates_with_ot) / len(candidates),
            "note": "contextual gene evidence only, not variant-level validation",
        },
    ]

    write_tsv(OUT_RESOLUTION, resolution_rows, RESOLUTION_FIELDS)
    write_tsv(OUT_ASSOCIATIONS, association_rows, ASSOCIATION_FIELDS)
    write_tsv(OUT_DATASOURCES, datasource_rows, DATASOURCE_FIELDS)
    write_tsv(OUT_CONTEXT, context_rows, CONTEXT_FIELDS)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["metric", "n", "denominator", "fraction", "note"],
    )
    input_records = {
        "section4_candidates": file_record(candidate_input),
        "hgnc_complete_set": file_record(HGNC),
        "open_targets_disease": file_record(
            OPEN_TARGETS / "disease/disease.parquet"
        ),
    }
    if TARGET_SUMMARY.exists():
        input_records["section4_provisional_target_summary"] = file_record(
            TARGET_SUMMARY
        )
    if TARGET_EVIDENCE.exists():
        input_records["section4_provisional_target_evidence"] = file_record(
            TARGET_EVIDENCE
        )
    write_manifest(
        OUT_MANIFEST,
        {
            "result_id": RESULT_ID,
            "analysis": "HGNC resolution and exact COPD Open Targets gene evidence",
            "open_targets_release": '"latest" mirror accessed 2026-09-21',
            "disease_scope": {
                "id": DISEASE_ID,
                "name": disease["name"],
                "code": disease["code"],
                "descendants_not_included": True,
            },
            "candidate_count": len(candidates),
            "unique_gene_tokens": len(unique_tokens),
            "resolved_gene_tokens": resolved_count,
            "resolved_unique_ensembl_targets": len(association_rows),
            "open_targets_associated_genes": associated_genes,
            "target_summary": target_summary_status,
            "software": {
                "python": __import__("platform").python_version(),
                "pyarrow": pyarrow.__version__,
            },
            "inputs": input_records,
            "open_targets_input_partitions": {
                "overall": [
                    file_record(path)
                    for path in sorted(
                        (OPEN_TARGETS / "association_overall_direct").glob(
                            "*.parquet"
                        )
                    )
                ],
                "by_datasource": [
                    file_record(path)
                    for path in sorted(
                        (OPEN_TARGETS / "association_by_datasource_direct").glob(
                            "*.parquet"
                        )
                    )
                ],
            },
            "outputs": {
                "gene_resolution": file_record(OUT_RESOLUTION),
                "open_targets_gene_associations": file_record(OUT_ASSOCIATIONS),
                "open_targets_datasources": file_record(OUT_DATASOURCES),
                "candidate_gene_context": file_record(OUT_CONTEXT),
                "summary": file_record(OUT_SUMMARY),
            },
            "interpretation_guardrail": (
                "An Open Targets gene-disease association supports COPD relevance of a "
                "mapped gene. It neither validates the candidate variant nor establishes "
                "that the candidate regulatory element controls that gene."
            ),
        },
    )
    LOG.write_text(
        f"result_id={RESULT_ID}\n"
        f"candidate_input={candidate_input}\n"
        f"candidate_count={len(candidates)}\n"
        f"unique_gene_tokens={len(unique_tokens)}\n"
        f"resolved_gene_tokens={resolved_count}\n"
        f"open_targets_associated_genes={associated_genes}\n"
        f"candidates_with_OpenTargets_gene_context={len(candidates_with_ot)}\n"
    )
    print(
        f"{RESULT_ID}: resolved {resolved_count:,}/{len(unique_tokens):,} gene tokens; "
        f"{associated_genes:,}/{len(association_rows):,} resolved genes had exact COPD "
        "Open Targets associations."
    )


if __name__ == "__main__":
    main()
