#!/usr/bin/env python3
"""Match predicted-causal COPD variants to GTEx v10 Lung significant cis-eQTLs.

Matches are exact on GRCh38 chromosome, 1-based position, REF, and ALT. A
position-only or allele-swapped row is deliberately not accepted.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path

import pyarrow
import pyarrow.parquet as pq

from common import (
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
    write_manifest,
    write_tsv,
)


RESULT_ID = "COPD-S5-R001"
GTEX_TAR = ROOT / "data/gtex/GTEx_Analysis_v10_eQTL.tar"
PAIR_MEMBER = (
    "GTEx_Analysis_v10_eQTL_updated/Lung.v10.eQTLs.signif_pairs.parquet"
)
EGENE_MEMBER = "GTEx_Analysis_v10_eQTL_updated/Lung.v10.eGenes.txt.gz"

OUT_PAIRS = RESULTS / f"{RESULT_ID}_GTEx_v10_Lung_exact_significant_pairs.tsv.gz"
OUT_AUDIT = RESULTS / f"{RESULT_ID}_GTEx_v10_Lung_candidate_audit.tsv"
OUT_SUMMARY = RESULTS / f"{RESULT_ID}_GTEx_v10_Lung_summary.tsv"
OUT_MANIFEST = RESULTS / f"{RESULT_ID}_manifest.json"
LOG = LOGS / "01_gtex_lung_eqtl.log"

PAIR_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "variant_class",
    "candidate_rsids",
    "gtex_variant_id",
    "gtex_tissue",
    "gene_id",
    "gene_id_versionless",
    "gene_name",
    "gene_biotype",
    "tss_distance",
    "alternate_allele_frequency",
    "minor_allele_samples",
    "minor_allele_count",
    "pval_nominal",
    "slope_per_gtex_alt_allele",
    "slope_se",
    "pval_nominal_threshold",
    "gene_min_pval_nominal",
    "gene_pval_beta",
    "match_definition",
]

AUDIT_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "candidate_rsids",
    "expected_gtex_variant_id",
    "query_status",
    "exact_significant_pair_count",
    "exact_significant_gene_count",
    "significant_genes",
    "minimum_nominal_p",
]


def load_egene_metadata(archive: tarfile.TarFile) -> dict[str, dict[str, str]]:
    member = archive.extractfile(EGENE_MEMBER)
    if member is None:
        raise FileNotFoundError(f"missing tar member: {EGENE_MEMBER}")
    result: dict[str, dict[str, str]] = {}
    with gzip.open(member, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["gene_id"]] = {
                "gene_name": row.get("gene_name", ""),
                "gene_biotype": row.get("biotype", ""),
            }
    return result


def expected_variant_id(row: dict[str, str]) -> str:
    chrom = candidate_core(row)["chromosome_grch38"]
    pos = integer(row.get("position_grch38", ""))
    ref = str(row.get("ref", "")).upper()
    alt = str(row.get("alt", "")).upper()
    if not chrom or pos is None or not ref or not alt:
        return ""
    if any(base not in "ACGT" for base in ref + alt):
        return ""
    return f"chr{chrom}_{pos}_{ref}_{alt}_b38"


def main() -> None:
    ensure_directories()
    candidate_input, candidates = load_candidates()
    candidates.sort(key=candidate_sort_key)
    wanted: dict[str, list[dict[str, str]]] = defaultdict(list)
    for candidate in candidates:
        variant_id = expected_variant_id(candidate)
        if variant_id:
            wanted[variant_id].append(candidate)

    pair_rows: list[dict[str, object]] = []
    scanned_pairs = 0
    with tarfile.open(GTEX_TAR, "r") as archive:
        gene_metadata = load_egene_metadata(archive)
        member = archive.extractfile(PAIR_MEMBER)
        if member is None:
            raise FileNotFoundError(f"missing tar member: {PAIR_MEMBER}")
        with tempfile.NamedTemporaryFile(suffix=".parquet") as temporary:
            shutil.copyfileobj(member, temporary, length=8 * 1024 * 1024)
            temporary.flush()
            parquet = pq.ParquetFile(temporary.name)
            columns = [
                "gene_id",
                "variant_id",
                "tss_distance",
                "af",
                "ma_samples",
                "ma_count",
                "pval_nominal",
                "slope",
                "slope_se",
                "pval_nominal_threshold",
                "min_pval_nominal",
                "pval_beta",
            ]
            for batch in parquet.iter_batches(batch_size=131_072, columns=columns):
                data = batch.to_pydict()
                scanned_pairs += batch.num_rows
                for index, variant_id in enumerate(data["variant_id"]):
                    candidate_matches = wanted.get(variant_id)
                    if not candidate_matches:
                        continue
                    gene_id = str(data["gene_id"][index])
                    gene = gene_metadata.get(gene_id, {})
                    for candidate in candidate_matches:
                        core = candidate_core(candidate)
                        pair_rows.append(
                            {
                                **core,
                                "gtex_variant_id": variant_id,
                                "gtex_tissue": "Lung",
                                "gene_id": gene_id,
                                "gene_id_versionless": gene_id.split(".")[0],
                                "gene_name": gene.get("gene_name", ""),
                                "gene_biotype": gene.get("gene_biotype", ""),
                                "tss_distance": data["tss_distance"][index],
                                "alternate_allele_frequency": data["af"][index],
                                "minor_allele_samples": data["ma_samples"][index],
                                "minor_allele_count": data["ma_count"][index],
                                "pval_nominal": data["pval_nominal"][index],
                                "slope_per_gtex_alt_allele": data["slope"][index],
                                "slope_se": data["slope_se"][index],
                                "pval_nominal_threshold": data[
                                    "pval_nominal_threshold"
                                ][index],
                                "gene_min_pval_nominal": data[
                                    "min_pval_nominal"
                                ][index],
                                "gene_pval_beta": data["pval_beta"][index],
                                "match_definition": (
                                    "exact GRCh38 chromosome, 1-based position, REF, and ALT; "
                                    "row is present in the GTEx v10 Lung significant-pairs file"
                                ),
                            }
                        )

    pair_rows.sort(
        key=lambda row: (
            integer(row["predicted_causal_priority_rank"]) or 10**12,
            str(row["candidate_record_id"]),
            numeric(row["pval_nominal"]) or 1.0,
            str(row["gene_id"]),
        )
    )
    by_candidate: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        by_candidate[str(row["candidate_record_id"])].append(row)

    audit_rows = []
    for candidate in candidates:
        core = candidate_core(candidate)
        candidate_id = str(core["candidate_record_id"])
        expected = expected_variant_id(candidate)
        matches = by_candidate.get(candidate_id, [])
        genes = sorted(
            {
                str(row["gene_name"] or row["gene_id_versionless"])
                for row in matches
            }
        )
        if not expected:
            status = "not_queryable_missing_or_noncanonical_alleles_or_coordinate"
        elif matches:
            status = "exact_significant_lung_cis_eqtl"
        else:
            status = "no_exact_significant_lung_cis_eqtl"
        audit_rows.append(
            {
                **core,
                "expected_gtex_variant_id": expected,
                "query_status": status,
                "exact_significant_pair_count": len(matches),
                "exact_significant_gene_count": len(genes),
                "significant_genes": ";".join(genes),
                "minimum_nominal_p": min(
                    (float(row["pval_nominal"]) for row in matches), default=None
                ),
            }
        )

    matched_candidates = sum(
        row["query_status"] == "exact_significant_lung_cis_eqtl"
        for row in audit_rows
    )
    queryable = sum(bool(row["expected_gtex_variant_id"]) for row in audit_rows)
    unique_genes = {
        str(row["gene_id_versionless"]) for row in pair_rows if row["gene_id_versionless"]
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
            "metric": "queryable_exact_GRCh38_REF_ALT",
            "n": queryable,
            "denominator": len(candidates),
            "fraction": queryable / len(candidates),
            "note": "canonical sequence alleles and coordinate available",
        },
        {
            "metric": "candidates_with_exact_significant_GTEx_v10_Lung_cis_eqtl",
            "n": matched_candidates,
            "denominator": queryable,
            "fraction": matched_candidates / queryable if queryable else None,
            "note": "exact GRCh38 chromosome-position-REF-ALT match",
        },
        {
            "metric": "exact_significant_variant_gene_pairs",
            "n": len(pair_rows),
            "denominator": "",
            "fraction": "",
            "note": "a candidate can be associated with multiple significant eGenes",
        },
        {
            "metric": "unique_significant_eGenes",
            "n": len(unique_genes),
            "denominator": "",
            "fraction": "",
            "note": "versionless Ensembl gene IDs",
        },
    ]

    write_tsv(OUT_PAIRS, pair_rows, PAIR_FIELDS)
    write_tsv(OUT_AUDIT, audit_rows, AUDIT_FIELDS)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["metric", "n", "denominator", "fraction", "note"],
    )
    write_manifest(
        OUT_MANIFEST,
        {
            "result_id": RESULT_ID,
            "analysis": "exact GTEx v10 Lung significant cis-eQTL lookup",
            "genome_build": "GRCh38",
            "tissue": "Lung",
            "match_definition": "exact chromosome, 1-based position, REF, and ALT",
            "candidate_count": len(candidates),
            "queryable_candidate_count": queryable,
            "matched_candidate_count": matched_candidates,
            "exact_pair_count": len(pair_rows),
            "source_pair_rows_scanned": scanned_pairs,
            "software": {
                "python": __import__("platform").python_version(),
                "pyarrow": pyarrow.__version__,
            },
            "inputs": {
                "section4_candidates": file_record(candidate_input),
                "gtex_v10_eqtl_tar": file_record(GTEX_TAR),
                "gtex_pair_member": PAIR_MEMBER,
                "gtex_egene_member": EGENE_MEMBER,
            },
            "outputs": {
                "exact_pairs": file_record(OUT_PAIRS),
                "candidate_audit": file_record(OUT_AUDIT),
                "summary": file_record(OUT_SUMMARY),
            },
            "interpretation_guardrail": (
                "GTEx Lung is bulk non-diseased tissue; an exact significant eQTL is "
                "molecular association evidence, not proof of COPD causality or cell-type specificity."
            ),
        },
    )
    LOG.write_text(
        f"result_id={RESULT_ID}\n"
        f"candidate_input={candidate_input}\n"
        f"candidate_count={len(candidates)}\n"
        f"source_pair_rows_scanned={scanned_pairs}\n"
        f"matched_candidates={matched_candidates}\n"
        f"exact_pairs={len(pair_rows)}\n"
    )
    print(
        f"{RESULT_ID}: {matched_candidates:,}/{queryable:,} queryable candidates "
        f"matched {len(pair_rows):,} exact significant Lung cis-eQTL pairs."
    )


if __name__ == "__main__":
    main()

