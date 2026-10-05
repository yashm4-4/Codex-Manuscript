#!/usr/bin/env python3
"""Map predicted causal COPD regulatory variants to provisional target genes.

This produces auditable candidate-to-gene hypotheses from GENCODE v50 TSS
proximity, GWAS Catalog mapped-gene links, selected COPD locus membership, and
direct coding overlap.  These are not treated as enhancer-to-gene proof.  Lung
eQTL evidence is added independently in Section 5.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import json
import platform
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases/COPD/04_modeling"
RESULTS = SECTION / "results"
INPUT = RESULTS / "COPD-S4-R005_predicted_causal_population_genetics.tsv"
GTF = ROOT / "data/gencode/gencode.v50.annotation.gtf.gz"

OUT_EVIDENCE = RESULTS / "COPD-S4-R006_candidate_target_evidence.tsv"
OUT_CANDIDATES = RESULTS / "COPD-S4-R006_candidate_target_summary.tsv"
OUT_GENES = RESULTS / "COPD-S4-R006_target_gene_ranking.tsv"
OUT_MANIFEST = RESULTS / "COPD-S4-R006_analysis_manifest.json"

PROXIMAL_WINDOW_BP = 100_000
CANONICAL = {f"chr{value}" for value in range(1, 23)} | {"chrX", "chrY"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).lower() in {"true", "1", "yes"}


def tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return sorted({part.strip() for part in str(value).split(";") if part.strip()})


def parse_attributes(attributes: str) -> dict[str, str]:
    return {
        key: value
        for key, value in re.findall(r'(\S+) "([^"]+)";', attributes)
    }


def load_genes() -> dict[str, list[dict[str, object]]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in CANONICAL:
                continue
            attrs = parse_attributes(fields[8])
            start, end = int(fields[3]), int(fields[4])
            tss = start if fields[6] == "+" else end
            genes[fields[0]].append(
                {
                    "gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
                    "gene_type": attrs.get("gene_type", ""),
                    "tss_1based": tss,
                    "strand": fields[6],
                }
            )
    for chrom in genes:
        genes[chrom].sort(
            key=lambda row: (int(row["tss_1based"]), str(row["gene_id"]))
        )
    return dict(genes)


def chrom_name(value: object) -> str:
    result = str(value).removeprefix("chr")
    if result.endswith(".0") and result[:-2].isdigit():
        result = result[:-2]
    return f"chr{result}"


def nearest(entries: list[dict[str, object]], position: int, protein_only: bool = False):
    selected = (
        [entry for entry in entries if entry["gene_type"] == "protein_coding"]
        if protein_only
        else entries
    )
    if not selected:
        return None
    positions = [int(entry["tss_1based"]) for entry in selected]
    index = bisect.bisect_left(positions, position)
    neighbors = selected[max(0, index - 1) : min(len(selected), index + 1)]
    return min(
        neighbors,
        key=lambda entry: (
            abs(int(entry["tss_1based"]) - position), str(entry["gene_id"])
        ),
    )


def add_evidence(
    rows: list[dict[str, object]],
    seen: set[tuple[str, str, str]],
    *,
    candidate: str,
    gene_name: str,
    gene_id: str = "",
    gene_type: str = "",
    method: str,
    distance_bp: float | int | None = None,
    interpretation: str,
    model_context: str,
) -> None:
    if not gene_name:
        return
    key = (candidate, gene_name, method)
    if key in seen:
        return
    seen.add(key)
    rows.append(
        {
            "candidate_record_id": candidate,
            "gene_name": gene_name,
            "gencode_gene_id": gene_id,
            "gene_type": gene_type,
            "mapping_method": method,
            "absolute_tss_distance_bp": distance_bp,
            "model_context": model_context,
            "evidence_interpretation": interpretation,
        }
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(INPUT, sep="\t", low_memory=False)
    genes = load_genes()
    by_name: dict[str, list[dict[str, object]]] = defaultdict(list)
    for entries in genes.values():
        for entry in entries:
            by_name[str(entry["gene_name"])].append(entry)

    evidence_rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    candidate_rows = []
    for variant in variants.itertuples(index=False):
        candidate = str(variant.candidate_record_id)
        chrom = chrom_name(variant.chromosome_grch38)
        position = int(variant.position_grch38)
        entries = genes.get(chrom, [])
        enhancer = as_bool(variant.predicted_causal_enhancer)
        silencer = as_bool(variant.predicted_causal_silencer)
        model_context = "both" if enhancer and silencer else "enhancer" if enhancer else "silencer"

        nearest_any = nearest(entries, position)
        nearest_protein = nearest(entries, position, protein_only=True)
        for method, entry in (
            ("nearest_TSS_any_gene", nearest_any),
            ("nearest_TSS_protein_coding", nearest_protein),
        ):
            if entry is not None:
                add_evidence(
                    evidence_rows,
                    seen,
                    candidate=candidate,
                    gene_name=str(entry["gene_name"]),
                    gene_id=str(entry["gene_id"]),
                    gene_type=str(entry["gene_type"]),
                    method=method,
                    distance_bp=abs(int(entry["tss_1based"]) - position),
                    interpretation="proximity hypothesis only; distance does not establish regulatory targeting",
                    model_context=model_context,
                )

        positions = [int(entry["tss_1based"]) for entry in entries]
        lower = bisect.bisect_left(positions, position - PROXIMAL_WINDOW_BP)
        upper = bisect.bisect_right(positions, position + PROXIMAL_WINDOW_BP)
        within_100kb = []
        for entry in entries[lower:upper]:
            distance = abs(int(entry["tss_1based"]) - position)
            within_100kb.append(str(entry["gene_name"]))
            add_evidence(
                evidence_rows,
                seen,
                candidate=candidate,
                gene_name=str(entry["gene_name"]),
                gene_id=str(entry["gene_id"]),
                gene_type=str(entry["gene_type"]),
                method="TSS_within_100kb",
                distance_bp=distance,
                interpretation="proximal candidate target; requires molecular linkage validation",
                model_context=model_context,
            )

        for gene_name in tokens(variant.linked_gwas_genes):
            matches = by_name.get(gene_name, [])
            match = matches[0] if len(matches) == 1 else {}
            add_evidence(
                evidence_rows,
                seen,
                candidate=candidate,
                gene_name=gene_name,
                gene_id=str(match.get("gene_id", "")),
                gene_type=str(match.get("gene_type", "")),
                method="GWAS_Catalog_mapped_gene_via_linked_tag",
                interpretation="association-based gene link propagated through r2>=0.8 LD; not enhancer-to-gene proof",
                model_context=model_context,
            )
        for gene_name in tokens(variant.assigned_selected_loci):
            matches = by_name.get(gene_name, [])
            match = matches[0] if len(matches) == 1 else {}
            add_evidence(
                evidence_rows,
                seen,
                candidate=candidate,
                gene_name=gene_name,
                gene_id=str(match.get("gene_id", "")),
                gene_type=str(match.get("gene_type", "")),
                method="selected_COPD_GWAS_locus_membership",
                interpretation="tag-mapped or physical gene-plus-100-kb locus assignment; not regulatory target proof",
                model_context=model_context,
            )
        for gene_name in tokens(variant.coding_CDS_gene_names):
            matches = by_name.get(gene_name, [])
            match = matches[0] if len(matches) == 1 else {}
            add_evidence(
                evidence_rows,
                seen,
                candidate=candidate,
                gene_name=gene_name,
                gene_id=str(match.get("gene_id", "")),
                gene_type=str(match.get("gene_type", "")),
                method="coding_CDS_overlap",
                distance_bp=0,
                interpretation="direct CDS overlap; this is a coding-gene assignment, not a distal regulatory target",
                model_context=model_context,
            )

        candidate_rows.append(
            {
                # R004 assigns the causal candidates ranks 1..337 in the full
                # candidate table; R005 preserves that source column name.
                "predicted_causal_priority_rank": variant.priority_rank_all_candidates,
                "candidate_record_id": candidate,
                "chromosome_grch38": chrom.removeprefix("chr"),
                "position_grch38": position,
                "ref": variant.ref,
                "alt": variant.alt,
                "model_context": model_context,
                "nearest_gene": "" if nearest_any is None else nearest_any["gene_name"],
                "nearest_gene_tss_distance_bp": np.nan
                if nearest_any is None
                else abs(int(nearest_any["tss_1based"]) - position),
                "nearest_protein_coding_gene": ""
                if nearest_protein is None
                else nearest_protein["gene_name"],
                "nearest_protein_coding_tss_distance_bp": np.nan
                if nearest_protein is None
                else abs(int(nearest_protein["tss_1based"]) - position),
                "genes_with_tss_within_100kb": ";".join(sorted(set(within_100kb))),
                "linked_gwas_genes": variant.linked_gwas_genes,
                "assigned_selected_loci": variant.assigned_selected_loci,
                "coding_CDS_gene_names": variant.coding_CDS_gene_names,
                "target_mapping_limit": "no matched severe-emphysema donor re2g or HiChIP map was available; Section 5 adds GTEx Lung eQTL evidence",
            }
        )

    evidence = pd.DataFrame(evidence_rows)
    candidate_summary = pd.DataFrame(candidate_rows).sort_values(
        "predicted_causal_priority_rank"
    )
    if evidence.empty:
        gene_ranking = pd.DataFrame(
            columns=[
                "target_gene_rank",
                "gene_name",
                "n_predicted_causal_variants",
                "n_enhancer_context",
                "n_silencer_context",
                "mapping_methods",
            ]
        )
    else:
        grouped_rows = []
        for gene_name, frame in evidence.groupby("gene_name", sort=False):
            contexts = frame.drop_duplicates("candidate_record_id").set_index(
                "candidate_record_id"
            )["model_context"]
            grouped_rows.append(
                {
                    "gene_name": gene_name,
                    "n_predicted_causal_variants": frame[
                        "candidate_record_id"
                    ].nunique(),
                    "n_enhancer_context": int(contexts.isin(["enhancer", "both"]).sum()),
                    "n_silencer_context": int(contexts.isin(["silencer", "both"]).sum()),
                    "mapping_methods": ";".join(sorted(set(frame["mapping_method"]))),
                    "minimum_TSS_distance_bp": frame[
                        "absolute_tss_distance_bp"
                    ].min(),
                }
            )
        gene_ranking = pd.DataFrame(grouped_rows).sort_values(
            ["n_predicted_causal_variants", "minimum_TSS_distance_bp", "gene_name"],
            ascending=[False, True, True],
            na_position="last",
        )
        gene_ranking.insert(0, "target_gene_rank", np.arange(1, len(gene_ranking) + 1))

    evidence.to_csv(OUT_EVIDENCE, sep="\t", index=False)
    candidate_summary.to_csv(OUT_CANDIDATES, sep="\t", index=False)
    gene_ranking.to_csv(OUT_GENES, sep="\t", index=False)

    outputs = {
        "candidate_target_evidence": OUT_EVIDENCE,
        "candidate_target_summary": OUT_CANDIDATES,
        "target_gene_ranking": OUT_GENES,
    }
    manifest = {
        "result_id": "COPD-S4-R006",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "gencode_release": "GENCODE v50 (Ensembl 116)",
        "proximity_window_bp": PROXIMAL_WINDOW_BP,
        "candidate_count": len(variants),
        "mapping_methods": [
            "nearest TSS for any GENCODE gene",
            "nearest protein-coding TSS",
            "TSS within 100 kb",
            "GWAS Catalog mapped gene propagated through linked tag",
            "selected COPD GWAS locus membership",
            "coding CDS overlap",
        ],
        "limitation": "all target assignments are provisional; matched-donor chromatin-contact or rE2G maps were unavailable and Lung eQTL evidence is evaluated separately",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {
            "predicted_causal_population_table": {
                "path": str(INPUT),
                "bytes": INPUT.stat().st_size,
                "sha256": sha256(INPUT),
            },
            "gencode_gtf": {
                "path": str(GTF),
                "bytes": GTF.stat().st_size,
                "sha256": sha256(GTF),
            },
        },
        "outputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in outputs.items()
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"Mapped {len(variants):,} predicted causal variants to "
        f"{gene_ranking.shape[0]:,} provisional genes through {len(evidence):,} evidence rows."
    )


if __name__ == "__main__":
    main()
