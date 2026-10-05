#!/usr/bin/env python3
"""Integrate Section 5 evidence without changing the Section 4 priority order."""

from __future__ import annotations

from collections import defaultdict

from common import (
    LOGS,
    RESULTS,
    SECTION,
    candidate_core,
    candidate_sort_key,
    ensure_directories,
    file_record,
    integer,
    is_true,
    load_candidates,
    numeric,
    read_tsv,
    write_manifest,
    write_tsv,
)


RESULT_ID = "COPD-S5-R005"
GTEX_AUDIT = RESULTS / "COPD-S5-R001_GTEx_v10_Lung_candidate_audit.tsv"
GTEX_PAIRS = RESULTS / "COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz"
MPRA_AUDIT = RESULTS / "COPD-S5-R002_MPRAbase_candidate_audit.tsv"
MPRA_ELEMENTS = RESULTS / "COPD-S5-R002_MPRAbase_v4_9_3_candidate_elements.tsv.gz"
GENE_CONTEXT = RESULTS / "COPD-S5-R003_candidate_gene_context.tsv.gz"
ACCESS_SUMMARY = RESULTS / "COPD-S5-R004_access_summary.tsv"

OUT_MATRIX = RESULTS / f"{RESULT_ID}_integrated_candidate_validation.tsv"
OUT_SUMMARY = RESULTS / f"{RESULT_ID}_integrated_summary.tsv"
OUT_REPORT = RESULTS / "section_5_computational_validation.md"
OUT_MANIFEST = RESULTS / f"{RESULT_ID}_manifest.json"
LOG = LOGS / "05_integrate_validation.log"

MATRIX_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "variant_class",
    "candidate_rsids",
    "priority_tier",
    "predicted_causal_enhancer",
    "predicted_causal_silencer",
    "linked_gwas_tag_ids",
    "linked_gwas_genes",
    "assigned_selected_loci",
    "gtex_lung_exact_significant_eqtl",
    "gtex_lung_significant_pair_count",
    "gtex_lung_significant_gene_count",
    "gtex_lung_significant_genes",
    "gtex_lung_minimum_nominal_p",
    "gtex_lung_gene_matches_candidate_mapping",
    "gtex_lung_matching_candidate_gene_count",
    "gtex_lung_matching_candidate_genes",
    "gtex_lung_gene_matches_provisional_target",
    "gtex_lung_matching_provisional_targets",
    "mprabase_any_element_evidence",
    "mprabase_coordinate_element_evidence",
    "mprabase_exact_rsid_metadata_evidence",
    "mprabase_matched_element_count",
    "mprabase_elements_with_scores",
    "mprabase_matched_dataset_count",
    "mprabase_matched_datasets",
    "candidate_gene_mapping_count",
    "hgnc_resolved_gene_count",
    "open_targets_COPD_associated_gene_count",
    "open_targets_COPD_associated_genes",
    "provisional_section4_target_count",
    "provisional_section4_targets",
    "variant_or_regional_evidence_channels",
    "gene_context_channels",
    "integrated_validation_status",
    "interpretation_scope",
]


def require_inputs() -> None:
    missing = [
        path
        for path in (
            GTEX_AUDIT,
            GTEX_PAIRS,
            MPRA_AUDIT,
            MPRA_ELEMENTS,
            GENE_CONTEXT,
            ACCESS_SUMMARY,
        )
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "run Section 5 component scripts before integration; missing:\n"
            + "\n".join(f"  - {path}" for path in missing)
        )


def main() -> None:
    ensure_directories()
    require_inputs()
    candidate_input, candidates = load_candidates()
    candidates.sort(key=candidate_sort_key)
    candidate_ids = {row["candidate_record_id"] for row in candidates}

    gtex_rows = read_tsv(GTEX_AUDIT)
    mpra_rows = read_tsv(MPRA_AUDIT)
    mpra_element_rows = read_tsv(MPRA_ELEMENTS)
    context_rows = read_tsv(GENE_CONTEXT)
    gtex = {row["candidate_record_id"]: row for row in gtex_rows}
    mpra = {row["candidate_record_id"]: row for row in mpra_rows}
    if set(gtex) != candidate_ids:
        raise ValueError("GTEx candidate audit does not match current Section 4 input")
    if set(mpra) != candidate_ids:
        raise ValueError("MPRAbase candidate audit does not match current Section 4 input")

    context_by_candidate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in context_rows:
        if row["candidate_record_id"] in candidate_ids:
            context_by_candidate[row["candidate_record_id"]].append(row)

    matrix_rows = []
    for candidate in candidates:
        core = candidate_core(candidate)
        candidate_id = str(core["candidate_record_id"])
        gtex_row = gtex[candidate_id]
        mpra_row = mpra[candidate_id]
        context = context_by_candidate.get(candidate_id, [])

        gtex_present = (
            gtex_row["query_status"] == "exact_significant_lung_cis_eqtl"
        )
        mpra_coordinate = integer(
            mpra_row["coordinate_matched_element_count"]
        ) not in {None, 0}
        mpra_id = integer(mpra_row["exact_rsid_matched_element_count"]) not in {
            None,
            0,
        }
        mpra_any = integer(mpra_row["matched_element_count"]) not in {None, 0}
        resolved_genes = {
            row["approved_symbol"]
            for row in context
            if row["resolution_status"] == "resolved" and row["approved_symbol"]
        }
        ot_genes = {
            row["approved_symbol"]
            for row in context
            if is_true(row["open_targets_exact_COPD_association_present"])
            and row["approved_symbol"]
        }
        provisional_targets = {
            row["approved_symbol"] or row["input_gene_token"]
            for row in context
            if "provisional_Section4_target:" in row["mapping_origins"]
        }
        mapped_context_genes = {
            row["approved_symbol"]
            for row in context
            if row["approved_symbol"]
            and any(
                origin in {
                    "GWAS_Catalog_mapped_gene",
                    "selected_GWAS_gene_locus",
                }
                for origin in row["mapping_origins"].split(";")
            )
        }
        gtex_genes = set(
            token
            for token in gtex_row["significant_genes"].split(";")
            if token
        )
        matching_mapped_genes = gtex_genes & mapped_context_genes
        matching_provisional_targets = gtex_genes & provisional_targets

        variant_channels = []
        if gtex_present:
            variant_channels.append("GTEx_v10_Lung_exact_significant_cis_eQTL")
        if mpra_id:
            variant_channels.append("MPRAbase_exact_candidate_rsid_metadata")
        if mpra_coordinate:
            variant_channels.append(
                "MPRAbase_roundtrip_exact_hg19_assayed_element_containment"
            )
        gene_channels = []
        if ot_genes:
            gene_channels.append("OpenTargets_exact_COPD_gene_association")
        if provisional_targets:
            gene_channels.append("Section4_provisional_target_mapping")

        if gtex_present and mpra_any:
            status = "exact_lung_eQTL_plus_MPRAbase_element_evidence"
        elif gtex_present:
            status = "exact_lung_eQTL_evidence"
        elif mpra_id:
            status = "exact_rsid_MPRAbase_metadata_evidence"
        elif mpra_coordinate:
            status = "regional_MPRAbase_element_evidence"
        elif ot_genes:
            status = "gene_context_only"
        else:
            status = "no_support_in_queried_primary_public_resources"

        matrix_rows.append(
            {
                **core,
                "gtex_lung_exact_significant_eqtl": gtex_present,
                "gtex_lung_significant_pair_count": gtex_row[
                    "exact_significant_pair_count"
                ],
                "gtex_lung_significant_gene_count": gtex_row[
                    "exact_significant_gene_count"
                ],
                "gtex_lung_significant_genes": gtex_row["significant_genes"],
                "gtex_lung_minimum_nominal_p": gtex_row["minimum_nominal_p"],
                "gtex_lung_gene_matches_candidate_mapping": bool(
                    matching_mapped_genes
                ),
                "gtex_lung_matching_candidate_gene_count": len(
                    matching_mapped_genes
                ),
                "gtex_lung_matching_candidate_genes": ";".join(
                    sorted(matching_mapped_genes)
                ),
                "gtex_lung_gene_matches_provisional_target": bool(
                    matching_provisional_targets
                ),
                "gtex_lung_matching_provisional_targets": ";".join(
                    sorted(matching_provisional_targets)
                ),
                "mprabase_any_element_evidence": mpra_any,
                "mprabase_coordinate_element_evidence": mpra_coordinate,
                "mprabase_exact_rsid_metadata_evidence": mpra_id,
                "mprabase_matched_element_count": mpra_row[
                    "matched_element_count"
                ],
                "mprabase_elements_with_scores": mpra_row[
                    "elements_with_reported_scores"
                ],
                "mprabase_matched_dataset_count": mpra_row[
                    "matched_dataset_count"
                ],
                "mprabase_matched_datasets": mpra_row["matched_datasets"],
                "candidate_gene_mapping_count": len(context),
                "hgnc_resolved_gene_count": len(resolved_genes),
                "open_targets_COPD_associated_gene_count": len(ot_genes),
                "open_targets_COPD_associated_genes": ";".join(sorted(ot_genes)),
                "provisional_section4_target_count": len(provisional_targets),
                "provisional_section4_targets": ";".join(
                    sorted(provisional_targets)
                ),
                "variant_or_regional_evidence_channels": ";".join(variant_channels),
                "gene_context_channels": ";".join(gene_channels),
                "integrated_validation_status": status,
                "interpretation_scope": (
                    "GTEx is exact variant-gene molecular association; MPRAbase coordinate "
                    "containment is regional experimental evidence unless exact allelic testing "
                    "is independently established; Open Targets is gene-level context"
                ),
            }
        )

    n = len(matrix_rows)
    gtex_count = sum(bool(row["gtex_lung_exact_significant_eqtl"]) for row in matrix_rows)
    mpra_any_count = sum(bool(row["mprabase_any_element_evidence"]) for row in matrix_rows)
    mpra_coordinate_count = sum(
        bool(row["mprabase_coordinate_element_evidence"]) for row in matrix_rows
    )
    mpra_id_count = sum(
        bool(row["mprabase_exact_rsid_metadata_evidence"]) for row in matrix_rows
    )
    both_count = sum(
        bool(row["gtex_lung_exact_significant_eqtl"])
        and bool(row["mprabase_any_element_evidence"])
        for row in matrix_rows
    )
    any_candidate_evidence = sum(
        bool(row["gtex_lung_exact_significant_eqtl"])
        or bool(row["mprabase_any_element_evidence"])
        for row in matrix_rows
    )
    gene_context_count = sum(
        int(row["open_targets_COPD_associated_gene_count"]) > 0
        for row in matrix_rows
    )
    gtex_mapped_gene_count = sum(
        bool(row["gtex_lung_gene_matches_candidate_mapping"])
        for row in matrix_rows
    )
    gtex_provisional_target_count = sum(
        bool(row["gtex_lung_gene_matches_provisional_target"])
        for row in matrix_rows
    )
    gtex_with_provisional_targets = sum(
        bool(row["gtex_lung_exact_significant_eqtl"])
        and int(row["provisional_section4_target_count"]) > 0
        for row in matrix_rows
    )
    none_count = n - any_candidate_evidence
    summary_rows = [
        {
            "metric": "predicted_causal_candidates",
            "n": n,
            "denominator": n,
            "fraction": 1.0,
            "note": "Section 4 predicted-causal union; priority order unchanged",
        },
        {
            "metric": "exact_significant_GTEx_v10_Lung_cis_eqtl",
            "n": gtex_count,
            "denominator": n,
            "fraction": gtex_count / n,
            "note": "exact GRCh38 chromosome-position-REF-ALT",
        },
        {
            "metric": "any_MPRAbase_element_evidence",
            "n": mpra_any_count,
            "denominator": n,
            "fraction": mpra_any_count / n,
            "note": "exact candidate-rsID metadata or round-trip-validated hg19 coordinate containment",
        },
        {
            "metric": "GTEx_eGene_matches_candidate_mapped_gene_or_selected_locus",
            "n": gtex_mapped_gene_count,
            "denominator": gtex_count,
            "fraction": gtex_mapped_gene_count / gtex_count if gtex_count else None,
            "note": "HGNC-approved symbol concordance; gene mapping remains contextual",
        },
        {
            "metric": "GTEx_eGene_matches_provisional_Section4_target",
            "n": gtex_provisional_target_count,
            "denominator": gtex_with_provisional_targets,
            "fraction": (
                gtex_provisional_target_count / gtex_with_provisional_targets
                if gtex_with_provisional_targets
                else None
            ),
            "note": "denominator is exact GTEx-supported candidates with an R006 provisional target set",
        },
        {
            "metric": "MPRAbase_coordinate_containment",
            "n": mpra_coordinate_count,
            "denominator": n,
            "fraction": mpra_coordinate_count / n,
            "note": "regional assay evidence; not assumed allele-specific",
        },
        {
            "metric": "MPRAbase_exact_candidate_rsid_metadata",
            "n": mpra_id_count,
            "denominator": n,
            "fraction": mpra_id_count / n,
            "note": "source focal-tag rsIDs excluded",
        },
        {
            "metric": "GTEx_and_MPRAbase_evidence",
            "n": both_count,
            "denominator": n,
            "fraction": both_count / n,
            "note": "intersection of exact Lung eQTL and any MPRAbase evidence",
        },
        {
            "metric": "any_variant_or_regional_public_evidence",
            "n": any_candidate_evidence,
            "denominator": n,
            "fraction": any_candidate_evidence / n,
            "note": "union of GTEx and MPRAbase; evidence scopes differ",
        },
        {
            "metric": "no_variant_or_regional_evidence_in_primary_queries",
            "n": none_count,
            "denominator": n,
            "fraction": none_count / n,
            "note": "not evidence against causality; resources and contexts are incomplete",
        },
        {
            "metric": "candidate_with_OpenTargets_COPD_gene_context",
            "n": gene_context_count,
            "denominator": n,
            "fraction": gene_context_count / n,
            "note": "gene-level context, excluded from variant-validation union",
        },
    ]

    write_tsv(OUT_MATRIX, matrix_rows, MATRIX_FIELDS)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["metric", "n", "denominator", "fraction", "note"],
    )

    mpra_candidate_word = "candidate" if mpra_any_count == 1 else "candidates"
    mpra_coordinate_word = (
        "candidate" if mpra_coordinate_count == 1 else "candidates"
    )
    if mpra_any_count == 1:
        matched = next(
            row for row in matrix_rows if row["mprabase_any_element_evidence"]
        )
        matched_element_rows = [
            row
            for row in mpra_element_rows
            if row["candidate_record_id"] == matched["candidate_record_id"]
        ]
        contexts = sorted(
            {
                value
                for row in matched_element_rows
                for value in row["cell_lines_or_tissues"].split(";")
                if value
            }
        )
        pmids = sorted(
            {row["dataset_pmid"] for row in matched_element_rows if row["dataset_pmid"]}
        )
        mpra_detail = (
            f" The sole match was {matched['candidate_record_id']} "
            f"({matched['candidate_rsids'] or 'no stable rsID'}), covered by "
            f"{len(matched_element_rows)} elements in "
            f"{matched['mprabase_matched_datasets'] or 'an MPRAbase dataset'}"
            f" ({';'.join(contexts) or 'context not reported'}; PMID "
            f"{';'.join(pmids) or 'not reported'})."
        )
    else:
        mpra_detail = ""
    report = f"""# Section 5: Computational Validation

## Scope and inputs

The analysis evaluated all {n:,} Section 4 predicted-causal regulatory candidate records without changing their pre-existing priority order. Three evidence types were kept separate: exact variant-level GTEx v10 Lung significant cis-eQTLs (COPD-S5-R001), MPRAbase v4.9.3 assayed-element evidence after audited genome-build conversion (COPD-S5-R002), and HGNC-resolved, exact-node Open Targets COPD gene associations (COPD-S5-R003). Resource and biobank availability was audited in COPD-S5-R004.

## Exact GTEx v10 Lung cis-eQTL evidence

Exact matching required the same GRCh38 chromosome, 1-based position, REF, and ALT. {gtex_count:,} of {n:,} predicted-causal candidates ({100 * gtex_count / n:.2f}%) matched at least one row in the GTEx v10 Lung significant-pairs file. For {gtex_mapped_gene_count:,} candidates, at least one significant eGene also matched an HGNC-resolved GWAS mapped gene or selected gene locus. A significant eGene matched at least one broad R006 provisional target for {gtex_provisional_target_count:,} of {gtex_with_provisional_targets:,} GTEx-supported candidates with target hypotheses ({100 * gtex_provisional_target_count / gtex_with_provisional_targets if gtex_with_provisional_targets else 0:.2f}%). R006 includes nearest-TSS and all-TSS-within-100-kb hypotheses, so this concordance is not a chromatin-contact assignment. These are molecular association results from bulk, non-diseased lung and do not by themselves establish COPD causality or a specific causal lung cell type.

## MPRAbase evidence

Candidate sites were represented as one-base GRCh38 BED intervals, lifted to hg19 with the UCSC chain, and required to map uniquely and return to the original GRCh38 base on reverse liftover. {mpra_any_count:,} {mpra_candidate_word} ({100 * mpra_any_count / n:.2f}%) had at least one MPRAbase element match. Coordinate containment was observed for {mpra_coordinate_count:,} {mpra_coordinate_word} and an exact candidate rsID in element metadata was observed for {mpra_id_count:,}.{mpra_detail} Coordinate containment shows that an assayed sequence covered the candidate base. It does not demonstrate that the REF and ALT alleles were both tested, and the observed HepG2 context is not a COPD lung-cell assay. Raw score scales cannot be compared across heterogeneous MPRA studies without study-specific calibration.

## Gene catalog context

At least one HGNC-resolved candidate-linked gene had an exact `MONDO_0005002` Open Targets association for {gene_context_count:,} candidates ({100 * gene_context_count / n:.2f}%). The tested gene universe includes broad R006 proximity hypotheses, including all TSSs within 100 kb, so this high coverage is descriptive context rather than evidence of enrichment. It is not independent evidence for the candidate variant or proof that a candidate element regulates that gene. All Section 4 target assignments remain explicitly labeled provisional.

## Integrated result

{any_candidate_evidence:,} candidates ({100 * any_candidate_evidence / n:.2f}%) had exact GTEx or MPRAbase regional/identifier evidence; {both_count:,} had both evidence types. The remaining {none_count:,} candidates had no support in these two primary public queries. That category is not a negative causal call because available resources do not cover all lung cell types, disease states, perturbations, alleles, or regulatory assay designs.

QTLbase was not used as primary evidence. Its public endpoint is reserved for exact candidate-rsID requests with returned-count auditing; GWAS source focal-tag rsIDs are not substituted for LD-proxy identifiers. No framework-named biobank had local exact-candidate COPD statistics. Controlled-access and missing-summary resources are therefore recorded as validation gaps rather than null replications.

## Reproducible outputs

- `COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz`: exact significant variant-gene pairs.
- `COPD-S5-R002_MPRAbase_v4_9_3_candidate_elements.tsv.gz`: candidate-element matches with assay metadata and evidence scope.
- `COPD-S5-R003_candidate_gene_context.tsv.gz`: HGNC resolution and exact COPD Open Targets context.
- `COPD-S5-R004_resource_access_audit.tsv` and `COPD-S5-R004_biobank_access_audit.tsv`: public and controlled-access audit.
- `COPD-S5-R005_integrated_candidate_validation.tsv`: one row per predicted-causal candidate, in Section 4 priority order.
- `COPD-S5-R005_integrated_summary.tsv`: denominators and evidence-overlap counts used above.
"""
    OUT_REPORT.write_text(report)

    write_manifest(
        OUT_MANIFEST,
        {
            "result_id": RESULT_ID,
            "analysis": "integrated computational validation matrix",
            "priority_policy": "preserve Section 4 predicted-causal priority order; do not rerank",
            "candidate_count": n,
            "exact_lung_eqtl_candidates": gtex_count,
            "exact_lung_eqtl_matching_mapped_gene_candidates": gtex_mapped_gene_count,
            "exact_lung_eqtl_matching_provisional_target_candidates": gtex_provisional_target_count,
            "mprabase_candidates": mpra_any_count,
            "both_primary_evidence_types": both_count,
            "any_variant_or_regional_public_evidence": any_candidate_evidence,
            "inputs": {
                "section4_candidates": file_record(candidate_input),
                "gtex_candidate_audit": file_record(GTEX_AUDIT),
                "gtex_exact_pairs": file_record(GTEX_PAIRS),
                "mprabase_candidate_audit": file_record(MPRA_AUDIT),
                "mprabase_elements": file_record(MPRA_ELEMENTS),
                "gene_context": file_record(GENE_CONTEXT),
                "access_summary": file_record(ACCESS_SUMMARY),
            },
            "outputs": {
                "integrated_candidate_matrix": file_record(OUT_MATRIX),
                "integrated_summary": file_record(OUT_SUMMARY),
                "section_report": file_record(OUT_REPORT),
            },
            "evidence_hierarchy": {
                "variant_level": "exact GTEx v10 Lung significant cis-eQTL; exact candidate rsID in MPRAbase metadata",
                "regional_experimental": "round-trip-validated candidate-base containment in an MPRAbase element",
                "gene_context": "HGNC-resolved exact-node Open Targets COPD association",
            },
        },
    )
    LOG.write_text(
        f"result_id={RESULT_ID}\n"
        f"candidate_input={candidate_input}\n"
        f"candidate_count={n}\n"
        f"exact_lung_eqtl_candidates={gtex_count}\n"
        f"mprabase_candidates={mpra_any_count}\n"
        f"both_primary_evidence_types={both_count}\n"
        f"any_variant_or_regional_public_evidence={any_candidate_evidence}\n"
    )
    print(
        f"{RESULT_ID}: integrated {n:,} candidates; {any_candidate_evidence:,} had "
        "GTEx or MPRAbase evidence."
    )


if __name__ == "__main__":
    main()
