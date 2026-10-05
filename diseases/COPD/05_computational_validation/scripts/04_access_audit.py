#!/usr/bin/env python3
"""Create an explicit public/controlled-resource and biobank access audit."""

from __future__ import annotations

from pathlib import Path

from common import (
    LOGS,
    RESULTS,
    ROOT,
    ensure_directories,
    file_record,
    write_manifest,
    write_tsv,
)


RESULT_ID = "COPD-S5-R004"
INVENTORY = ROOT / "data/validation_resources.tsv"
OUT_RESOURCES = RESULTS / f"{RESULT_ID}_resource_access_audit.tsv"
OUT_BIOBANKS = RESULTS / f"{RESULT_ID}_biobank_access_audit.tsv"
OUT_SUMMARY = RESULTS / f"{RESULT_ID}_access_summary.tsv"
OUT_MANIFEST = RESULTS / f"{RESULT_ID}_manifest.json"
LOG = LOGS / "04_access_audit.log"

RESOURCE_FIELDS = [
    "resource",
    "category",
    "COPD_relevance",
    "evidence_unit",
    "local_path",
    "local_status",
    "access_class",
    "section5_use_status",
    "reason_or_limitation",
]

BIOBANK_FIELDS = [
    "biobank",
    "population_scope",
    "local_COPD_data",
    "summary_level_access",
    "participant_level_access",
    "section5_use_status",
    "required_next_step",
    "interpretation",
]


def local_status(relative: str) -> str:
    if not relative or relative == "-":
        return "not_local"
    path = ROOT / "data" / relative
    if not path.exists():
        return "declared_path_missing"
    if path.is_file():
        return f"local_file:{path.stat().st_size}_bytes"
    files = sum(1 for item in path.rglob("*") if item.is_file())
    return f"local_directory:{files}_files"


def main() -> None:
    ensure_directories()
    resources = [
        {
            "resource": "GTEx v10",
            "category": "molecular QTL",
            "COPD_relevance": "high; Lung tissue",
            "evidence_unit": "exact GRCh38 variant-gene significant cis-eQTL pair",
            "local_path": "gtex/GTEx_Analysis_v10_eQTL.tar",
            "access_class": "public bulk summary data",
            "section5_use_status": "analyzed in COPD-S5-R001",
            "reason_or_limitation": "bulk non-diseased lung; not COPD-specific or cell-type-specific",
        },
        {
            "resource": "QTLbase2",
            "category": "multi-omic QTL",
            "COPD_relevance": "conditional on exact candidate rsID and tissue/QTL type",
            "evidence_unit": "public exact-rsID endpoint association",
            "local_path": "-",
            "access_class": "public web endpoint; no local bulk mirror",
            "section5_use_status": "not used as primary evidence",
            "reason_or_limitation": "query only exact candidate rsIDs when endpoint stability and returned-count audit can be verified; linked focal-tag rsIDs are not candidate IDs",
        },
        {
            "resource": "xQTL Atlas (NIAGADS)",
            "category": "multi-omic QTL",
            "COPD_relevance": "potentially relevant after lung track selection",
            "evidence_unit": "track-specific association",
            "local_path": "xqtl_atlas",
            "access_class": "public metadata; bulk associations via NIAGADS DSS",
            "section5_use_status": "metadata audited; no disease-specific bulk query",
            "reason_or_limitation": "association files are not locally mirrored and no lung-specific track was prespecified",
        },
        {
            "resource": "xQTL Serve",
            "category": "molecular QTL",
            "COPD_relevance": "superseded",
            "evidence_unit": "legacy portal query",
            "local_path": "-",
            "access_class": "legacy public portal unavailable",
            "section5_use_status": "not used",
            "reason_or_limitation": "functionally superseded by xQTL Atlas in the project inventory",
        },
        {
            "resource": "ProteomeVariation",
            "category": "protein QTL",
            "COPD_relevance": "potential candidate-rsID pQTL evidence",
            "evidence_unit": "portal significant pQTL record",
            "local_path": "-",
            "access_class": "public portal; no local bulk mirror",
            "section5_use_status": "not queried",
            "reason_or_limitation": "GTEx Lung exact cis-eQTL was prioritized; portal results require a separate exact-rsID count audit",
        },
        {
            "resource": "MetaBrain",
            "category": "molecular QTL",
            "COPD_relevance": "low tissue relevance; brain",
            "evidence_unit": "GRCh38 brain cis-eQTL association",
            "local_path": "metabrain/2021-07-23-release",
            "access_class": "local research use with redistribution restriction",
            "section5_use_status": "not used for COPD primary validation",
            "reason_or_limitation": "brain tissue does not match the lung-centered COPD model; reuse terms prohibit redistribution",
        },
        {
            "resource": "SingleBrain v2",
            "category": "molecular QTL",
            "COPD_relevance": "low tissue relevance; brain",
            "evidence_unit": "top brain association",
            "local_path": "singlebrain/top_associations",
            "access_class": "public compact subset; full statistics not local",
            "section5_use_status": "not used",
            "reason_or_limitation": "brain-specific and only top-association subset is locally available",
        },
        {
            "resource": "PsychENCODE",
            "category": "molecular QTL/regulatory",
            "COPD_relevance": "low tissue relevance; brain",
            "evidence_unit": "brain molecular association",
            "local_path": "-",
            "access_class": "public processed products; primary data controlled",
            "section5_use_status": "not used",
            "reason_or_limitation": "tissue mismatch and primary data require controlled access",
        },
        {
            "resource": "MPRAbase v4.9.3",
            "category": "MPRA",
            "COPD_relevance": "cross-context experimental regulatory assay evidence",
            "evidence_unit": "assayed library element containing candidate base or exact candidate rsID metadata",
            "local_path": "mprabase/mprabase_v4_9.3.db",
            "access_class": "public curated SQLite release",
            "section5_use_status": "analyzed in COPD-S5-R002",
            "reason_or_limitation": "most records are hg19 and assay scores are heterogeneous; coordinate containment is not automatically allele-specific",
        },
        {
            "resource": "BrainTF",
            "category": "TF ChIP-seq",
            "COPD_relevance": "low tissue relevance; brain",
            "evidence_unit": "brain-region and cell-type TF peak overlap",
            "local_path": "braintf",
            "access_class": "public author-provided consolidated peak collections",
            "section5_use_status": "not used",
            "reason_or_limitation": "brain-specific TF binding does not match the lung-centered COPD models",
        },
        {
            "resource": "Open Targets Platform",
            "category": "gene-disease catalog",
            "COPD_relevance": "exact MONDO_0005002 gene-disease context",
            "evidence_unit": "direct target-disease association and datasource aggregation",
            "local_path": "open_targets",
            "access_class": "public bulk parquet mirror",
            "section5_use_status": "analyzed in COPD-S5-R003",
            "reason_or_limitation": "gene-level context does not validate a candidate variant or enhancer-gene link",
        },
        {
            "resource": "HGNC",
            "category": "gene nomenclature",
            "COPD_relevance": "resolves candidate-linked gene identifiers",
            "evidence_unit": "approved symbol, alias, and stable gene identifiers",
            "local_path": "hgnc/hgnc_complete_set.txt",
            "access_class": "public bulk table",
            "section5_use_status": "analyzed in COPD-S5-R003",
            "reason_or_limitation": "nomenclature evidence only; no disease association implied",
        },
        {
            "resource": "OMIM",
            "category": "gene-disease catalog",
            "COPD_relevance": "potential Mendelian gene context",
            "evidence_unit": "licensed gene-phenotype entry",
            "local_path": "-",
            "access_class": "license or API key required",
            "section5_use_status": "not used",
            "reason_or_limitation": "licensed data are absent from the shared store",
        },
        {
            "resource": "DisGeNET",
            "category": "gene-disease catalog",
            "COPD_relevance": "potential aggregated gene-disease context",
            "evidence_unit": "registered/licensed gene-disease association",
            "local_path": "-",
            "access_class": "registration or license required",
            "section5_use_status": "not used",
            "reason_or_limitation": "release-specific licensed data are absent; Open Targets is locally reproducible",
        },
        {
            "resource": "Neale Lab UK Biobank LDSC",
            "category": "population genetics summary",
            "COPD_relevance": "trait-level heritability context only",
            "evidence_unit": "UK Biobank trait SNP-heritability topline",
            "local_path": "neale_ukbb_ldsc/ukb31063_h2_topline.02Oct2019.tsv.gz",
            "access_class": "public summary table",
            "section5_use_status": "not used for exact candidate validation",
            "reason_or_limitation": "not a candidate-variant association resource",
        },
        {
            "resource": "HuGE Calculator",
            "category": "variant-disease portal",
            "COPD_relevance": "scope mismatch",
            "evidence_unit": "gene-by-metabolic-phenotype query",
            "local_path": "-",
            "access_class": "public portal",
            "section5_use_status": "not used",
            "reason_or_limitation": "current project inventory characterizes the tool as metabolic-disease focused rather than an exact COPD variant resource",
        },
        {
            "resource": "MAtCH",
            "category": "twin-study catalog",
            "COPD_relevance": "heritability context",
            "evidence_unit": "study-level twin record",
            "local_path": "-",
            "access_class": "public portal",
            "section5_use_status": "not used for candidate validation",
            "reason_or_limitation": "not variant-level evidence",
        },
    ]
    for row in resources:
        row["local_status"] = local_status(row["local_path"])

    biobanks = [
        {
            "biobank": "All of Us",
            "population_scope": "United States, diverse ancestry",
            "local_COPD_data": "none",
            "summary_level_access": "no COPD candidate summary set mirrored in this project",
            "participant_level_access": "controlled research environment/authorization",
            "section5_use_status": "not analyzed",
            "required_next_step": "approved workspace plus prespecified COPD phenotype and candidate analysis",
            "interpretation": "absence here is an access/data-availability gap, not negative replication",
        },
        {
            "biobank": "Million Veteran Program",
            "population_scope": "United States veterans, multi-ancestry",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate summary set mirrored",
            "participant_level_access": "controlled collaborative/application access",
            "section5_use_status": "not analyzed",
            "required_next_step": "authorized collaboration and harmonized COPD phenotype",
            "interpretation": "absence here is not negative replication",
        },
        {
            "biobank": "UK Biobank",
            "population_scope": "United Kingdom population cohort",
            "local_COPD_data": "public Neale LDSC topline only; no exact candidate COPD statistics",
            "summary_level_access": "some public derived results; analysis-specific summaries vary",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not used for exact candidate replication",
            "required_next_step": "approved project and prespecified spirometry/diagnosis phenotype",
            "interpretation": "trait-level public files cannot substitute for candidate association testing",
        },
        {
            "biobank": "Biobank Japan",
            "population_scope": "Japanese clinical cohort",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate COPD release mirrored",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "locate an appropriate COPD summary release or obtain approved access",
            "interpretation": "important ancestry replication gap",
        },
        {
            "biobank": "Estonian Biobank",
            "population_scope": "Estonian population cohort",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate COPD release mirrored",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "approved access and phenotype harmonization",
            "interpretation": "absence here is not negative replication",
        },
        {
            "biobank": "FinnGen",
            "population_scope": "Finnish founder population",
            "local_COPD_data": "none",
            "summary_level_access": "public releases may provide phenotype summary statistics; none mirrored for this analysis",
            "participant_level_access": "controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "select release/endpoint and perform exact allele/build-harmonized lookup",
            "interpretation": "founder-population replication remains a planned analysis",
        },
        {
            "biobank": "Qatar Biobank",
            "population_scope": "Qatari and regional population cohort",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate COPD release mirrored",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "approved access and phenotype harmonization",
            "interpretation": "important regional ancestry gap",
        },
        {
            "biobank": "Genes & Health",
            "population_scope": "British Bangladeshi and Pakistani communities",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate COPD release mirrored",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "approved collaboration and candidate/phenotype plan",
            "interpretation": "important ancestry and rare-variant follow-up gap",
        },
        {
            "biobank": "China Kadoorie Biobank",
            "population_scope": "Chinese population cohort",
            "local_COPD_data": "none",
            "summary_level_access": "no exact candidate COPD release mirrored",
            "participant_level_access": "application-controlled",
            "section5_use_status": "not analyzed",
            "required_next_step": "approved access or appropriate public COPD summary release",
            "interpretation": "important East Asian replication gap",
        },
    ]

    summary_rows = [
        {
            "metric": "framework_named_or_replacement_resources_audited",
            "n": len(resources),
            "denominator": len(resources),
            "fraction": 1.0,
            "note": "includes local replacements documented in validation_resources.tsv",
        },
        {
            "metric": "resources_used_in_primary_section5_analysis",
            "n": sum("analyzed in" in row["section5_use_status"] for row in resources),
            "denominator": len(resources),
            "fraction": sum("analyzed in" in row["section5_use_status"] for row in resources)
            / len(resources),
            "note": "GTEx, MPRAbase, Open Targets, and HGNC",
        },
        {
            "metric": "biobanks_audited",
            "n": len(biobanks),
            "denominator": len(biobanks),
            "fraction": 1.0,
            "note": "all biobanks named by Framework Section 5",
        },
        {
            "metric": "biobanks_with_local_exact_candidate_COPD_statistics",
            "n": 0,
            "denominator": len(biobanks),
            "fraction": 0.0,
            "note": "do not interpret unavailable replication as a null result",
        },
    ]

    write_tsv(OUT_RESOURCES, resources, RESOURCE_FIELDS)
    write_tsv(OUT_BIOBANKS, biobanks, BIOBANK_FIELDS)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["metric", "n", "denominator", "fraction", "note"],
    )
    write_manifest(
        OUT_MANIFEST,
        {
            "result_id": RESULT_ID,
            "analysis": "Section 5 public, controlled-access, and biobank resource audit",
            "scope": (
                "local availability and project access state; not a legal determination "
                "or a claim that external access policies will remain unchanged"
            ),
            "inventory_date": "2026-09-25",
            "audit_date": "2026-10-01",
            "resources_audited": len(resources),
            "biobanks_audited": len(biobanks),
            "inputs": {"shared_validation_inventory": file_record(INVENTORY)},
            "outputs": {
                "resource_access_audit": file_record(OUT_RESOURCES),
                "biobank_access_audit": file_record(OUT_BIOBANKS),
                "summary": file_record(OUT_SUMMARY),
            },
        },
    )
    LOG.write_text(
        f"result_id={RESULT_ID}\n"
        f"resources_audited={len(resources)}\n"
        f"biobanks_audited={len(biobanks)}\n"
        "biobanks_with_local_exact_candidate_COPD_statistics=0\n"
    )
    print(
        f"{RESULT_ID}: audited {len(resources)} resources and {len(biobanks)} biobanks."
    )


if __name__ == "__main__":
    main()
