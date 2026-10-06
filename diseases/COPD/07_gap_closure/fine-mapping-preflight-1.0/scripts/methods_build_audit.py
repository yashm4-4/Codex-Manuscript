#!/usr/bin/env python3
"""Assemble documentation/metadata audit; no LD values or scientific fits.

Reads only this subtask's saved public documentation and resource metadata.
Checks byte integrity, summarizes metadata and writes a new audit JSON once.
"""
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

STAGE = Path(__file__).resolve().parents[1]
AUDIT = STAGE / "audits/methods_ld"


def record_file(path):
    data = path.read_bytes()
    return {"path": str(path.relative_to(STAGE)), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def main():
    output = AUDIT / "methods_ld_audit.json"
    if output.exists():
        raise RuntimeError("Refuse to overwrite methods audit")
    records = []
    for path in sorted((AUDIT / "sources").glob("*.json")):
        record = json.loads(path.read_text())
        body = STAGE / record["saved_path"]
        actual = record_file(body)
        assert actual["bytes"] == record["saved_bytes"], path
        assert actual["sha256"] == record["saved_sha256"], path
        record["metadata_file"] = record_file(path)
        record["byte_integrity"] = "PASS"
        record["semantic_status"] = (
            "ACCESS_FAILURE_RETAINED" if record["status"] != "HTTP_SUCCESS" else
            "BOUNDED_PREFIX_NOT_FULL_DOCUMENT" if record.get("truncated") else
            "PUBLIC_DOCUMENT_OR_METADATA_RETRIEVED_NOT_LOCUS_QC"
        )
        if record["source_id"] == "susie_get_cs_reference":
            record["semantic_status"] = "HTML_META_REDIRECT_STUB_NOT_METHOD_DOCUMENT"
        if record["source_id"] == "bbj_ld_github":
            record["semantic_status"] = "SEARCH_ZERO_RESULTS_NOT_PROOF_OF_ABSENCE"
        if record["source_id"] == "gnomad_ld_docs":
            record["semantic_status"] = "GENERAL_RELEASE_DOCUMENT_NOT_SPECIFIC_LD_CONTRACT"
        if record["source_id"] == "finemap_official_http":
            record["transport_limitation"] = "HTTP_ONLY_AFTER_TWO_HTTPS_CERTIFICATE_FAILURES; public documentation, not executable download"
        records.append(record)

    population_ns = {"AFR": 6636, "AMR": 980, "CSA": 8876, "EAS": 2709, "EUR": 420531, "MID": 1599}
    pan_populations = []
    for pop in population_ns:
        text = (AUDIT / f"sources/panukb_{pop}_blockmatrix_metadata.body").read_text()
        header = {key: int(re.search(r'"' + key + r'":(\d+)', text).group(1))
                  for key in ("blockSize", "nRows", "nCols")}
        assert header["nRows"] == header["nCols"]
        population = {
            "ancestry": pop, "documented_overall_cohort_n": population_ns[pop],
            "ld_matrix_global_n": None,
            "ld_matrix_global_n_status": "GLOBAL_FIELD_EXISTS_BUT_VALUE_NOT_DECODED; overall cohort N is not asserted as exact LD or phenotype N",
            "matrix_metadata": header,
            "matrix_url": f"s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.bm",
            "variant_index_url": f"s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.variant.ht",
            "matrix_values_accessed": False,
            "evidence_ids": [f"panukb_{pop}_blockmatrix_metadata", "panukb_overview_docs"],
        }
        if pop in ("EUR", "EAS"):
            md = json.loads(gzip.decompress((AUDIT / f"sources/panukb_{pop}_variant_metadata.body").read_bytes()))
            population["variant_index_schema"] = md["table_type"]
            population["variant_index_writer_version"] = md["hail_version"]
            population["variant_index_row_count"] = sum(md["components"]["partition_counts"]["counts"])
            assert population["variant_index_row_count"] == header["nRows"]
            population["evidence_ids"].append(f"panukb_{pop}_variant_metadata")
        pan_populations.append(population)

    rows = list(csv.DictReader((AUDIT / "sources/onekg_sample_panel.body").read_text().splitlines(), delimiter="\t"))
    one_kg_n = {"total": len(rows), "superpopulation": dict(Counter(row["super_pop"] for row in rows)),
                "population": dict(Counter(row["pop"] for row in rows))}
    assert one_kg_n["total"] == 2504

    resources = [
        {
            "resource_id": "PANUKB_LD_RELEASE", "status": "PUBLIC_METADATA_VERIFIED_LOCUS_QC_NOT_RUN",
            "native_build": "GRCh37", "matrix_type": "signed, covariate-adjusted dosage Gram/correlation; triangular, 10 Mb banded",
            "allele_coding": "Intended ALT dosage: compute_ld_matrix.py -> get_filtered_mt -> hail.import_bgen(dosage); Hail defines expected ALT count; index alleles[0] REF",
            "documentation_filters": "INFO >0.8; MAC >20 (resource documentation, not a new fine-mapping readiness rule)",
            "release_provenance_discrepancy": "Archived current source uses min_mac=19 with helper >= and no explicit INFO filter in that call chain; cannot assume current master reproduces historical release; resolve actual release QC/index before inference",
            "scaling_caveat": "Source normalizes dosage before residualization, then writes Z Z^T/n without explicit post-residualization standardization; future released-diagonal/statistic-scaling validation required",
            "covariates": ["intercept", "age", "sex", "age^2", "age*sex", "age^2*sex", "PC1-PC10"],
            "phenotype_sample_match": "NOT_ESTABLISHED; same-project population LD need not be exact phenotype/subgroup sample LD",
            "population_metadata": pan_populations,
            "documented_all_population_size_tb": 43.3, "documented_eur_size_tb": 14.1,
            "future_access": "Hail-compatible locus-only block extraction plus exact variant index; no full matrix download authorized",
            "evidence_ids": ["panukb_ld_docs", "panukb_hail_docs", "panukb_flat_docs", "panukb_ld_code", "panukb_genotype_resource_code", "hail_bgen_import_contract"],
        },
        {
            "resource_id": "BBJ_MATCHED_DOSAGE_LD", "status": "PUBLIC_MATCHED_SIGNED_MATRIX_NOT_VERIFIED",
            "native_build": "GRCh37/hg19 for audited original GWAS; any future LD release must independently match",
            "ld_sample_n": None, "release_accession": None, "matrix_type": "Original in-sample dosage LD reported in method paper; downloadable signed R/index not established",
            "available_not_equivalent": ["PheWeb/hum0197 GWAS files", "Public fine-mapping output archives", "Controlled-access genotype routes"],
            "proposed_action": "Obtain/verify appropriate BBJ dosage LD and index or a scientifically justified Japanese reference under a new authorized plan; no generic EUR substitution",
            "summary_statistics_access_update": "Alternative-source audit reports original Ishigaki COPD hum0014.v17.COPD.v1.zip signed schema/access verified; this does not resolve LD access",
            "summary_statistics_url": "https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.COPD.v1.zip",
            "source_update_provenance": "Coordination from /root/external_v2_model_contract; no GWAS data read by methods/LD subtask",
            "evidence_ids": ["bbj_pheweb_downloads", "bbj_hum0197_availability", "bbj_finemapping_availability", "bbj_finemapping_methods"],
        },
        {
            "resource_id": "GBMI_DIAGNOSTIC_GNOMAD_LD", "status": "DIAGNOSTIC_PRECEDENT_NOT_MATCHED_FINE_MAPPING_R",
            "native_build": "gnomAD v2 LD GRCh37; published workflow lifted identity to GRCh38 for GBMI",
            "ld_sample_n": None, "matrix_type": "SLALOM uses external ancestry correlations and weighted r2, not a demonstrated joint signed meta-analysis R",
            "ancestries_in_published_diagnostic": ["African", "Admixed American", "East Asian", "Finnish", "non-Finnish European"],
            "missing_ancestries_in_published_diagnostic": ["Central/South Asian", "Middle Eastern"],
            "blocker": "Per-variant cohort/sample heterogeneity, overlapping cohorts and mixed ancestries require an explicit compatible covariance/likelihood; ancestry-specific meta-analysis remains multi-cohort",
            "prohibited_substitution": "Single EUR matrix for mixed-ancestry meta-analysis or a weighted r2 diagnostic approximation as signed-R fine-mapping input",
            "evidence_ids": ["gbmi_slalom_europepmc_xml", "gbmi_resources_current"],
        },
        {
            "resource_id": "FINNGEN_PUBLIC_SISU_LD", "status": "R2_DPRIME_BROWSER_NOT_SIGNED_COHORT_R",
            "native_build": "Require exact release/reference confirmation for future dense signed matrix",
            "ld_sample_n": None, "panels_documented": {"sisu3": "through DF7", "sisu4": "DF8-DF9", "sisu42": "DF10 onward"},
            "matrix_type": "Public ld_server parse_ld returns query-to-window r2 and d_prime, not signed Pearson r",
            "access_route": "Official FAQ says all-pairwise bulk download can be arranged on request; no contact made",
            "blocker": "Exact R11 phenotype/cohort sample match, signed dense matrix, index, allele coding and N not verified; SISu imputation panel is not automatically study-sample LD",
            "evidence_ids": ["finngen_ld_browser", "finngen_ld_server_readme", "finngen_ld_server_source", "finngen_pairwise_ld_access", "finngen_methods_markdown"],
        },
        {
            "resource_id": "ONEKG_PHASE3", "status": "PUBLIC_SAMPLE_METADATA_VERIFIED_NOT_PRIMARY_LD_DEFAULT",
            "release": "20130502 integrated_call_samples_v3 panel; phase 3", "native_build": "GRCh37; GRCh38 reanalysis/high-coverage releases are distinct",
            "sample_counts_from_public_panel": one_kg_n,
            "matrix_type": "No matrix obtained/computed; genotype-derived signed LD would need a separately frozen construction contract",
            "limitations": ["Small ancestry/subpopulation reference compared with biobank GWAS", "Reference size requirements depend on study size/LD/variant frequency", "JPT is only 104 samples in this release", "Hard-call versus dosage/imputation differences", "New 3202 high-coverage release includes additional related samples, not an interchangeable 2504 panel"],
            "proposed_use": "At most an explicitly justified future sensitivity; no default primary substitute for absent matched LD",
            "evidence_ids": ["onekg_sample_panel", "thousand_genomes_igsr", "reference_size_pmc_html"],
        },
    ]
    mapping = [
        {"dataset": "Sakornsakolpat 2019 / GCST007692", "preferred_ld": "Original cohort/meta-analysis-compatible LD", "status": "MATCHED_LD_NOT_ESTABLISHED", "reason": "PanUKB EUR is not automatically matched to the original heterogeneous cohort/meta-analysis"},
        *[{"dataset": "Kim UKB marginal COPD " + stratum, "accession": accession,
           "source_native_build": "GRCh37", "cases": cases, "controls": controls,
           "preferred_ld": "Exact subgroup/analysis-sample and covariate-compatible UKB dosage LD",
           "status": "STRATUM_MATCHED_LD_NOT_VERIFIED",
           "reason": "Full PanUKB EUR cohort LD cannot be labeled exact in-sample for smoking-selected subgroups or another association test; possible reference use needs new justification/QC",
           "metadata_provenance": "Root's all-104 source audit; no source GWAS rows read here",
           "source_doi": "10.1093/aje/kwaa227"}
          for accession, stratum, cases, controls in (
              ("GCST90016588", "ever-smoker", 12446, 59145),
              ("GCST90016589", "never-smoker", 8631, 120544),
              ("GCST90016593", "current-smoker", 4589, 10001),
              ("GCST90016594", "noncurrent-smoker", 16488, 169688))],
        *[{"dataset": "Ishigaki BBJ direct COPD / hum0014.v17.COPD.v1 " + sex,
           "accession": accession, "source_native_build": "GRCh37",
           "preferred_ld": "Exact eligible Japanese BBJ dosage LD and index for matching sex/phenotype sample",
           "status": "SUMMARY_ACCESS_VERIFIED_BY_ALTERNATIVES_AUDIT_LD_PENDING",
           "reason": "Public signed GWAS access does not prove public original LD availability or matching sex/subgroup sample",
           "metadata_provenance": "Root and alternatives source audits; no source GWAS rows read here"}
          for accession, sex in (("GCST90013709", "bothsex"), ("GCST90013746", "male"), ("GCST90013781", "female"))],
        {"dataset": "PanUKB ICD10 J44 EUR", "preferred_ld": "PanUKB EUR ldadj matrix and indexed alleles", "status": "BEST_DOCUMENTED_PUBLIC_SAME_PROJECT_ROUTE_CONDITIONAL", "reason": "Secondary EHR phenotype; release/scaling/sample/statistic/variant/locus QC still required"},
        {"dataset": "PanUKB J44 AFR/CSA or other ancestry strata", "preferred_ld": "Corresponding PanUKB ancestry matrix", "status": "PUBLIC_METADATA_AVAILABLE_PHENOTYPE_POWER_AND_LOCUS_QC_PENDING", "reason": "Do not replace with EUR; population-level N is not phenotype case/control N"},
        {"dataset": "GBMI COPD EUR/AFR/EAS/AMR/mixed", "preferred_ld": "Cohort/ancestry-compatible signed sufficient-statistic covariance", "status": "META_ANALYSIS_LIKELIHOOD_AND_LD_PENDING", "reason": "gnomAD weighted-r2 SLALOM precedent is diagnostic, not a matched likelihood"},
        {"dataset": "FinnGen R11 J10_COPD", "preferred_ld": "Release/phenotype-matched Finnish signed dosage LD", "status": "PUBLIC_SIGNED_DENSE_LD_NOT_VERIFIED", "reason": "SISu browser yields r2/d_prime, not required matrix"},
    ]
    gates = [
        ("phenotype", "Exact disease/trait definition, direction, controls, tier, ancestry, sex/subgroup, source release; no weaker-phenotype relabeling"),
        ("full_file", "Genome-wide signed statistics include null measurements; complete download/hash/schema/row validation, not header/range proof"),
        ("signed_statistic", "Verified effect allele and scale; finite beta/valid positive OR; positive SE or valid signed statistic; P encoding/test and N resolved"),
        ("identity", "Native assembly/reference and exact allelic contrast verified; deterministic duplicate/palindromic/indel rules; no rsID-only join"),
        ("ld_origin", "LD release, cohort/subgroup, ancestry, sample N, imputation, dosages/covariates and signed counted allele documented"),
        ("density", "All retained variants share exact ordered signed LD; missingness/attrition and lead/strong-variant losses audited; no universal coverage percent invented"),
        ("ld_numerics", "Triangle reconstruction, scaling/diagonal, symmetry, range, PSD/rank/condition and distance bands verified under frozen precision tolerances; no silent repairs"),
        ("sample_compatibility", "Per-variant case/control/nominal/effective N, test and cohort weights match approved likelihood; no automatic Neff substitution"),
        ("consistency", "Future z/R conditional diagnostics reviewed; no universal lambda threshold or fit-improving allele flip"),
        ("overlap", "Cohort overlap/repeated releases/controls cataloged; no unsupported independent replication or meta-combination"),
        ("locus_boundary", "GWAS-only deterministic intervals, MHC/complex regions excluded/deferred, boundary/long-range coverage checked"),
        ("execution_freeze", "User authorization plus pinned packages/binaries/environment/options/seeds and full signed input/LD freeze"),
        ("postfit", "Future convergence, cap sensitivity, CS purity, search/ELBO and instability assessed; method agreement does not fix input errors"),
    ]
    config = {
        "primary_method": "SuSiE-RSS multi-signal, conditional on all required input/LD gates",
        "sensitivity_method": "FINEMAP same-input prespecified sensitivity, not automatic consensus or LD repair",
        "prior": "Uniform variant prior; no functional/model/benchmark/candidate-derived priors",
        "proposed_primary_max_signals": 10,
        "signal_cap_interpretation": "Capacity/design precedent, not known number of causal variants; unresolved saturation requires separate action",
        "possible_prospective_cap_sensitivity": [5, 20],
        "cap_sensitivity_status": "DESIGN_EXAMPLE_NOT_AUTHORIZED_EXECUTION",
        "credible_set_coverage": 0.95,
        "documented_susie_cs_min_abs_r": 0.5,
        "purity_interpretation": "CS reporting filter, not panel eligibility threshold; r, not r2",
        "residual_variance": "Off for reference LD; exact in-sample/statistic match required before considering estimation",
        "observed_documentation_versions": {"susieR": "0.16.6", "FINEMAP_official_site": "1.4.2"},
        "execution_versions": "NOT_INSTALLED_OR_PINNED_FOR_EXECUTION; must freeze actual package/source or binary hashes in future",
        "new_mismatch_corrections": "Not silently adopted from changing documentation",
        "sources": ["susie_rss_paper", "susie_rss_reference", "susie_diagnostics", "finemap_official_http", "bbj_finemapping_methods"],
    }
    loci = {
        "selection": "GWAS own validated association test only; no candidates/genes/models/benchmarks",
        "proposed_genomewide_p_threshold": 5e-8,
        "proposed_initial_window_bp": 3000000,
        "window_rule": "plus/minus 1.5Mb, chromosome clipped, transitive overlap merge; dense nonsignificant rows retained",
        "tie_rule": "association significance; chromosome; position; exact allelic key; source row ID",
        "mhc_mask_proposal": {"chrom": "6", "start1": 25000000, "end1": 36000000,
                              "builds": ["GRCh37", "GRCh38"], "interpretation": "Deliberately conservative same-number native-build operational mask; NOT asserted exact liftover; exclude whole intersecting merged interval",
                              "source_grch38_finngen_mask": "chr6:25-34Mb inclusive", "wider_precedent": "BBJ/FinnGen/UKBB methods chr6:25-36Mb"},
        "complex_loci": "Defer if long-range/multisignal/structural/boundary/size conditions exceed frozen validated plan; no memory-motivated trimming",
        "initial_chromosomes": "autosomes only; sex-chromosome contract separate",
        "expected_eligible_loci": None,
        "expected_loci_status": "NOT_DETERMINED; no locus scan permitted/performed",
    }
    result = {
        "schema_version": "1.0", "stage": STAGE.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "audit_status": "DOCUMENTATION_AND_METADATA_AUDIT_COMPLETE_WITH_EXPLICIT_UNRESOLVED_GATES",
        "scope": {"LD_matrix_values_downloaded": False, "LD_matrices_computed": False,
                  "fine_mapping_executed": False, "scientific_locus_qc_executed": False,
                  "GWAS_candidate_rows_or_model_benchmark_files_opened": False,
                  "prior_stage_artifacts_modified": False, "external_contact_sent": False,
                  "public_documentation_code_metadata_and_sample_manifest_only": True},
        "readiness_warning": "Access to a schema or metadata is not complete input validation and is not inference readiness",
        "resources": resources, "dataset_ld_recommendations": mapping,
        "gates": [{"gate_id": k, "requirement": v, "scientific_status": "NOT_RUN_PREFLIGHT_ONLY"} for k, v in gates],
        "method_design": config, "locus_design": loci,
        "numeric_rule_provenance": [
            {"rule": "P<5e-8", "role": "Conventional GWAS-only selection proposal", "source_id": "finngen_methods_markdown"},
            {"rule": "3Mb total window", "role": "Initial design precedent, not proof of LD containment", "source_id": "finngen_methods_markdown"},
            {"rule": "10 maximum effects;95%CS", "role": "Published pipeline design precedent and conditional posterior coverage", "source_ids": ["finngen_methods_markdown", "susie_rss_reference"]},
            {"rule": "CS minimum abs(r)=0.5", "role": "Documented method reporting filter, not input eligibility", "source_id": "susie_rss_reference"},
            {"rule": "INFO/MAC PanUKB", "role": "Resource filters; current source/documentation discrepancy unresolved", "source_id": "panukb_ld_docs"},
            {"rule": "LD overlap %, lambda, frequency delta, PSD epsilon", "role": "NO universal hard readiness cutoff established or fabricated", "source_id": None},
        ],
        "compute_storage_plan": {
            "gpu_required": False, "runtime_measured": False, "expected_locus_count": None,
            "dense_float64_matrix_bytes_formula": "8*p*p",
            "illustrative_dense_matrix_bytes": {str(p): 8*p*p for p in (5000,10000,20000,50000)},
            "planning_workspace_matrix_copies": "3-6 provisional engineering headroom, not measured guarantee",
            "future_job_constraints": "CPU-only locus-bounded jobs, memory-limited concurrency; no whole 43.3TB matrix download or silent locus truncation",
            "cloud_egress_request_cost": "UNKNOWN; public access is not proof of zero compute/egress cost",
            "total_cpu_hours_or_wall_time": "NOT_ESTIMATED_WITHOUT_AUTHORIZED_LOCI_AND_RESOURCE_BENCHMARK",
        },
        "source_access_counts": dict(Counter(r["status"] for r in records)),
        "source_body_bytes": sum(r["saved_bytes"] for r in records),
        "source_count": len(records), "source_records": records,
        "source_byte_integrity": "PASS",
        "collection_receipts": [record_file(x) for x in sorted(AUDIT.glob("source_collection_*.json"))],
        "design_documents": [record_file(AUDIT / x) for x in ("RISK_DIRECTION_CONTRACT.md", "METHOD_AND_ELIGIBILITY_DESIGN.md")],
        "scripts": [record_file(x) for x in sorted((STAGE / "scripts").glob("methods*.py"))],
        "collection_commands": [
            "python scripts/methods_collect_sources.py --run-id initial",
            "python scripts/methods_collect_followup_sources.py --run-id followup",
            "python scripts/methods_collect_contract_sources.py --run-id contract",
            "python scripts/methods_collect_primary_fallback.py --run-id primary_fallback",
            "python scripts/methods_build_audit.py",
        ],
        "command_interpretation": "Paths are stage-relative; collection is public metadata/documentation GET only; this build performs byte checks and metadata counts, not statistical analysis",
    }
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["audit_status"], "source_count": len(records),
                      "source_body_bytes": result["source_body_bytes"],
                      "source_integrity": "PASS", "output": str(output.relative_to(STAGE)),
                      "output_sha256": record_file(output)["sha256"]}, indent=2))


if __name__ == "__main__":
    main()
