# COPD V2 fine-mapping and risk-direction preflight

Stage: `fine-mapping-preflight-1.0`. Audit date: 2026-10-06 UTC. Baseline: `1f41df67dcf18d626a3862a38877919d2d3114a4`.

## Decision

**Seven direct-COPD accessions have verified publicly accessible signed summary-statistic schemas, but no dataset or locus is yet execution-cleared for fine-mapping.** These are four European UK Biobank smoking-stratum GWAS from Kim and three Japanese BBJ sex strata from Ishigaki. Availability of the full-release endpoint and representative records does not establish complete-file integrity, locus density, exact signed-LD matching or statistical compatibility. The expected number of eligible loci remains **unknown**; zero execution-cleared loci is not a finding of zero COPD association loci.

There is no currently fully ready primary dataset. The best practical **conditional European primary track** is Kim's directly spirometric COPD in ever smokers, **GCST90016588**, with **GCST90016589** retained as a separate never-smoker stratum. This recommendation is based on direct ascertainment, actual signed-file access and the prospect of resolving UKB-derived LD—not on candidate overlap, functional benchmarks, model performance or observed fine-mapping results. It is explicitly an ever-smoker/never-smoker analysis, not a replacement labelled as the broader Sakornsakolpat meta-analysis. Current/noncurrent strata overlap the first partition and are separate sensitivities.

Ishigaki's original physician-diagnosed Japanese COPD **GCST90013709** is a conditional East Asian primary track if compatible dense Japanese signed LD can be obtained. Male/female releases are sex-stratum sensitivities, not independent replications of the combined cohort. Original BBJ clinical ascertainment remains the frozen moderate-confidence primary classification; it is not silently equated with later EHR/text-mined BBJ phenotypes.

The strongest documented public **secondary** route is **Pan-UKB ICD10 J44 EUR**, with 11,536 cases and 408,995 controls and publicly documented same-project dosage LD. It remains an EHR/ICD sensitivity, not direct clinically/spirometrically adjudicated COPD. LD release/scaling, sample and test compatibility still need validation. No technically convenient source has been promoted into the primary phenotype class.

The preferred broad direct-COPD anchor, **Sakornsakolpat 2019 / GCST007692**, is **not executable from the public material retrieved**. An author/consortium request for original complete signed statistics, preferably ancestry-specific with compatible LD information, is required. Its small public dbGaP display export is incomplete and unsigned. No investigator request was sent.

Machine-readable conclusions are in [recommendations](provenance/recommendations.json) and the [104-accession final eligibility matrix](tables/study_final_eligibility_matrix.tsv).

## Audit population and evidence coverage

All 104 frozen core COPD accessions were audited, not only studies supporting V1 candidates. Their classes remain unchanged: 27 direct COPD, 41 EHR COPD, one ML-liability, six chronic bronchitis, two smoking-interaction coefficients, two mixed EHR/self-report, two exacerbation/infection endpoints, 16 gene/burden/CNV designs, one progression/decline, one within-case severity and five unresolved. Non-COPD lung-function results are contextual evidence only; they are not added to a COPD fine-mapping denominator.

The current Catalog v2 API returned metadata for all 104. Standard or advertised summary-statistics directories returned HTTP200 for 69 and HTTP404 for 35. **68 accessions had statistic-file headers and representative records inspected, across 138 raw/formatted/harmonized files.** These counts include aggregate and nonprimary designs; they are not counts of fine-mappable GWAS. The remaining HTTP200 directory was CKB's suspension notice. A Catalog flag, directory response or sample record is not a usability pass.

The [Catalog technical matrix](tables/study_phenotype_summary_stat_eligibility.tsv) retains all requested fields: publication, phenotype/controls, cohorts, ancestry and current N descriptions, genotyping context, native/harmonized metadata builds, identifiers, effect/other alleles, effect/SE/P/AF/N/QC columns, source URLs, exact inspected prefixes and hashes, and explicit unknowns. [Header audit](tables/header_schema_audit.tsv) and [file inventory](tables/source_accession_file_inventory.tsv) preserve each file rather than collapsing raw and harmonized representations. The **final** matrix integrates provider evidence, notably BBJ; a Catalog-only “no file found” is not its final availability decision.

Independent source-specific reviews cover all 27 direct accessions and 25 provider dataset/ancestry records. These are overlapping views, not 129 independent studies. Their scope limitations are recorded row by row: some historical accessions received current API/FTP checks and frozen publication adjudication, not a newly exhaustive worldwide file search. Failure to verify a source is not proof that it is globally unavailable.

Public acquisition receipts retain URL, UTC, HTTP/error state, byte range, local path, size and SHA-256. The [unified ledger](provenance/source_access_ledger.tsv) contains **802 saved responses, 64,623,735 bytes**. Prefix hashes cover only saved bytes, not the entire advertised remote file. Four small dbGaP exports were downloaded completely for release-scope/schema auditing; no complete genome-wide GWAS or LD matrix values were downloaded. Source HTML/error bodies are distinguished from scientific payloads. The source manifest and final freeze give authoritative counts/hashes.

## Direct-COPD signed sources

| Accession | Exact analysis | Cases / controls | Verified signed fields | Remaining critical requirement |
|---|---|---:|---|---|
| GCST90016588 | UKB ever-smoker spirometric COPD | 12,446 / 59,145 | Effect/other alleles, OR, SE, P; native GRCh37 | Stratum-compatible LD, SE/test lineage, whole-file and locus QC |
| GCST90016589 | UKB never-smoker spirometric COPD | 8,631 / 120,544 | Same native schema | Separate estimand; matched LD and QC |
| GCST90016593 | UKB current-smoker spirometric COPD | 4,589 / 10,001 | Same native schema | Overlapping-partition sensitivity; matched LD and QC |
| GCST90016594 | UKB noncurrent-smoker spirometric COPD | 16,488 / 169,688 | Same native schema | Overlapping-partition sensitivity; matched LD and QC |
| GCST90013709 | Original BBJ physician-diagnosed COPD, both sexes | 3,315 / 201,592 | Allele2=ALT effect; BETA/SE/P, AF/N/quality; GRCh37 | Dense Japanese cohort-compatible signed LD, full-file and locus QC |
| GCST90013746 | Original BBJ male COPD | 2,855 / 103,089 | Male BETA.x/SE.x/p.value.x, Allele2 effect | Sex-specific LD/sample contract; no per-variant N in sex file |
| GCST90013781 | Original BBJ female COPD | 460 / 98,503 | Female BETA.y/SE.y/p.value.y, Allele2 effect | Small case count, sex-specific LD/sample contract; no per-variant N |

The [seven-source table](tables/accessible_signed_direct_COPD_sources.tsv), [direct review](audits/direct/DIRECT_COPD_SOURCE_AUDIT.md) and [provider review](audits/alternatives/ALTERNATIVES_AUDIT.md) retain exact evidence and limitations.

Kim's direct case phenotype is pre-bronchodilator FEV1/FVC <0.7 and FEV1 <80% predicted. Ordinary smoking-stratum logistic GWAS is distinct from the two interaction-model main-effect accessions. Raw files omit AF, per-variant N and imputation-quality values. Placeholder columns in harmonized files do not cure that missingness. The paper/PLINK documentation supports a log-OR SE interpretation, but the Catalog conversion lineage and future numerical test consistency remain prerequisites. Raw and harmonized effect alleles can be reversed with reciprocal OR; their columns must never be mixed. [Kim methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8096488/), [actual ever-smoker directory](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016588/).

The original BBJ ZIP is **2,244,311,450 bytes** and was sampled by bounded HTTP206 ranges, not fully downloaded. Archive directory offsets allowed actual all-sex and sex-stratified autosomal/chrX headers to be decoded. The provider README resolves allele and male/female suffix semantics; it also distinguishes a sex-meta coefficient from either sex-specific coefficient. The initial proposed fine-mapping design is autosomal only; checking a chrX header did not authorize chrX inference. JENGER timeouts and an absent Catalog mirror do not negate the successful original-provider route. [BBJ README](https://humandbs.dbcls.jp/files/hum0014/hum0014_v17_v18_v21_README.txt), [sample counts](https://humandbs.dbcls.jp/files/hum0014/hum0014_v17_v18_v21_sample_size.xlsx), [archive schema evidence](audits/alternatives/bbj_archive_schema_evidence.json).

## Why Sakornsakolpat and other important anchors remain blocked

The Sakornsakolpat paper reports 35,735 cases, 222,076 controls and 6,224,355 tested variants. Its publicly exposed `pha004766` file has **19,373 records** and explicitly absolute coefficients. It is not a complete signed meta-analysis release. Effect-allele display inconsistencies, blank sample-size fields, and display-versus-submitted build differences add further obstacles. The retained original file and its dictionary, not a Catalog availability flag, establish these limitations. [Primary study](https://www.nature.com/articles/s41588-018-0342-2), [public dbGaP export](https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000179/analyses/phs000179.pha004766.txt).

The main discovery combines heterogeneous ancestry cohorts, albeit predominantly European; this is supported by the primary methods and the specific dbGaP analysis page, not inferred solely from aggregate ancestry labels. The paper's use of a 10,000-person UKB LD reference does not provide an openly verified matched matrix for this preflight. A future single-European substitution would require scientific justification and cannot stand in for ancestry-specific/meta-analysis-compatible LD. [dbGaP analysis provenance](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?study_id=phs000179.v6.p2&pha=4766).

Hobbs 2017 exports similarly contain only 25,412 / 25,338 / 25,350 records and explicitly absolute coefficients. Selected Stage-2 replication must not be treated as a dense final genome-wide analysis. Joo's male/female lead tables do not replace complete signed files. Moll's exonic design does not cover dense regulatory loci; its figshare appendix is a DOCX, not a full common-variant release. Older Cho and other unresolved sources require original-release discovery/request and build/LD clarification. These are source-specific limitations, not negative COPD findings. [Direct-source audit and exact export hashes](audits/direct/DIRECT_COPD_SOURCE_AUDIT.md).

Published locus totals are contextual only: Sakornsakolpat reports 82; Hobbs reports 13 Stage-1 and 22 after selected follow-up; Joo reports 17 male and 14 female; Kim's 48/55 are joint interaction-test loci, **not** its four marginal-stratum counts. None is substituted for a prospective eligible-locus count. No Catalog-hit reconstruction or candidate-selected locus enumeration was performed.

## Secondary resources and ancestry-specific cautions

| Resource | Verified public evidence | Disposition |
|---|---|---|
| Pan-UKB J44 EUR | ALT beta/SE/negative-log10-P, exact phenotype manifest, public same-project LD metadata | Best conditional secondary; preserve EUR columns, not cross-ancestry meta columns |
| Pan-UKB J44 AFR / CSA | Corresponding ancestry-specific fields and sample counts | Separate small-case strata; power/QC limitations and matching ancestry LD remain explicit |
| Pan-UKB 496 / 496.21 | Actual files and ICD/PheCode mapping | Broad airway obstruction/bronchitis-related sensitivity; not silently pooled with J44 or direct COPD |
| GBMI AFR/AMR/EAS/EUR/mixed | Official source workbook, actual signed file prefixes, per-variant cohort/sample fields | Secondary; ancestry-specific meta-analysis still needs appropriate covariance; no single-EUR LD for mixed results |
| MVP PheCode 496 | Released ancestry-specific/mixed Catalog prefixes | Secondary; sampled EUR SE values are NA despite the column; OR/CI alone are not automatically beta/SE-ready |
| FinnGen R11 J10_COPD | Signed release-specific prefix and manifest: 21,617 / 372,627 | Finnish registry sensitivity; exact phenotype code freeze and dense signed release-matched LD unresolved |
| UKB WGS J44 GCST90473707 | Signed beta/SE and N/INFO fields | Secondary; older imputation LD cannot be presumed to cover WGS rare variants |
| GEMINI v1 COPD | COPD-specific README and signed release | UKB+FinnGen EHR meta-analysis, not direct-COPD consortium data; overlap/covariance unresolved |
| CKB GCST90246122 | Catalog suspension notice plus current provider access policy | Encrypted downloads require decryption-key application; no keys requested or restriction bypassed |

Provider definitions, builds, sizes, sample counts, allele/P encodings and source-specific failures are fully retained in the [25-row inventory](tables/provider_dataset_ancestry_inventory.tsv). In particular, negative-log10 P is not P; auxiliary ALT is not necessarily the harmonized effect allele; a standard-error column containing NA is not usable SE; and public documentation or published fine-mapping results do not prove access to the original matching LD. [Pan-UKB schema](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html), [CKB access terms](https://pheweb.ckbiobank.org/about).

## LD feasibility and risk-direction contract

No study has passed an exact matching signed-LD and locus-QC contract in this preflight. The [resource inventory](tables/LD_resource_inventory.tsv), [source compatibility matrix](tables/LD_source_compatibility_matrix.tsv) and [104-study LD requirements](tables/study_LD_requirements.tsv) separate preferred LD from actually verified access and compatibility.

Pan-UKB is the best-documented public same-project route for its own ancestry-specific EHR GWAS. Its matrices are signed dosage-based, covariate-adjusted, triangular and 10-Mb banded. Metadata and the keyed variant index were inspected, not matrix values. Exact matrix N, phenotype-subset/test matching, released-source provenance and scaling remain open; the archived calculation's normalization order requires diagonal/scale reconciliation before treating the matrix as Pearson R. Full-population Pan-UKB LD is not exact in-sample LD for Kim's smoking subsets merely because both are UKB. [Pan-UKB LD documentation](https://pan.ukbb.broadinstitute.org/docs/ld/index.html), [detailed audited contract](audits/methods_ld/METHOD_AND_ELIGIBILITY_DESIGN.md).

No openly accessible dense original BBJ matrix with an adequate matching contract was verified. The original study also drew Japanese controls from ToMMo/IMM/JPHC/J-MICC; “BBJ-matched” must therefore resolve this pooled sample/control contract, not assume every participant belongs to BBJ. GBMI's weighted-r² diagnostic precedent is not the required signed meta-analysis likelihood. FinnGen's public browser/API exposes r²/D-prime, not a verified dense signed Pearson matrix. The V1 1000 Genomes proxy resource was neither reused nor recalculated; small external ancestry panels are at most a separately justified sensitivity, not a default workaround.

The [risk-direction contract](audits/methods_ld/RISK_DIRECTION_CONTRACT.md) requires signed, correctly oriented disease association: positive beta/log(OR), or OR >1, means the encoded effect allele is disease-increasing; negative beta, or 0<OR<1, reverses the biallelic contrast. Zero/tie, unknown outcome coding, unsigned coefficients or ambiguous allele identity yield no increasing-allele assignment. Trait-increasing ML/lung-function/severity effects stay in their own estimand. No per-variant risk alleles were assigned here.

Risk cannot come from REF/ALT, major/minor status, frequency, ancestral state, LD alone, model direction, GTEx NES, nearest gene or biological expectation. Exact build/reference/allele harmonization must handle swaps, strand complements, palindromes, multiallelic contrasts, indels and duplicate conflicts. A verified allele swap transforms both signed z and both LD axes consistently; r² cannot supply missing signs. No sign flip may be chosen solely because it improves summary-statistic/LD agreement.

## Prospective method, locus rules and gates

The proposed method is **multi-signal SuSiE-RSS**, with **FINEMAP** as a same-input prespecified sensitivity. Uniform variant priors preserve the candidate/model/benchmark firewall. A high PIP is statistical prioritization conditional on data and model, not proof of biological causality. [SuSiE-RSS primary method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299).

The [13 prospective gates](tables/proposed_fine_mapping_eligibility_gates.tsv) require phenotype validity, locus-complete signed statistics, verified effect/test/SE/N, exact variant/build/allele consistency, sufficiently dense overlap with appropriate signed LD, matrix integrity, and diagnostics for summary-statistic/LD inconsistencies. Missingness and excluded-variant ledgers must expose loss of strong or boundary signals. No arbitrary percentage-overlap or mismatch cutoff is invented to make a source pass. Actual numerical tolerances and software/configuration hashes must be frozen before future execution.

The candidate-blind proposal selects each dataset's own valid **P <5×10⁻⁸** associations, initially uses **±1.5-Mb windows**, and transitively merges overlaps while retaining eligible nonsignificant variants. It proposes autosomes first, a conservative native-build chr6:25–36-Mb complex-MHC exclusion, deterministic ties and explicit boundary/long-range-LD checks. This operational mask is not claimed to be an exact cross-build liftover. No locus has been selected under these rules. [Rules](tables/proposed_locus_definition_rules.tsv), [numerical-rule evidence](tables/numerical_design_rule_evidence.tsv).

The design proposes up to ten signals and nominal 95% credible sets, with prespecified cap/convergence/purity/boundary sensitivity—not outcome-driven choice. Detailed source justifications, the distinction between absolute r and r², reference-LD residual-variance restrictions, software/version pinning and failure handling are in the [method design](audits/methods_ld/METHOD_AND_ELIGIBILITY_DESIGN.md). These are recommendations frozen as a preflight result, not approval to execute or permissive replacements for missing inputs.

## Overlap, compute and storage

[Cohort-overlap considerations](STUDY_OVERLAP_CONSIDERATIONS.md) preserve UKB, BBJ and consortium reuse; smoking partitions, sex strata and meta-analysis components must not be counted as independent replication. Repeated Catalog YAML sample entries must not be summed. N_cases/N_controls, nominal N, per-variant N and any likelihood-specific effective N remain distinct.

Future fine-mapping is CPU-only. A dense float64 matrix needs `8p²` bytes: 10,000 variants require 0.8 GB and 50,000 require 20 GB for one matrix, before multiple working copies and fit state. These are engineering estimates, not measured runtimes. Exact CPU hours and total storage cannot be responsibly estimated without eligible loci and the final resource contract. Source-file sizes and formulas bound planning: BBJ's original archive is 2.24 GB; Pan-UKB J44 is 1.83 GB; the full Pan-UKB LD release is documented as 43.3 TB and must not be copied wholesale. Future access should be approved locus-block extraction, never silent region truncation to fit memory. [Compute/storage plan](provenance/compute_storage_plan.json).

## Deliverables, validation and stop

The bundle includes all source/audit scripts and receipts; the 104-study final and Catalog-detail matrices; source/header/risk-direction/ancestry tables; 27 direct-source and 25 provider reviews; LD inventory/compatibility; source-specific limitations; proposed gates/locus/method rules; overlap and compute plans; [priority sources for future validation](tables/priority_datasets_for_future_input_validation.tsv); [sources requiring access or signed-file resolution](tables/datasets_requiring_access_or_signed_source_resolution.tsv); and [unsuitable or separate-estimand designs](tables/datasets_unsuitable_or_separate_estimand.tsv).

The [immediately executable fine-mapping table](tables/immediately_executable_fine_mapping_datasets.tsv) intentionally has no data rows. That empty list is not an assertion that no usable public route exists: it distinguishes promising accessible inputs from completed LD/locus prerequisites. All seven direct signed releases and the conditional Pan-UKB secondary route are explicit priorities for investigator review.

Independent final validation is recorded in [its receipt](provenance/independent_validation.json); the [freeze](provenance/freeze.json) and [artifact checksum ledger](provenance/artifact_checksums.tsv) bind the completed bundle. Earlier tracked artifacts must match the baseline, with only append-only V2 activity/decision/result updates. The validation checks integrity and evidence/contract consistency, not unexecuted locus scientific adequacy.

**Stop for investigator review.** No statistical fine-mapping, full-GWAS association/locus scan, V2 scoring, 337-candidate reranking, benchmark re-evaluation, LD calculation, retraining, retuning, colocalization, target-gene analysis, manuscript revision, external contact, commit or push was performed. No earlier scientific stage was revised.
