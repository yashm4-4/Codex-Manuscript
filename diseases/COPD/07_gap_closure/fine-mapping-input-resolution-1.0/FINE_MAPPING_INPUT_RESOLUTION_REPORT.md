# Fine-mapping input resolution / execution readiness 1.0

**All three tracks are NOT CLEARED. Zero loci are execution-cleared, and future SuSiE-RSS/FINEMAP execution is not authorized by this evidence.**

This completed stage resolves full-file inputs and performs diagnostics only. It contains no PIPs, credible sets, posterior effects, candidate/model scores, functional prioritization, target-gene assignments, colocalization or manuscript revision. It stops for investigator review. Publication of this new stage is not authorized; no commit or push was performed.

The governing design is the frozen `../fine-mapping-preflight-1.0/`. The exact new authorization is `provenance/user_authorization.txt`. Baseline Git commit is `54dfd85d457ea4370a497af88bb90ebe10cf26d1`. All previous stages are immutable; only the three shared registers receive new appended rows.

## Readiness counts

| Track | Complete GWAS rows | Merged loci / non-MHC | All summary gates | All LD gates | Execution-cleared | Final state |
|---|---:|---:|---:|---:|---:|---|
| A | 9,886,853 | 11 / 10 | 0 | 0 | 0 | NOT CLEARED |
| B | 8,678,470 | 5 / 5 | 0 | 0 | 0 | NOT CLEARED |
| C | 28,987,534 | 7 / 6 | 0 | 0 | 0 | NOT CLEARED |

There are **23 track-specific merged GWAS-defined intervals: 21 non-MHC and two deferred MHC intervals**; overlapping intervals across tracks are not unique genomic regions. “All summary gates” means a complete method-ready statistic, identity and sample-size contract. It does not mean absence of valid signed associations. “All LD gates” includes origin, density, numerical behavior, sample relationship and summary-statistic covariance compatibility. A basic numerical pass alone does not close this gate. Known unmet prerequisites justify NOT CLEARED for every track; no track is labeled execution-cleared by default or by technical convenience.

Track A remains the primary practical direct-COPD track; Track B is the independent Japanese ancestry track; Track C remains a secondary EHR/ICD sensitivity. A and C share UK Biobank participants and cannot be presented as independent replication. BBJ’s pooled Japanese case/control composition also requires an appropriate LD/sample contract. See `COHORT_OVERLAP_AND_SCOPE.md`.

## Track A: Kim ever-smoker COPD, GCST90016588

The complete native GRCh37 file contains **9,886,853 rows**, with the published provider MD5 verified. Complete formatted and harmonized Catalog releases were audited only to resolve transport/effect semantics; the native file remains authoritative. Complete-file acquisition, checksums and schema/missingness evidence are indexed in `tables/full_file_acquisition.tsv`, `tables/complete_file_integrity.tsv` and `tables/whole_file_audit_inventory.tsv`.

The source OR refers to the encoded effect allele. `beta = log(OR)` and supplied log-OR SE are supported by source methods/software semantics and the whole-file audit. All **1,623 significant rows** pass the frozen decimal-rounding consistency contract. Across the full file, 9,861,764 pass and **25,089 nonsignificant rows fail** that prospective precision gate, including 436 locus rows. Some failure strings display floating-point serialization tails; the numerical cause is not established for every failure. All failures are retained, not repaired, and the tolerance was not relaxed. Raw/harmonized allele reversal and reciprocal-OR/sign behavior were checked; unmatched transport records were not joined by rsID alone.

The study’s 12,446 cases and 59,145 controls give nominal N=71,591. Per-variant N, AF and INFO are absent in the authoritative native schema. The nominal total was recorded separately, not silently substituted for per-variant or likelihood-specific effective N.

Eleven merged intervals include one deferred MHC interval, leaving ten non-MHC loci. All-P locus tables preserve 170,904 rows, of which 108,870 are non-MHC. Reference checks, proven allele swaps/complements and uniquely normalized indels are recorded. There are 23,664 prospective palindromic rows with unresolved Kim strand/frequency evidence, including the chr2 lead. Duplicates and ambiguous indels remain explicit.

For the ten non-MHC loci, exact matching yields **91,193 ordered diagnostic variants** from 92,494 eligible summary rows; **1,301 eligible rows are absent from the LD index**. Across the source significant rows, **224/1,522 are missing or excluded**, and **9/10 leads are retained**. Nine loci lose at least one significant row. The chr4:105,262,244–108,341,962 locus has no significant identity/LD loss but still lacks a closed sample/statistic/reference-LD contract. All ten retained matrices were diagnosed; none is promoted to an execution input.

The LD source is `UKBB.EUR.ldadj.bm` from Pan-UKB. For Kim this is **external full-European UKB reference LD**, not matched in-sample ever-smoker LD. No exact smoking-subgroup covariance was obtained. Descriptive residual and spectral diagnostics cannot establish that smoking-subset versus full-UKB mismatch is acceptable by themselves. Principal evidence: `tracks/A/TRACK_A_REPORT.md`, `tracks/A/audit.json`, `ld/A/gwas_ld_overlap_summary.json` and `ld/LD_SOURCE_CONTRACT.md`.

## Track B: original combined BBJ COPD, GCST90013709

The complete combined-autosomal member contains **8,678,470 rows**. The original source ZIP, combined member, size/range provenance, ZIP-member CRCs and gzip integrity are preserved. A premature transfer was resumed against consistent byte-range/ETag evidence; the failed attempt remains recorded. Local SHA256 hashes bind the bytes; no nonexistent independent provider checksum is claimed.

Allele2 is the ALT effect allele and BETA supplies the signed approximate log-OR contrast. Per-variant N ranges from **204,905 to 204,907**; 11,542 rows differ from nominal total N. An initial implementation incorrectly excluded nonnominal N despite the already-frozen positive-N criterion. That implementation was corrected under the unchanged gate, and its prior artifacts are retained. The corrected final result is **five loci and 462 significant rows**, with no MHC locus. Sex-specific and chrX archive members were not analyzed as additional tracks.

All 52,013 prospective locus rows are retained and their REF alleles match GRCh37. **51,989 signed directions are verified**; 24 symbolic-allele rows remain unresolved. A conflicting four-row duplicate identity is outside these loci and remains documented.

Printed SE is related to the SPA-adjusted test rather than an independently validated ordinary Wald likelihood. Agreement of beta/SE with SPA P therefore does not close the Gaussian/RSS covariance contract. A contemporaneous source-score reconstruction fails 338/9,141 deterministic sampled rows, leaving exact software/score-column lineage unresolved. No score, P value, SE or allele orientation was repaired.

No legitimate, currently accessible **dense signed Japanese LD matrix** closes the original BBJ plus pooled Japanese controls’ sample and variant contract. Controlled BBJ/JGA genotype/reference access is a future route, not completed access. ToMMo recombination maps and frequency-only resources are not signed LD matrices. No 1000 Genomes EAS substitute was used; r-squared was not converted to signed correlation. Consequently all five LD and execution gates remain unmet. Principal evidence: `tracks/B/TRACK_B_REPORT.md`, `tracks/B/results/audit.json` and `tracks/B/results/ld_source_contract.json`.

## Track C: Pan-UKB ICD10 J44 EUR

The complete public file contains **28,987,534 rows**, with provider MD5 and gzip integrity verified. Only the EUR J44 both-sex association was analyzed. There are 23,861,813 finite valid EUR signed-statistic rows and 18,170,886 autosomal source-QC-passing rows. EUR-missing and low-confidence rows remain accounted for. One nonmissing invalid row has SE=0 and infinite negative-log10-P, was source-low-confidence, and remains preserved as an exclusion.

The seven merged intervals contain **334 significant seeds**, six non-MHC loci and one deferred MHC interval. The chr20 interval is clipped at the GRCh37 chromosome end. All 259,252 locus rows are retained and all source REF alleles match the reference. Twenty-four normalized-identity collision rows remain unresolved; none is significant. The provider `high_quality` annotation was not a prospectively required gate and was not used to remove significant rows. No functional/gene columns from the contextual variant manifest were used.

The source association uses ALT effect, signed beta, supplied SE and **negative-log10-P**. Exact phenotype metadata identifies SAIGE 0.36.3 and 11,536 cases / 408,995 controls. The pinned source code derives SE from beta and SPA P. Whole-file beta/SE/P agreement is therefore arithmetic consistency, not independent evidence for the covariance assumed by RSS. The current phenotype metadata container and older flat-file release are distinguished; their precise execution lineage is not assumed.

The exact same-project signed LD index contains **420,542 samples**, while the J44 phenotype count is **420,531**. This small difference is documented; the number 11 is not itself an invented rejection threshold. What remains unresolved is the exact sample/release/test covariance relationship and appropriate per-variant/effective-N contract. Broader variant-manifest AN/2 is not substituted for phenotype or LD N.

All eligible summary rows and all **333 non-MHC significant seeds** match the native signed LD index. The six diagnostic matrices contain 17,801; 27,382; 20,055; 18,706; 22,575; and 19,448 variants, respectively. No significant signal is lost in these joins. This strong technical coverage does not close the sampling and likelihood gates. Principal evidence: `tracks/C/REPORT.md`, `tracks/C/results/audit.json`, `ld/C/gwas_ld_overlap_summary.json` and `ld/panukb_J44_EUR_exact_phenotype_metadata.json`.

## LD representation and diagnostics

Only blocks intersecting the GWAS-defined A/C non-MHC loci were obtained from `https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.ldadj.bm/`, with the accompanying `UKBB.EUR.ldadj.variant.ht` index. The 43-TB resource was not downloaded wholesale. Original block bytes, response receipts, index partitions and variant order are retained. B has no acquired matrix; both MHC intervals are deferred.

The provider stores a signed, upper-triangular, 10-Mb-banded, covariate-residualized Gram representation. Historical source corroborates standardization before covariate projection and division by sample N. The exact executed release configuration remains unbound: historical code and the actual chrX-containing index do not establish one exact code snapshot. This limitation is explicit in `ld/source_contracts.json`.

Each original upper triangle is preserved as `source_upper_triangle.npy`. A separate diagnostic correlation is computed explicitly as `Gij / sqrt(Gii Gjj)` after checking finite positive diagonals and the absent lower triangle. This documented representation conversion is not a nearest-PSD, ridge, sign, or off-band repair. All retained pairs lie within their locus and the 10-Mb stored radius. Missing source blocks were checked explicitly before treating sparse representation as legitimate. Exact ordered-input hashes and dimensions bind matrices to allele-indexed z values.

Float64 numerical tolerances were frozen before matrix results in `config/ld_diagnostic_contract.json`; descriptive residual specifications were separately frozen before evaluation in `config/ld_residual_contract.json`. Full eigenvalue/rank/condition checks preserve negative eigenvalues and singular directions. The null covariance mismatch scalar s and lead-conditional residuals are descriptive diagnostics only. No result-derived s/residual/coverage threshold was invented; multiple association signals can affect residuals. No matrix was fitted to improve allele signs or repaired for execution.

Basic numerical status counts: **A: 10/10; C: 6/6**. These are distinct from the zero complete scientific LD gates. The full numeric evidence, including raw diagonals, eigenvalue tolerances, rank deficiency, null-space z energy and residual distributions, is in `tables/LD_numerical_QC.tsv`, per-locus `numeric_diagnostics.json`, and `ld/extraction_validation.json`.

All 16 diagnosed matrices have a numerical rank below their dimension at the frozen threshold, and 16 exceed the positive-subspace condition warning. The descriptive mismatch s estimates range from 2.42872e-15 to 1.87188e-05. These small estimates do not demonstrate that the sample/statistic contract is closed. Conversely, this audit does not claim a gross PSD, symmetry or correlation-range failure where those tests passed. Unknown provenance/sampling requirements and the recorded variant losses remain the reasons execution is not cleared.

| Track | Locus | Matrix variants | Basic numerical status | Numerical rank | Condition warning | Descriptive s |
|---|---|---:|---|---:|---|---:|
| A | A_chr2_228001745_231092415 | 10919 | PASS_NUMERICS_ONLY | 10651 | True | 8.09353771649385e-07 |
| A | A_chr3_167160231_170450432 | 8604 | PASS_NUMERICS_ONLY | 8325 | True | 3.737695902981158e-09 |
| A | A_chr4_88300892_91560531 | 10162 | PASS_NUMERICS_ONLY | 9826 | True | 2.8489991148455975e-15 |
| A | A_chr4_105262244_108341962 | 7998 | PASS_NUMERICS_ONLY | 7709 | True | 3.1770738706075364e-06 |
| A | A_chr4_143727600_147016737 | 8008 | PASS_NUMERICS_ONLY | 7564 | True | 1.1441159906564816e-08 |
| A | A_chr5_146267580_149356522 | 8888 | PASS_NUMERICS_ONLY | 8693 | True | 3.737692033100373e-09 |
| A | A_chr6_141148056_144366387 | 8125 | PASS_NUMERICS_ONLY | 7889 | True | 2.478215106617483e-06 |
| A | A_chr15_70078560_73206624 | 7732 | PASS_NUMERICS_ONLY | 7639 | True | 3.7365686594321126e-15 |
| A | A_chr15_77212117_80575335 | 10454 | PASS_NUMERICS_ONLY | 10215 | True | 3.737689295633784e-09 |
| A | A_chr16_73811812_77016739 | 10303 | PASS_NUMERICS_ONLY | 10195 | True | 2.4287159878294812e-15 |
| C | C_L001 | 17801 | PASS_NUMERICS_ONLY | 16446 | True | 6.5365805479491615e-06 |
| C | C_L003 | 27382 | PASS_NUMERICS_ONLY | 27052 | True | 1.1147485892286572e-05 |
| C | C_L004 | 20055 | PASS_NUMERICS_ONLY | 19316 | True | 1.3060777055198132e-07 |
| C | C_L005 | 18706 | PASS_NUMERICS_ONLY | 18060 | True | 1.8718797870498033e-05 |
| C | C_L006 | 22575 | PASS_NUMERICS_ONLY | 22038 | True | 5.509437591970134e-07 |
| C | C_L007 | 19448 | PASS_NUMERICS_ONLY | 19240 | True | 5.759526484790782e-15 |

Rank deficiency or low-support z energy requires scientific review; a positive-subspace condition number must not be interpreted as invertibility of the entire matrix. The saved numerical-null-space energy describes the low-eigenvalue subspace under the frozen rank threshold, not proof of an exact algebraic null space. No singular matrix was silently inverted. All diagnostic matrices and failures remain in the frozen payload. Runtime provenance also records a resource-only correction: later queued workers use four threads within the actual eight-physical-core Slurm allocation; earlier running solvers were left intact, with the original script preserved. No algorithm, float64 precision, tolerance, input or readiness gate changed.

## Signed disease-risk direction

`tables/gwas_only_risk_direction.tsv.gz` preserves **482,169 rows**, including unresolved and deferred-MHC rows. Assigned directions are A **140,886**, B **51,989**, C **175,749**, totaling **368,624**. The table contains source and analysis identities, effect/other alleles, beta/log-OR or OR, SE, z, P semantics, disease-increasing allele, direction status and harmonization status. Native versus normalized allele labels are documented in `tables/variant_ledgers_metadata.json`; no unsupported normalized risk allele is inferred.

These are directions of signed GWAS contrasts, including nonsignificant variants. They are not causal claims, individual disease predictions, PIPs or candidate priorities. No direction comes from REF/ALT labels alone, AF, ancestry, ancestral state, regulatory scores or a functional model. Unresolved direction rows have no assigned disease-increasing allele.

## Future method and remaining blockers

The frozen prospective design remains **multi-signal SuSiE-RSS with uniform variant priors**, followed, if separately authorized and technically resolved, by **same-input FINEMAP sensitivity**. `config/future_finemapping_spec.json` pins R 4.4.3, susieR 0.16.6, FINEMAP 1.4.2, Reference-LAPACK/BLAS tag v3.12.1, float64, L=10, convergence settings, 95% credible-set design and purity threshold 0.5 absolute correlation. These are future settings; no credible sets were constructed. Source archives/hashes are retained. A complete compiled environment and dependency lock are not claimed.

Before execution, the following remain necessary:

- **A:** resolve missing/unresolved strong identities and missing LD variants where applicable; obtain or justify smoking-subgroup-compatible signed covariance or external-reference approximation; close the variant/effective-N and statistical covariance contract; review retained precision failures and diagnostic residual/rank evidence without relaxing gates.
- **B:** legitimate dense Japanese genotype/reference access with exact allele/order/dosage/sample provenance; resolve original pooled-cohort compatibility and score/software lineage; establish a defensible SPA/statistic covariance and effective-N design; then run the currently unavailable LD diagnostics.
- **C:** resolve exact GWAS/LD sample and release linkage, per-variant/effective N, SPA-calibrated statistic versus RSS likelihood/covariance, and the diagnostic rank/residual implications. Same-project availability and complete strong-signal overlap alone do not establish this.
- **Software:** build and validate the pinned R/dependency/BLAS environment and record actual binary/library hashes; resolve FINEMAP’s undocumented seed control and provider-binary authenticity issue before its stochastic sensitivity. No fictional FINEMAP `--seed` flag was introduced.
- **Authorization:** investigator review and new explicit permission to run statistical fine-mapping after eligible inputs are genuinely cleared.

`execution_inputs/manifest.json` therefore contains **zero packages**. Existing ordered summary/LD files are labeled diagnostic materials, not approved execution inputs. There is no implication that every track must succeed before a genuinely cleared future track can proceed; each track/locus requires its own applicable gates and subsequent authorization.

Sakornsakolpat 2019 / GCST007692 remains separate. The public 19,373 absolute-coefficient records were not reconstructed into signed statistics. `SAKORN2019_DATA_REQUEST_DRAFT.md` is an **unsent** author/consortium request for complete signed statistics and ancestry/cohort/LD information. No external message was sent.

## Artifact map and freeze

- Complete-file acquisition/integrity/schema/effect audit: `tables/full_file_acquisition.tsv`, `tables/complete_file_integrity.tsv`, `tables/whole_file_audit_inventory.tsv`, `tables/effect_scale_summary.tsv` and authoritative `tracks/{A,B,C}/` outputs.
- Exact candidate-blind locus boundaries: `tables/prospective_loci.tsv`; complete nonsignificant rows remain in track locus-variant tables.
- Harmonization and risk direction: `tables/variant_harmonization_inventory.tsv`, `tables/gwas_only_risk_direction.tsv.gz` and `tables/variant_ledgers_metadata.json`.
- Signed LD/source/ordering/exclusion evidence: `ld/source_contracts.json`, `ld/{A,C}/{locus_id}/`, `tables/GWAS_LD_compatibility.tsv` and `ld/extraction_validation.json`.
- Track/locus/gate decisions: `tables/track_readiness.tsv`, `tables/locus_readiness.tsv`, `tables/readiness_gates.tsv` and `tables/readiness_integration.json`.
- Software design and no-cleared-input record: `config/future_finemapping_spec.json`, `execution_inputs/manifest.json`.
- Provenance, synthetic tests, preserved implementation/access failures and independent reviews: `provenance/`, track result validation and receipt directories.

The exact report path is `diseases/COPD/07_gap_closure/fine-mapping-input-resolution-1.0/FINE_MAPPING_INPUT_RESOLUTION_REPORT.md`. Final freeze status is recorded in `provenance/freeze.json`; the SHA256 ledger is `provenance/artifact_checksums.tsv` and the independent verification receipt is `provenance/freeze_verification.json`. The ledger binds scientific inputs, outputs, matrices, original source blocks, diagnostics, reports and provenance. Only inventoried interpreter bytecode is excluded; no acquired scientific payload is silently omitted. Large local artifacts are inventoried for future storage/publication review.

`provenance/stage_integrity_validation.json` verifies the earlier preflight payload, unchanged baseline history, authorized register-only tracked changes and append-only original register prefixes. The old preflight register hashes now identify preserved prefixes; they are not falsely asserted to hash the expanded whole registers. The current complete registers are separately bound by the new freeze. HEAD remains at the baseline; this stage stays local.

**Stop for investigator review. No statistical fine-mapping module was started.**
