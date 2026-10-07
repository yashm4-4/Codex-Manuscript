# Track B — original BBJ combined COPD input resolution

**NOT CLEARED.** The complete signed summary file is acquired and audited, but no defensible dense signed Japanese LD source was obtained. No posterior inference, statistical fine-mapping, candidate scoring or model scoring was run.

## Complete file

Original release: [hum0014.v17.COPD.v1.zip](https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.COPD.v1.zip); combined-sex autosomal member `COPD.auto.rsq07.mac10.txt.gz`. Container 2,244,311,450 bytes, SHA-256 `a2473c13b81713a42427cbebbcefebb744ac3acb377dcd92412830c83fee3054`. The first full response ended early; exact resume ranges, ETags and the initial failure are preserved. The completed ZIP size and every member's CRC32 pass. The combined gzip is 901,221,471 bytes, SHA-256 `3c19eb1c4eb618ce18f62cf979caefb70abacae6a2418047e5feebf3fd6bca5e`; both streaming passes reached the gzip end and validated compression. Other archived chromosome-X and sex-specific files were only covered by container CRC, not scientifically analyzed.

The provider N field is variant-specific (204,905–204,907), with 11,542 rows below nominal N. The frozen source-row gate requires positive N. A coding error that initially treated nominal-N deviation as failure was corrected to that unchanged rule; the original attempt and rationale are retained in `implementation_correction.json` and its prior-attempt folder. This correction restored the chr12 locus; only the final top-level tables are authoritative.

The complete file contains **8,678,470 autosomal rows** across all 22 autosomes, 20 space-delimited fields, provider hg19 coordinates. Original row order is chromosome-lexicographic, preserved by source-row identifiers; locus ordering is deterministic numeric chromosome order. Per-variant N, AF_Allele2, AF.Cases, AF.Controls, MAC, Rsq, score, SPA convergence and normal-versus-SPA P are audited. See `whole_file_schema_missingness.tsv`, `whole_file_qc.tsv` and `whole_file_audit.json` for exact missingness/ranges and all failures.

There are 3 duplicate excess rows and 3,313 multiple-identity source positions. The duplicated chr17:1144632 C/CT source identity has four differing rows, outside all selected loci; it is not silently collapsed. No rsID column is present: identity is build/chromosome/position/REF/ALT, never an rsID-only join.

## Signed effect and test contract

The provider defines Allele1=REF, Allele2=ALT/effect allele, BETA as effect of ALT, SE as its uncertainty, and `p.value` as the SAIGE saddlepoint P; `p.value.NA` is the normal-approximation P. The original physician-diagnosed COPD analysis includes 3,315 cases and 201,592 controls, modeled jointly using SAIGE 0.29.4.2 with age, sex and five PCs. Controls include BBJ unrelated diagnoses and ToMMo/IMM/JPHC/J-MICC population cohorts. The original study used SHAPEIT2/minimac3 with 1000 Genomes Phase3 imputation, Rsq>=0.7 and MAC>=10; that imputation source does not authorize substituting a 1KG LD matrix.

`z=BETA/SE` and both P columns are compared for every source row, with tolerances frozen before evaluation. 0 rows exceed the frozen P-versus-z tolerance for SPA P and 399,244 for normal-approximation P. The deterministic evidence table additionally preserves every 1,000th source row plus all qualifying significant rows, spanning every chromosome and P range.

SAIGE source lineage shows that the SPA branch derives SE from its approximate log-OR and SPA P. Exact compiled v0.29.4.2 metadata is retained; adjacent and contemporaneous original source snapshots are retained with bounded version claims (`sources/saige_SOURCE_LINEAGE_NOTE.md`). Consequently, beta/SE versus SPA P agreement is expected by construction; it is not independent proof of an ordinary Gaussian likelihood. No SE, P, effect or allele is repaired to force concordance. The naive beta=Tstat/varT check is descriptive only because the SPA branch internally scales by sqrt(MAC); a source-derived branch expression is separately tested in `score_scale_sample_audit.json` and fails338/9141 sampled rows, including chr12 per-variant-N rows. This score-column lineage remains unresolved; no branch is selected merely to improve agreement. The ALT signed contrast supports estimated disease direction independently of fine-mapping readiness.

## Prospective loci and identity

The frozen candidate-blind P<5e-8, +/-1.5-Mb, transitive-merge rules produce **5 prospective loci**, from 462 source-valid significant rows. No resulting merged interval intersects the conservative chr6:25–36Mb mask; no interval reaches a chromosome edge. Exact intervals and deterministic lead/source-row tie resolution are in `loci.tsv`. All nonsignificant source variants inside the intervals are retained.

Across these loci, 52,013 rows were checked against the independently acquired GRCh37 reference, with 52,013 reference matches and 0 mismatches. There are 0 normalized identity changes, 7,280 palindromic SNPs, and 0 conflicting normalized identities. No allele swap or strand flip was applied. Indels/MNVs receive explicit normalized identities while originals remain retained. Full row ledgers, excluded rows and unassessed LD coverage remain explicit.

`risk_direction.tsv.gz` assigns estimated disease-increasing alleles for 51,989 verified signed source contrasts. This includes nonsignificant estimates and is not a claim of causal direction or statistically proven risk. No regulatory predictions are joined.

## LD and readiness

Official BBJ/NBDC resources expose controlled-access genotype routes (including original BBJ arrays and Japanese WGS imputation references). The GRCh37 BBJ7K reference has 7,472 BBJ WGS participants plus 2,504 1KG participants and is controlled access; it is not the original full COPD analysis sample. BBJ1K/2K references likewise require access and independent external-reference justification. The public ToMMo LD map is recombination distance from 96 individuals, not signed Pearson correlation; jMorp marginal frequencies cannot supply LD signs. Later PheWeb resources do not establish the original release's dense signed matrix. Exact source contracts and receipts are in `ld_source_contract.json` and `sources/`.

All **5 prospective loci contain verified signed source rows**; 0 have every original row fully resolved (symbolic/nonsequence rows remain explicit), **0 passing all summary-statistic likelihood gates**, **0 passing LD gates**, and **0 execution-cleared**. No matrix numerics or summary-LD residuals can be truthfully reported without an acquired matrix. The original mixed-model/SPA statistic covariance, sample-size design, signed matrix provenance/coverage/numerics and locus boundary compatibility remain open. All five loci retain NOT CLEARED. The independent Japanese recruitment offers a distinct ancestry track from UKB, but exact participant overlap is not inferred from labels; male/female analyses were not expanded or treated as replications.

No future SuSiE-RSS execution is authorized by this evidence. No completed execution input package is created. Investigator review is required before any statistical fine-mapping.
