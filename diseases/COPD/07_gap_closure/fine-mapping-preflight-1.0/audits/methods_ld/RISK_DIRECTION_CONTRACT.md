# Prospective signed association / allele / LD contract

Stage: `fine-mapping-preflight-1.0`. Design only. No risk alleles, variants, loci, posterior probabilities or credible sets have been computed here. This contract is candidate-blind and must precede any later ingestion, harmonization or inference authorization.

## What a risk direction means

A disease-risk-increasing allele is identifiable only from a valid signed association for a correctly oriented **COPD case-versus-control outcome**, with a verified effect allele and its corresponding other allele. Positive log-odds beta, or OR > 1, makes the effect allele risk-increasing; negative beta, or 0 < OR < 1, makes the other allele risk-increasing for that exact biallelic contrast. Beta = 0 or OR = 1 supplies no increasing direction. An estimate's sign is not proof of association or biological causality. Uncertainty, significance and phenotype class must remain separate fields.

This is an algebraic interpretation of the reported association, not a rule that ALT is risky. Explicit case/control coding, transformation, effect scale, study accession, ancestry, sex stratum and release are mandatory. A quantitative lung-function, ML/liability, within-case severity/progression, smoking-interaction or bronchitis effect remains **trait-increasing/decreasing in its own phenotype**; it must not silently acquire a direct COPD-susceptibility risk label. Reversed case coding requires an explicitly documented outcome reversal. Unknown coding means unresolved direction.

Never infer risk from REF/ALT status, minor/major allele, allele frequency, ancestral state, LD magnitude/sign alone, eQTL NES, nearest gene, published functional outcomes, or V1/V2 model output. Do not combine signs across overlapping or phenotypically different studies to manufacture consensus.

## Required future row contract

Preserve source strings alongside normalized representations:

| Layer | Required fields and acceptance rule |
|---|---|
| Source identity | URL, release/accession, exact object hash, row identifier and source schema; full genome-wide file rather than a hit list or header sample |
| Association | Effect allele, other allele, signed beta or positive finite OR with declared scale, positive finite SE or valid signed test statistic, P encoding/test definition, per-variant N or documented constant N, case/control totals and applicable QC flags |
| Phenotype | Exact case definition, control exclusions, outcome direction, quantitative transformation, direct/EHR/ML/severity/subphenotype/interaction/non-COPD tier |
| Variant identity | Native assembly, chromosome, 1-based position, exact REF/ALT, explicit contrast at multiallelic sites, normalization/reference provenance, rsID only as an auxiliary identifier |
| LD identity | Release, ancestry/cohort/sample definition, dosage-versus-hard-call input, imputation and covariate processing, matrix row index, indexed alleles and counted allele, LD sample N and exact local source objects |
| Harmonization | Source-to-LD allele mapping, swap/complement action, signed orientation multiplier, retained/excluded reason, duplicate/conflict status; retain the original values |

P-only or unsigned-Z files cannot establish direction. Do not reconstruct SE from an arbitrary P column: a saddlepoint-adjusted, mixed-model score, meta-analysis or rounded P need not equal the Wald P of beta/SE. The analysis-statistic and sample-size contract must be resolved first. If beta/SE is a valid Wald statistic, preserve that signed z; if the source supplies a score statistic, use a method configuration explicitly compatible with it. Current [SuSiE-RSS documentation](https://stephenslab.github.io/susieR/reference/susie_rss.html) distinguishes Wald and score inputs; this is a version-sensitive contract, not permission to recode a source test to improve fit.

## Signed LD transformation

For future standardized genotypes in the source LD counted-allele orientation, let `d_j` be +1 for the chosen analysis allele and -1 for a verified biallelic swap. With `D = diag(d)`, the consistent transformation is `z_analysis = D z_source` and `R_analysis = D R_source D`. Swapping only beta/z or only an LD row is invalid. Reorder both matrix axes by the **same exact variant key** as the statistics. Complementing nucleotide labels without swapping which physical allele is counted does not itself negate z. Strand, build and allele identity must be reconciled before deciding any sign change.

Keep `R` as signed correlation, not r², |r|, D-prime, LD scores or a proxy list. The sign cannot be recovered from r² by square rooting. [FINEMAP's official input specification](http://www.christianbenner.com/) requires same-ordered Pearson correlations and explicitly discusses effect-allele consistency; [SuSiE's diagnostic vignette](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html) illustrates false signals from allele mismatches. A diagnostic suggestion to flip is not independent allele evidence: never flip an allele solely to improve z/LD agreement.

Exact duplicate records with identical declared provenance can be deduplicated under a frozen rule; conflicting duplicates block that contrast. Keep multiallelic contrasts separate and do not join on position or rsID alone. For indels, require the reference assembly and a reproducible left-normalized representation with verified reference bases before any cross-file join. Unresolved complements and palindromic A/T or C/G alleles are excluded with reasons; trustworthy strand/reference and matching allele-frequency evidence may resolve them, but frequency never determines risk direction. No universal frequency-difference rescue threshold is established by this preflight.

## Resource-specific cautions

PanUKB flat-file documentation explicitly identifies ALT as the effect allele in GRCh37. Archived LD source calls `get_filtered_mt(..., entry_fields=['dosage'])`, ultimately `hail.import_bgen`; [Hail's official import contract](https://hail.is/docs/0.2/methods/impex.html#hail.methods.import_bgen) defines dosage as expected ALT count and the first indexed allele as REF. This documents the intended sign chain, not a successful local harmonization. Verify the released index and every retained row in the future. [PanUKB allele/schema documentation](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html)

The public PanUKB matrix is triangular and distance-banded. Its archived calculation standardizes dosages **before** covariate residualization and then writes `Z Zᵀ/n`, without an explicit second variance normalization. Therefore local diagonal/scaling must be inspected and reconciled with the documented correlation contract; one cannot assume a unit diagonal from the word “LD.” A justified conversion from a residualized Gram matrix to correlation would require positive diagonals and matched statistic scaling, recorded as a prospective derived artifact, never a silent repair. No matrix entries were accessed or transformed in this stage. [PanUKB computation code](https://github.com/atgu/ukbb_pan_ancestry/blob/master/compute_ld_matrix.py)

## Prospective tests, not executed scientific calculations

Before later inference, synthetic unit tests must demonstrate: one verified allele swap negates the corresponding z and both LD off-diagonal row/column signs; a pure strand complement preserves physical orientation; a permutation moves z and both R axes together; identical REF/ALT labels do not imply the same effect allele; ambiguous, conflicting or missing inputs are rejected. Tests must include reversed case coding, OR = 1, negative/nonfinite OR, nonpositive SE, mixed builds, multiallelic collisions, indels and exact duplicate provenance. This document defines those tests; it does not execute harmonization or assign direction to any COPD variant.
