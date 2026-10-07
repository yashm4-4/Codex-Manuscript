# Independent pre-freeze review

Review status: **PROVISIONAL; STAGE COMPLETION AND FREEZE PENDING**.

This review inspected the existing Track A/C audit products, Track B products, root LD acquisition/join/extraction/diagnostic scripts, numerical contracts, recorded validations, and selected provenance. It did not download GWAS data, rerun a scientific audit, calculate matrix diagnostics, execute statistical fine-mapping, or change scientific outputs. The reviewer did not implement Track A, Track C, or root matrix processing. The reviewer did implement Track B and the future-software specification earlier; review of those components is an internal consistency check, not an independent replication. Checks below describe a snapshot while root matrix diagnostics were still running. Final report, readiness matrix, historical/register validation, and stage checksum freeze must receive a subsequent review.

## Canonical counts and interpretation

| Track | Complete source rows | Merged GWAS loci | Non-MHC loci | Deferred MHC loci | Significant seed rows | Signed direction rows in prospective loci |
|---|---:|---:|---:|---:|---:|---:|
| A, Kim ever-smokers | 9,886,853 | 11 | 10 | 1 | 1,623 | 140,886 including MHC; 92,487 non-MHC |
| B, BBJ combined COPD | 8,678,470 | 5 | 5 | 0 | 462 | 51,989 |
| C, Pan-UKB J44 EUR | 28,987,534 source union rows | 7 | 6 | 1 | 334 | 175,749 including MHC; 125,967 non-MHC |

The source-union row count for C must not be described as the number of valid EUR statistics. C reports 23,861,813 valid signed EUR rows and 18,170,886 autosomal rows passing source confidence. Counts of source-derived signed contrasts are not counts of execution-ready SuSiE/FINEMAP inputs. The 21 non-MHC prospective loci across the three tracks are track-specific intervals, not 21 independent genomic regions.

The current zero execution-clearance conclusions are defensible. A's supplied OR/SE contract and significant statistics pass the frozen precision check, but nine of ten non-MHC loci have unresolved significant identities; the remaining locus also lacks a closed sample/N/reference-LD contract. B has an explicit signed ALT-effect source contract and 51,989 reference-verified locus contrasts; its exact score/test/software-branch lineage, RSS-compatible sampling model, and dense signed Japanese LD are unresolved. C has source-signed ALT effects and dense same-project signed LD, but SPA-calibrated SE arithmetic does not independently establish a Gaussian/Wald likelihood, and sample/release/statistic covariance compatibility remains unresolved.

One label needs care in the final integrated report: B's `source_schema_signed_identity_pass_nonMHC_loci: 0` is a conservative complete-input gate and must not be paraphrased as no verified BBJ effect direction or no verified source identities. The same audit explicitly reports five loci with verified signed rows. Report complete-summary-likelihood clearance separately from signed-contrast validity.

## Preserved corrections and data exclusions

The final B count is five loci and 462 significant source-valid rows. The earlier four-locus/164-row attempt wrongly made nominal-N disagreement a row failure, although the frozen rule required positive N. `tracks/B/results/implementation_correction.json`, the archived first outputs, and `scripts/initial_nominal_N_audit.py` preserve that error. The correction restores chr12 under the original rule; it is not a relaxed threshold. Original lower-N source values remain unchanged.

A preserves the failed first reference-output serialization and validation under `results/attempt1_reference_serialization/`. Canonical source-string restoration retains the original P/OR/SE strings, without changing scientific decisions. `tracks/A/results/artifact_validation.json` records all listed checks passing. Its 25,089 whole-file precision failures are nonsignificant and remain quarantined; 436 lie in prospective intervals. Palindromic/duplicate/unproven identities remain explicit, including unresolved significant variants.

C preserves its invalid SE=0/negative-log-P=infinity source row outside the locus seeds, and its normalized duplicate pairs remain excluded without selecting the more significant representation. Provider `high_quality` includes ancestry-frequency/gnomAD criteria and is correctly kept as an annotation rather than automatically removing strong GWAS rows. The canonical diagnostic join therefore retains strong high-quality-false rows, including the chr20 lead. C's 51 recorded product checks pass. These are existing product-validation results, not an independently repeated full-file audit.

## Frozen contracts

Direct SHA-256 recomputation matched all seven inspected contract receipts:

- A effect contract: `067d27a9f0ff391fce0d1496a64487e278ae452666e1d8a8dc7a41aef1c8bf79`.
- A formatted transport contract: `f54c68a73db552560e2f0dc3bdfa34adbfeb173ee15a9ade2fe198fd7a499c7d`.
- A reference contract: `272bd53cc3512c80643886786649b3625efc104a6044e440f4668dfaaa240a93`.
- B tolerances: `315534b227e9538116734f6894bde18e991ba18cf986aa73718d956cc27aaa2e`.
- C GWAS contract: `81394dd6bf940ca7b380bf8e4da61db9816f5c9bbb851376058e3e0c42898628`.
- Root LD diagnostic contract: `5f87da29e70eafca3dcc6f9707ae8a0ecf0ca7847b4ae02c6d9414be8587a058`.
- Root LD residual contract: `a55d0d9a874128a797af9260730d2775bc60ffba587500ac3d8446c4535fc289`.

A's receipt transparently states it is retrospective, based on unchanged configuration contents and pre-audit mtimes. It must not be described as a contemporaneously signed external attestation. Root numerical contracts predate numerical results; no data-derived clearance threshold is introduced. Floating-point symmetry, range, PSD, rank, and conditioning rules are machine-precision/dimension based. The descriptive mismatch parameter has no invented acceptance cutoff. Statistical/test/sample provenance gates cannot be overridden by a numerical pass.

## Root LD processing review

The archived provider code constructs a covariate-residualized dosage Gram matrix, then stores only an upper triangle within its genomic radius. Root preserves that source matrix. The separately named deterministic conversion `G_ij / sqrt(G_ii G_jj)` is explicitly a diagnostic correlation representation, not proof that the matrix is the covariance required by the released SAIGE statistics. Raw and derived matrix hashes must remain distinct in the final freeze.

C joins by exact native chromosome/position/REF/ALT keys, requires valid signed statistics and reference/identity checks, asserts ALT effect orientation, rejects repeated LD indices, and orders inputs by provider index. A's join additionally verifies the archived reference normalization chain and exact provider-index hash; the source effect-to-LD-ALT multiplier comes from independently established physical allele coding. It transforms z and beta into the native LD-ALT basis, which is consistent with leaving the source matrix in that basis. It does not infer signs from r-squared, P values, residual improvement, or rsID alone. A preserves original scalar strings and writes transformed values only in named derived columns.

The matrix extraction script takes exactly the sorted unique indices from the ordered input and records its hash, matrix hash, shape, and float64 type. At this review snapshot, all six C matrix-extraction records existed, with sizes 17,801; 27,382; 20,055; 18,706; 22,575; and 19,448. Direct recomputation confirmed all six ordered-input hashes still matched their extraction records. Matrix payload hashes themselves were not independently recomputed in this bounded review. A extraction/join completion remained pending.

Both block plans report zero absent source-grid parts: A requests 346 unique blocks and C 203. Final acquisition validation must confirm every planned object has an intact successful receipt and payload before interpreting any sparse matrix zero. The reviewed acquisition code records errors but does not itself reject a nonempty `absent_parts` list. All current intervals are inside the provider radius; final validation should explicitly preserve this condition.

Root diagnostic code preserves the source triangle, checks finite positive diagonals and the unused lower triangle, derives and saves the diagnostic correlation, evaluates symmetry/range and the full eigenspectrum, reports positive-subspace conditioning/null-space z energy, and computes descriptive lead-conditional and null-covariance mismatch quantities. There are no variant-effect posterior fits or statistical fine-mapping calls in these scripts. No ridge, nearest-PSD projection, eigenvalue clipping, arbitrary sign repair, or silent off-band fill is implemented. The synthetic numerical test log reports three tests passing. No real-locus `numeric_diagnostics.json` existed yet at this snapshot, so no real-locus numerical outcome is independently endorsed here.

## Items required before final freeze

1. Complete the intended root A/C acquisition, overlap, extraction, and numerical products, or explicitly record any uncompleted diagnostic as a blocker. Confirm exact lead/strong-variant coverage and preserve failed attempts.
2. Independently verify source-block hashes/receipts, raw and derived matrix hashes, ordered-input hashes, index uniqueness/order, and matrix dimensions. The numerical script does not itself reassert the extraction-input hash; the final validator must do so.
3. Reconcile final track/locus readiness to the canonical corrected counts. Keep signed source validity, complete method-specific summary-likelihood gates, LD numerical checks, scientific LD compatibility, and execution clearance distinct. A numerical pass or C's same-project label alone cannot clear execution.
4. Check the final report, machine-readable readiness, future-software design limitations, and all exact artifact paths against the completed products. Proposed software source versions are not an installed or validated future runtime; FINEMAP seed control and complete binary/dependency lineage remain future blockers.
5. Verify historical stage payload/checksum preservation and the complete append-only register-prefix integrity after root appends the new entries. At the review snapshot, `git status --short` showed only the untracked new stage; this supports, but does not replace, the final historical hash comparison.
6. Generate and independently verify the final stage inventory, payload checksums, freeze record, and final validation. No stage-level freeze or publication claim is supported by this provisional review.

The inspected processing scripts and audit products contain no evidence of candidate/model scoring, use of the 337 candidates to define loci, per-variant posterior/credible-set inference, functional annotations as priors, target-gene assignment, or colocalization. Archived third-party source and future method configuration are documentation, not evidence of execution. Historical integrity and absence-of-prohibited-execution claims remain subject to root's final provenance validation.
