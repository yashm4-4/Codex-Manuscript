# BBJ effect/test source lineage

The provider README and original Ishigaki methods both identify SAIGE 0.29.4.2. The exact archived compiled package was downloaded solely for source-lineage inspection, not installed or executed. Its DESCRIPTION confirms Version 0.29.4.2, Date 2018-09-21, and ALT-oriented BETA and Tstat. The original 0.29.4.2-tag URL was 404 and its receipt is preserved.

The immutable source snapshot 75173102920f717a3f82cc2e4a7c35aa1f6d0abc is dated 2018-09-21, the same date as the compiled archive, but its repository DESCRIPTION says 0.29.4 without the patch suffix. The immutable adjacent tag 0.29.4.4 (440b99e1153863868c5818f53b9103ded07a08ec) and the earlier 0.29 source (42bf1c3b5955d158a3b46afd14170a780e47e5db) provide corroborating algorithmic lineage. Exact binary-to-source build equivalence is not asserted.

In these snapshots scoreTest_SPAGMMAT_binaryTrait estimates an approximate log-OR from null-model score/variance, then sets SE = abs(logOR / qnorm(SPA_P / 2)). The fast normal branch uses the unadjusted score P in the analogous conversion. The branch internally converts high-frequency ALT to minor-allele dosage, and explicitly reverses BETA and Tstat on output to restore ALT orientation. This documented internal operation is not authorization to flip source rows to improve diagnostics.

Thus beta/SE versus SPA P agreement is expected and does not independently establish ordinary Gaussian/Wald likelihood suitability. The provider's p.value.NA retains the normal-approximation P; discrepancies between it and the SPA P are preserved in the whole-file diagnostic. Source direction can be interpreted separately from RSS execution readiness. No new SE, P, effect estimate or allele orientation is substituted.

Sources: provider README https://humandbs.dbcls.jp/files/hum0014/hum0014_v17_v18_v21_README.txt ; original study https://pmc.ncbi.nlm.nih.gov/articles/PMC7968075/ ; SAIGE method https://pmc.ncbi.nlm.nih.gov/articles/PMC6119127/ ; archived software/source URLs are in local JSON receipts.
