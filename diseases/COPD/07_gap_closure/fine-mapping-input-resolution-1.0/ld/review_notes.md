**Independent review of root LD diagnostic implementation**

Reviewed `scripts/match_pan_gwas_ld.py`, `scripts/extract_ld_matrices.py`, `scripts/ld_numeric_diagnostics.py` and `config/ld*_contract.json`. The exact reviewed hashes and review time are in `tracks/C/results/root_ld_code_review.json`. The reviewer changed no root diagnostic scripts or matrices.

The code uses exact native chromosome/position/REF/ALT keys, independently checked source identities, unique sorted LD indices and ALT-effect orientation. Index order is passed unchanged to Hail's submatrix extraction and bound to the ordered-table hash. Reconstruction adds the missing opposite triangle and restores the original diagonal before explicit Gram-to-correlation normalization. Source triangles and derived matrices are distinct artifacts. The code does not flip alleles to improve a residual, repair PSD, infer variant effects, or compute posterior probabilities. These choices are appropriate for the stated diagnostic scope.

Two contract-traceability improvements were requested before freezing:

- Bind a same-chromosome and `max(position)-min(position)<=10,000,000` check to each extracted matrix, including source index order and locus edges. The six non-MHC C intervals are each below 3.25 Mb, so no observed C locus violates the source band; an explicit artifact is still needed to verify rather than assume coverage.
- Record rank deficiency and whether positive-subspace conditioning exceeds the frozen warning threshold separately. The running script's `PASS_NUMERICS_ONLY` label checks PSD, symmetry and correlation range. It must not be read as approval of singular-support compatibility, good conditioning, correct source-statistic covariance or final eligibility. Null-space z energy remains a diagnostic requiring interpretation.

Root chose a separate post-extraction band/condition audit to bind these facts without changing the already-running diagnostic source. Its final artifact should be consulted with the numeric reports and source contracts; this note does not assert that the future artifact has already passed.

The bounded scalar mismatch optimizer includes domain checks and endpoint comparisons and records convergence. Its reported objective/parameter is descriptive; a bounded scalar search does not independently certify a global optimum for every possible multimodal objective and supplies no universal scientific clearance threshold. Pairwise residuals likewise do not distinguish an input problem from additional true signals by themselves.

The exact C diagnostic join Boolean is `qc_status == SOURCE_STATISTIC_QC_PASS` plus a VERIFIED harmonization status. Root's independently computed C counts agree with that Boolean for every non-MHC locus. Provider `high_quality=false` remains an annotation rather than an automatic rejection. All source-significant C records remain represented. The 24 normalization-collision rows have P>5e-8 and remain unresolved; none was silently chosen over its alternative representation.

No source-code review or matrix numerical result overrides the separate unresolved release/test/sample contract described in `LD_SOURCE_CONTRACT.md`. Final readiness and software/input freezing belong to the root stage; no fine-mapping execution is approved here.
