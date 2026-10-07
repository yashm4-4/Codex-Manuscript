# Target-gene evidence execution contract 1.0

Inputs are the frozen R010 order, R006 mapping, GTEx v10 significant-pair
archive, phenotype R005, signed A/B/C GWAS ledger, and the three preflight
contracts. Their hashes must equal the preflight receipt before execution.
The previously audited 337-site bidirectional liftOver and the GRCh37 FASTA
are read-only local references; V1 2,001-base reference windows validate the
GRCh38 candidate REF. Source files live in ignored `sources/` and can be
reacquired from the URLs in `source_acquisition_receipts.tsv`.

Crosswalks use chromosome, 1-based coordinate, REF and ALT. An rsID alone
cannot establish identity. Signed effects require a verified source contrast,
reference-consistent allele pair, unique point liftOver and round trip. An
indel additionally requires matching reference-anchored alleles and equal
left-normalized representations. Palindromic SNVs are withheld because the
available point liftOver does not independently resolve strand orientation.
Missing candidates outside the prior GWAS slices mean absent from the *local
ledger*, not absent from the complete original GWAS. A missing significant
GTEx row does not establish absence of a cis-eQTL.

Published Moloc is imported at the reported GRCh38 two-megabase window,
gene/feature and QTL class. Table 5 gives one highest-PPA record per window;
its best-colocalized SNP is metadata, not a causal candidate claim. A
candidate overlapping a published window receives locus context only. rE2G
element overlap is predictive. ABC and rE2G in ENCSR528UQX share a lineage;
only rE2G is counted. V1 GENCODE/Catalog proximity and model scores do not
increase the target tier. The frozen hierarchy is applied literally, with
coding overlap separate. No de novo colocalization is executed.
