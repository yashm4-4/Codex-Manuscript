# Section 3 traceability handoff

The following rows are ready to append to the project-level registers. They
were kept out of the shared files while Section 2 and Section 4 agents were
writing concurrently.

## `sources.tsv`

```tsv
COPD-SRC-038	3	ENCODE Project. SCREEN Registry of candidate cis-regulatory elements, Registry V3, GRCh38.	https://downloads.wenglab.org/V3/GRCh38-cCREs.bed	2026-09-21	database release	verified	Local immutable source under data/screen_ccres; 1,063,878 canonical intervals
COPD-SRC-039	3	Ensembl. Regulatory Build for Homo sapiens GRCh38, release 116.	https://ftp.ensembl.org/pub/release-116/regulation/homo_sapiens/GRCh38/annotation/	2026-09-21	database release	verified	Local immutable regulatory-features GFF3 under data/ensembl_regulation; 380,818 canonical intervals
COPD-SRC-040	3	FANTOM Consortium. FANTOM5 CAGE enhancer atlas, hg38_latest.	https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/	2026-09-21	database release	verified	Local immutable BED under data/fantom5; 63,285 canonical enhancer intervals
COPD-SRC-041	3	UCSC Genome Browser. RepeatMasker annotations for the GRCh38/hg38 assembly.	https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz	2026-09-21	database table	verified	Local immutable table under data/repeats; 5,317,286 canonical intervals retained
COPD-SRC-042	3	ENCODE Project and Boyle Lab. ENCODE hg38 blacklist version 2.	https://github.com/Boyle-Lab/Blacklist	2026-09-21	reference dataset	verified	Local immutable BED under data/blacklist; 636 canonical intervals
COPD-SRC-043	3	GENCODE Project. GENCODE human release 50 comprehensive gene annotation for GRCh38.	https://www.gencodegenes.org/human/release_50.html	2026-09-21	database release	verified	Local immutable GTF under data/gencode; 387,360 unique canonical gene-aware CDS intervals after normalization
```

## `results_register.tsv`

```tsv
COPD-S3-R001	3	Normalized GRCh38 regulatory reference inventory	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R001_regulatory_reference_inventory.tsv	complete	COPD-SRC-038 through COPD-SRC-042	2026-10-01	Five shared sources normalized to six-column canonical-chromosome BED with source and output SHA-256 values
COPD-S3-R002	3	Severe-emphysema lung enhancer and silencer burden in selected COPD GWAS-gene loci	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R002_lobe_union_summary.tsv	complete	COPD-SRC-018;COPD-SRC-035;COPD-SRC-043	2026-10-01	Three donor-matched lung lobes; 152 selected and 140 replicated unambiguous GENCODE v50 loci; full limitations in computational report
COPD-S3-R003	3	Frozen tag and ancestry-matched LD-proxy regulatory classification	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz	complete	COPD-SRC-018;COPD-SRC-035;COPD-SRC-037 through COPD-SRC-043	2026-10-01	15,389 candidate records retained; 15,386 BED-eligible; GENCODE CDS, four donor element sets, SCREEN, Ensembl, FANTOM5, RepeatMasker, and blacklist annotated with nonexclusive and explicit exclusive summaries
```

## `activity_log.tsv`

```tsv
2026-10-01	3	COPD-S3-A001	Download and validate selected severe-emphysema lung peak files	ENCODE REST API; Python requests	Donor ENCDO520EJG; upper-right, lower-right, lower-left lung; released pseudoreplicated GRCh38 ATAC, H3K27ac, H3K27me3 narrowPeak	diseases/COPD/03_regulatory_landscape/data/encode_peaks/	complete	Nine files; 26,659,866 compressed bytes; SHA-256 and donor metadata recorded
2026-10-01	3	COPD-S3-A002	Normalize shared regulatory reference resources	Python 3.9.15; pandas 2.0.3	GRCh38; 0-based half-open BED; chr1-22,X,Y; source files immutable	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R001_regulatory_reference_inventory.tsv	complete	Five resources; 6,825,903 retained canonical intervals; zero invalid intervals
2026-10-01	3	COPD-S3-A003	Build lobe-specific and donor-union regulatory elements and quantify GWAS-gene locus burden	Python 3.9.15; pandas 2.0.3	Same-lobe ATAC-histone overlap at least 1 bp; GENCODE v50 gene body plus or minus 100 kb; donor union merges overlapping or book-ended intervals	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R002_lobe_union_summary.tsv	complete	16 element sets; 32 summary rows; 2,432 per-gene burden rows
2026-10-01	3	COPD-S3-A004	Independently validate regulatory-element and locus-overlap counts	bedtools 2.31.1	Six same-lobe refined counts and eight donor-union locus-overlap counts	diseases/COPD/03_regulatory_landscape/logs/04_validate_regulatory_burden.log	complete	14 of 14 exact count checks passed
2026-10-01	3	COPD-S3-A005	Classify the frozen GWS-tag and ancestry-matched LD-proxy candidate set	Python 3.9.15; pandas 2.0.3; bedtools 2.31.1	GRCh38; at least 1 bp overlap; all 15,389 source records retained; three unlocalized tags explicit; role and ancestry-panel provenance preserved	diseases/COPD/03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz	complete	15,386 BED-eligible records intersected with 10 feature sets; 660-tag-identifier view and nonexclusive and exclusive summaries generated
2026-10-01	3	COPD-S3-A006	Independently validate frozen candidate classification	bedtools 2.31.1; Python standard library	Direct bedtools recomputation for 10 feature sets plus input preservation, explicit missingness, tag expansion, partition, source-ID, and manifest-hash audits	diseases/COPD/03_regulatory_landscape/logs/06_validate_candidate_classification.log	complete	All 10 independent overlap counts matched; all 81 stratum-by-hierarchy partitions summed; complete validation passed
```

## `decisions.tsv`

```tsv
COPD-DEC-008	2026-10-01	3	Operational regulatory-element definitions	Preliminary enhancer is H3K27ac; preliminary silencer is H3K27me3; refined elements are same-lobe ATAC peaks overlapping the corresponding mark by at least 1 bp	Implements both framework definitions while preserving donor and anatomical matching	confirmed
COPD-DEC-009	2026-10-01	3	GWAS-gene regulatory-locus coordinate rule	Use one exact, unambiguous GENCODE v50 gene record plus 100,000 bp on each side; exclude Y_RNA, GUSBP5, and CYP2B7P because exact symbols map to multiple records	Prevents arbitrary coordinate assignment and makes the operational flank explicit	confirmed
COPD-DEC-010	2026-10-01	3	Three-lobe donor union	Merge overlapping or book-ended same-definition intervals across the three lobes	Counts distinct covered regulatory regions rather than summing repeated lobe-specific peak calls	confirmed
COPD-DEC-011	2026-10-01	3	Descriptive peak filtering	Do not remove promoter-overlapping or ENCODE-blacklisted peaks from the Section 3 descriptive burden counts	Preserves released ENCODE calls; downstream model controls and candidate QC apply their own exclusions	confirmed
COPD-DEC-012	2026-10-01	3	Candidate overlap and ineligible-coordinate handling	Require at least 1 bp of GRCh38 BED overlap; retain all source candidate records; assign NA_not_bed_eligible rather than false to feature flags for the three unlocalized tags	Prevents loss of unresolved Catalog tags and distinguishes untested coordinates from tested non-overlap	confirmed
COPD-DEC-013	2026-10-01	3	Exclusive candidate reporting hierarchy	Report preliminary, refined, and comprehensive precedence schemes with bed-ineligible records separated; in the comprehensive scheme use coding CDS, refined enhancer, refined silencer, preliminary enhancer, preliminary silencer, other known regulatory, RepeatMasker, then other	Produces auditable complete partitions while retaining all mutually nonexclusive source flags; precedence is a reporting convention rather than biological rank	confirmed
COPD-DEC-014	2026-10-01	3	Candidate reporting units and ancestry provenance	Use unique reference records as the primary unit, add a 660-tag-identifier view, report GWS-tag and LD-proxy roles separately, and retain AFR/AMR/EAS/EUR/SAS proxy-panel strata as overlapping memberships	Avoids conflating eight shared-reference tag pairs or treating overlapping tag/proxy and panel strata as additive independent groups	confirmed
```
