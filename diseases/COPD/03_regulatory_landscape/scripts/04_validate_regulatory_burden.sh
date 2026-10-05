#!/usr/bin/env bash
set -euo pipefail

# Independent bedtools checks for COPD-S3-R002. The production analysis is
# implemented in Python, so this script deliberately uses a different engine.

section_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bedtools_bin=${BEDTOOLS_BIN:-/usr/local/apps/bedtools/2.31.1/bin/bedtools}
if [[ ! -x "$bedtools_bin" ]]; then
    bedtools_bin=$(command -v bedtools)
fi
"$bedtools_bin" --version

validation_tmp=$(mktemp -d)
trap 'rm -rf "$validation_tmp"' EXIT

peak_dir="$section_dir/data/encode_peaks"
locus_table="$section_dir/results/COPD-S3-R002_gene_loci_grch38.tsv"
element_dir="$section_dir/data/copd_donor_regulatory_elements"

awk -F '\t' 'BEGIN {OFS="\t"} NR > 1 {print $6,$9,$10,$1}' \
    "$locus_table" > "$validation_tmp/selected.bed"
awk -F '\t' 'BEGIN {OFS="\t"} NR > 1 && $5 == "True" {print $6,$9,$10,$1}' \
    "$locus_table" > "$validation_tmp/replicated.bed"

check_value() {
    local label=$1
    local observed=$2
    local expected=$3
    if [[ "$observed" != "$expected" ]]; then
        printf 'FAIL\t%s\tobserved=%s\texpected=%s\n' "$label" "$observed" "$expected" >&2
        exit 1
    fi
    printf 'PASS\t%s\tobserved=%s\texpected=%s\n' "$label" "$observed" "$expected"
}

# ATAC peaks with at least one same-lobe histone-peak overlap.
while read -r lobe element atac histone expected; do
    observed=$(
        "$bedtools_bin" intersect -u \
            -a "$peak_dir/$atac.bed.gz" \
            -b "$peak_dir/$histone.bed.gz" |
        awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {n++} END {print n+0}'
    )
    check_value "same_lobe_refined:${lobe}:${element}" "$observed" "$expected"
done <<'EOF'
upper_right enhancer ENCFF906HOT ENCFF149QUM 91110
upper_right silencer ENCFF906HOT ENCFF685CNW 13981
lower_right enhancer ENCFF189LZA ENCFF299FWI 99601
lower_right silencer ENCFF189LZA ENCFF591PDC 19142
lower_left enhancer ENCFF899TQV ENCFF014OZD 113178
lower_left silencer ENCFF899TQV ENCFF725NIC 22432
EOF

# Distinct donor-union regions overlapping either locus union.
while read -r locus_set element definition expected; do
    if [[ "$locus_set" == "replicated" ]]; then
        loci="$validation_tmp/replicated.bed"
    else
        loci="$validation_tmp/selected.bed"
    fi
    regions="$element_dir/donor_union.$element.$definition.grch38.bed.gz"
    observed=$("$bedtools_bin" intersect -u -a "$regions" -b "$loci" | wc -l | tr -d ' ')
    check_value "donor_union_locus_overlap:${locus_set}:${definition}:${element}" "$observed" "$expected"
done <<'EOF'
replicated enhancer preliminary_histone_peak 1892
replicated silencer preliminary_histone_peak 1201
replicated enhancer refined_accessibility_histone_overlap 1815
replicated silencer refined_accessibility_histone_overlap 282
selected enhancer preliminary_histone_peak 2027
selected silencer preliminary_histone_peak 1241
selected enhancer refined_accessibility_histone_overlap 1953
selected silencer refined_accessibility_histone_overlap 298
EOF

printf 'PASS\tall_checks\t14_of_14\n'
