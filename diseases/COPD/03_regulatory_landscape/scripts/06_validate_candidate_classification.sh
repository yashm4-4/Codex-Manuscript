#!/usr/bin/env bash
set -euo pipefail

# Independent checks for COPD-S3-R003. Production annotation is implemented in
# Python; the primary overlap checks below deliberately use bedtools directly.

section_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
project_dir=$(cd "$section_dir/../../.." && pwd)
candidate_bed="$project_dir/diseases/COPD/02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.bed"
candidate_tsv="$project_dir/diseases/COPD/02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.tsv.gz"
classification="$section_dir/results/COPD-S3-R003_candidate_classification.tsv.gz"
tag_classification="$section_dir/results/COPD-S3-R003_gws_tag_identifier_classification.tsv.gz"
nonexclusive="$section_dir/results/COPD-S3-R003_nonexclusive_summary.tsv"
exclusive="$section_dir/results/COPD-S3-R003_exclusive_summary.tsv"
feature_definitions="$section_dir/results/COPD-S3-R003_feature_definitions.tsv"
manifest="$section_dir/results/COPD-S3-R003_analysis_manifest.json"

bedtools_bin=${BEDTOOLS_BIN:-/usr/local/apps/bedtools/2.31.1/bin/bedtools}
if [[ ! -x "$bedtools_bin" ]]; then
    bedtools_bin=$(command -v bedtools)
fi
"$bedtools_bin" --version

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

check_value "candidate_bed_rows" "$(wc -l < "$candidate_bed" | tr -d ' ')" "15386"

donor_dir="$section_dir/data/copd_donor_regulatory_elements"
reference_dir="$section_dir/data/regulatory_references"
cds_bed="$section_dir/data/candidate_classification/gencode_v50_cds.grch38.canonical.bed.gz"

# Each value is the number of distinct candidate records with at least 1 bp of
# overlap. These calls do not consume the production intersection tables.
while IFS=$'\t' read -r feature feature_bed expected; do
    observed=$(
        "$bedtools_bin" intersect -u -a "$candidate_bed" -b "$feature_bed" |
            wc -l |
            tr -d ' '
    )
    check_value "independent_overlap:${feature}" "$observed" "$expected"
done <<EOF
coding_CDS	$cds_bed	257
donor_preliminary_enhancer	$donor_dir/donor_union.enhancer.preliminary_histone_peak.grch38.bed.gz	1377
donor_preliminary_silencer	$donor_dir/donor_union.silencer.preliminary_histone_peak.grch38.bed.gz	612
donor_refined_enhancer	$donor_dir/donor_union.enhancer.refined_accessibility_histone_overlap.grch38.bed.gz	588
donor_refined_silencer	$donor_dir/donor_union.silencer.refined_accessibility_histone_overlap.grch38.bed.gz	86
screen_ccre	$reference_dir/screen_ccre_v3.grch38.canonical.bed.gz	2247
ensembl_regulatory	$reference_dir/ensembl_regulatory_v116.grch38.canonical.bed.gz	1276
fantom5_enhancer	$reference_dir/fantom5_cage_enhancers.grch38.canonical.bed.gz	162
repeatmasker	$reference_dir/ucsc_repeatmasker.grch38.canonical.bed.gz	8131
encode_blacklist	$reference_dir/encode_blacklist_v2.grch38.canonical.bed.gz	21
EOF

python - "$candidate_tsv" "$candidate_bed" "$classification" "$tag_classification" \
    "$nonexclusive" "$exclusive" "$feature_definitions" "$manifest" <<'PY'
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import sys

(
    candidate_path,
    bed_path,
    classification_path,
    tag_path,
    nonexclusive_path,
    exclusive_path,
    feature_definitions_path,
    manifest_path,
) = map(Path, sys.argv[1:])


def fail(label, observed, expected):
    raise SystemExit(f"FAIL\t{label}\tobserved={observed}\texpected={expected}")


def check(label, observed, expected):
    if observed != expected:
        fail(label, observed, expected)
    print(f"PASS\t{label}\tobserved={observed}\texpected={expected}")


def read_tsv(path, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


inputs = read_tsv(candidate_path, compressed=True)
outputs = read_tsv(classification_path, compressed=True)
tags = read_tsv(tag_path, compressed=True)
check("candidate_input_rows", len(inputs), 15389)
check("classification_rows", len(outputs), 15389)
check("tag_identifier_rows", len(tags), 660)

input_fields = list(inputs[0])
for row_number, (source, result) in enumerate(zip(inputs, outputs), start=2):
    for field in input_fields:
        if source[field] != result[field]:
            fail(
                f"input_preservation:line_{row_number}:{field}",
                result[field],
                source[field],
            )
print("PASS\tinput_columns_preserved\t15389_of_15389")

input_ids = [row["candidate_record_id"] for row in inputs]
output_ids = [row["candidate_record_id"] for row in outputs]
check("unique_candidate_input_ids", len(set(input_ids)), 15389)
check("unique_candidate_output_ids", len(set(output_ids)), 15389)
if input_ids != output_ids:
    fail("candidate_order_and_ids", "changed", "identical")
print("PASS\tcandidate_order_and_ids\tidentical")

with bed_path.open() as handle:
    bed_ids = [line.rstrip("\n").split("\t")[3] for line in handle]
eligible_ids = {row["candidate_record_id"] for row in inputs if row["bed_eligible"] == "True"}
check("unique_candidate_bed_ids", len(set(bed_ids)), 15386)
if set(bed_ids) != eligible_ids:
    fail("bed_id_set", "mismatch", "eligible_candidate_ids")
print("PASS\tbed_id_set\t15386_eligible_ids")

expected_ineligible = {
    "unmatched_tag:15q25.1",
    "unmatched_tag:rs139284640",
    "unmatched_tag:rs751872749",
}
ineligible = [row for row in outputs if row["bed_eligible"] == "False"]
if {row["candidate_record_id"] for row in ineligible} != expected_ineligible:
    fail("bed_ineligible_ids", "unexpected", sorted(expected_ineligible))
print("PASS\tbed_ineligible_ids\t3_expected_tags")

direct_features = {
    "coding_CDS": 257,
    "donor_preliminary_enhancer": 1377,
    "donor_preliminary_silencer": 612,
    "donor_refined_enhancer": 588,
    "donor_refined_silencer": 86,
    "screen_ccre": 2247,
    "ensembl_regulatory": 1276,
    "fantom5_enhancer": 162,
    "repeatmasker": 8131,
    "encode_blacklist": 21,
}
for feature, expected in direct_features.items():
    check(
        f"classification_overlap:{feature}",
        sum(row[feature] == "True" for row in outputs),
        expected,
    )
    if any(row[feature] != "NA_not_bed_eligible" for row in ineligible):
        fail(f"ineligible_state:{feature}", "not_explicit", "NA_not_bed_eligible")
print("PASS\tineligible_feature_states\t10_of_10_explicit")

derived_features = {
    "donor_preliminary_any": 1956,
    "donor_refined_any": 639,
    "known_regulatory_any": 2638,
    "any_tested_biological_annotation": 10649,
    "unannotated_across_tested_biological_features": 4737,
}
for feature, expected in derived_features.items():
    check(
        f"classification_derived:{feature}",
        sum(row[feature] == "True" for row in outputs),
        expected,
    )
    if any(row[feature] != "NA_not_bed_eligible" for row in ineligible):
        fail(f"ineligible_state:{feature}", "not_explicit", "NA_not_bed_eligible")

expected_exclusive = {
    "exclusive_class_framework_preliminary": {
        "bed_ineligible": 3,
        "coding_CDS": 257,
        "preliminary_enhancer": 1329,
        "preliminary_silencer": 556,
        "repetitive_only_in_hierarchy": 7504,
        "other": 5740,
    },
    "exclusive_class_framework_refined": {
        "bed_ineligible": 3,
        "coding_CDS": 257,
        "refined_enhancer": 561,
        "refined_silencer": 46,
        "repetitive_only_in_hierarchy": 8000,
        "other": 6522,
    },
    "exclusive_class_comprehensive": {
        "bed_ineligible": 3,
        "coding_CDS": 257,
        "refined_enhancer": 561,
        "refined_silencer": 46,
        "preliminary_enhancer": 835,
        "preliminary_silencer": 518,
        "other_known_regulatory": 1515,
        "repetitive_only_in_hierarchy": 6917,
        "other": 4737,
    },
}
for column, expected_counts in expected_exclusive.items():
    observed = defaultdict(int)
    for row in outputs:
        observed[row[column]] += 1
    if dict(observed) != expected_counts:
        fail(f"overall_hierarchy:{column}", dict(observed), expected_counts)
    print(f"PASS\toverall_hierarchy:{column}\t{sum(observed.values())}_records")

expected_tag_ids = {
    tag
    for row in inputs
    if row["is_gws_tag"] == "True"
    for tag in row["gws_tag_ids"].split(";")
    if tag
}
observed_tag_ids = [row["gws_tag_id"] for row in tags]
check("unique_tag_identifier_ids", len(set(observed_tag_ids)), 660)
if set(observed_tag_ids) != expected_tag_ids:
    fail("tag_identifier_set", "mismatch", "expanded_input_tag_ids")
check("tag_identifier_bed_eligible", sum(row["bed_eligible"] == "True" for row in tags), 657)

nonexclusive = read_tsv(nonexclusive_path)
overall = {
    row["feature"]: int(row["n_overlapping"])
    for row in nonexclusive
    if row["stratum_type"] == "overall" and row["stratum"] == "all_candidate_records"
}
expected_overall = {**direct_features, **derived_features}
for feature, expected in expected_overall.items():
    check(f"nonexclusive_summary:{feature}", overall.get(feature), expected)

exclusive = read_tsv(exclusive_path)
groups = defaultdict(list)
for row in exclusive:
    groups[(row["stratum_type"], row["stratum"], row["hierarchy"])].append(row)
check("exclusive_summary_groups", len(groups), 81)
for key, rows in groups.items():
    totals = {int(row["n_total_records"]) for row in rows}
    if len(totals) != 1:
        fail(f"exclusive_denominator:{key}", sorted(totals), "one_value")
    observed = sum(int(row["n_records"]) for row in rows)
    expected = totals.pop()
    if observed != expected:
        fail(f"exclusive_sum:{key}", observed, expected)
print("PASS\texclusive_summary_partition_sums\t81_of_81")

feature_definitions = read_tsv(feature_definitions_path)
expected_source_ids = {
    "coding_CDS": "COPD-SRC-043",
    "donor_preliminary_enhancer": "COPD-SRC-018",
    "donor_preliminary_silencer": "COPD-SRC-018",
    "donor_refined_enhancer": "COPD-SRC-018",
    "donor_refined_silencer": "COPD-SRC-018",
    "screen_ccre": "COPD-SRC-038",
    "ensembl_regulatory": "COPD-SRC-039",
    "fantom5_enhancer": "COPD-SRC-040",
    "repeatmasker": "COPD-SRC-041",
    "encode_blacklist": "COPD-SRC-042",
}
observed_source_ids = {row["feature"]: row["source_id"] for row in feature_definitions}
if observed_source_ids != expected_source_ids:
    fail("feature_source_ids", observed_source_ids, expected_source_ids)
print("PASS\tfeature_source_ids\t10_of_10")

manifest = json.loads(manifest_path.read_text())
check("manifest_candidate_rows", manifest["candidate_input_rows"], 15389)
check("manifest_candidate_bed_rows", manifest["candidate_bed_rows"], 15386)
manifest_outputs = {Path(item["path"]).name: item for item in manifest["outputs"]}
for path in (
    classification_path,
    tag_path,
    nonexclusive_path,
    exclusive_path,
    feature_definitions_path,
):
    entry = manifest_outputs.get(path.name)
    if entry is None:
        fail(f"manifest_output:{path.name}", "missing", "present")
    if entry["sha256"] != sha256(path):
        fail(f"manifest_hash:{path.name}", entry["sha256"], sha256(path))
print("PASS\tmanifest_output_hashes\t5_of_5_checked")
print("PASS\tall_candidate_classification_checks\tcomplete")
PY
