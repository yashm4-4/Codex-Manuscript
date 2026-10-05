#!/usr/bin/env python3
"""Annotate prioritized COPD variants with frequency and ancestral state.

Allele frequencies come from the local 1000 Genomes Phase 3 GRCh38 remap.
Ancestral states are queried only for predicted causal variants through the
official Ensembl REST VEP and Variation endpoints and cached verbatim enough to
make every assignment auditable.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam
import requests


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases/COPD/04_modeling"
RESULTS = SECTION / "results"
VCF_DIR = ROOT / "data/1000genomes"
INPUT = RESULTS / "COPD-S4-R004_prioritized_candidates.tsv.gz"

OUT_ALL = RESULTS / "COPD-S4-R005_population_annotated_candidates.tsv.gz"
OUT_CAUSAL = RESULTS / "COPD-S4-R005_predicted_causal_population_genetics.tsv"
OUT_SUMMARY = RESULTS / "COPD-S4-R005_population_summary.tsv"
OUT_ANCESTRAL = RESULTS / "COPD-S4-R005_ancestral_allele_audit.tsv"
OUT_RAW = SECTION / "data/COPD-S4-R005_ensembl_ancestral_cache.json.gz"
OUT_MANIFEST = RESULTS / "COPD-S4-R005_analysis_manifest.json"

SUPERPOPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")
SAMPLE_LIST_DIR = ROOT / "diseases/COPD/02_gwas/data/ld_work"
ENSEMBL_REST = "https://rest.ensembl.org"
USER_AGENT = "COPD-regulatory-genomics-workflow/1.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def chromosome(value: object) -> str:
    if pd.isna(value):
        return ""
    result = str(value).removeprefix("chr")
    if result.endswith(".0") and result[:-2].isdigit():
        result = result[:-2]
    return result


def info_allele_value(record: pysam.VariantRecord, field: str, index: int) -> float:
    value = record.info.get(field)
    if value is None:
        return np.nan
    if isinstance(value, tuple):
        if index >= len(value) or value[index] is None:
            return np.nan
        return float(value[index])
    return float(value)


def info_allele_count(record: pysam.VariantRecord, field: str, index: int) -> int | None:
    value = record.info.get(field)
    if value is None:
        return None
    if isinstance(value, tuple):
        if index >= len(value) or value[index] is None:
            return None
        return int(value[index])
    return int(value)


def population_sample_lists() -> dict[str, list[str]]:
    result = {}
    for population in SUPERPOPULATIONS:
        path = SAMPLE_LIST_DIR / f"samples_{population}_reference.txt"
        result[population] = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return result


def genotype_frequency(
    record: pysam.VariantRecord, alternate_index: int, samples: list[str]
) -> tuple[int, int, float]:
    allele_number = 0
    allele_count = 0
    target = alternate_index + 1
    for sample in samples:
        genotype = record.samples[sample].get("GT")
        if genotype is None:
            continue
        for allele in genotype:
            if allele is None or allele < 0:
                continue
            allele_number += 1
            allele_count += int(allele == target)
    frequency = allele_count / allele_number if allele_number else np.nan
    return allele_count, allele_number, frequency


def annotate_frequencies(candidates: pd.DataFrame) -> pd.DataFrame:
    handles: dict[str, pysam.VariantFile] = {}
    population_samples = population_sample_lists()
    rows: list[dict[str, Any]] = []
    for row in candidates.itertuples(index=False):
        chrom = chromosome(row.chromosome_grch38)
        ref = "" if pd.isna(row.ref) else str(row.ref).upper()
        alt = "" if pd.isna(row.alt) else str(row.alt).upper()
        if not chrom or pd.isna(row.position_grch38) or not ref or not alt:
            rows.append(
                {
                    "candidate_record_id": row.candidate_record_id,
                    "frequency_match_status": "unavailable_candidate_alleles_or_coordinate",
                }
            )
            continue
        vcf_path = VCF_DIR / (
            f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        )
        if not vcf_path.exists():
            rows.append(
                {
                    "candidate_record_id": row.candidate_record_id,
                    "frequency_match_status": "chromosome_vcf_unavailable",
                }
            )
            continue
        if chrom not in handles:
            handles[chrom] = pysam.VariantFile(str(vcf_path))
        position = int(row.position_grch38)
        exact = None
        exact_index = -1
        for record in handles[chrom].fetch(chrom, position - 1, position):
            if record.pos != position or record.ref.upper() != ref:
                continue
            for index, observed_alt in enumerate(record.alts or ()):
                if observed_alt.upper() == alt:
                    exact = record
                    exact_index = index
                    break
            if exact is not None:
                break
        if exact is None:
            rows.append(
                {
                    "candidate_record_id": row.candidate_record_id,
                    "frequency_match_status": "not_exactly_matched_in_1000g",
                }
            )
            continue
        global_ac = info_allele_count(exact, "AC", exact_index)
        global_an_raw = exact.info.get("AN")
        global_an = int(global_an_raw) if global_an_raw is not None else None
        annotation: dict[str, Any] = {
            "candidate_record_id": row.candidate_record_id,
            "frequency_match_status": "exact_ref_alt_match",
            "global_alt_ac": global_ac,
            "global_allele_number": global_an,
            "global_alt_af": (
                global_ac / global_an
                if global_ac is not None and global_an not in (None, 0)
                else info_allele_value(exact, "AF", exact_index)
            ),
        }
        for population in SUPERPOPULATIONS:
            allele_count, allele_number, frequency = genotype_frequency(
                exact, exact_index, population_samples[population]
            )
            annotation[f"{population}_alt_ac"] = allele_count
            annotation[f"{population}_allele_number"] = allele_number
            annotation[f"{population}_alt_af"] = frequency
        rows.append(annotation)
    for handle in handles.values():
        handle.close()

    annotations = pd.DataFrame(rows)
    result = candidates.merge(
        annotations, on="candidate_record_id", how="left", validate="one_to_one"
    )
    af_columns = [f"{population}_alt_af" for population in SUPERPOPULATIONS]
    result["max_superpopulation_alt_af"] = result[af_columns].max(axis=1, skipna=True)
    result["min_superpopulation_alt_af"] = result[af_columns].min(axis=1, skipna=True)
    maximum = result["max_superpopulation_alt_af"]
    result["frequency_class"] = np.select(
        [maximum.ge(0.05), maximum.ge(0.01) & maximum.lt(0.05), maximum.lt(0.01)],
        ["common", "low_frequency", "rare"],
        default="unavailable",
    )
    af_matrix = result[af_columns]
    result["strict_population_specific"] = (
        af_matrix.ge(0.01).sum(axis=1).eq(1) & af_matrix.lt(0.001).sum(axis=1).eq(4)
    )
    result.loc[af_matrix.notna().sum(axis=1).lt(5), "strict_population_specific"] = False
    result["population_specific_definition"] = (
        "ALT AF >=0.01 in exactly one 1000G superpopulation and <0.001 in each of the other four"
    )
    return result


def minimal_alleles(ref: str, alt: str) -> tuple[str, str]:
    left, right = ref.upper(), alt.upper()
    # Event representation used by VEP (for example C/CT becomes -/T).
    while left and right and left[-1] == right[-1]:
        left, right = left[:-1], right[:-1]
    while left and right and left[0] == right[0]:
        left, right = left[1:], right[1:]
    return left or "-", right or "-"


def allele_pair_matches(candidate_ref: str, candidate_alt: str, allele_string: str) -> bool:
    observed = {item.upper() or "-" for item in allele_string.split("/")}
    direct = {candidate_ref.upper(), candidate_alt.upper()}
    minimal = set(minimal_alleles(candidate_ref, candidate_alt))
    return direct.issubset(observed) or minimal.issubset(observed)


def request_json(
    session: requests.Session,
    method: str,
    endpoint: str,
    *,
    payload: object | None = None,
    params: dict[str, object] | None = None,
) -> object:
    url = f"{ENSEMBL_REST}{endpoint}"
    for attempt in range(6):
        response = session.request(
            method,
            url,
            json=payload,
            params=params,
            timeout=120,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        if response.status_code == 429:
            delay = float(response.headers.get("Retry-After", 2 ** attempt))
            time.sleep(min(delay, 30.0))
            continue
        if response.status_code >= 500:
            time.sleep(min(2 ** attempt, 30))
            continue
        response.raise_for_status()
        return response.json()
    raise RuntimeError(f"Ensembl request failed after retries: {method} {url}")


def causal_input_hash(causal: pd.DataFrame) -> str:
    values = "\n".join(sorted(causal["candidate_record_id"].astype(str))) + "\n"
    return hashlib.sha256(values.encode()).hexdigest()


def write_json_gzip(path: Path, value: object) -> None:
    serialized = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0
    ) as compressed:
        compressed.write(serialized)


def read_json_gzip(path: Path) -> object:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def query_ancestral(causal: pd.DataFrame) -> dict[str, object]:
    expected_hash = causal_input_hash(causal)
    if OUT_RAW.exists():
        cached = read_json_gzip(OUT_RAW)
        if isinstance(cached, dict) and cached.get("candidate_input_sha256") == expected_hash:
            return cached

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    release_information = request_json(session, "GET", "/info/data")
    candidate_to_variants: dict[str, list[str]] = {}
    valid = causal[
        causal["chromosome_grch38"].notna()
        & causal["position_grch38"].notna()
        & causal["ref"].notna()
        & causal["alt"].notna()
    ].copy()

    records = list(valid.itertuples(index=False))
    for start in range(0, len(records), 150):
        batch = records[start : start + 150]
        token_to_candidate: dict[str, str] = {}
        vcf_rows = []
        candidate_lookup = {}
        for offset, row in enumerate(batch):
            token = f"COPD{start + offset:06d}"
            candidate = str(row.candidate_record_id)
            token_to_candidate[token] = candidate
            candidate_lookup[candidate] = row
            vcf_rows.append(
                f"{chromosome(row.chromosome_grch38)} {int(row.position_grch38)} "
                f"{token} {str(row.ref).upper()} {str(row.alt).upper()} . . ."
            )
        response = request_json(
            session,
            "POST",
            "/vep/human/region",
            payload={"variants": vcf_rows},
            params={"pick": 1, "variant_class": 1},
        )
        if not isinstance(response, list):
            raise ValueError("unexpected Ensembl VEP response")
        for item in response:
            token = str(item.get("id", ""))
            candidate = token_to_candidate.get(token)
            if candidate is None:
                # The third VCF field is also recoverable from the input string.
                fields = str(item.get("input", "")).split()
                if len(fields) >= 3:
                    candidate = token_to_candidate.get(fields[2])
            if candidate is None:
                continue
            row = candidate_lookup[candidate]
            identifiers = []
            normalized_input_matches = allele_pair_matches(
                str(row.ref), str(row.alt), str(item.get("allele_string", ""))
            )
            for colocated in item.get("colocated_variants", []) or []:
                identifier = str(colocated.get("id", ""))
                allele_string = str(colocated.get("allele_string", ""))
                # In repetitive indels, the colocated dbSNP mapping is often
                # left-normalized across the repeat (for example a poly-T
                # allele string), while VEP's input allele remains -/T.  A
                # normalized input match therefore identifies the colocated
                # record even when literal alleles differ.
                if identifier and (
                    normalized_input_matches
                    or allele_pair_matches(str(row.ref), str(row.alt), allele_string)
                ):
                    identifiers.append(identifier)
            candidate_to_variants[candidate] = sorted(
                set(identifiers), key=lambda value: (not value.startswith("rs"), value)
            )

    identifiers = sorted(
        {identifier for values in candidate_to_variants.values() for identifier in values}
    )
    variation_records: dict[str, object] = {}
    for start in range(0, len(identifiers), 150):
        batch = identifiers[start : start + 150]
        response = request_json(
            session,
            "POST",
            "/variation/human",
            payload={"ids": batch},
        )
        if not isinstance(response, dict):
            raise ValueError("unexpected Ensembl Variation response")
        variation_records.update(response)

    cache: dict[str, object] = {
        "source": "Ensembl REST VEP and Variation endpoints",
        "base_url": ENSEMBL_REST,
        "accessed_utc": datetime.now(timezone.utc).isoformat(),
        "release_information": release_information,
        "candidate_input_sha256": expected_hash,
        "candidate_count": len(causal),
        "candidate_to_variation_ids": candidate_to_variants,
        "variation_records": variation_records,
    }
    write_json_gzip(OUT_RAW, cache)
    return cache


def ancestral_annotations(causal: pd.DataFrame, cache: dict[str, object]) -> pd.DataFrame:
    candidate_to_ids = cache.get("candidate_to_variation_ids", {})
    variation_records = cache.get("variation_records", {})
    rows = []
    for row in causal.itertuples(index=False):
        candidate = str(row.candidate_record_id)
        ref = "" if pd.isna(row.ref) else str(row.ref).upper()
        alt = "" if pd.isna(row.alt) else str(row.alt).upper()
        chrom = chromosome(row.chromosome_grch38)
        identifiers = candidate_to_ids.get(candidate, []) if isinstance(candidate_to_ids, dict) else []
        selected_id = ""
        ancestral = ""
        mapping_alleles = ""
        mapping_reference_equivalent = ""
        mapping_alternate_equivalent = ""
        for identifier in identifiers:
            record = variation_records.get(identifier, {}) if isinstance(variation_records, dict) else {}
            for mapping in record.get("mappings", []) if isinstance(record, dict) else []:
                if (
                    str(mapping.get("assembly_name", "")) != "GRCh38"
                    or str(mapping.get("seq_region_name", "")) != chrom
                ):
                    continue
                value = mapping.get("ancestral_allele")
                if value not in (None, "", "-") or "-" in minimal_alleles(ref, alt):
                    selected_id = str(identifier)
                    ancestral = "" if value is None else str(value).upper()
                    mapping_alleles = str(mapping.get("allele_string", ""))
                    observed_mapping_alleles = [
                        value.upper() or "-" for value in mapping_alleles.split("/")
                    ]
                    if observed_mapping_alleles:
                        mapping_reference_equivalent = observed_mapping_alleles[0]
                        length_delta = len(alt) - len(ref)
                        equivalent = [
                            allele
                            for allele in observed_mapping_alleles[1:]
                            if len(allele.replace("-", ""))
                            - len(mapping_reference_equivalent.replace("-", ""))
                            == length_delta
                        ]
                        if len(equivalent) == 1:
                            mapping_alternate_equivalent = equivalent[0]
                    break
            if selected_id:
                break

        normalized_ref, normalized_alt = minimal_alleles(ref, alt) if ref and alt else ("", "")
        if ancestral in {ref, normalized_ref}:
            status = "ancestral_matches_reference"
            derived = "ALT"
        elif ancestral in {alt, normalized_alt}:
            status = "ancestral_matches_alternate"
            derived = "REF"
        elif ancestral and ancestral == mapping_reference_equivalent:
            status = "ancestral_matches_normalized_mapping_reference"
            derived = "ALT"
        elif ancestral and ancestral == mapping_alternate_equivalent:
            status = "ancestral_matches_normalized_mapping_alternate"
            derived = "REF"
        elif identifiers and ancestral:
            status = "ancestral_not_mappable_to_candidate_alleles"
            derived = ""
        elif identifiers:
            status = "ancestral_allele_unavailable"
            derived = ""
        else:
            status = "no_exact_ensembl_variation_identifier"
            derived = ""
        rows.append(
            {
                "candidate_record_id": candidate,
                "ensembl_variation_ids": ";".join(identifiers),
                "selected_ensembl_variation_id": selected_id,
                "ensembl_mapping_alleles": mapping_alleles,
                "mapping_reference_equivalent": mapping_reference_equivalent,
                "mapping_alternate_equivalent": mapping_alternate_equivalent,
                "ancestral_allele": ancestral,
                "ancestral_call_status": status,
                "derived_candidate_allele": derived,
            }
        )
    return pd.DataFrame(rows)


def add_derived_frequencies(data: pd.DataFrame) -> pd.DataFrame:
    data["global_derived_allele_frequency"] = np.select(
        [data["derived_candidate_allele"].eq("ALT"), data["derived_candidate_allele"].eq("REF")],
        [data["global_alt_af"], 1.0 - data["global_alt_af"]],
        default=np.nan,
    )
    derived_columns = []
    for population in SUPERPOPULATIONS:
        alt_column = f"{population}_alt_af"
        derived_column = f"{population}_derived_allele_frequency"
        data[derived_column] = np.select(
            [data["derived_candidate_allele"].eq("ALT"), data["derived_candidate_allele"].eq("REF")],
            [data[alt_column], 1.0 - data[alt_column]],
            default=np.nan,
        )
        derived_columns.append(derived_column)
    data["max_superpopulation_derived_allele_frequency"] = data[derived_columns].max(
        axis=1, skipna=True
    )
    data["global_DAF_gt_0_5"] = data["global_derived_allele_frequency"].gt(0.5).astype(
        "boolean"
    )
    data["any_superpopulation_DAF_gt_0_5"] = data[
        "max_superpopulation_derived_allele_frequency"
    ].gt(0.5).astype("boolean")
    unresolved = data["derived_candidate_allele"].eq("")
    data.loc[unresolved, ["global_DAF_gt_0_5", "any_superpopulation_DAF_gt_0_5"]] = pd.NA
    return data


def make_summary(causal: pd.DataFrame) -> pd.DataFrame:
    n = len(causal)
    rows = []

    def add(metric: str, mask: pd.Series, denominator: int = n, note: str = "") -> None:
        count = int(mask.fillna(False).sum())
        rows.append(
            {
                "metric": metric,
                "n": count,
                "denominator": denominator,
                "fraction": count / denominator if denominator else np.nan,
                "note": note,
            }
        )

    add("predicted_causal_variants", pd.Series(True, index=causal.index))
    matched = causal["frequency_match_status"].eq("exact_ref_alt_match")
    add("frequency_resolved", matched)
    for category in ("common", "low_frequency", "rare", "unavailable"):
        add(f"frequency_class_{category}", causal["frequency_class"].eq(category))
    add(
        "strict_population_specific",
        causal["strict_population_specific"],
        note="ALT AF >=0.01 in exactly one superpopulation and <0.001 in the other four",
    )
    ancestral_resolved = causal["derived_candidate_allele"].isin(["REF", "ALT"])
    add("ancestral_state_resolved", ancestral_resolved)
    resolved_n = int(ancestral_resolved.sum())
    add(
        "global_DAF_gt_0_5",
        causal.loc[ancestral_resolved, "global_DAF_gt_0_5"],
        denominator=resolved_n,
        note="denominator is predicted causal variants with an Ensembl ancestral allele mappable to REF or ALT",
    )
    add(
        "any_superpopulation_DAF_gt_0_5",
        causal.loc[ancestral_resolved, "any_superpopulation_DAF_gt_0_5"],
        denominator=resolved_n,
        note="denominator is predicted causal variants with an Ensembl ancestral allele mappable to REF or ALT",
    )
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(INPUT, sep="\t", low_memory=False)
    candidates["predicted_causal_regulatory"] = as_bool(
        candidates["predicted_causal_regulatory"]
    )
    annotated = annotate_frequencies(candidates)
    causal = annotated[annotated["predicted_causal_regulatory"]].copy()

    try:
        cache = query_ancestral(causal)
        ancestral = ancestral_annotations(causal, cache)
        ancestral_query_status = "completed"
        ancestral_query_error = ""
    except Exception as exc:
        ancestral_query_status = "failed"
        ancestral_query_error = f"{type(exc).__name__}: {exc}"
        ancestral = pd.DataFrame(
            {
                "candidate_record_id": causal["candidate_record_id"],
                "ensembl_variation_ids": "",
                "selected_ensembl_variation_id": "",
                "ensembl_mapping_alleles": "",
                "mapping_reference_equivalent": "",
                "mapping_alternate_equivalent": "",
                "ancestral_allele": "",
                "ancestral_call_status": "ensembl_query_failed",
                "derived_candidate_allele": "",
            }
        )
        cache = {
            "source": "Ensembl REST VEP and Variation endpoints",
            "base_url": ENSEMBL_REST,
            "accessed_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_input_sha256": causal_input_hash(causal),
            "candidate_count": len(causal),
            "status": ancestral_query_status,
            "error": ancestral_query_error,
        }
        write_json_gzip(OUT_RAW, cache)

    annotated = annotated.merge(
        ancestral, on="candidate_record_id", how="left", validate="one_to_one"
    )
    for column in (
        "ensembl_variation_ids",
        "selected_ensembl_variation_id",
        "ensembl_mapping_alleles",
        "mapping_reference_equivalent",
        "mapping_alternate_equivalent",
        "ancestral_allele",
        "ancestral_call_status",
        "derived_candidate_allele",
    ):
        annotated[column] = annotated[column].fillna("")
    annotated = add_derived_frequencies(annotated)
    causal = annotated[annotated["predicted_causal_regulatory"]].copy()
    summary = make_summary(causal)

    annotated.to_csv(
        OUT_ALL,
        sep="\t",
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    causal.to_csv(OUT_CAUSAL, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    ancestral.to_csv(OUT_ANCESTRAL, sep="\t", index=False)

    outputs = {
        "population_annotated_candidates": OUT_ALL,
        "predicted_causal_population_genetics": OUT_CAUSAL,
        "population_summary": OUT_SUMMARY,
        "ancestral_allele_audit": OUT_ANCESTRAL,
        "ensembl_ancestral_cache": OUT_RAW,
    }
    manifest = {
        "result_id": "COPD-S4-R005",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "frequency_source": "1000 Genomes Phase 3 20190312 GRCh38 remap",
        "frequency_calculation": "global ALT AF recalculated from exact AC/AN; superpopulation ALT AF recalculated from called GT alleles in the five audited effective sample lists rather than rounded INFO AF fields",
        "superpopulations": list(SUPERPOPULATIONS),
        "frequency_class_definition": {
            "common": "maximum alternate-allele frequency across five superpopulations >=0.05",
            "low_frequency": "maximum alternate-allele frequency >=0.01 and <0.05",
            "rare": "maximum alternate-allele frequency <0.01",
            "population_specific": "alternate-allele frequency >=0.01 in exactly one superpopulation and <0.001 in each of the other four",
        },
        "ancestral_source": "Ensembl REST VEP and Variation endpoints",
        "ancestral_query_status": ancestral_query_status,
        "ancestral_query_error": ancestral_query_error,
        "candidate_records": len(annotated),
        "predicted_causal_records": len(causal),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pysam": pysam.__version__,
            "requests": requests.__version__,
        },
        "inputs": {
            "prioritized_candidates": {
                "path": str(INPUT),
                "bytes": INPUT.stat().st_size,
                "sha256": sha256(INPUT),
            }
        },
        "outputs": {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in outputs.items()
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        f"Annotated {len(annotated):,} candidates and {len(causal):,} predicted causal variants; "
        f"ancestral query status: {ancestral_query_status}."
    )


if __name__ == "__main__":
    main()
