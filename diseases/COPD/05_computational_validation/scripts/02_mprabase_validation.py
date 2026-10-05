#!/usr/bin/env python3
"""Audit COPD candidate overlap with MPRAbase v4.9.3.

MPRAbase human genomic coordinates are hg19. Candidate positions are lifted as
one-base BED intervals from GRCh38 to hg19 and required to map uniquely and
round-trip exactly before coordinate comparison. Coordinate containment means
that an assayed library element covered the base; it is not treated as proof
that both candidate alleles were experimentally compared.
"""

from __future__ import annotations

import bisect
import re
import shutil
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path

from common import (
    DATA,
    LOGS,
    RESULTS,
    ROOT,
    candidate_core,
    candidate_sort_key,
    clean_chromosome,
    ensure_directories,
    extract_rsids,
    file_record,
    integer,
    load_candidates,
    numeric,
    sha256,
    tokens,
    write_manifest,
    write_tsv,
)


RESULT_ID = "COPD-S5-R002"
DATABASE = ROOT / "data/mprabase/mprabase_v4_9.3.db"
HG38_TO_HG19 = ROOT / "data/liftover/hg38ToHg19.over.chain.gz"
HG19_TO_HG38 = ROOT / "data/liftover/hg19ToHg38.over.chain.gz"

OUT_ELEMENTS = RESULTS / f"{RESULT_ID}_MPRAbase_v4_9_3_candidate_elements.tsv.gz"
OUT_LIFT = RESULTS / f"{RESULT_ID}_GRCh38_to_hg19_liftover_audit.tsv"
OUT_AUDIT = RESULTS / f"{RESULT_ID}_MPRAbase_candidate_audit.tsv"
OUT_SUMMARY = RESULTS / f"{RESULT_ID}_MPRAbase_summary.tsv"
OUT_MANIFEST = RESULTS / f"{RESULT_ID}_manifest.json"
LOG = LOGS / "02_mprabase_validation.log"

LIFT_INPUT = DATA / f"{RESULT_ID}_candidate_sites_GRCh38.bed"
LIFT_MAPPED = DATA / f"{RESULT_ID}_candidate_sites_hg19.bed"
LIFT_UNMAPPED = DATA / f"{RESULT_ID}_candidate_sites_hg19_unmapped.bed"
ROUNDTRIP_INPUT = DATA / f"{RESULT_ID}_candidate_sites_hg19_for_roundtrip.bed"
ROUNDTRIP_MAPPED = DATA / f"{RESULT_ID}_candidate_sites_GRCh38_roundtrip.bed"
ROUNDTRIP_UNMAPPED = DATA / f"{RESULT_ID}_candidate_sites_GRCh38_roundtrip_unmapped.bed"

COORDINATE_RE = re.compile(r"^chr([^:]+):(\d+)-(\d+)$")
CANONICAL = {str(value) for value in range(1, 23)} | {"X", "Y"}

LIFT_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "candidate_rsids",
    "liftover_token",
    "liftover_status",
    "hg19_mapping_count",
    "chromosome_hg19",
    "position_hg19",
    "roundtrip_status",
    "coordinate_match_eligible",
]

ELEMENT_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "variant_class",
    "candidate_rsids",
    "chromosome_hg19",
    "position_hg19",
    "match_channels",
    "coordinate_relation",
    "library_element_id",
    "library_element_name",
    "library_id",
    "library_name",
    "dataset_id",
    "dataset_name",
    "dataset_pmid",
    "dataset_reference",
    "mprabase_genome_build",
    "mprabase_element_coordinate",
    "sequence_length",
    "exact_candidate_rsids_in_element_metadata",
    "n_element_sample_scores",
    "minimum_reported_score",
    "maximum_reported_score",
    "sample_ids",
    "sample_names",
    "cell_lines_or_tissues",
    "library_strategies",
    "evidence_scope",
]

AUDIT_FIELDS = [
    "candidate_record_id",
    "predicted_causal_priority_rank",
    "chromosome_grch38",
    "position_grch38",
    "ref",
    "alt",
    "candidate_rsids",
    "liftover_status",
    "chromosome_hg19",
    "position_hg19",
    "roundtrip_status",
    "coordinate_match_eligible",
    "mprabase_query_status",
    "matched_element_count",
    "coordinate_matched_element_count",
    "exact_rsid_matched_element_count",
    "elements_with_reported_scores",
    "matched_dataset_count",
    "matched_datasets",
]


def locate_liftover() -> Path:
    in_path = shutil.which("liftOver")
    candidates = [
        Path(in_path) if in_path else None,
        Path("/usr/local/apps/ucsc/503/bin/x86_64/liftOver"),
    ]
    for candidate in candidates:
        if candidate is not None and candidate.exists() and candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "UCSC liftOver not found; load ucsc/503 or update locate_liftover()"
    )


def read_bed(path: Path) -> list[tuple[str, int, int, str]]:
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return rows


def run_liftover(
    binary: Path,
    input_bed: Path,
    chain: Path,
    mapped: Path,
    unmapped: Path,
) -> str:
    process = subprocess.run(
        [
            str(binary),
            "-minMatch=0.95",
            str(input_bed),
            str(chain),
            str(mapped),
            str(unmapped),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return (process.stdout + process.stderr).strip()


def liftover_candidates(
    candidates: list[dict[str, str]], binary: Path
) -> tuple[list[dict[str, object]], str]:
    token_to_candidate: dict[str, dict[str, str]] = {}
    original_bed: dict[str, tuple[str, int, int]] = {}
    input_lines = []
    for index, candidate in enumerate(candidates):
        token = f"COPDS5V{index:06d}"
        candidate["_liftover_token"] = token
        token_to_candidate[token] = candidate
        chrom = clean_chromosome(candidate.get("chromosome_grch38", ""))
        position = integer(candidate.get("position_grch38", ""))
        if chrom in CANONICAL and position is not None and position > 0:
            interval = (f"chr{chrom}", position - 1, position)
            original_bed[token] = interval
            input_lines.append(f"{interval[0]}\t{interval[1]}\t{interval[2]}\t{token}\n")
    LIFT_INPUT.write_text("".join(input_lines))

    forward_log = run_liftover(
        binary, LIFT_INPUT, HG38_TO_HG19, LIFT_MAPPED, LIFT_UNMAPPED
    )
    mapped_by_token: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for chrom, start, end, token in read_bed(LIFT_MAPPED):
        mapped_by_token[token].append((chrom, start, end))

    roundtrip_input_lines = []
    for token, mappings in sorted(mapped_by_token.items()):
        for mapping_index, (chrom, start, end) in enumerate(mappings):
            roundtrip_token = f"{token}__M{mapping_index}"
            roundtrip_input_lines.append(
                f"{chrom}\t{start}\t{end}\t{roundtrip_token}\n"
            )
    ROUNDTRIP_INPUT.write_text("".join(roundtrip_input_lines))
    reverse_log = run_liftover(
        binary,
        ROUNDTRIP_INPUT,
        HG19_TO_HG38,
        ROUNDTRIP_MAPPED,
        ROUNDTRIP_UNMAPPED,
    )
    reverse_by_token: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for chrom, start, end, token in read_bed(ROUNDTRIP_MAPPED):
        original_token = token.split("__M", 1)[0]
        reverse_by_token[original_token].append((chrom, start, end))

    audit_rows = []
    for candidate in candidates:
        token = candidate["_liftover_token"]
        core = candidate_core(candidate)
        mappings = mapped_by_token.get(token, [])
        if token not in original_bed:
            status = "not_queryable_noncanonical_or_missing_coordinate"
        elif not mappings:
            status = "unmapped"
        elif len(mappings) == 1:
            status = "unique_one_base_mapping"
        else:
            status = "multiple_mappings"

        hg19_chrom = clean_chromosome(mappings[0][0]) if len(mappings) == 1 else ""
        hg19_position = mappings[0][1] + 1 if len(mappings) == 1 else ""
        reverse = reverse_by_token.get(token, [])
        if len(mappings) != 1:
            roundtrip_status = "not_tested_without_unique_forward_mapping"
        elif not reverse:
            roundtrip_status = "roundtrip_unmapped"
        elif len(reverse) != 1:
            roundtrip_status = "roundtrip_multiple_mappings"
        elif reverse[0] == original_bed[token]:
            roundtrip_status = "roundtrip_exact"
        else:
            roundtrip_status = "roundtrip_nonexact"
        eligible = status == "unique_one_base_mapping" and roundtrip_status == "roundtrip_exact"
        audit_rows.append(
            {
                **core,
                "liftover_token": token,
                "liftover_status": status,
                "hg19_mapping_count": len(mappings),
                "chromosome_hg19": hg19_chrom,
                "position_hg19": hg19_position,
                "roundtrip_status": roundtrip_status,
                "coordinate_match_eligible": eligible,
            }
        )
    return audit_rows, f"forward: {forward_log}\nreverse: {reverse_log}".strip()


def load_database_metadata(connection: sqlite3.Connection):
    libraries = {
        row[0]: {"library_name": row[1] or "", "dataset_id": row[2] or ""}
        for row in connection.execute(
            "select library_id, library_name, datasets_id from designed_library"
        )
    }
    datasets = {
        row[0]: {
            "dataset_name": row[1] or "",
            "dataset_pmid": row[2] or "",
            "dataset_reference": row[3] or "",
        }
        for row in connection.execute(
            "select datasets_id, datasets_name, PMID, Reference from datasets"
        )
    }
    return libraries, datasets


def aggregate_scores(
    connection: sqlite3.Connection, element_ids: list[str]
) -> dict[str, dict[str, object]]:
    aggregation: dict[str, dict[str, object]] = {}
    for start in range(0, len(element_ids), 500):
        batch = element_ids[start : start + 500]
        placeholders = ",".join("?" for _ in batch)
        query = f"""
            select es.library_element_id, es.score, s.sample_id, s.sample_name,
                   s.Cell_line_tissue, s.Library_strategy
            from element_score es
            left join sample s on es.sample_id = s.sample_id
            where es.library_element_id in ({placeholders})
        """
        for element_id, score, sample_id, sample_name, tissue, strategy in connection.execute(
            query, batch
        ):
            item = aggregation.setdefault(
                element_id,
                {
                    "scores": [],
                    "sample_ids": set(),
                    "sample_names": set(),
                    "tissues": set(),
                    "strategies": set(),
                },
            )
            value = numeric(score)
            if value is not None:
                item["scores"].append(value)
            if sample_id:
                item["sample_ids"].add(str(sample_id).strip())
            if sample_name:
                item["sample_names"].add(str(sample_name).strip())
            if tissue:
                item["tissues"].add(str(tissue).strip())
            if strategy:
                item["strategies"].add(str(strategy).strip())
    return aggregation


def match_mprabase(
    connection: sqlite3.Connection,
    candidates: list[dict[str, str]],
    lift_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    candidate_lookup = {row["candidate_record_id"]: row for row in candidates}
    eligible_by_chrom: dict[str, list[tuple[int, str]]] = defaultdict(list)
    candidate_rsid_map: dict[str, set[str]] = defaultdict(set)
    for lift in lift_rows:
        candidate_id = str(lift["candidate_record_id"])
        if lift["coordinate_match_eligible"]:
            eligible_by_chrom[str(lift["chromosome_hg19"])].append(
                (int(lift["position_hg19"]), candidate_id)
            )
        for rsid in tokens(candidate_lookup[candidate_id].get("candidate_rsids", "")):
            candidate_rsid_map[rsid.lower()].add(candidate_id)
    for chrom in eligible_by_chrom:
        eligible_by_chrom[chrom].sort()

    matches: dict[tuple[str, str], dict[str, object]] = {}
    source_rows = 0
    hg19_rows = 0
    parseable_hg19_rows = 0
    exact_rsid_metadata_rows = 0
    for element_id, library_id, build, coordinate, name, sequence in connection.execute(
        """select library_element_id, library_id, genome_build,
                  element_coordinate, library_element_name, sequence
           from library_sequence"""
    ):
        source_rows += 1
        build = (build or "").strip()
        coordinate = (coordinate or "").strip()
        name = (name or "").strip()
        element_rsids = set(extract_rsids(f"{coordinate} {name}"))
        id_candidates: set[str] = set()
        for rsid in element_rsids:
            id_candidates.update(candidate_rsid_map.get(rsid, set()))
        if id_candidates:
            exact_rsid_metadata_rows += 1
        for candidate_id in id_candidates:
            key = (candidate_id, element_id)
            item = matches.setdefault(
                key,
                {
                    "candidate_record_id": candidate_id,
                    "library_element_id": element_id,
                    "library_id": library_id,
                    "mprabase_genome_build": build,
                    "mprabase_element_coordinate": coordinate,
                    "library_element_name": name,
                    "sequence_length": len(sequence or ""),
                    "match_channels": set(),
                    "coordinate_relations": set(),
                    "exact_rsids": set(),
                },
            )
            item["match_channels"].add("exact_candidate_rsid_in_element_metadata")
            item["exact_rsids"].update(
                element_rsids & set(tokens(candidate_lookup[candidate_id]["candidate_rsids"]))
            )

        if build != "hg19":
            continue
        hg19_rows += 1
        match = COORDINATE_RE.fullmatch(coordinate)
        if not match:
            continue
        parseable_hg19_rows += 1
        chrom = clean_chromosome(match.group(1))
        interval_start = int(match.group(2))
        interval_end = int(match.group(3))
        if interval_end < interval_start:
            interval_start, interval_end = interval_end, interval_start
        positioned = eligible_by_chrom.get(chrom, [])
        if not positioned:
            continue
        positions = [item[0] for item in positioned]
        left = bisect.bisect_left(positions, interval_start)
        right = bisect.bisect_right(positions, interval_end)
        for position, candidate_id in positioned[left:right]:
            if interval_start == interval_end == position:
                relation = "single_base_exact_coordinate"
            elif position in {interval_start, interval_end}:
                relation = "interval_boundary_coordinate"
            else:
                relation = "interval_interior_coordinate"
            key = (candidate_id, element_id)
            item = matches.setdefault(
                key,
                {
                    "candidate_record_id": candidate_id,
                    "library_element_id": element_id,
                    "library_id": library_id,
                    "mprabase_genome_build": build,
                    "mprabase_element_coordinate": coordinate,
                    "library_element_name": name,
                    "sequence_length": len(sequence or ""),
                    "match_channels": set(),
                    "coordinate_relations": set(),
                    "exact_rsids": set(),
                },
            )
            item["match_channels"].add("roundtrip_exact_hg19_coordinate_containment")
            item["coordinate_relations"].add(relation)

    libraries, datasets = load_database_metadata(connection)
    element_ids = sorted({key[1] for key in matches})
    score_metadata = aggregate_scores(connection, element_ids)
    lift_by_candidate = {str(row["candidate_record_id"]): row for row in lift_rows}
    rows: list[dict[str, object]] = []
    for (candidate_id, _), match in matches.items():
        candidate = candidate_lookup[candidate_id]
        core = candidate_core(candidate)
        lift = lift_by_candidate[candidate_id]
        library = libraries.get(str(match["library_id"]), {})
        dataset = datasets.get(str(library.get("dataset_id", "")), {})
        scores = score_metadata.get(str(match["library_element_id"]), {})
        values = scores.get("scores", [])
        channels = sorted(match["match_channels"])
        coordinate_channel = "roundtrip_exact_hg19_coordinate_containment" in channels
        exact_id_channel = "exact_candidate_rsid_in_element_metadata" in channels
        if exact_id_channel and coordinate_channel:
            scope = "exact candidate rsID metadata plus assayed-element coordinate containment"
        elif exact_id_channel:
            scope = "exact candidate rsID in element metadata"
        else:
            scope = (
                "candidate base lies within an hg19-indexed assayed element; "
                "candidate allele comparison is not established"
            )
        rows.append(
            {
                **core,
                "chromosome_hg19": lift["chromosome_hg19"],
                "position_hg19": lift["position_hg19"],
                "match_channels": ";".join(channels),
                "coordinate_relation": ";".join(sorted(match["coordinate_relations"])),
                "library_element_id": match["library_element_id"],
                "library_element_name": match["library_element_name"],
                "library_id": match["library_id"],
                "library_name": library.get("library_name", ""),
                "dataset_id": library.get("dataset_id", ""),
                "dataset_name": dataset.get("dataset_name", ""),
                "dataset_pmid": dataset.get("dataset_pmid", ""),
                "dataset_reference": dataset.get("dataset_reference", ""),
                "mprabase_genome_build": match["mprabase_genome_build"],
                "mprabase_element_coordinate": match["mprabase_element_coordinate"],
                "sequence_length": match["sequence_length"],
                "exact_candidate_rsids_in_element_metadata": ";".join(
                    sorted(match["exact_rsids"], key=lambda item: int(item[2:]))
                ),
                "n_element_sample_scores": len(values),
                "minimum_reported_score": min(values) if values else None,
                "maximum_reported_score": max(values) if values else None,
                "sample_ids": ";".join(sorted(scores.get("sample_ids", set()))),
                "sample_names": ";".join(sorted(scores.get("sample_names", set()))),
                "cell_lines_or_tissues": ";".join(sorted(scores.get("tissues", set()))),
                "library_strategies": ";".join(sorted(scores.get("strategies", set()))),
                "evidence_scope": scope,
            }
        )
    rows.sort(
        key=lambda row: (
            integer(row["predicted_causal_priority_rank"]) or 10**12,
            str(row["candidate_record_id"]),
            str(row["dataset_id"]),
            str(row["library_element_id"]),
        )
    )
    source_counts = {
        "library_sequence_rows": source_rows,
        "hg19_rows": hg19_rows,
        "parseable_hg19_coordinate_rows": parseable_hg19_rows,
        "rows_with_exact_candidate_rsid_metadata": exact_rsid_metadata_rows,
    }
    return rows, source_counts


def main() -> None:
    ensure_directories()
    candidate_input, candidates = load_candidates()
    candidates.sort(key=candidate_sort_key)
    liftover_binary = locate_liftover()
    lift_rows, liftover_log = liftover_candidates(candidates, liftover_binary)

    connection = sqlite3.connect(f"file:{DATABASE}?mode=ro", uri=True)
    try:
        element_rows, source_counts = match_mprabase(
            connection, candidates, lift_rows
        )
        database_counts = {
            table: connection.execute(f"select count(*) from {table}").fetchone()[0]
            for table in (
                "datasets",
                "designed_library",
                "library_sequence",
                "sample",
                "element_score",
            )
        }
    finally:
        connection.close()

    elements_by_candidate: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in element_rows:
        elements_by_candidate[str(row["candidate_record_id"])].append(row)
    lift_by_candidate = {str(row["candidate_record_id"]): row for row in lift_rows}
    audit_rows = []
    for candidate in candidates:
        core = candidate_core(candidate)
        candidate_id = str(core["candidate_record_id"])
        lift = lift_by_candidate[candidate_id]
        matches = elements_by_candidate.get(candidate_id, [])
        coordinate_matches = [
            row
            for row in matches
            if "roundtrip_exact_hg19_coordinate_containment"
            in str(row["match_channels"])
        ]
        id_matches = [
            row
            for row in matches
            if "exact_candidate_rsid_in_element_metadata"
            in str(row["match_channels"])
        ]
        scored = [row for row in matches if int(row["n_element_sample_scores"]) > 0]
        datasets = sorted(
            {
                str(row["dataset_name"] or row["dataset_id"])
                for row in matches
                if row["dataset_name"] or row["dataset_id"]
            }
        )
        if matches:
            status = "mprabase_element_evidence_found"
        elif not lift["coordinate_match_eligible"] and not core["candidate_rsids"]:
            status = "not_queryable_no_roundtrip_coordinate_or_exact_rsid"
        else:
            status = "no_mprabase_element_evidence"
        audit_rows.append(
            {
                **core,
                "liftover_status": lift["liftover_status"],
                "chromosome_hg19": lift["chromosome_hg19"],
                "position_hg19": lift["position_hg19"],
                "roundtrip_status": lift["roundtrip_status"],
                "coordinate_match_eligible": lift["coordinate_match_eligible"],
                "mprabase_query_status": status,
                "matched_element_count": len(matches),
                "coordinate_matched_element_count": len(coordinate_matches),
                "exact_rsid_matched_element_count": len(id_matches),
                "elements_with_reported_scores": len(scored),
                "matched_dataset_count": len(datasets),
                "matched_datasets": ";".join(datasets),
            }
        )

    candidate_count = len(candidates)
    eligible_count = sum(bool(row["coordinate_match_eligible"]) for row in lift_rows)
    any_match_count = sum(int(row["matched_element_count"]) > 0 for row in audit_rows)
    coordinate_match_count = sum(
        int(row["coordinate_matched_element_count"]) > 0 for row in audit_rows
    )
    exact_id_count = sum(
        int(row["exact_rsid_matched_element_count"]) > 0 for row in audit_rows
    )
    summary_rows = [
        {
            "metric": "predicted_causal_candidates",
            "n": candidate_count,
            "denominator": candidate_count,
            "fraction": 1.0,
            "note": "Section 4 predicted-causal union",
        },
        {
            "metric": "unique_roundtrip_exact_GRCh38_to_hg19_sites",
            "n": eligible_count,
            "denominator": candidate_count,
            "fraction": eligible_count / candidate_count,
            "note": "one-base UCSC liftOver, minMatch=0.95, exact hg19-to-GRCh38 round trip",
        },
        {
            "metric": "candidates_in_any_MPRAbase_element",
            "n": any_match_count,
            "denominator": candidate_count,
            "fraction": any_match_count / candidate_count,
            "note": "union of exact candidate-rsID metadata and coordinate containment",
        },
        {
            "metric": "candidates_in_hg19_indexed_MPRAbase_elements",
            "n": coordinate_match_count,
            "denominator": eligible_count,
            "fraction": coordinate_match_count / eligible_count if eligible_count else None,
            "note": "regional assay evidence; not necessarily an allelic MPRA",
        },
        {
            "metric": "candidates_with_exact_rsid_in_MPRAbase_element_metadata",
            "n": exact_id_count,
            "denominator": sum(bool(row["candidate_rsids"]) for row in audit_rows),
            "fraction": (
                exact_id_count / sum(bool(row["candidate_rsids"]) for row in audit_rows)
                if any(bool(row["candidate_rsids"]) for row in audit_rows)
                else None
            ),
            "note": "source focal-tag rsIDs are excluded from candidate identifiers",
        },
        {
            "metric": "matched_candidate_element_records",
            "n": len(element_rows),
            "denominator": "",
            "fraction": "",
            "note": "one row per candidate and MPRAbase library element",
        },
        {
            "metric": "matched_unique_MPRAbase_datasets",
            "n": len({row["dataset_id"] for row in element_rows}),
            "denominator": database_counts["datasets"],
            "fraction": len({row["dataset_id"] for row in element_rows})
            / database_counts["datasets"],
            "note": "scores are not comparable across studies without study-specific calibration",
        },
    ]

    write_tsv(OUT_LIFT, lift_rows, LIFT_FIELDS)
    write_tsv(OUT_ELEMENTS, element_rows, ELEMENT_FIELDS)
    write_tsv(OUT_AUDIT, audit_rows, AUDIT_FIELDS)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        ["metric", "n", "denominator", "fraction", "note"],
    )
    write_manifest(
        OUT_MANIFEST,
        {
            "result_id": RESULT_ID,
            "analysis": "MPRAbase v4.9.3 candidate element validation",
            "candidate_genome_build": "GRCh38",
            "mprabase_coordinate_build": "hg19",
            "coordinate_conversion": {
                "representation": "one-base BED interval [GRCh38 position-1, position)",
                "forward_chain": file_record(HG38_TO_HG19),
                "reverse_chain": file_record(HG19_TO_HG38),
                "minimum_match": 0.95,
                "eligibility": "unique forward mapping and exact one-base round trip",
                "binary": file_record(liftover_binary),
            },
            "mprabase_coordinate_interpretation": (
                "chr:start-end strings are compared as closed coordinate spans because "
                "the database includes point elements encoded with start=end; boundary "
                "and interior matches are retained as separate audit categories"
            ),
            "database_counts": database_counts,
            "source_scan_counts": source_counts,
            "candidate_count": candidate_count,
            "liftover_eligible_count": eligible_count,
            "matched_candidate_count": any_match_count,
            "candidate_element_record_count": len(element_rows),
            "inputs": {
                "section4_candidates": file_record(candidate_input),
                "mprabase_sqlite": file_record(DATABASE),
            },
            "outputs": {
                "candidate_elements": file_record(OUT_ELEMENTS),
                "liftover_audit": file_record(OUT_LIFT),
                "candidate_audit": file_record(OUT_AUDIT),
                "summary": file_record(OUT_SUMMARY),
                "forward_input_bed": file_record(LIFT_INPUT),
                "forward_mapped_bed": file_record(LIFT_MAPPED),
                "forward_unmapped_bed": file_record(LIFT_UNMAPPED),
                "roundtrip_input_bed": file_record(ROUNDTRIP_INPUT),
                "roundtrip_mapped_bed": file_record(ROUNDTRIP_MAPPED),
                "roundtrip_unmapped_bed": file_record(ROUNDTRIP_UNMAPPED),
            },
            "interpretation_guardrail": (
                "Coordinate containment establishes that a library element covering the "
                "candidate base was assayed. It does not establish that REF and ALT were "
                "both tested, that the reported score is allele-specific, or that the "
                "assay context is disease-relevant."
            ),
        },
    )
    LOG.write_text(
        f"result_id={RESULT_ID}\n"
        f"candidate_input={candidate_input}\n"
        f"candidate_count={candidate_count}\n"
        f"liftover_eligible={eligible_count}\n"
        f"matched_candidates={any_match_count}\n"
        f"candidate_element_records={len(element_rows)}\n"
        f"liftover_binary_sha256={sha256(liftover_binary)}\n"
        f"{liftover_log}\n"
    )
    print(
        f"{RESULT_ID}: {any_match_count:,}/{candidate_count:,} candidates had "
        f"MPRAbase element evidence ({len(element_rows):,} candidate-element records)."
    )


if __name__ == "__main__":
    main()

