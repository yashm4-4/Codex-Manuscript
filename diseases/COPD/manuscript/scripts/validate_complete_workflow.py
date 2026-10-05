#!/usr/bin/env python3
"""Read-only audit of the COPD analysis outputs and assembled manuscript.

The only files written by this program are the requested validation table and
its log in ``diseases/COPD/manuscript/results``.  All analysis outputs and the
four global traceability registers are opened read-only.

The default mode is stage-aware: missing Section 6 or manuscript deliverables
are warnings while those stages are still being assembled.  Use
``--strict-final`` for the final gate, where planned-but-missing deliverables
and traceability handoffs are failures.  Existing malformed or inconsistent
files fail in both modes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional, Sequence
from urllib.parse import unquote


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[4]
COPD = SCRIPT.parents[2]
MANUSCRIPT = COPD / "manuscript"
OUTPUT = MANUSCRIPT / "results" / "COPD-workflow-validation.tsv"
LOG = MANUSCRIPT / "results" / "COPD-workflow-validation.log"

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
RESULT_ID_RE = re.compile(r"COPD-S[1-6]-R\d{3}")
FORBIDDEN_RE = re.compile(r"\b(?:TODO|TBD)\b|\u2014", re.IGNORECASE)
LOCAL_IMAGE_RE = re.compile(
    r"!\[[^\]]*\]\(\s*(?P<target><[^>]+>|[^\s)]+)(?:\s+[^)]*)?\s*\)"
)


def _raise_csv_limit() -> None:
    """Set the largest CSV field size supported by the current interpreter."""

    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


_raise_csv_limit()


def rel(path: Path) -> str:
    """Return a portable project-relative label when possible."""

    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except (OSError, ValueError):
        return str(path)


def clean_cell(value: object) -> str:
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Audit:
    def __init__(self, strict_final: bool, fast: bool) -> None:
        self.strict_final = strict_final
        self.fast = fast
        self.rows: list[dict[str, str]] = []
        self.hash_cache: dict[Path, str] = {}
        self._counter = 0

    def add(
        self,
        category: str,
        scope: str,
        status: str,
        observed: object,
        expected: object,
        details: object,
    ) -> None:
        if status not in {"PASS", "WARN", "FAIL"}:
            raise ValueError(f"invalid audit status: {status}")
        self._counter += 1
        self.rows.append(
            {
                "check_id": f"WFV-{self._counter:04d}",
                "category": clean_cell(category),
                "scope": clean_cell(scope),
                "status": status,
                "observed": clean_cell(observed),
                "expected": clean_cell(expected),
                "details": clean_cell(details),
            }
        )

    def planned_status(self) -> str:
        return "FAIL" if self.strict_final else "WARN"

    def hash(self, path: Path) -> Optional[str]:
        resolved = path.resolve()
        if resolved in self.hash_cache:
            return self.hash_cache[resolved]
        size = resolved.stat().st_size
        if self.fast and size >= 1_000_000_000:
            return None
        if size >= 1_000_000_000:
            print(
                f"Hashing declared large file ({size / 1_000_000_000:.2f} GB): {rel(path)}",
                file=sys.stderr,
                flush=True,
            )
        value = sha256_file(resolved)
        self.hash_cache[resolved] = value
        return value


EXPECTED_CORE: dict[str, list[str]] = {
    "Section 1": [
        "01_background/results/section_1_background.md",
        "01_background/results/epidemiology.tsv",
        "01_background/results/heritability.tsv",
        "01_background/results/genes_pathway_tissue.tsv",
        "01_background/results/cell_tissue_resources.tsv",
        "01_background/results/treatments.tsv",
    ],
    "Section 2": [
        "02_gwas/results/section_2_gwas.md",
        "02_gwas/results/COPD-S2-R001_studies_core.tsv",
        "02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv",
        "02_gwas/results/COPD-S2-R003A_genes_multistudy.tsv",
        "02_gwas/results/COPD-S2-R004_coding_noncoding_summary.tsv",
        "02_gwas/results/COPD-S2-R005_locus_table_gencode_v50.tsv",
        "02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.tsv.gz",
        "02_gwas/results/COPD-S2-R006F_summary.tsv",
        "02_gwas/results/COPD-S2-R007A_length_tests.tsv",
        "02_gwas/results/COPD-S2-R008E_summary.tsv",
    ],
    "Section 3": [
        "03_regulatory_landscape/results/section_3_regulatory_landscape.md",
        "03_regulatory_landscape/results/section_3_evidence_synthesis.md",
        "03_regulatory_landscape/results/COPD-S3-R001_regulatory_reference_inventory.tsv",
        "03_regulatory_landscape/results/COPD-S3-R002_gene_loci_grch38.tsv",
        "03_regulatory_landscape/results/COPD-S3-R002_lobe_union_summary.tsv",
        "03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz",
        "03_regulatory_landscape/results/COPD-S3-R004_expression_evidence.tsv",
        "03_regulatory_landscape/results/COPD-S3-R004_functional_variant_evidence.tsv",
        "03_regulatory_landscape/results/COPD-S3-R004_regulatory_element_evidence.tsv",
    ],
    "Section 4": [
        "04_modeling/results/section_4_modeling.md",
        "04_modeling/results/COPD-S4-R001_model_performance.tsv",
        "04_modeling/results/COPD-S4-R001_evaluation_manifest.json",
        "04_modeling/results/COPD-S4-R002_candidate_sequence_audit.tsv.gz",
        "04_modeling/results/COPD-S4-R002_candidate_sequence_manifest.json",
        "04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz",
        "04_modeling/results/COPD-S4-R003_scoring_manifest.json",
        "04_modeling/results/COPD-S4-R004_predicted_causal_regulatory_variants.tsv",
        "04_modeling/results/COPD-S4-R004_prioritized_candidates.tsv.gz",
        "04_modeling/results/COPD-S4-R004_analysis_manifest.json",
        "04_modeling/results/COPD-S4-R005_predicted_causal_population_genetics.tsv",
        "04_modeling/results/COPD-S4-R005_analysis_manifest.json",
        "04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv",
        "04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv",
        "04_modeling/results/COPD-S4-R006_analysis_manifest.json",
        "04_modeling/results/COPD-S4-R007_summary.tsv",
        "04_modeling/results/COPD-S4-R007_interpretation.md",
        "04_modeling/results/COPD-S4-R007_analysis_manifest.json",
        "04_modeling/results/COPD-S4-R008_candidate_tfbs_summary.tsv.gz",
        "04_modeling/results/COPD-S4-R008_motif_rankings.tsv",
        "04_modeling/results/COPD-S4-R008_tf_rankings.tsv",
        "04_modeling/results/COPD-S4-R008_methods.md",
        "04_modeling/results/COPD-S4-R008_analysis_manifest.json",
        "04_modeling/results/COPD-S4-R009_literature_functional_variant_recovery.tsv",
        "04_modeling/results/COPD-S4-R009_summary.tsv",
        "04_modeling/results/COPD-S4-R010_THE_LIST.tsv",
        "04_modeling/results/COPD-S4-R010_analysis_manifest.json",
    ],
    "Section 5": [
        "05_computational_validation/results/section_5_computational_validation.md",
        "05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_candidate_audit.tsv",
        "05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz",
        "05_computational_validation/results/COPD-S5-R001_manifest.json",
        "05_computational_validation/results/COPD-S5-R002_GRCh38_to_hg19_liftover_audit.tsv",
        "05_computational_validation/results/COPD-S5-R002_MPRAbase_candidate_audit.tsv",
        "05_computational_validation/results/COPD-S5-R002_MPRAbase_v4_9_3_candidate_elements.tsv.gz",
        "05_computational_validation/results/COPD-S5-R002_manifest.json",
        "05_computational_validation/results/COPD-S5-R003_HGNC_gene_resolution.tsv",
        "05_computational_validation/results/COPD-S5-R003_OpenTargets_COPD_gene_associations.tsv",
        "05_computational_validation/results/COPD-S5-R003_candidate_gene_context.tsv.gz",
        "05_computational_validation/results/COPD-S5-R003_manifest.json",
        "05_computational_validation/results/COPD-S5-R004_resource_access_audit.tsv",
        "05_computational_validation/results/COPD-S5-R004_biobank_access_audit.tsv",
        "05_computational_validation/results/COPD-S5-R004_manifest.json",
        "05_computational_validation/results/COPD-S5-R005_integrated_candidate_validation.tsv",
        "05_computational_validation/results/COPD-S5-R005_integrated_summary.tsv",
        "05_computational_validation/results/COPD-S5-R005_manifest.json",
        "05_computational_validation/results/COPD-S5_validation_checks.tsv",
    ],
}


EXPECTED_SECTION6 = [
    "06_experimental_validation/results/section_6_experimental_validation.md",
    "06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv",
    "06_experimental_validation/results/COPD-S6-R002_selection_audit.tsv.gz",
    "06_experimental_validation/results/COPD-S6-R003_MPRA_constructs.tsv",
    "06_experimental_validation/results/COPD-S6-R004_candidate_validation_plan.tsv",
    "06_experimental_validation/results/COPD-S6-R004_TF_first_plan.tsv",
    "06_experimental_validation/results/COPD-S6-R004_cell_context_controls.tsv",
    "06_experimental_validation/results/COPD-S6-R005_collaborator_README.md",
    "06_experimental_validation/results/COPD-S6-R005_collaborator_manifest.tsv",
    "06_experimental_validation/results/COPD-S6-R005_collaborator_package.tar.gz",
    "06_experimental_validation/results/COPD-S6-R005_package_manifest.json",
    "06_experimental_validation/results/COPD-S6-R006_analysis_manifest.json",
    "06_experimental_validation/results/COPD-S6_validation_checks.tsv",
    "06_experimental_validation/scripts/01_build_experimental_design.py",
    "06_experimental_validation/scripts/02_validate_design.py",
    "06_experimental_validation/scripts/03_build_collaborator_package.py",
    "06_experimental_validation/scripts/run_section6.py",
]


UNIQUE_TABLES: dict[str, str] = {
    "02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.tsv.gz": "candidate_record_id",
    "03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz": "candidate_record_id",
    "04_modeling/results/COPD-S4-R002_candidate_sequence_audit.tsv.gz": "candidate_record_id",
    "04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz": "candidate_record_id",
    "04_modeling/results/COPD-S4-R004_prioritized_candidates.tsv.gz": "candidate_record_id",
    "04_modeling/results/COPD-S4-R004_predicted_causal_regulatory_variants.tsv": "candidate_record_id",
    "04_modeling/results/COPD-S4-R005_population_annotated_candidates.tsv.gz": "candidate_record_id",
    "04_modeling/results/COPD-S4-R005_predicted_causal_population_genetics.tsv": "candidate_record_id",
    "04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv": "candidate_record_id",
    "04_modeling/results/COPD-S4-R010_THE_LIST.tsv": "candidate_record_id",
    "05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_candidate_audit.tsv": "candidate_record_id",
    "05_computational_validation/results/COPD-S5-R002_MPRAbase_candidate_audit.tsv": "candidate_record_id",
    "05_computational_validation/results/COPD-S5-R005_integrated_candidate_validation.tsv": "candidate_record_id",
}


CANONICAL_TABLES: list[tuple[str, str]] = [
    (
        "04_modeling/results/COPD-S4-R005_predicted_causal_population_genetics.tsv",
        "priority_rank_all_candidates",
    ),
    (
        "04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv",
        "predicted_causal_priority_rank",
    ),
    (
        "04_modeling/results/COPD-S4-R010_THE_LIST.tsv",
        "predicted_causal_priority_rank",
    ),
    (
        "05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_candidate_audit.tsv",
        "predicted_causal_priority_rank",
    ),
    (
        "05_computational_validation/results/COPD-S5-R002_MPRAbase_candidate_audit.tsv",
        "predicted_causal_priority_rank",
    ),
    (
        "05_computational_validation/results/COPD-S5-R005_integrated_candidate_validation.tsv",
        "predicted_causal_priority_rank",
    ),
]


def check_expected_files(audit: Audit) -> None:
    for section, paths in EXPECTED_CORE.items():
        missing = [path for path in paths if not (COPD / path).is_file()]
        audit.add(
            "expected_files",
            section,
            "FAIL" if missing else "PASS",
            f"{len(paths) - len(missing)}/{len(paths)} present",
            "all listed core outputs present",
            "missing: " + "; ".join(missing) if missing else "all core files found",
        )

    section6 = COPD / "06_experimental_validation"
    s6_results = sorted(path for path in (section6 / "results").glob("*") if path.is_file())
    missing = [path for path in EXPECTED_SECTION6 if not (COPD / path).is_file()]
    if not s6_results:
        status = audit.planned_status()
        details = "Section 6 is not yet populated"
    elif missing:
        status = "FAIL"
        details = "missing: " + "; ".join(missing)
    else:
        status = "PASS"
        details = "all final Section 6 files found"
    audit.add(
        "expected_files",
        "Section 6",
        status,
        f"{len(EXPECTED_SECTION6) - len(missing)}/{len(EXPECTED_SECTION6)} expected files present",
        "complete shortlist, construct, validation-plan, collaborator-package, manifest, and script set",
        details,
    )


def scan_one_tsv(path: Path) -> tuple[list[str], int, list[str]]:
    problems: list[str] = []
    rows = 0
    try:
        with open_text(path) as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            try:
                header = next(reader)
            except StopIteration:
                return [], 0, ["empty file"]
            if header:
                header[0] = header[0].lstrip("\ufeff")
            if not header or all(value == "" for value in header):
                problems.append("empty header")
            width = len(header)
            for line_number, row in enumerate(reader, start=2):
                rows += 1
                if len(row) != width:
                    problems.append(
                        f"record ending near physical line {line_number} has {len(row)} fields; expected {width}"
                    )
                    if len(problems) >= 10:
                        break
        return header, rows, problems
    except (OSError, UnicodeError, csv.Error) as exc:
        return [], rows, [f"read error: {type(exc).__name__}: {exc}"]


def check_tsv_rectangularity(audit: Audit) -> None:
    paths = sorted(COPD.rglob("*.tsv")) + sorted(COPD.rglob("*.tsv.gz"))
    paths = [path for path in paths if path.resolve() != OUTPUT.resolve()]
    failures: list[str] = []
    total_rows = 0
    for path in paths:
        _header, rows, problems = scan_one_tsv(path)
        total_rows += rows
        if problems:
            failures.append(f"{rel(path)} ({'; '.join(problems)})")
    audit.add(
        "tsv_rectangularity",
        "all COPD TSV/TSV.GZ files",
        "FAIL" if failures else "PASS",
        f"{len(paths) - len(failures)}/{len(paths)} valid; {total_rows} data rows scanned",
        "every existing TSV readable with one stable field count and unique headers",
        " | ".join(failures[:20]) if failures else "no malformed TSV records or headers found",
    )


def read_dict_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("missing header")
        rows = list(reader)
        return list(reader.fieldnames), rows


def check_unique_column(
    audit: Audit, path: Path, column: str, category: str = "id_uniqueness"
) -> None:
    if not path.is_file():
        audit.add(category, rel(path), "FAIL", "file missing", f"unique nonempty {column}", "")
        return
    try:
        header, rows = read_dict_rows(path)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        audit.add(category, rel(path), "FAIL", "unreadable", f"unique nonempty {column}", exc)
        return
    if column not in header:
        audit.add(category, rel(path), "FAIL", "column missing", column, "")
        return
    values = [row.get(column, "").strip() for row in rows]
    empty = sum(not value for value in values)
    counts = Counter(value for value in values if value)
    duplicates = sorted(value for value, count in counts.items() if count > 1)
    status = "FAIL" if empty or duplicates else "PASS"
    details = []
    if empty:
        details.append(f"{empty} empty IDs")
    if duplicates:
        details.append("duplicates: " + ", ".join(duplicates[:20]))
    audit.add(
        category,
        rel(path),
        status,
        f"{len(rows)} rows; {len(counts)} distinct nonempty IDs",
        f"{len(rows)} unique nonempty {column} values",
        "; ".join(details) if details else "identifier column is unique and complete",
    )


def check_ids(audit: Audit) -> None:
    register_specs = [
        (COPD / "sources.tsv", "source_id"),
        (COPD / "results_register.tsv", "result_id"),
        (COPD / "decisions.tsv", "decision_id"),
        (COPD / "activity_log.tsv", "action_id"),
    ]
    for path, column in register_specs:
        check_unique_column(audit, path, column)
    for relative, column in UNIQUE_TABLES.items():
        check_unique_column(audit, COPD / relative, column)


def parse_int(value: str) -> Optional[int]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    integer = int(number)
    return integer if number == integer else None


def candidate_rank_map(
    path: Path, rank_column: str
) -> tuple[list[str], dict[str, int], list[str]]:
    header, rows = read_dict_rows(path)
    problems: list[str] = []
    if "candidate_record_id" not in header:
        return [], {}, ["candidate_record_id missing"]
    if rank_column not in header:
        return [], {}, [f"{rank_column} missing"]
    ordered_ids: list[str] = []
    ranks: dict[str, int] = {}
    for row_number, row in enumerate(rows, start=2):
        candidate = row.get("candidate_record_id", "").strip()
        rank = parse_int(row.get(rank_column, ""))
        if not candidate:
            problems.append(f"row {row_number}: empty candidate_record_id")
            continue
        if rank is None:
            problems.append(f"row {row_number}: invalid {rank_column}")
            continue
        ordered_ids.append(candidate)
        if candidate in ranks:
            problems.append(f"duplicate candidate: {candidate}")
        ranks[candidate] = rank
    return ordered_ids, ranks, problems


def check_candidate_consistency(audit: Audit) -> None:
    canonical_path = (
        COPD / "04_modeling/results/COPD-S4-R004_predicted_causal_regulatory_variants.tsv"
    )
    try:
        canonical_order, canonical_ranks, canonical_problems = candidate_rank_map(
            canonical_path, "predicted_causal_priority_rank"
        )
    except Exception as exc:  # malformed input is also reported by TSV audit
        audit.add(
            "candidate_consistency",
            rel(canonical_path),
            "FAIL",
            "unreadable",
            "canonical candidate table",
            exc,
        )
        return
    expected_ranks = list(range(1, len(canonical_order) + 1))
    observed_ranks = [canonical_ranks.get(candidate) for candidate in canonical_order]
    if observed_ranks != expected_ranks:
        canonical_problems.append("canonical ranks are not contiguous 1..N in file order")
    audit.add(
        "candidate_consistency",
        rel(canonical_path),
        "FAIL" if canonical_problems else "PASS",
        f"{len(canonical_order)} rows; {len(canonical_ranks)} candidate IDs",
        "unique candidates ranked contiguously 1..N",
        "; ".join(canonical_problems[:20]) if canonical_problems else "canonical S4 candidate universe",
    )
    if canonical_problems:
        return

    canonical_set = set(canonical_order)
    for relative, rank_column in CANONICAL_TABLES:
        path = COPD / relative
        if not path.is_file():
            audit.add(
                "candidate_consistency",
                rel(path),
                "FAIL",
                "missing",
                f"same {len(canonical_set)} candidates and canonical ranks",
                "",
            )
            continue
        try:
            ids, ranks, problems = candidate_rank_map(path, rank_column)
        except Exception as exc:
            audit.add(
                "candidate_consistency", rel(path), "FAIL", "unreadable", "canonical set/ranks", exc
            )
            continue
        observed_set = set(ids)
        missing = sorted(canonical_set - observed_set)
        extra = sorted(observed_set - canonical_set)
        rank_mismatch = sorted(
            candidate
            for candidate in canonical_set & observed_set
            if ranks.get(candidate) != canonical_ranks[candidate]
        )
        order_mismatch = ids != canonical_order
        if missing:
            problems.append(f"missing {len(missing)} canonical candidates")
        if extra:
            problems.append(f"contains {len(extra)} noncanonical candidates")
        if rank_mismatch:
            problems.append(f"{len(rank_mismatch)} canonical-rank mismatches")
        if order_mismatch:
            problems.append("physical row order differs from canonical rank order")
        audit.add(
            "candidate_consistency",
            rel(path),
            "FAIL" if problems else "PASS",
            f"{len(ids)} rows; {len(observed_set)} IDs; {len(rank_mismatch)} rank mismatches",
            f"same {len(canonical_set)} IDs, ranks, and order as S4-R004",
            "; ".join(problems[:20]) if problems else "candidate count, membership, rank, and order agree",
        )

    check_section6_candidates(audit, canonical_ranks)


def check_section6_candidates(audit: Audit, canonical_ranks: dict[str, int]) -> None:
    results = COPD / "06_experimental_validation" / "results"
    candidate_tables: list[Path] = []
    for path in sorted(results.glob("*.tsv")) + sorted(results.glob("*.tsv.gz")):
        header, _rows, problems = scan_one_tsv(path)
        if not problems and "candidate_record_id" in header:
            candidate_tables.append(path)
    if not candidate_tables:
        audit.add(
            "candidate_consistency",
            "Section 6 candidate handoff",
            audit.planned_status(),
            "0 candidate-bearing result tables",
            "at least one Section 6 table keyed by candidate_record_id",
            "Section 6 is not yet populated" if not audit.strict_final else "final candidate handoff missing",
        )
        return

    canonical_set = set(canonical_ranks)
    for path in candidate_tables:
        header, rows = read_dict_rows(path)
        ids = [row.get("candidate_record_id", "").strip() for row in rows]
        unique_ids = {candidate for candidate in ids if candidate}
        extra = sorted(unique_ids - canonical_set)
        problems: list[str] = []
        if not unique_ids:
            problems.append("no nonempty candidate IDs")
        if extra:
            problems.append(f"{len(extra)} IDs absent from the S4 canonical set")
        canonical_rank_column = next(
            (
                column
                for column in ("predicted_causal_priority_rank", "priority_rank_all_candidates")
                if column in header
            ),
            None,
        )
        if canonical_rank_column:
            mismatches = 0
            invalid = 0
            for row in rows:
                candidate = row.get("candidate_record_id", "").strip()
                rank = parse_int(row.get(canonical_rank_column, ""))
                if rank is None:
                    invalid += 1
                elif candidate in canonical_ranks and rank != canonical_ranks[candidate]:
                    mismatches += 1
            if invalid:
                problems.append(f"{invalid} invalid canonical ranks")
            if mismatches:
                problems.append(f"{mismatches} canonical-rank mismatches")

        local_rank_columns = [
            column
            for column in header
            if "rank" in column.lower() and column != canonical_rank_column
        ]
        for column in local_rank_columns:
            values = [parse_int(row.get(column, "")) for row in rows]
            nonempty = [value for value in values if value is not None]
            if nonempty and len(set(nonempty)) == len(nonempty):
                if sorted(nonempty) != list(range(1, len(nonempty) + 1)):
                    problems.append(f"{column} is unique but not contiguous 1..N")
        audit.add(
            "candidate_consistency",
            rel(path),
            "FAIL" if problems else "PASS",
            f"{len(rows)} rows; {len(unique_ids)} canonical candidates represented",
            "candidate IDs are an S4 subset and carried canonical ranks agree",
            "; ".join(problems) if problems else "Section 6 candidate membership/ranks are coherent",
        )


def check_registered_outputs(audit: Audit) -> None:
    register = COPD / "results_register.tsv"
    if not register.is_file():
        audit.add(
            "traceability",
            rel(register),
            "FAIL",
            "missing",
            "readable result register",
            "",
        )
        return
    header, rows = read_dict_rows(register)
    required = {"result_id", "output_path", "status"}
    if not required.issubset(header):
        audit.add(
            "traceability",
            rel(register),
            "FAIL",
            "missing columns",
            ", ".join(sorted(required)),
            "",
        )
        return
    complete = [
        row
        for row in rows
        if row.get("status", "").strip().lower().startswith("complete")
    ]
    missing: list[str] = []
    for row in complete:
        raw = row.get("output_path", "").strip()
        path = Path(raw)
        if not path.is_absolute():
            path = ROOT / path
        if not path.is_file():
            missing.append(f"{row.get('result_id', '')}: {raw}")
    audit.add(
        "traceability",
        "registered complete outputs",
        "FAIL" if missing else "PASS",
        f"{len(complete) - len(missing)}/{len(complete)} registered outputs found",
        "every complete result-register path exists",
        "; ".join(missing[:20]) if missing else "all registered complete output paths resolve",
    )

    registered = {row.get("result_id", "").strip() for row in rows}
    discovered: set[str] = set()
    for directory in COPD.glob("0[1-6]_*/results"):
        for path in directory.glob("*"):
            match = RESULT_ID_RE.search(path.name)
            if match:
                discovered.add(match.group(0))
    unregistered = sorted(discovered - registered)
    audit.add(
        "traceability",
        "result-ID register coverage",
        audit.planned_status() if unregistered else "PASS",
        f"{len(discovered) - len(unregistered)}/{len(discovered)} discovered result IDs registered",
        "each result ID represented in the global result register by final handoff",
        "register handoff required for: " + ", ".join(unregistered)
        if unregistered
        else "all discovered result IDs registered",
    )


def _looks_like_path(value: object) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False
    value = value.strip()
    if value.startswith(("http://", "https://")):
        return False
    return (
        "/" in value
        or value.startswith(".")
        or bool(
            re.search(
                r"\.(?:tsv(?:\.gz)?|csv|json|md|py|sh|bed(?:\.gz)?|fa(?:sta)?|gz|h5|hdf5|keras|parquet|txt|bw|meme)$",
                value,
                re.IGNORECASE,
            )
        )
    )


def resolve_declared_path(
    raw: str, manifest: Path, basename_index: dict[str, list[Path]]
) -> Optional[Path]:
    raw = raw.strip()
    if not raw or raw.startswith(("http://", "https://")):
        return None
    candidate = Path(raw)
    candidates: list[Path] = []
    if candidate.is_absolute():
        candidates.append(candidate)
        markers = ["/disease-regulatory-genomics-workflow/", "/Writing_Scientific_Manuscript_Codex/"]
        for marker in markers:
            if marker in raw:
                suffix = raw.split(marker, 1)[1]
                candidates.append(ROOT / suffix)
    else:
        candidates.extend([ROOT / candidate, manifest.parent / candidate])
        for ancestor in manifest.parents:
            candidates.append(ancestor / candidate)
            if ancestor == ROOT:
                break
    for path in candidates:
        try:
            if path.is_file():
                return path
        except OSError:
            continue
    matches = basename_index.get(candidate.name, [])
    if len(matches) == 1:
        return matches[0]
    return None


def _preferred_path_value(mapping: dict[str, Any], base: Optional[str]) -> Optional[str]:
    keys: list[str] = []
    if base:
        trimmed = re.sub(r"_at_run$", "", base)
        keys.extend(
            [
                base,
                f"{base}_path",
                trimmed,
                f"{trimmed}_path",
            ]
        )
    keys.extend(["path", "output_path", "resolved_project_path", "file", "filename"])
    for key in keys:
        value = mapping.get(key)
        if _looks_like_path(value):
            return str(value)
    if base:
        normalized_base = re.sub(r"_(?:path|file|filename)$", "", base)
        for key, value in mapping.items():
            normalized_key = re.sub(r"_(?:path|file|filename)$", "", key)
            if normalized_key == normalized_base and _looks_like_path(value):
                return str(value)
        # A large flat key/value manifest may contain one unrelated path (most
        # often the script path).  Do not use that as a fallback for every
        # unpaired checksum declaration.
        return None
    path_values = [
        str(value)
        for key, value in mapping.items()
        if ("path" in key.lower() or key.lower() in {"file", "script", "input", "output"})
        and _looks_like_path(value)
    ]
    return path_values[0] if len(path_values) == 1 else None


ManifestRef = tuple[str, Optional[str], str, Optional[int]]


def json_manifest_refs(value: Any, prefix: str = "") -> list[ManifestRef]:
    refs: list[ManifestRef] = []
    if isinstance(value, dict):
        direct_sha = value.get("sha256")
        if isinstance(direct_sha, str) and SHA256_RE.fullmatch(direct_sha.strip()):
            raw = _preferred_path_value(value, None)
            bytes_value = value.get("bytes", value.get("size_bytes"))
            expected_bytes = int(bytes_value) if isinstance(bytes_value, int) else None
            refs.append((prefix or "sha256", raw, direct_sha.lower(), expected_bytes))
        for key, sha in value.items():
            match = re.match(r"^(.+)_sha256(?:_at_run)?$", key)
            if not match or not isinstance(sha, str) or not SHA256_RE.fullmatch(sha.strip()):
                continue
            base = match.group(1)
            raw = _preferred_path_value(value, base)
            expected_bytes = None
            for size_key in (f"{base}_bytes", f"{base}_size_bytes", "bytes", "size_bytes"):
                size_value = value.get(size_key)
                if isinstance(size_value, int):
                    expected_bytes = size_value
                    break
            label = f"{prefix}.{key}" if prefix else key
            refs.append((label, raw, sha.lower(), expected_bytes))
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else key
            refs.extend(json_manifest_refs(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            refs.extend(json_manifest_refs(child, f"{prefix}[{index}]"))
    return refs


def tsv_manifest_refs(
    path: Path,
    basename_index: dict[str, list[Path]],
    path_hints: dict[str, list[str]],
) -> list[ManifestRef]:
    header, rows = read_dict_rows(path)
    refs: list[ManifestRef] = []
    lower = {name.lower(): name for name in header}
    copd_relative_key = lower.get("path_relative_to_copd_root")
    path_key = lower.get("path", copd_relative_key)
    if "sha256" in lower and path_key is not None:
        sha_key = lower["sha256"]
        for index, row in enumerate(rows, start=2):
            sha = row.get(sha_key, "").strip()
            if not SHA256_RE.fullmatch(sha):
                continue
            size = None
            size_key = lower.get("bytes", lower.get("size_bytes"))
            if size_key is not None:
                raw_size = row.get(size_key, "").strip()
                size = int(raw_size) if raw_size.isdigit() else None
            raw_path = row.get(path_key, "").strip()
            if raw_path and path_key == copd_relative_key:
                raw_path = f"diseases/COPD/{raw_path}"
            refs.append((f"row {index}", raw_path or None, sha.lower(), size))
        return refs

    if set(lower) & {"key", "field"} and "value" in lower:
        key_name = lower.get("key", lower.get("field", ""))
        value_name = lower["value"]
        flat = {row.get(key_name, ""): row.get(value_name, "") for row in rows}
        for key, sha in flat.items():
            match = re.match(r"^(.+)_sha256(?:_at_run)?$", key)
            if not match or not SHA256_RE.fullmatch(sha.strip()):
                continue
            base = match.group(1)
            raw = _preferred_path_value(flat, base)
            hinted = path_hints.get(base, [])
            if raw is None and len(set(hinted)) == 1:
                raw = hinted[0]
            if raw is None and base.startswith("output_"):
                filename = base[len("output_") :]
                matches = basename_index.get(filename, [])
                if len(matches) == 1:
                    raw = str(matches[0])
            if raw is None:
                # Some early key/value manifests recorded a hash but omitted
                # its path. Resolve only an unambiguous filename/token match;
                # otherwise leave the declaration explicitly unresolvable.
                stem = re.sub(r"_(?:input|output|file)$", "", base)
                exact_matches: list[Path] = []
                for suffix in (".tsv", ".tsv.gz", ".bed", ".bed.gz"):
                    exact_matches.extend(basename_index.get(stem + suffix, []))
                if len(set(exact_matches)) == 1:
                    raw = str(exact_matches[0])
                elif not exact_matches:
                    tokens = [
                        token
                        for token in re.split(r"_+", stem.lower())
                        if token not in {"input", "output", "file"}
                    ]
                    token_matches = [
                        candidate
                        for paths in basename_index.values()
                        for candidate in paths
                        if tokens
                        and all(token in candidate.name.lower() for token in tokens)
                        and candidate.name.endswith((".tsv", ".tsv.gz", ".bed", ".bed.gz"))
                    ]
                    if len(set(token_matches)) == 1:
                        raw = str(token_matches[0])
            size = None
            for size_key in (f"{base}_bytes", f"{base}_size_bytes"):
                raw_size = flat.get(size_key, "")
                if raw_size.isdigit():
                    size = int(raw_size)
                    break
            refs.append((key, raw, sha.lower(), size))
    return refs


def check_manifests(audit: Audit) -> None:
    manifests = sorted(COPD.rglob("*manifest*.json")) + sorted(COPD.rglob("*manifest*.tsv"))
    basename_index: dict[str, list[Path]] = defaultdict(list)
    for path in COPD.rglob("*"):
        if path.is_file():
            basename_index[path.name].append(path)

    # Build cross-manifest hints before validation. Several Section 2 summary
    # manifests omit an input path that is declared by the same key in the
    # preceding execution manifest.
    path_hints: dict[str, list[str]] = defaultdict(list)
    for manifest in [path for path in manifests if path.suffix.lower() == ".tsv"]:
        try:
            header, rows = read_dict_rows(manifest)
        except (OSError, UnicodeError, csv.Error, ValueError):
            continue
        lowered = {column.lower(): column for column in header}
        if "value" not in lowered or not (set(lowered) & {"key", "field"}):
            continue
        key_column = lowered.get("key", lowered.get("field", ""))
        value_column = lowered["value"]
        for row in rows:
            key = row.get(key_column, "").strip()
            value = row.get(value_column, "").strip()
            if not _looks_like_path(value):
                continue
            normalized = re.sub(r"_path$", "", key)
            path_hints[normalized].append(value)

    for manifest in manifests:
        try:
            if manifest.suffix.lower() == ".json":
                with manifest.open("r", encoding="utf-8") as handle:
                    refs = json_manifest_refs(json.load(handle))
            else:
                refs = tsv_manifest_refs(manifest, basename_index, path_hints)
        except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, ValueError) as exc:
            audit.add(
                "manifest_hashes", rel(manifest), "FAIL", "manifest unreadable", "valid manifest", exc
            )
            continue

        deduplicated: list[ManifestRef] = []
        seen: set[tuple[Optional[str], str]] = set()
        for ref in refs:
            key = (ref[1], ref[2])
            if key not in seen:
                seen.add(key)
                deduplicated.append(ref)
        refs = deduplicated

        verified = 0
        unresolved: list[str] = []
        missing: list[str] = []
        mismatches: list[str] = []
        byte_mismatches: list[str] = []
        skipped: list[str] = []
        for label, raw, declared, expected_bytes in refs:
            if raw is None:
                unresolved.append(label)
                continue
            resolved = resolve_declared_path(raw, manifest, basename_index)
            if resolved is None:
                missing.append(f"{label} -> {raw}")
                continue
            if expected_bytes is not None and resolved.stat().st_size != expected_bytes:
                byte_mismatches.append(
                    f"{label}: {resolved.stat().st_size} bytes != {expected_bytes}"
                )
            observed = audit.hash(resolved)
            if observed is None:
                skipped.append(label)
            elif observed != declared:
                mismatches.append(f"{label} -> {rel(resolved)}")
            else:
                verified += 1

        hard = missing + mismatches + byte_mismatches
        if hard:
            status = "FAIL"
        elif unresolved or skipped:
            status = audit.planned_status() if audit.strict_final and unresolved else "WARN"
        else:
            status = "PASS"
        details: list[str] = []
        if missing:
            details.append("missing: " + ", ".join(missing[:8]))
        if mismatches:
            details.append("SHA mismatch: " + ", ".join(mismatches[:8]))
        if byte_mismatches:
            details.append("byte mismatch: " + ", ".join(byte_mismatches[:8]))
        if unresolved:
            details.append("unresolvable SHA declarations: " + ", ".join(unresolved[:12]))
        if skipped:
            details.append("large hashes skipped by --fast: " + ", ".join(skipped[:8]))
        if not details:
            details.append("all resolvable declared SHA-256 values agree")
        audit.add(
            "manifest_hashes",
            rel(manifest),
            status,
            f"{verified} verified; {len(unresolved)} unresolvable; {len(skipped)} skipped; {len(hard)} errors",
            "all path-resolvable declared hashes and byte counts agree",
            "; ".join(details),
        )


def choose_main_manuscript() -> tuple[Optional[Path], list[Path]]:
    direct = sorted(
        path
        for path in MANUSCRIPT.glob("*.md")
        if path.name.lower() != "readme.md" and path.stat().st_size > 0
    )
    if direct:
        candidates = direct
    else:
        candidates = sorted(
            path
            for path in MANUSCRIPT.rglob("*.md")
            if path.name.lower() != "readme.md"
            and "results" not in path.relative_to(MANUSCRIPT).parts
            and "scripts" not in path.relative_to(MANUSCRIPT).parts
            and path.stat().st_size > 0
        )
    if not candidates:
        return None, []
    preferred = [
        path
        for path in candidates
        if re.search(r"(?:manuscript|copd)", path.stem, re.IGNORECASE)
    ]
    pool = preferred or candidates
    main = max(pool, key=lambda path: path.stat().st_size)
    return main, candidates


def check_manuscript_files(audit: Audit) -> Optional[Path]:
    main, candidates = choose_main_manuscript()
    if main is None:
        audit.add(
            "manuscript_files",
            "single-file manuscript",
            audit.planned_status(),
            "not present",
            "one nonempty main Markdown manuscript",
            "manuscript assembly is still pending",
        )
    else:
        audit.add(
            "manuscript_files",
            "single-file manuscript",
            "FAIL" if len(candidates) != 1 else "PASS",
            f"{len(candidates)} candidate Markdown files; selected {rel(main)}",
            "exactly one main Markdown manuscript",
            "candidates: " + ", ".join(rel(path) for path in candidates),
        )

    maps = sorted(MANUSCRIPT.rglob("*.tsv"))
    maps = [
        path
        for path in maps
        if path.resolve() != OUTPUT.resolve()
        and "result" in path.name.lower()
        and ("map" in path.name.lower() or "manuscript" in path.name.lower())
    ]
    audit.add(
        "manuscript_files",
        "result-to-manuscript map",
        audit.planned_status() if not maps else "PASS",
        f"{len(maps)} matching map files",
        "at least one TSV result-to-manuscript map",
        ", ".join(rel(path) for path in maps) if maps else "readiness map not yet present",
    )
    return main


def check_figure_links(audit: Audit, manuscript: Optional[Path]) -> None:
    if manuscript is None:
        audit.add(
            "figure_links",
            "main manuscript",
            audit.planned_status(),
            "manuscript absent",
            "all embedded local figure links resolve",
            "deferred until manuscript assembly",
        )
        return
    text = manuscript.read_text(encoding="utf-8")
    targets = [match.group("target").strip("<>") for match in LOCAL_IMAGE_RE.finditer(text)]
    problems: list[str] = []
    for target in targets:
        decoded = unquote(target.split("#", 1)[0].split("?", 1)[0])
        if target.startswith(("http://", "https://", "data:")):
            problems.append(f"nonlocal figure link: {target}")
            continue
        target_path = Path(decoded)
        if target_path.is_absolute():
            problems.append(f"absolute figure link: {target}")
            continue
        resolved = manuscript.parent / target_path
        if not resolved.is_file():
            problems.append(f"missing figure: {target}")
    if not targets:
        status = audit.planned_status()
        problems.append("no embedded figures found")
    else:
        status = "FAIL" if problems else "PASS"
    audit.add(
        "figure_links",
        rel(manuscript),
        status,
        f"{len(targets)} image links; {len(problems)} problems",
        "one or more relative local figure links and every target exists",
        "; ".join(problems) if problems else "all embedded figure assets resolve",
    )

    callouts = [
        ("S" if match.group(1) else "", int(match.group(2)))
        for match in re.finditer(r"\bFig(?:ure)?\.?\s*(S?)(\d+)\b", text, re.IGNORECASE)
    ]
    main_numbers = []
    for prefix, number in callouts:
        if not prefix and number not in main_numbers:
            main_numbers.append(number)
    expected = list(range(1, max(main_numbers) + 1)) if main_numbers else []
    order_ok = main_numbers == expected
    audit.add(
        "figure_links",
        "figure callout numbering",
        "PASS" if order_ok else "FAIL",
        ",".join(str(value) for value in main_numbers) or "none",
        "first main-figure callouts appear as contiguous 1..N",
        "callout order is coherent" if order_ok else "out-of-order or skipped main-figure number",
    )


def expand_citation_group(group: str) -> list[int]:
    values: list[int] = []
    normalized = group.replace("\u2013", "-").replace("\u2212", "-")
    for token in re.split(r"\s*,\s*", normalized):
        if re.fullmatch(r"\d+", token):
            values.append(int(token))
        elif re.fullmatch(r"\d+\s*-\s*\d+", token):
            start, end = [int(value) for value in re.split(r"\s*-\s*", token)]
            if start <= end and end - start <= 1000:
                values.extend(range(start, end + 1))
    return values


def check_citations(audit: Audit, manuscript: Optional[Path]) -> None:
    if manuscript is None:
        audit.add(
            "citation_numbering",
            "main manuscript",
            audit.planned_status(),
            "manuscript absent",
            "bidirectional consecutive numbered citations",
            "deferred until manuscript assembly",
        )
        return
    text = manuscript.read_text(encoding="utf-8")
    heading = re.search(r"(?im)^#{1,6}\s+References\s*$", text)
    if not heading:
        audit.add(
            "citation_numbering",
            rel(manuscript),
            "FAIL",
            "References heading missing",
            "numbered reference list",
            "",
        )
        return
    body = text[: heading.start()]
    reference_text = text[heading.end() :]
    reference_numbers: list[int] = []
    for match in re.finditer(r"(?m)^\s*(?:\[(\d+)\]|(\d+)[.)])\s+", reference_text):
        reference_numbers.append(int(match.group(1) or match.group(2)))
    if not reference_numbers:
        audit.add(
            "citation_numbering",
            rel(manuscript),
            "FAIL",
            "0 numbered references",
            "numbered reference entries",
            "References section has no recognized numbered entries",
        )
        return
    max_reference = max(reference_numbers)
    cited: list[int] = []
    matches: list[tuple[int, str]] = []
    for match in re.finditer(r"\[([0-9][0-9,\s\-\u2013]*)\](?!\()", body):
        matches.append((match.start(), match.group(1)))
    for match in re.finditer(r"\(([0-9][0-9,\s\-\u2013]*)\)", body):
        matches.append((match.start(), match.group(1)))
    for _position, group in sorted(matches):
        for number in expand_citation_group(group):
            if number <= max_reference:
                cited.append(number)
    unique_first: list[int] = []
    for number in cited:
        if number not in unique_first:
            unique_first.append(number)

    ref_counts = Counter(reference_numbers)
    duplicate_refs = sorted(number for number, count in ref_counts.items() if count > 1)
    expected_refs = list(range(1, max_reference + 1))
    missing_entries = sorted(set(expected_refs) - set(reference_numbers))
    out_of_range = sorted(set(cited) - set(reference_numbers))
    uncited = sorted(set(reference_numbers) - set(cited))
    introduction_ok = unique_first == expected_refs
    problems: list[str] = []
    if duplicate_refs:
        problems.append("duplicate reference numbers: " + ",".join(map(str, duplicate_refs)))
    if missing_entries:
        problems.append("missing reference numbers: " + ",".join(map(str, missing_entries)))
    if out_of_range:
        problems.append("citations without entries: " + ",".join(map(str, out_of_range)))
    if uncited:
        problems.append("uncited reference entries: " + ",".join(map(str, uncited)))
    if not introduction_ok:
        problems.append(
            "first-appearance order is "
            + (",".join(map(str, unique_first[:30])) if unique_first else "empty")
        )
    audit.add(
        "citation_numbering",
        rel(manuscript),
        "FAIL" if problems else "PASS",
        f"{len(reference_numbers)} entries; {len(set(cited))} cited numbers",
        "unique consecutive entries, first cited in 1..N order, no orphan/uncited entries",
        "; ".join(problems) if problems else "bidirectional citation numbering is coherent",
    )


def check_forbidden_strings(audit: Audit) -> None:
    paths: set[Path] = set(COPD.rglob("*.md"))
    paths.update(MANUSCRIPT.rglob("*.tex"))
    paths.update(MANUSCRIPT.rglob("*.tsv"))
    excluded = {SCRIPT.resolve(), OUTPUT.resolve(), LOG.resolve()}
    hits: list[str] = []
    scanned = 0
    for path in sorted(paths):
        if path.resolve() in excluded or "scripts" in path.relative_to(COPD).parts:
            continue
        scanned += 1
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            hits.append(f"{rel(path)} unreadable: {exc}")
            continue
        for line_number, line in enumerate(text.splitlines(), start=1):
            matches = list(FORBIDDEN_RE.finditer(line))
            if matches:
                tokens = ",".join(match.group(0) for match in matches)
                hits.append(f"{rel(path)}:{line_number} [{tokens}]")
                if len(hits) >= 100:
                    break
    audit.add(
        "forbidden_strings",
        "COPD reports and manuscript deliverables",
        "FAIL" if hits else "PASS",
        f"{scanned} files scanned; {len(hits)} matching lines",
        "no TODO, TBD, or Unicode em dash",
        "; ".join(hits[:30]) if hits else "no forbidden strings found",
    )


def check_section_self_audits(audit: Audit) -> None:
    specs = [
        (
            COPD / "05_computational_validation/results/COPD-S5_validation_checks.tsv",
            "Section 5",
        ),
        (
            COPD / "06_experimental_validation/results/COPD-S6_validation_checks.tsv",
            "Section 6",
        ),
    ]
    for path, label in specs:
        if not path.is_file():
            status = audit.planned_status() if label == "Section 6" else "FAIL"
            audit.add(
                "section_self_audit", rel(path), status, "missing", f"all {label} checks PASS", ""
            )
            continue
        header, rows = read_dict_rows(path)
        if "status" not in header:
            audit.add(
                "section_self_audit", rel(path), "FAIL", "status column missing", "all checks PASS", ""
            )
            continue
        bad = [row for row in rows if row.get("status", "").strip().upper() != "PASS"]
        audit.add(
            "section_self_audit",
            rel(path),
            "FAIL" if bad else "PASS",
            f"{len(rows) - len(bad)}/{len(rows)} PASS",
            f"all {label} production checks PASS",
            "; ".join(
                row.get("check", row.get("check_name", row.get("check_id", "unnamed")))
                for row in bad[:20]
            )
            if bad
            else f"{label} self-audit is clean",
        )


def write_outputs(audit: Audit, started: datetime) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["check_id", "category", "scope", "status", "observed", "expected", "details"]
    temporary = OUTPUT.with_name(OUTPUT.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit.rows)
    os.replace(temporary, OUTPUT)

    counts = Counter(row["status"] for row in audit.rows)
    findings = [row for row in audit.rows if row["status"] != "PASS"]
    ended = datetime.now(timezone.utc)
    lines = [
        "COPD complete-workflow validation",
        f"started_utc: {started.isoformat()}",
        f"finished_utc: {ended.isoformat()}",
        f"mode: {'strict-final' if audit.strict_final else 'stage-aware'}",
        f"hash_mode: {'fast (>=1 GB skipped)' if audit.fast else 'full'}",
        f"checks: {len(audit.rows)}",
        f"PASS: {counts['PASS']}",
        f"WARN: {counts['WARN']}",
        f"FAIL: {counts['FAIL']}",
        f"table: {rel(OUTPUT)}",
    ]
    if findings:
        lines.append("findings:")
        for row in findings[:30]:
            lines.append(
                f"  {row['check_id']} {row['status']} {row['category']} | {row['scope']} | {row['details']}"
            )
        if len(findings) > 30:
            lines.append(f"  ... {len(findings) - 30} additional findings in the TSV")
    else:
        lines.append("findings: none")
    LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(
        f"COPD workflow validation: {counts['PASS']} PASS, {counts['WARN']} WARN, "
        f"{counts['FAIL']} FAIL ({len(audit.rows)} checks)"
    )
    print(f"Table: {rel(OUTPUT)}")
    print(f"Log:   {rel(LOG)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict-final",
        action="store_true",
        help="treat missing Section 6/manuscript deliverables and traceability handoffs as failures",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="skip recomputing declared hashes for individual files >=1 GB (reported as WARN)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now(timezone.utc)
    audit = Audit(strict_final=args.strict_final, fast=args.fast)
    check_expected_files(audit)
    check_tsv_rectangularity(audit)
    check_ids(audit)
    check_candidate_consistency(audit)
    check_registered_outputs(audit)
    check_manifests(audit)
    manuscript = check_manuscript_files(audit)
    check_figure_links(audit, manuscript)
    check_citations(audit, manuscript)
    check_forbidden_strings(audit)
    check_section_self_audits(audit)
    write_outputs(audit, started)
    return 1 if any(row["status"] == "FAIL" for row in audit.rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
