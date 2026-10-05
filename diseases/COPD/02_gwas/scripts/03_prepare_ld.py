#!/usr/bin/env python3
"""Prepare ancestry-matched inputs for COPD genome-wide-significant LD expansion.

Each tag variant is assigned only to 1000 Genomes super-populations represented
in the *initial* (discovery) cohorts of the GWAS Catalog association(s) in which
that tag reaches P <= 5e-8.  The script deliberately has no default population:
an ancestry that cannot be mapped to a 1000 Genomes super-population remains in
the audit table instead of silently being treated as European.

Stable result files from the LD run use the COPD-S2-R006 prefix.  This step
writes reproducible working inputs under data/ld_work and a preparation log.
All paths are derived from this file unless overridden on the command line.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


SCRIPT = Path(__file__).resolve()
SECTION_ROOT = SCRIPT.parents[1]
REPO_ROOT = SECTION_ROOT.parents[2]
DEFAULT_SHARED_DATA = REPO_ROOT / "data"
WORK = SECTION_ROOT / "data" / "ld_work"
LOGS = SECTION_ROOT / "logs"

PANELS = ("AFR", "AMR", "EAS", "EUR", "SAS")

# Ordered, explicit mappings.  Composite Catalog labels are handled because
# each phrase is searched independently in the complete label.
ANCESTRY_PHRASES: tuple[tuple[str, str], ...] = (
    ("African American or Afro-Caribbean", "AFR"),
    ("African unspecified", "AFR"),
    ("Hispanic or Latin American", "AMR"),
    ("East Asian", "EAS"),
    ("European", "EUR"),
    ("South Asian", "SAS"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gws-associations",
        type=Path,
        default=SECTION_ROOT / "results" / "COPD-S2-R002_gws_associations.tsv",
    )
    parser.add_argument(
        "--gws-tags",
        type=Path,
        default=SECTION_ROOT / "results" / "COPD-S2-R002_gws_unique_tag_variants.tsv",
    )
    parser.add_argument(
        "--ancestries",
        type=Path,
        default=SECTION_ROOT / "data" / "copd_ancestries.tsv",
    )
    parser.add_argument(
        "--sample-panel",
        type=Path,
        default=(
            DEFAULT_SHARED_DATA
            / "1000genomes"
            / "integrated_call_samples_v3.20130502.ALL.panel"
        ),
    )
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument(
        "--log", type=Path, default=LOGS / "COPD-S2-ld_prep.log"
    )
    return parser.parse_args()


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("COPD-S2-LD-prep")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ\t%(levelname)s\t%(message)s")
    formatter.converter = lambda *args: datetime.now(timezone.utc).timetuple()
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ancestry_panels(label: object) -> set[str]:
    if pd.isna(label):
        return set()
    text = str(label)
    return {panel for phrase, panel in ANCESTRY_PHRASES if phrase in text}


def main() -> int:
    args = parse_args()
    logger = setup_logging(args.log)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    associations = pd.read_csv(args.gws_associations, sep="\t", low_memory=False)
    tags = pd.read_csv(args.gws_tags, sep="\t", low_memory=False)
    ancestries = pd.read_csv(args.ancestries, sep="\t", low_memory=False)
    panel = pd.read_csv(args.sample_panel, sep="\t", low_memory=False)

    require_columns(
        associations,
        [
            "normalized_variant_id",
            "STUDY ACCESSION",
            "P",
            "genome_wide_significant",
            "phenotype_scope",
        ],
        "GWS associations",
    )
    require_columns(
        tags,
        [
            "normalized_variant_id",
            "chromosome",
            "position_grch38_catalog",
            "reported_risk_alleles",
            "genome_wide_significant",
            "phenotype_scope",
        ],
        "GWS tags",
    )
    require_columns(
        ancestries,
        [
            "STUDY ACCESSION",
            "STAGE",
            "BROAD ANCESTRAL CATEGORY",
            "copd_category",
        ],
        "ancestries",
    )
    require_columns(panel, ["sample", "super_pop"], "1000 Genomes panel")

    if not associations["phenotype_scope"].eq("core_copd").all():
        raise AssertionError("non-core phenotype leaked into GWS associations")
    if not tags["phenotype_scope"].eq("core_copd").all():
        raise AssertionError("non-core phenotype leaked into GWS tags")
    if not associations["genome_wide_significant"].astype(bool).all():
        raise AssertionError("non-significant row leaked into GWS associations")
    if not tags["genome_wide_significant"].astype(bool).all():
        raise AssertionError("non-significant row leaked into GWS tags")
    if tags["normalized_variant_id"].duplicated().any():
        raise AssertionError("GWS tag table is not unique by normalized_variant_id")

    discovery = ancestries.loc[
        ancestries["copd_category"].eq("core_copd")
        & ancestries["STAGE"].eq("initial")
    ].copy()
    discovery["mapped_panels"] = discovery["BROAD ANCESTRAL CATEGORY"].map(
        ancestry_panels
    )

    # The association table, not the all-significance tag aggregate, controls
    # ancestry assignment so a non-significant reuse of the same tag cannot add
    # an unrelated panel.
    variant_studies = (
        associations.groupby("normalized_variant_id")["STUDY ACCESSION"]
        .agg(lambda values: sorted(set(values.dropna().astype(str))))
        .to_dict()
    )
    study_to_rows: dict[str, pd.DataFrame] = {
        study: group.copy()
        for study, group in discovery.groupby("STUDY ACCESSION", sort=False)
    }

    focal_rows: list[dict[str, object]] = []
    assignment_rows: list[dict[str, object]] = []
    unmapped_rows: list[dict[str, object]] = []

    for row in tags.itertuples(index=False):
        variant = str(row.normalized_variant_id)
        studies = variant_studies.get(variant, [])
        support: dict[str, set[str]] = defaultdict(set)
        categories: dict[str, set[str]] = defaultdict(set)
        all_categories: set[str] = set()

        for study in studies:
            study_rows = study_to_rows.get(study)
            if study_rows is None or study_rows.empty:
                unmapped_rows.append(
                    {
                        "normalized_variant_id": variant,
                        "study_accession": study,
                        "ancestry_category": "",
                        "reason": "no_core_initial_ancestry_row",
                    }
                )
                continue
            # Use vectorized iteration over the named Series; this remains
            # stable for Catalog columns containing spaces.
            for category_value in study_rows["BROAD ANCESTRAL CATEGORY"]:
                category = "" if pd.isna(category_value) else str(category_value)
                all_categories.add(category)
                mapped = ancestry_panels(category)
                if not mapped:
                    unmapped_rows.append(
                        {
                            "normalized_variant_id": variant,
                            "study_accession": study,
                            "ancestry_category": category,
                            "reason": "no_1000g_superpopulation_mapping",
                        }
                    )
                for population in mapped:
                    support[population].add(study)
                    categories[population].add(category)

        supported = sorted(support)
        focal_rows.append(
            {
                "phenotype_scope": "core_copd",
                "normalized_variant_id": variant,
                "identifier_type": row.identifier_type,
                "chromosome": row.chromosome,
                "position_grch38_catalog": row.position_grch38_catalog,
                "reported_risk_alleles": row.reported_risk_alleles,
                "best_p": row.best_p,
                "gws_study_accessions": ";".join(studies),
                "initial_ancestry_categories": ";".join(sorted(all_categories)),
                "supported_panels": ";".join(supported),
                "n_supported_panels": len(supported),
                "ancestry_mapping_basis": (
                    "GWAS Catalog core-COPD initial-stage cohorts for GWS "
                    "association rows; no default panel"
                ),
            }
        )
        for population in supported:
            assignment_rows.append(
                {
                    "phenotype_scope": "core_copd",
                    "normalized_variant_id": variant,
                    "panel": population,
                    "supporting_studies": ";".join(sorted(support[population])),
                    "supporting_ancestry_categories": ";".join(
                        sorted(categories[population])
                    ),
                    "n_supporting_studies": len(support[population]),
                }
            )

    focal = pd.DataFrame(focal_rows).sort_values(
        ["chromosome", "position_grch38_catalog", "normalized_variant_id"],
        na_position="last",
        kind="stable",
    )
    assignments = pd.DataFrame(assignment_rows).sort_values(
        ["panel", "normalized_variant_id"], kind="stable"
    )
    unmapped = pd.DataFrame(
        unmapped_rows,
        columns=[
            "normalized_variant_id",
            "study_accession",
            "ancestry_category",
            "reason",
        ],
    ).drop_duplicates()

    focal.to_csv(args.work_dir / "focal_variants.tsv", sep="\t", index=False)
    assignments.to_csv(
        args.work_dir / "focal_panel_assignments.tsv", sep="\t", index=False
    )
    unmapped.to_csv(
        args.work_dir / "unmapped_ancestry_assignments.tsv", sep="\t", index=False
    )

    panel_counts: list[dict[str, object]] = []
    for population in PANELS:
        samples = panel.loc[panel["super_pop"].eq(population), "sample"].dropna()
        if samples.duplicated().any():
            raise AssertionError(f"duplicate sample IDs in {population} panel")
        samples.to_csv(
            args.work_dir / f"samples_{population}.txt", index=False, header=False
        )
        panel_counts.append(
            {
                "panel": population,
                "n_reference_samples": len(samples),
                "n_assigned_focal_variants": int(
                    assignments["panel"].eq(population).sum()
                ),
            }
        )
    pd.DataFrame(panel_counts).to_csv(
        args.work_dir / "panel_counts.tsv", sep="\t", index=False
    )

    no_panel = focal.loc[focal["n_supported_panels"].eq(0)]
    if not no_panel.empty:
        logger.warning("%d focal variants have no supported panel", len(no_panel))

    manifest_rows = [
        ("run_utc", datetime.now(timezone.utc).isoformat()),
        ("python", platform.python_version()),
        ("pandas", pd.__version__),
        ("script", str(SCRIPT)),
        ("script_sha256", sha256_file(SCRIPT)),
        ("gws_associations", str(args.gws_associations.resolve())),
        ("gws_associations_sha256", sha256_file(args.gws_associations)),
        ("gws_tags", str(args.gws_tags.resolve())),
        ("gws_tags_sha256", sha256_file(args.gws_tags)),
        ("ancestries", str(args.ancestries.resolve())),
        ("ancestries_sha256", sha256_file(args.ancestries)),
        ("sample_panel", str(args.sample_panel.resolve())),
        ("sample_panel_sha256", sha256_file(args.sample_panel)),
        ("n_gws_tags", len(focal)),
        ("n_focal_panel_assignments", len(assignments)),
        ("n_tags_without_supported_panel", len(no_panel)),
        ("population_order", ",".join(PANELS)),
        ("ancestry_stage", "initial"),
        ("genome_build", "GRCh38"),
    ]
    pd.DataFrame(manifest_rows, columns=["field", "value"]).to_csv(
        args.work_dir / "prep_manifest.tsv", sep="\t", index=False
    )

    logger.info("GWS focal tags: %d", len(focal))
    logger.info("Focal-panel assignments: %d", len(assignments))
    for item in panel_counts:
        logger.info(
            "%s: %d samples, %d assigned tags",
            item["panel"],
            item["n_reference_samples"],
            item["n_assigned_focal_variants"],
        )
    logger.info("Tags without supported panel: %d", len(no_panel))
    logger.info("Wrote LD preparation inputs to %s", args.work_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
