#!/usr/bin/env python3
"""Create and verify a deterministic COPD Section 6 collaborator package."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases/COPD/06_experimental_validation"
RESULTS = SECTION / "results"
DATA = SECTION / "data"
SCRIPTS = SECTION / "scripts"

MEMBER_MANIFEST = RESULTS / "COPD-S6-R005_collaborator_manifest.tsv"
PACKAGE = RESULTS / "COPD-S6-R005_collaborator_package.tar.gz"
PACKAGE_MANIFEST = RESULTS / "COPD-S6-R005_package_manifest.json"
ARCHIVE_ROOT = "COPD_Section6_handoff"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def payload() -> List[Tuple[Path, str]]:
    return [
        (SECTION / "README.md", "section overview and reproducible run command"),
        (RESULTS / "COPD-S6-R001_candidate_shortlist.tsv", "shortlist with exact upstream evidence"),
        (RESULTS / "COPD-S6-R002_selection_audit.tsv.gz", "all-candidate selection and LD-redundancy audit"),
        (RESULTS / "COPD-S6-R003_MPRA_constructs.tsv", "allele-specific MPRA insert designs"),
        (RESULTS / "COPD-S6-R004_candidate_validation_plan.tsv", "candidate-specific staged validation plan"),
        (RESULTS / "COPD-S6-R004_TF_first_plan.tsv", "candidate-by-model TF perturbation hypotheses"),
        (RESULTS / "COPD-S6-R004_cell_context_controls.tsv", "cell-context and control matrix"),
        (DATA / "COPD-S6-R005_upstream_inputs.tsv", "upstream provenance and checksums"),
        (RESULTS / "COPD-S6-R005_collaborator_README.md", "collaborator-facing interpretation and handoff instructions"),
        (RESULTS / "COPD-S6-R006_analysis_manifest.json", "analysis parameters, inputs, outputs, and guardrails"),
        (RESULTS / "section_6_experimental_validation.md", "integrated Section 6 report"),
        (RESULTS / "COPD-S6_validation_checks.tsv", "machine-readable validation checks"),
        (SCRIPTS / "01_build_experimental_design.py", "reproducible design builder"),
        (SCRIPTS / "02_validate_design.py", "reproducible design validator"),
        (SCRIPTS / "03_build_collaborator_package.py", "deterministic package builder"),
        (SCRIPTS / "run_section6.py", "Section 6 workflow runner"),
    ]


def archive_name(path: Path) -> str:
    return f"{ARCHIVE_ROOT}/{path.resolve().relative_to(SECTION.resolve())}"


def add_bytes(tar: tarfile.TarFile, data: bytes, name: str) -> None:
    info = tarfile.TarInfo(name=name)
    info.size = len(data)
    info.mtime = 0
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mode = 0o644
    tar.addfile(info, io.BytesIO(data))


def main() -> None:
    members = payload()
    missing = [relative(path) for path, _ in members if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing package payload:\n" + "\n".join(missing))

    manifest_rows: List[Dict[str, object]] = []
    for path, role in members:
        manifest_rows.append(
            {
                "archive_path": archive_name(path),
                "workspace_path": relative(path),
                "role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "included_in_package": True,
            }
        )
    manifest_frame = pd.DataFrame(manifest_rows).sort_values("archive_path")
    manifest_frame.to_csv(MEMBER_MANIFEST, sep="\t", index=False)
    manifest_bytes = MEMBER_MANIFEST.read_bytes()
    manifest_archive_name = f"{ARCHIVE_ROOT}/results/{MEMBER_MANIFEST.name}"

    with PACKAGE.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle:
            with tarfile.open(fileobj=gzip_handle, mode="w") as tar:
                for row in manifest_frame.to_dict(orient="records"):
                    path = ROOT / str(row["workspace_path"])
                    add_bytes(tar, path.read_bytes(), str(row["archive_path"]))
                add_bytes(tar, manifest_bytes, manifest_archive_name)

    expected_hashes = {
        str(row["archive_path"]): str(row["sha256"])
        for row in manifest_frame.to_dict(orient="records")
    }
    expected_hashes[manifest_archive_name] = sha256_bytes(manifest_bytes)
    observed_hashes: Dict[str, str] = {}
    with tarfile.open(PACKAGE, "r:gz") as tar:
        for member in tar.getmembers():
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            observed_hashes[member.name] = sha256_bytes(extracted.read())
    if observed_hashes != expected_hashes:
        missing_names = sorted(set(expected_hashes) - set(observed_hashes))
        extra_names = sorted(set(observed_hashes) - set(expected_hashes))
        mismatched = sorted(
            name
            for name in set(expected_hashes) & set(observed_hashes)
            if expected_hashes[name] != observed_hashes[name]
        )
        raise ValueError(
            "Package verification failed: "
            f"missing={missing_names}, extra={extra_names}, checksum_mismatch={mismatched}"
        )

    package_manifest = {
        "result_id": "COPD-S6-R005",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "collaborator_package_built_and_checksum_verified",
        "experimental_status": "all_experiments_proposed_not_performed",
        "archive": {
            "path": relative(PACKAGE),
            "bytes": PACKAGE.stat().st_size,
            "sha256": sha256(PACKAGE),
            "format": "gzip-compressed POSIX tar",
            "deterministic_metadata": "gzip and tar mtimes zero; uid/gid zero; fixed member order",
        },
        "member_manifest": {
            "path": relative(MEMBER_MANIFEST),
            "bytes": MEMBER_MANIFEST.stat().st_size,
            "sha256": sha256(MEMBER_MANIFEST),
        },
        "payload_file_count_excluding_member_manifest": len(manifest_frame),
        "archive_regular_file_count": len(observed_hashes),
        "archive_checksums_verified": True,
        "excluded_by_design": [
            "full paired 2001-bp upstream FASTAs; their paths and SHA-256 hashes are recorded in the upstream manifest",
            "MPRA vector, adapters, promoter, and barcode library; these require collaborator/platform decisions",
            "CRISPR guide and pegRNA sequences; feasibility and off-target design have not been performed",
            "wet-lab results; no Section 6 experiment has been performed",
        ],
    }
    PACKAGE_MANIFEST.write_text(json.dumps(package_manifest, indent=2, sort_keys=True) + "\n")

    print(f"Payload files: {len(manifest_frame)}")
    print(f"Archive regular files: {len(observed_hashes)}")
    print(f"Archive bytes: {PACKAGE.stat().st_size}")
    print(f"Archive SHA-256: {sha256(PACKAGE)}")
    print(f"Archive: {PACKAGE}")


if __name__ == "__main__":
    main()
