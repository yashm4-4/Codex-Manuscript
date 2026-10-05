#!/usr/bin/env python3
"""Run the reproducible COPD Section 5 pipeline in dependency order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parent


def run(script: str) -> None:
    command = [sys.executable, str(SCRIPTS / script)]
    print(f"Running {' '.join(command)}", flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="write only the public/controlled-access audit (does not require Section 4)",
    )
    args = parser.parse_args()
    run("04_access_audit.py")
    if args.audit_only:
        return
    for script in (
        "01_gtex_lung_eqtl.py",
        "02_mprabase_validation.py",
        "03_gene_catalog_evidence.py",
        "05_integrate_validation.py",
        "99_validate_section5.py",
    ):
        run(script)


if __name__ == "__main__":
    main()

