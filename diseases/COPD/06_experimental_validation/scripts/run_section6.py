#!/usr/bin/env python3
"""Run COPD Section 6 design, validation, and packaging in dependency order."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parent


def main() -> None:
    for script in [
        "01_build_experimental_design.py",
        "02_validate_design.py",
        "03_build_collaborator_package.py",
    ]:
        path = SCRIPTS / script
        print(f"Running {path}", flush=True)
        subprocess.run([sys.executable, str(path)], check=True)


if __name__ == "__main__":
    main()
