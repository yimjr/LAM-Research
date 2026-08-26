"""Run the public-data pilot in the prescribed order."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(script: str) -> None:
    command = [sys.executable, str(ROOT / "scripts" / script)]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    run("download_geo.py")
    run("prepare_matrix.py")
    run("qc_and_preprocess.py")
    run("analyze_lam_states.py")
    run("robustness_tests.py")


if __name__ == "__main__":
    main()
