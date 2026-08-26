"""Run the implemented core-reproduction, parallel robustness and discovery path."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    steps = [
        ["scripts/download_auxiliary_geo.py"],
        ["scripts/prepare_external_data.py"],
        ["scripts/build_reproduction_baseline.py"],
        ["scripts/reproduce_core.py"],
        ["scripts/run_targeted_robustness.py"],
        ["scripts/explore_lam_hypotheses.py"],
    ]
    for script in steps:
        print(f"[reproduction-plan] running {' '.join(script)}", flush=True)
        subprocess.run([sys.executable, *script], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
