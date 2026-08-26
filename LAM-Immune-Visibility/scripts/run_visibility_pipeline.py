"""Run the isolated LAMCORE immune-visibility pipeline."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from common import PROJECT_ROOT, project_relative, write_json


STAGES = {
    "score": "score_visibility_modules.py",
    "summarize": "summarize_patient_states.py",
    "associations": "analyze_state_associations.py",
    "immune": "analyze_immune_context.py",
    "spatial": "analyze_spatial_axes.py",
    "retention": "bridge_sirolimus_retention.py",
    "ranking": "rank_candidate_antigens.py",
    "report": "build_reports.py",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["all", *STAGES], default="all")
    args = parser.parse_args()
    stages = list(STAGES) if args.stage == "all" else [args.stage]
    completed = []
    for stage in stages:
        script = PROJECT_ROOT / "scripts" / STAGES[stage]
        print(f"\n=== {stage}: {script.name} ===", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=PROJECT_ROOT, check=True)
        completed.append(stage)
    write_json(PROJECT_ROOT / "manifests" / "pipeline_run_manifest.json", {
        "stages": completed,
        "python": project_relative(sys.executable),
        "source_write_policy": "read-only source projects; all outputs under this project",
    })
    print("\nCompleted:", ", ".join(completed))


if __name__ == "__main__":
    main()
