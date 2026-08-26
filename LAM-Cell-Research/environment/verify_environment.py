"""Verify that the project-local Python environment is usable."""

from __future__ import annotations

import importlib
import json
import os
import platform
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))


REQUIRED = [
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "matplotlib",
    "seaborn",
    "scrublet",
    "harmonypy",
    "yaml",
    "requests",
    "jupyterlab",
]


def main() -> None:
    if sys.version_info[:2] != (3, 12):
        raise SystemExit(
            f"Expected Python 3.12.x in the project venv, got {platform.python_version()}"
        )

    versions: dict[str, str] = {}
    missing: list[str] = []
    for name in REQUIRED:
        try:
            module = importlib.import_module(name)
        except Exception as exc:  # pragma: no cover - diagnostic path
            missing.append(f"{name}: {exc}")
            continue
        versions[name] = getattr(module, "__version__", "unknown")

    if missing:
        raise SystemExit("Missing or broken dependencies:\n" + "\n".join(missing))

    print(
        json.dumps(
            {
                "python": platform.python_version(),
                "executable": sys.executable,
                "platform": platform.platform(),
                "packages": versions,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
