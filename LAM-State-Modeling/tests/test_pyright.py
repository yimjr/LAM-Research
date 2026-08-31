from __future__ import annotations

import json
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def test_pyright_targets_the_project_venv_and_excludes_large_artifacts():
    config = json.loads((PROJECT_DIR / "pyrightconfig.json").read_text(encoding="utf-8"))
    assert config["pythonVersion"] == "3.12"
    assert config["venvPath"] == "/mnt/py-env/venvs"
    assert config["venv"] == "LAM-State-Modeling"
    assert config["typeCheckingMode"] == "basic"
    assert {"data", "results", "reports"}.issubset(set(config["exclude"]))


def test_pyright_dev_dependencies_are_pinned_for_this_environment():
    lock = (PROJECT_DIR / "environment/requirements-dev.lock").read_text(encoding="utf-8")
    assert "pyright==1.1.411" in lock
    assert "nodejs-wheel-binaries==24.19.0" in lock
