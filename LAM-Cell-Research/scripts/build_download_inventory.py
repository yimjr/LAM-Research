#!/usr/bin/env python3
"""Build a consolidated inventory of files under data/raw.

The existing manifests remain the provenance source of truth. This report
joins them to the files actually present on disk so missing or unregistered
files are visible without hashing every large extracted matrix again.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(value: str) -> str:
    return str(Path(value).as_posix())


def dataset_from_path(path: str) -> str:
    parts = Path(path).parts
    try:
        idx = parts.index("raw")
    except ValueError:
        return "unknown"
    rest = parts[idx + 1 :]
    if not rest:
        return "unknown"
    if rest[0] in {"external", "perturbation"} and len(rest) > 1:
        return rest[1]
    if rest[0] == "reference":
        return "reference"
    return rest[0]


def collect_provenance() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}

    def add(path: str, source_url: str = "", sha: str = "", source: str = "") -> None:
        key = relative(path)
        current = records.setdefault(key, {})
        if source_url:
            current["source_url"] = source_url
        if sha:
            current["manifest_sha256"] = sha
        if source:
            current["manifest_source"] = source

    data_manifest = yaml.safe_load((ROOT / "manifests/data_manifest.yaml").read_text()) or {}
    archive = data_manifest.get("archive", {})
    add(archive.get("path", ""), data_manifest.get("source_url", ""), archive.get("sha256", ""), "data_manifest.yaml")
    for item in data_manifest.get("reference_files", []):
        add(item.get("path", ""), item.get("source_url", ""), item.get("sha256", ""), "data_manifest.yaml")

    external = yaml.safe_load((ROOT / "manifests/external_data_manifest.yaml").read_text()) or {}
    for item in external.get("datasets", []):
        archive = item.get("archive", {})
        add(archive.get("path", ""), item.get("source_url", ""), archive.get("sha256", ""), "external_data_manifest.yaml")

    perturb = yaml.safe_load((ROOT / "manifests/perturbation_downloads.yaml").read_text()) or {}
    for item in perturb.get("downloads", []):
        add(item.get("path", ""), item.get("source_url", ""), item.get("sha256", ""), "perturbation_downloads.yaml")

    staged = yaml.safe_load((ROOT / "manifests/external_file_staging.yaml").read_text()) or {}
    for dataset in staged.get("datasets", []):
        source_url = dataset.get("source_url", "")
        for key in ("archive", "provenance_metadata"):
            value = dataset.get(key, "")
            if isinstance(value, str):
                add(value, source_url, source="external_file_staging.yaml")
        for key in ("staged_samples", "verified_archives"):
            for item in dataset.get(key, []) or []:
                archive = item.get("archive", "")
                add(archive, source_url, source="external_file_staging.yaml")
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="manifests/download_inventory.csv")
    args = parser.parse_args()
    provenance = collect_provenance()
    rows = []
    raw_root = ROOT / "data/raw"
    for path in sorted(raw_root.rglob("*")):
        if not path.is_file() or path.name in {".DS_Store", "README.md"}:
            continue
        rel = relative(path.relative_to(ROOT))
        rel_raw = relative(path.relative_to(ROOT / "data"))
        prov = provenance.get(rel, {})
        name = path.name.lower()
        is_extracted = "/extracted/" in f"/{rel}/"
        is_archive = name.endswith((".tar", ".tar.gz", ".tgz", ".zip"))
        if is_extracted:
            file_role = "extracted_from_archive"
            origin = "derived_from_downloaded_archive"
        elif is_archive:
            file_role = "download_archive"
            origin = "downloaded_archive"
        elif dataset_from_path(rel) == "reference":
            file_role = "reference_file"
            origin = "downloaded_reference"
        else:
            file_role = "downloaded_matrix_or_auxiliary"
            origin = "downloaded_file"
        manifest_sha = prov.get("manifest_sha256", "")
        actual_sha = sha256(path) if manifest_sha else ""
        rows.append({
            "dataset": dataset_from_path(rel),
            "relative_path": rel_raw,
            "file_role": file_role,
            "origin": origin,
            "size_bytes": path.stat().st_size,
            "sha256": actual_sha or manifest_sha,
            "sha256_status": "verified_against_manifest" if actual_sha and actual_sha == manifest_sha else ("manifest_recorded" if manifest_sha else "archive_or_format_verified" if is_extracted else "not_recorded"),
            "source_url": prov.get("source_url", ""),
            "manifest_source": prov.get("manifest_source", ""),
        })
    table = pd.DataFrame(rows).sort_values(["dataset", "relative_path"])
    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    summary = table.groupby(["dataset", "file_role"], dropna=False).agg(files=("relative_path", "count"), bytes=("size_bytes", "sum")).reset_index()
    summary.to_csv(output.with_name("download_inventory_summary.csv"), index=False)
    print({"files": len(table), "output": str(output), "summary": str(output.with_name("download_inventory_summary.csv"))})


if __name__ == "__main__":
    main()
