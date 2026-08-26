#!/usr/bin/env python3
"""Download and verify public processed perturbation matrices only."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, force: bool) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists():
        partial.unlink()
    with requests.get(url, stream=True, timeout=(30, 180)) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", "0"))
        with partial.open("wb") as handle, tqdm(total=total or None, unit="B", unit_scale=True, desc=destination.name) as progress:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    progress.update(len(chunk))
    partial.replace(destination)


def verify_gzip_text(path: Path) -> tuple[bool, str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        first_line = handle.readline()
    if not first_line.strip():
        raise RuntimeError(f"Compressed text file has no readable header: {path}")
    return True, first_line[:200].rstrip("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs="*", default=["GSE179044", "GSE84476"])
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--config", default="config/perturbation_datasets.yaml")
    args = parser.parse_args()
    cfg = yaml.safe_load((ROOT / args.config).read_text())
    records = []
    for dataset in args.datasets:
        if dataset not in cfg["datasets"]:
            raise ValueError(f"Unknown perturbation dataset: {dataset}")
        for item in cfg["datasets"][dataset]["processed_files"]:
            destination = ROOT / item["destination"]
            download(item["url"], destination, args.force)
            valid, header = verify_gzip_text(destination)
            records.append({
                "dataset": dataset,
                "label": item["label"],
                "source_url": cfg["datasets"][dataset]["source_url"],
                "download_url": item["url"],
                "path": str(destination.relative_to(ROOT)),
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                "size_bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                "gzip_text_verified": valid,
                "header_preview": header,
                "fastq_downloaded": False,
            })
    manifest = ROOT / "manifests" / "perturbation_downloads.yaml"
    manifest.write_text(yaml.safe_dump({"downloads": records}, sort_keys=False))
    print(json.dumps(records, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
