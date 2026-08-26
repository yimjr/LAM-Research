#!/usr/bin/env python3
"""Download public processed GEO archives for the discovery-validation phase.

This utility downloads only GEO supplementary archives. It refuses to extract
FASTQ/FQ members and writes a hash/size/source record before any downstream
conversion to AnnData. It does not claim that an archive is ready for analysis;
the archive inventory must still be inspected by the dataset-specific reader.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
URLS = {
    "GSE190260": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190260/suppl/GSE190260_RAW.tar",
    "GSE217108": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217108/suppl/GSE217108_RAW.tar",
    "GSE302356": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302356/suppl/GSE302356_RAW.tar",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, force: bool) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
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


def archive_inventory(path: Path) -> list[str]:
    if not tarfile.is_tarfile(path):
        raise RuntimeError(f"Downloaded file is not a valid tar archive: {path}")
    with tarfile.open(path, "r:*") as archive:
        members = archive.getnames()
    forbidden = [m for m in members if m.lower().endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
    if forbidden:
        raise RuntimeError("FASTQ/FQ members found; this project does not download raw sequencing data: " + ", ".join(forbidden[:5]))
    return members


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", nargs="+", choices=sorted(URLS))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    records = []
    for accession in args.accessions:
        raw_dir = ROOT / "data" / "raw" / "external" / accession
        raw_dir.mkdir(parents=True, exist_ok=True)
        url = URLS[accession]
        archive_path = raw_dir / f"{accession}_RAW.tar"
        download(url, archive_path, force=args.force)
        members = archive_inventory(archive_path)
        records.append(
            {
                "accession": accession,
                "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "download_url": url,
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                "archive_path": str(archive_path.relative_to(ROOT)),
                "size_bytes": archive_path.stat().st_size,
                "sha256": sha256(archive_path),
                "member_count": len(members),
                "fastq_members_rejected": True,
                "inventory_preview": members[:100],
                "status": "downloaded_inventory_checked_not_yet_converted",
            }
        )
    manifest = ROOT / "manifests" / "external_geo_downloads.yaml"
    manifest.write_text(yaml.safe_dump({"downloads": records}, sort_keys=False))
    print(json.dumps(records, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
