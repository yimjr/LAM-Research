"""Download and unpack processed GEO archives used by the reproduction plan."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tarfile
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


def safe_extract(archive: tarfile.TarFile, destination: Path) -> list[str]:
    destination = destination.resolve()
    extracted: list[str] = []
    for member in archive.getmembers():
        target = (destination / member.name).resolve()
        try:
            target.relative_to(destination)
        except ValueError as exc:
            raise RuntimeError(f"Unsafe archive member: {member.name}") from exc
        archive.extract(member, destination)
        extracted.append(member.name)
    return extracted


def unpack_nested_archives(root: Path) -> list[str]:
    unpacked: list[str] = []
    queue = sorted(root.rglob("*.tar")) + sorted(root.rglob("*.tar.gz"))
    seen: set[Path] = set()
    while queue:
        archive_path = queue.pop(0)
        if archive_path in seen or not archive_path.is_file():
            continue
        seen.add(archive_path)
        name = archive_path.name
        for suffix in (".tar.gz", ".tar"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        destination = archive_path.parent / name
        destination.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive_path, "r:*") as archive:
            members = safe_extract(archive, destination)
        unpacked.extend(f"{archive_path.name}:{member}" for member in members)
        queue.extend(sorted(destination.rglob("*.tar")) + sorted(destination.rglob("*.tar.gz")))
    return unpacked


def download(url: str, destination: Path, force: bool) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists():
        partial.unlink()
    with requests.get(url, stream=True, timeout=(30, 300)) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", "0"))
        with partial.open("wb") as handle, tqdm(
            total=total or None,
            unit="B",
            unit_scale=True,
            desc=destination.name,
        ) as progress:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    progress.update(len(chunk))
    partial.replace(destination)


def process_dataset(spec: dict, force: bool) -> dict:
    accession = str(spec["accession"])
    raw_dir = ROOT / "data" / "raw" / accession
    raw_dir.mkdir(parents=True, exist_ok=True)
    archive_path = raw_dir / f"{accession}_RAW.tar"
    extracted_dir = raw_dir / "extracted"
    download(str(spec["download_url"]), archive_path, force)
    if not tarfile.is_tarfile(archive_path):
        raise RuntimeError(f"Downloaded file is not a valid tar archive: {archive_path}")
    with tarfile.open(archive_path, "r:*") as archive:
        members = archive.getnames()
    requested = [str(sample["gsm"]) for sample in spec.get("include_samples", [])]
    missing = [gsm for gsm in requested if not any(gsm in member for member in members)]
    if missing:
        raise RuntimeError(f"Required samples missing from {accession} archive: {', '.join(missing)}")
    if extracted_dir.exists() and force:
        shutil.rmtree(extracted_dir)
    extracted_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:*") as archive:
        top_members = safe_extract(archive, extracted_dir)
    nested_members = unpack_nested_archives(extracted_dir)
    return {
        "accession": accession,
        "role": spec["role"],
        "source_url": spec["source_url"],
        "download_url": spec["download_url"],
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "archive": {
            "path": str(archive_path.relative_to(ROOT)),
            "size_bytes": archive_path.stat().st_size,
            "sha256": sha256(archive_path),
        },
        "requested_samples": requested,
        "archive_members": top_members,
        "nested_members": nested_members,
        "extracted_dir": str(extracted_dir.relative_to(ROOT)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", action="append", help="Only process this accession; repeatable")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "config" / "external_datasets.yaml").read_text())
    selected = set(args.dataset or [])
    specs = [spec for spec in config["external_datasets"] if not selected or spec["accession"] in selected]
    if not specs:
        raise ValueError("No external dataset matched --dataset")
    manifests = [process_dataset(spec, args.force) for spec in specs]
    output = ROOT / "manifests" / "external_data_manifest.yaml"
    previous = yaml.safe_load(output.read_text()) if output.exists() else {}
    previous_by_accession = {
        str(dataset["accession"]): dataset for dataset in previous.get("datasets", [])
    }
    previous_by_accession.update({str(dataset["accession"]): dataset for dataset in manifests})
    merged = [previous_by_accession[key] for key in sorted(previous_by_accession)]
    output.write_text(yaml.safe_dump({"status": "downloaded_and_extracted", "datasets": merged}, sort_keys=False))
    print(json.dumps(manifests, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"download_auxiliary_geo.py failed: {exc}", file=sys.stderr)
        raise
