"""Download and safely unpack the public processed GEO archive for GSE135851."""

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
    queue = sorted(root.rglob("*.tar")) + sorted(root.rglob("*.tar.gz")) + sorted(root.rglob("*.tgz"))
    seen: set[Path] = set()
    while queue:
        archive_path = queue.pop(0)
        if archive_path in seen or not archive_path.is_file():
            continue
        seen.add(archive_path)
        name = archive_path.name
        for suffix in (".tar.gz", ".tgz", ".tar"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        destination = archive_path.parent / name
        destination.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive_path, "r:*") as archive:
            members = safe_extract(archive, destination)
        unpacked.extend(f"{archive_path.name}:{member}" for member in members)
        queue.extend(
            sorted(destination.rglob("*.tar"))
            + sorted(destination.rglob("*.tar.gz"))
            + sorted(destination.rglob("*.tgz"))
        )
    return unpacked


def download(url: str, destination: Path, force: bool) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists():
        partial.unlink()
    with requests.get(url, stream=True, timeout=(30, 180)) as response:
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    accession = config["accession"]
    series_prefix = accession[:6] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}/{accession}/suppl/{accession}_RAW.tar"
    raw_dir = ROOT / config["paths"]["raw"] / accession
    raw_dir.mkdir(parents=True, exist_ok=True)
    archive_path = raw_dir / f"{accession}_RAW.tar"
    extracted_dir = raw_dir / "extracted"
    manifest_path = ROOT / config["paths"]["manifests"] / "data_manifest.yaml"

    download(url, archive_path, force=args.force)
    if not tarfile.is_tarfile(archive_path):
        raise RuntimeError(f"Downloaded file is not a valid tar archive: {archive_path}")

    with tarfile.open(archive_path, "r:*") as archive:
        member_names = archive.getnames()
    targets = config["targets"]
    missing = [
        target["gsm"]
        for target in targets
        if not any(target["gsm"] in member for member in member_names)
    ]
    if missing:
        raise RuntimeError(
            "Required target samples are missing from the GEO archive: " + ", ".join(missing)
        )

    if extracted_dir.exists() and args.force:
        shutil.rmtree(extracted_dir)
    extracted_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:*") as archive:
        top_members = safe_extract(archive, extracted_dir)
    nested_members = unpack_nested_archives(extracted_dir)

    manifest = {
        "status": "downloaded_and_extracted",
        "accession": accession,
        "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        "download_url": url,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "archive": {
            "path": str(archive_path.relative_to(ROOT)),
            "size_bytes": archive_path.stat().st_size,
            "sha256": sha256(archive_path),
        },
        "target_samples": [target["gsm"] for target in targets],
        "archive_members": top_members,
        "nested_members": nested_members,
        "extracted_dir": str(extracted_dir.relative_to(ROOT)),
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps(manifest["archive"], indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"download_geo.py failed: {exc}", file=sys.stderr)
        raise
