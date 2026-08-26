"""Download configured public files and record hashes.

This intentionally downloads only explicitly selected files. Large human single-cell
and spatial archives remain staged until a modality-specific input is selected.
"""

from __future__ import annotations

import argparse
import datetime as dt
import time
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

from common import ROOT, sha256_file, write_json


def download(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return {"status": "existing", "path": str(destination.relative_to(ROOT)), "sha256": sha256_file(destination), "size": destination.stat().st_size}
    partial = destination.with_suffix(destination.suffix + ".part")
    last_error = None
    for attempt in range(4):
        try:
            with requests.get(url, stream=True, timeout=(120, 1800)) as response:
                response.raise_for_status()
                total = int(response.headers.get("content-length", "0"))
                with partial.open("wb") as handle, tqdm(total=total or None, unit="B", unit_scale=True, desc=destination.name) as progress:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                            progress.update(len(chunk))
            break
        except requests.RequestException as exc:
            last_error = exc
            if partial.exists():
                partial.unlink()
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    if last_error is not None and not partial.exists():
        raise last_error
    partial.replace(destination)
    return {"status": "downloaded", "path": str(destination.relative_to(ROOT)), "sha256": sha256_file(destination), "size": destination.stat().st_size}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", action="append", required=True)
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "config" / "datasets.yaml").read_text())
    selected = set(args.dataset)
    records = []
    for accession, spec in config["datasets"].items():
        if accession not in selected:
            continue
        for file_spec in spec.get("files", []):
            destination = ROOT / file_spec["path"]
            record = {"accession": accession, "url": file_spec["url"], "downloaded_at": dt.datetime.now(dt.timezone.utc).isoformat()}
            record.update(download(file_spec["url"], destination))
            records.append(record)
    if not records:
        raise SystemExit("No selected dataset has configured downloadable files")
    write_json(ROOT / "manifests" / "download_manifest.json", {"records": records})
    for record in records:
        print(record)


if __name__ == "__main__":
    main()
