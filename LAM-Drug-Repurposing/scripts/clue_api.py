"""Small, auditable client for the CLUE/CMap L1000 API.

The project stores gene symbols, while the CLUE batch-query API expects Entrez
gene IDs.  This client therefore keeps three operations explicit:

1. ``prepare``: inspect the local signature table and write a query plan;
2. ``submit``: resolve symbols through CLUE's gene service and submit a JSON
   batch query to ``/api/jobs``;
3. ``poll`` / ``download``: retrieve job status and the completed archive.

Authentication is read only from ``CLUE_API_KEY`` (or the legacy-compatible
``CLUE_USER_KEY``) in the process environment.  The key is never written to
disk or printed.  Direction-reversal signatures are excluded unless explicitly
requested by the caller; the default signature file already follows the
project's scientific plan.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests

from common import ROOT, write_json


DEFAULT_BASE_URL = "https://api.clue.io/api"
DEFAULT_SIGNATURE = "results/signatures/GSE179044_cmap_query_signatures.csv"
DEFAULT_RESULT_DIR = "results/cmap"
DEFAULT_RAW_DIR = "data/raw/CLUE"
DEFAULT_CHUNK_SIZE = 200


class ClueApiError(RuntimeError):
    """An API or local-input error with a user-facing message."""


def _clean_base_url(value: str) -> str:
    return value.rstrip("/")


def _key_from_environment() -> str:
    key = os.environ.get("CLUE_API_KEY") or os.environ.get("CLUE_USER_KEY")
    if not key:
        raise ClueApiError(
            "No CLUE credential found. Set CLUE_API_KEY in the shell "
            "(do not paste it into project files or chat)."
        )
    return key.strip()


def _json_or_text(response: requests.Response) -> Any:
    try:
        return response.json()
    except ValueError:
        return response.text[:2000]


def _summarize_error(payload: Any) -> str:
    if isinstance(payload, dict):
        messages: list[str] = []
        for value in payload.values():
            if isinstance(value, dict):
                for field in ("errors", "warnings"):
                    entries = value.get(field, [])
                    if isinstance(entries, list):
                        messages.extend(
                            str(entry.get("text", entry)) if isinstance(entry, dict) else str(entry)
                            for entry in entries
                        )
        if messages:
            return "; ".join(messages[:8])
    return str(payload)[:1000]


class ClueClient:
    def __init__(self, api_key: str | None = None, base_url: str | None = None, timeout: int = 60):
        self.api_key = api_key or _key_from_environment()
        self.base_url = _clean_base_url(base_url or os.environ.get("CLUE_API_BASE_URL", DEFAULT_BASE_URL))
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json", "user_key": self.api_key})

    def request_json(self, method: str, path: str, **kwargs: Any) -> Any:
        url = f"{self.base_url}/{path.lstrip('/')}"
        try:
            response = self.session.request(method, url, timeout=self.timeout, **kwargs)
        except requests.RequestException as exc:
            raise ClueApiError(f"CLUE request failed ({method} {path}): {exc}") from exc
        payload = _json_or_text(response)
        if not response.ok:
            raise ClueApiError(
                f"CLUE returned HTTP {response.status_code} for {method} {path}: "
                f"{_summarize_error(payload)}"
            )
        return payload

    def metadata(self, endpoint: str) -> Any:
        return self.request_json("GET", endpoint)

    def resolve_symbols(self, symbols: Iterable[str], chunk_size: int = DEFAULT_CHUNK_SIZE) -> pd.DataFrame:
        unique = sorted({str(symbol).strip().upper() for symbol in symbols if str(symbol).strip()})
        rows: list[dict[str, Any]] = []
        for start in range(0, len(unique), chunk_size):
            chunk = unique[start : start + chunk_size]
            query_filter = {
                "where": {"gene_symbol": {"inq": chunk}},
                "fields": {"gene_id": 1, "gene_symbol": 1, "l1000_type": 1},
                "limit": len(chunk) * 2,
            }
            payload = self.request_json("GET", "genes", params={"filter": json.dumps(query_filter, separators=(",", ":"))})
            entries = payload.get("data", payload) if isinstance(payload, dict) else payload
            if isinstance(entries, dict):
                entries = entries.get("results", entries.get("rows", []))
            if not isinstance(entries, list):
                raise ClueApiError("Unexpected response shape from CLUE /genes endpoint")
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                symbol = str(entry.get("gene_symbol", "")).upper()
                entrez = entry.get("gene_id", entry.get("entrez_id", entry.get("gene_entrez_id")))
                if symbol and entrez not in (None, ""):
                    rows.append({
                        "gene_symbol": symbol,
                        "entrez_id": str(entrez),
                        "l1000_type": entry.get("l1000_type", ""),
                    })
        mapping = pd.DataFrame(rows, columns=["gene_symbol", "entrez_id", "l1000_type"])
        if mapping.empty:
            return mapping
        mapping = mapping.drop_duplicates(["gene_symbol", "entrez_id"])
        mapping["is_bing"] = mapping.l1000_type.astype(str).str.lower().isin({"landmark", "best inferred", "best_inferred", "bing"})
        return mapping.sort_values(["gene_symbol", "entrez_id"]).reset_index(drop=True)

    def submit_batch(self, name: str, up_gmt: str, down_gmt: str, tool_id: str = "sig_fastgutc_tool") -> Any:
        payload = {
            "tool_id": tool_id,
            "name": name,
            "uptag-cmapfile": up_gmt,
            "dntag-cmapfile": down_gmt,
            "data_type": "L1000",
            "dataset": "Touchstone",
            "ignoreWarnings": True,
        }
        return self.request_json("POST", "jobs", json=payload, headers={"Content-Type": "application/json"})

    def poll(self, job_id: str) -> Any:
        return self.request_json("GET", f"jobs/findByJobId/{job_id}")

    def download(self, url: str, destination: Path) -> None:
        if url.startswith("//"):
            url = "https:" + url
        elif url.startswith("/"):
            raise ClueApiError("Refusing to download an unexpected relative CLUE URL")
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            # The completed archive is normally on S3. Do not reuse the
            # authenticated CLUE session here, or the user_key could be sent
            # to a third-party download host.
            response = requests.get(url, headers={"Accept": "application/octet-stream"}, timeout=self.timeout, stream=True)
            response.raise_for_status()
            with destination.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        except requests.RequestException as exc:
            raise ClueApiError(f"CLUE result download failed: {exc}") from exc


def _signature_table(path: Path, include_reversal: bool = False) -> pd.DataFrame:
    if not path.exists():
        raise ClueApiError(f"Signature file does not exist: {path}")
    table = pd.read_csv(path)
    required = {"contrast", "gene", "direction"}
    missing = required - set(table.columns)
    if missing:
        raise ClueApiError(f"Signature file is missing columns: {sorted(missing)}")
    table = table.copy()
    table["gene"] = table["gene"].astype(str).str.upper().str.strip()
    if not include_reversal:
        excluded = {"direction_reversal", "reversal"}
        table = table.loc[~table.contrast.astype(str).str.lower().isin(excluded)]
    return table.loc[table.direction.isin(["up", "down"])].copy()


def _plan_from_symbols(table: pd.DataFrame, max_genes: int = 150) -> dict[str, Any]:
    contrasts: list[dict[str, Any]] = []
    for contrast, group in table.groupby("contrast", sort=True):
        up = sorted(group.loc[group.direction == "up", "gene"].dropna().unique())[:max_genes]
        down = sorted(group.loc[group.direction == "down", "gene"].dropna().unique())[:max_genes]
        contrasts.append({"contrast": str(contrast), "n_up": len(up), "n_down": len(down), "up_symbols": up, "down_symbols": down})
    return {"signature_file": str(DEFAULT_SIGNATURE), "max_genes_per_direction": max_genes, "contrasts": contrasts}


def _mapping_lookup(mapping: pd.DataFrame) -> dict[str, str]:
    if mapping.empty:
        return {}
    # Keep the first stable mapping; duplicate symbol→ID entries are retained in
    # the audit CSV but not duplicated in a query gene set.
    return mapping.drop_duplicates("gene_symbol").set_index("gene_symbol")["entrez_id"].to_dict()


def _gmt_rows(table: pd.DataFrame, mapping: pd.DataFrame, max_genes: int = 150) -> tuple[str, str, pd.DataFrame]:
    lookup = _mapping_lookup(mapping)
    up_rows: list[str] = []
    down_rows: list[str] = []
    audit_rows: list[dict[str, Any]] = []
    for contrast, group in table.groupby("contrast", sort=True):
        for direction, target in (("up", up_rows), ("down", down_rows)):
            genes = sorted(group.loc[group.direction == direction, "gene"].dropna().unique())[:max_genes]
            ids = []
            for symbol in genes:
                entrez = lookup.get(symbol)
                audit_rows.append({"contrast": contrast, "direction": direction, "gene_symbol": symbol, "entrez_id": entrez or "", "mapped": bool(entrez)})
                if entrez:
                    ids.append(str(entrez))
            if ids:
                target.append("\t".join([str(contrast), "", *ids]))
    return "\n".join(up_rows) + "\n", "\n".join(down_rows) + "\n", pd.DataFrame(audit_rows)


def _job_id(response: Any) -> str | None:
    if not isinstance(response, dict):
        return None
    candidates = [response.get("job_id")]
    result = response.get("result")
    if isinstance(result, dict):
        candidates.append(result.get("job_id"))
    for candidate in candidates:
        if candidate not in (None, ""):
            return str(candidate)
    return None


def command_prepare(args: argparse.Namespace) -> None:
    path = ROOT / args.signature
    table = _signature_table(path, include_reversal=args.include_reversal)
    plan = _plan_from_symbols(table, args.max_genes)
    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(out_dir / "clue_query_plan.json", plan)
    print(json.dumps({"status": "prepared", "contrasts": len(plan["contrasts"]), "out": str(out_dir / "clue_query_plan.json")}, ensure_ascii=False))


def command_doctor(args: argparse.Namespace) -> None:
    client = ClueClient(base_url=args.base_url)
    result = {}
    for endpoint in ("dataTypes", "datasets"):
        payload = client.metadata(endpoint)
        if isinstance(payload, list):
            result[endpoint] = {"status": "ok", "n_items": len(payload)}
        elif isinstance(payload, dict):
            result[endpoint] = {"status": "ok", "keys": sorted(payload.keys())[:20]}
        else:
            result[endpoint] = {"status": "ok", "type": type(payload).__name__}
    print(json.dumps(result, ensure_ascii=False))


def command_submit(args: argparse.Namespace) -> None:
    signature_path = ROOT / args.signature
    table = _signature_table(signature_path, include_reversal=args.include_reversal)
    client = ClueClient(base_url=args.base_url)
    symbols = table.gene.dropna().unique().tolist()
    mapping = client.resolve_symbols(symbols)
    up_gmt, down_gmt, audit = _gmt_rows(table, mapping, args.max_genes)
    if audit.empty:
        raise ClueApiError("No query genes were prepared")
    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(out_dir / "clue_gene_mapping.csv", index=False)
    audit.to_csv(out_dir / "clue_query_gene_audit.csv", index=False)
    (out_dir / "clue_uptag.gmt").write_text(up_gmt)
    (out_dir / "clue_dntag.gmt").write_text(down_gmt)
    mapped_fraction = float(audit.mapped.mean()) if len(audit) else 0.0
    if mapped_fraction < args.min_mapped_fraction:
        raise ClueApiError(
            f"Only {mapped_fraction:.1%} of signature genes mapped through CLUE; "
            "inspect clue_query_gene_audit.csv before submitting."
        )
    response = client.submit_batch(args.name, up_gmt, down_gmt)
    job_id = _job_id(response)
    record = {
        "status": "submitted" if job_id else "submitted_without_parsed_job_id",
        "job_id": job_id,
        "name": args.name,
        "signature": args.signature,
        "n_signature_genes": int(len(audit)),
        "mapped_fraction": mapped_fraction,
        "include_direction_reversal": bool(args.include_reversal),
        "api_base_url": client.base_url,
        "response": response,
    }
    filename = f"clue_job_{re.sub(r'[^A-Za-z0-9_.-]+', '_', job_id or 'submission')}.json"
    write_json(out_dir / filename, record)
    print(json.dumps({"status": record["status"], "job_id": job_id, "mapped_fraction": mapped_fraction, "record": str(out_dir / filename)}, ensure_ascii=False))


def command_poll(args: argparse.Namespace) -> None:
    client = ClueClient(base_url=args.base_url)
    payload = client.poll(args.job_id)
    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(out_dir / f"clue_job_{args.job_id}_status.json", payload)
    summary = {"job_id": args.job_id, "status": payload.get("status") if isinstance(payload, dict) else None, "download_status": payload.get("download_status") if isinstance(payload, dict) else None}
    print(json.dumps(summary, ensure_ascii=False))


def command_download(args: argparse.Namespace) -> None:
    client = ClueClient(base_url=args.base_url)
    payload = client.poll(args.job_id)
    url = payload.get("download_url") if isinstance(payload, dict) else None
    if not url:
        raise ClueApiError("No download_url in the current CLUE job response; poll again after completion")
    destination = ROOT / args.out_dir / f"clue_job_{args.job_id}.tar.gz"
    client.download(str(url), destination)
    print(json.dumps({"status": "downloaded", "job_id": args.job_id, "path": str(destination)}, ensure_ascii=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default=None, help=f"CLUE API base URL (default: {DEFAULT_BASE_URL})")
    sub = parser.add_subparsers(dest="command", required=True)

    prepare = sub.add_parser("prepare", help="write a local query plan without contacting CLUE")
    prepare.add_argument("--signature", default=DEFAULT_SIGNATURE)
    prepare.add_argument("--out-dir", default=DEFAULT_RESULT_DIR)
    prepare.add_argument("--max-genes", type=int, default=150)
    prepare.add_argument("--include-reversal", action="store_true")
    prepare.set_defaults(func=command_prepare)

    doctor = sub.add_parser("doctor", help="check authenticated CLUE metadata endpoints")
    doctor.set_defaults(func=command_doctor)

    submit = sub.add_parser("submit", help="resolve genes and submit a batch L1000 query")
    submit.add_argument("--signature", default=DEFAULT_SIGNATURE)
    submit.add_argument("--out-dir", default=DEFAULT_RESULT_DIR)
    submit.add_argument("--name", default="LAM residual and escape programs")
    submit.add_argument("--max-genes", type=int, default=150)
    submit.add_argument("--min-mapped-fraction", type=float, default=0.5)
    submit.add_argument("--include-reversal", action="store_true")
    submit.set_defaults(func=command_submit)

    poll = sub.add_parser("poll", help="poll a submitted CLUE job")
    poll.add_argument("job_id")
    poll.add_argument("--out-dir", default=DEFAULT_RESULT_DIR)
    poll.set_defaults(func=command_poll)

    download = sub.add_parser("download", help="download the completed CLUE job archive")
    download.add_argument("job_id")
    download.add_argument("--out-dir", default=DEFAULT_RAW_DIR)
    download.set_defaults(func=command_download)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except ClueApiError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
