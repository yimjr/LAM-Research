"""Resumable downloader that guards against servers ignoring HTTP Range requests."""

from __future__ import annotations

import argparse
from pathlib import Path

from urllib.error import HTTPError
from urllib.request import Request, urlopen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    args.destination.parent.mkdir(parents=True, exist_ok=True)
    offset = args.destination.stat().st_size if args.destination.exists() else 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    request = Request(args.url, headers=headers)
    try:
        response = urlopen(request, timeout=1800)
    except HTTPError as exc:
        raise SystemExit(f"HTTP {exc.code}: {exc.reason}") from exc
    with response:
        if offset and response.status != 206:
            # A 200 response means the server ignored Range. Restart rather than append.
            offset = 0
        mode = "ab" if offset else "wb"
        with args.destination.open(mode) as handle:
            while chunk := response.read(1024 * 1024):
                if chunk:
                    handle.write(chunk)


if __name__ == "__main__":
    main()
