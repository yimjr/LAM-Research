"""Build a transcript-to-gene map for a kallisto expression table from a GTF."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


ATTR = re.compile(r'(gene_id|transcript_id|gene_name) "([^"]+)"')


def stable_id(value: str) -> str:
    return value.split(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("expression", type=Path)
    parser.add_argument("gtf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with gzip.open(args.expression, "rt") as handle:
        wanted = {stable_id(line.split(None, 1)[0]) for index, line in enumerate(handle) if index}

    mapped: dict[str, tuple[str, str]] = {}
    with gzip.open(args.gtf, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = dict(ATTR.findall(fields[8]))
            transcript = stable_id(attrs.get("transcript_id", ""))
            if transcript in wanted:
                mapped[transcript] = (stable_id(attrs["gene_id"]), attrs.get("gene_name", ""))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target_id", "gene_id", "gene_symbol"])
        for transcript in sorted(wanted):
            gene_id, gene_symbol = mapped.get(transcript, ("", ""))
            writer.writerow([transcript, gene_id, gene_symbol])
    print({"requested": len(wanted), "mapped": len(mapped), "unmapped": len(wanted) - len(mapped)})


if __name__ == "__main__":
    main()
