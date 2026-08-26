"""Convert GEO 10x Matrix Market triplets into a combined AnnData object."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml


ROOT = Path(__file__).resolve().parents[1]


def open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else path.open("rt")


def read_vector(path: Path) -> list[str]:
    with open_text(path) as handle:
        values = [line.rstrip("\n\r").split("\t") for line in handle if line.strip()]
    if not values:
        raise RuntimeError(f"Empty annotation file: {path}")
    return [row[0] for row in values]


def read_genes(path: Path) -> tuple[list[str], list[str]]:
    with open_text(path) as handle:
        rows = [line.rstrip("\n\r").split("\t") for line in handle if line.strip()]
    if not rows:
        raise RuntimeError(f"Empty gene annotation file: {path}")
    ids = [row[0] for row in rows]
    symbols = [row[1] if len(row) > 1 and row[1] else row[0] for row in rows]
    return ids, symbols


def candidate_files(root: Path, tokens: list[str], kind: str) -> list[Path]:
    all_files = [path for path in root.rglob("*") if path.is_file()]
    token_files = [path for path in all_files if any(token.lower() in str(path).lower() for token in tokens)]
    if kind == "matrix":
        return [path for path in token_files if re.search(r"matrix.*mtx", path.name, re.I)]
    if kind == "barcodes":
        return [path for path in token_files if "barcode" in path.name.lower()]
    return [path for path in token_files if "gene" in path.name.lower() or "feature" in path.name.lower()]


def closest_file(matrix: Path, candidates: list[Path]) -> Path | None:
    if not candidates:
        return None
    same_parent = [path for path in candidates if path.parent == matrix.parent]
    return sorted(same_parent or candidates, key=lambda path: (len(path.parts), str(path)))[0]


def load_sample(root: Path, target: dict) -> ad.AnnData:
    tokens = [target["gsm"], target["sample_id"]]
    matrix = sorted(candidate_files(root, tokens, "matrix"), key=lambda path: (len(path.parts), str(path)))
    if not matrix:
        raise RuntimeError(f"No matrix.mtx file found for {target['sample_id']} ({target['gsm']})")
    matrix_path = matrix[0]
    barcode_path = closest_file(matrix_path, candidate_files(root, tokens, "barcodes"))
    gene_path = closest_file(matrix_path, candidate_files(root, tokens, "genes"))
    if barcode_path is None or gene_path is None:
        raise RuntimeError(
            f"Incomplete 10x triplet for {target['sample_id']}: "
            f"matrix={matrix_path}, barcodes={barcode_path}, genes={gene_path}"
        )

    with gzip.open(matrix_path, "rb") if matrix_path.name.endswith(".gz") else matrix_path.open("rb") as handle:
        matrix_data = sp.csr_matrix(scipy.io.mmread(handle))
    barcodes = read_vector(barcode_path)
    gene_ids, gene_symbols = read_genes(gene_path)
    if matrix_data.shape == (len(gene_ids), len(barcodes)):
        matrix_data = matrix_data.T.tocsr()
    elif matrix_data.shape != (len(barcodes), len(gene_ids)):
        raise RuntimeError(
            f"Matrix shape mismatch for {target['sample_id']}: {matrix_data.shape}; "
            f"genes={len(gene_ids)}, barcodes={len(barcodes)}"
        )

    cell_ids = [f"{target['sample_id']}:{barcode}" for barcode in barcodes]
    obs = pd.DataFrame(
        {
            "sample_id": target["sample_id"],
            "donor_id": target["donor_id"],
            "tissue": target["tissue"],
            "condition": target["condition"],
            "assay": target["assay"],
            "source_accession": target["gsm"],
        },
        index=cell_ids,
    )
    var = pd.DataFrame(
        {"gene_id": gene_ids, "gene_symbol": gene_symbols},
        index=pd.Index(gene_symbols, dtype=str),
    )
    var.index = var.index.astype(str)
    var.index = pd.Index(var.index).map(lambda name: name if name else "unknown_gene")
    var.index = pd.Index(var.index).astype(str)
    result = ad.AnnData(X=matrix_data, obs=obs, var=var)
    result.var_names = result.var_names.astype(str)
    result.obs_names = result.obs_names.astype(str)
    result.var_names_make_unique()
    result.uns["source_files"] = {
        "matrix": str(matrix_path.relative_to(root)),
        "barcodes": str(barcode_path.relative_to(root)),
        "genes": str(gene_path.relative_to(root)),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis.yaml")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / args.config).read_text())
    accession = config["accession"]
    extracted = ROOT / config["paths"]["raw"] / accession / "extracted"
    if not extracted.exists():
        raise FileNotFoundError(f"Missing extracted GEO data: {extracted}; run download_geo.py first")

    objects = [load_sample(extracted, target) for target in config["targets"]]
    combined = ad.concat(objects, join="outer", merge="same", label="source_sample", index_unique=None)
    combined.obs_names_make_unique()
    combined.var_names_make_unique()
    combined.uns["data_contract"] = {
        "counts_layer": "counts will be created during preprocessing",
        "source_accession": accession,
        "target_sample_count": len(objects),
    }
    output_dir = ROOT / config["paths"]["interim"]
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{accession}_combined_raw.h5ad"
    combined.write_h5ad(output, compression="gzip")
    print(f"Wrote {output}: {combined.n_obs} cells x {combined.n_vars} genes")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"prepare_matrix.py failed: {exc}", file=sys.stderr)
        raise
