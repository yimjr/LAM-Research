"""Convert the selected processed matrices from auxiliary GEO series."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))


def read_h5_sample(path: Path, sample: dict, accession: str) -> ad.AnnData:
    obj = sc.read_10x_h5(path, gex_only=True)
    obj.var_names = obj.var_names.astype(str)
    obj.var_names_make_unique()
    obj.obs_names = [f"{sample['sample_id']}:{name}" for name in obj.obs_names.astype(str)]
    obj.obs["sample_id"] = sample["sample_id"]
    obj.obs["donor_id"] = sample["donor_id"]
    obj.obs["tissue"] = "lung"
    obj.obs["assay"] = "scRNA"
    obj.obs["condition"] = "control"
    obj.obs["source_accession"] = accession
    obj.layers["counts"] = obj.X.copy()
    obj.var["gene_symbol"] = obj.var_names.astype(str)
    return obj


def read_dge(path: Path, sample: dict, accession: str) -> ad.AnnData:
    table = pd.read_csv(path, sep="\t", compression="gzip")
    if table.shape[1] < 2 or str(table.columns[0]).upper() != "GENE":
        raise ValueError(f"Unexpected DGE header: {path}")
    genes = table.iloc[:, 0].astype(str).to_numpy()
    barcodes = table.columns[1:].astype(str).to_numpy()
    values = table.iloc[:, 1:].to_numpy(dtype=np.int32, copy=True)
    matrix = sp.csr_matrix(values.T)
    obj = ad.AnnData(X=matrix, obs=pd.DataFrame(index=[f"{sample['sample_id']}:{x}" for x in barcodes]), var=pd.DataFrame(index=genes))
    obj.var_names_make_unique()
    obj.obs["sample_id"] = sample["sample_id"]
    obj.obs["donor_id"] = sample["donor_id"]
    obj.obs["tissue"] = "uterus"
    obj.obs["assay"] = "scRNA"
    obj.obs["condition"] = "control"
    obj.obs["source_accession"] = accession
    obj.var["gene_symbol"] = obj.var_names.astype(str)
    obj.layers["counts"] = obj.X.copy()
    return obj


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["GSE122960", "GSE118180", "all"], default="all")
    args = parser.parse_args()
    output_dir = ROOT / "data" / "processed" / "external"
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.dataset in ("GSE122960", "all"):
        spec = [
            ("GSM3489182", "GSE122960_Donor_01", "GSE122960_Donor_01"),
            ("GSM3489187", "GSE122960_Donor_03", "GSE122960_Donor_03"),
            ("GSM3489189", "GSE122960_Donor_04", "GSE122960_Donor_04"),
            ("GSM3489191", "GSE122960_Donor_05", "GSE122960_Donor_05"),
            ("GSM3489193", "GSE122960_Donor_06", "GSE122960_Donor_06"),
            ("GSM3489195", "GSE122960_Donor_07", "GSE122960_Donor_07"),
        ]
        objects = []
        for gsm, sample_id, donor_id in spec:
            file = sorted((ROOT / "data/raw/GSE122960/extracted").glob(f"{gsm}_*_filtered_gene_bc_matrices_h5.h5"))
            if len(file) != 1:
                raise FileNotFoundError(f"Expected one filtered H5 for {gsm}, found {file}")
            objects.append(read_h5_sample(file[0], {"sample_id": sample_id, "donor_id": donor_id}, "GSE122960"))
        combined = ad.concat(objects, join="outer", merge="same", index_unique=None)
        combined.write_h5ad(output_dir / "GSE122960_normal_lung.h5ad", compression="gzip")
        print("GSE122960", combined.shape)

    if args.dataset in ("GSE118180", "all"):
        file = ROOT / "data/raw/GSE118180/extracted/GSM3320143_WT_Uterus_out_gene_exon_tagged.dge.txt.gz"
        if not file.exists():
            raise FileNotFoundError(file)
        obj = read_dge(file, {"sample_id": "GSE118180_Wildtype_Uterus", "donor_id": "GSE118180_Wildtype_Uterus"}, "GSE118180")
        obj.write_h5ad(output_dir / "GSE118180_wildtype_uterus.h5ad", compression="gzip")
        print("GSE118180", obj.shape)


if __name__ == "__main__":
    main()
