#!/usr/bin/env python3
"""Convert staged external processed matrices into analysis-ready AnnData.

RNA/GEX is converted here. ATAC, spatial, and Xenium files remain separate
modalities and are registered in manifests/external_modalities.yaml; they are
not forced into an RNA expression matrix.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/tmp/lammpl")

MARKER_GENES = ["PMEL", "MLANA", "MITF", "ACTA2", "ESR1", "FIGF", "VEGFD", "CTSK"]


def read_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def read_matrix(path: Path) -> ad.AnnData:
    if path.suffix.lower() == ".h5":
        obj = sc.read_10x_h5(path, gex_only=True)
    else:
        obj = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True, cache=False)
    obj.var_names = obj.var_names.astype(str)
    original_symbols = obj.var_names.astype(str).to_numpy()
    obj.var["gene_symbol"] = original_symbols
    obj.var_names_make_unique()
    obj.var["gene_symbol_upper"] = obj.var["gene_symbol"].astype(str).str.upper().to_numpy()
    return obj


def column_values(obj: ad.AnnData, gene: str) -> np.ndarray:
    symbols = obj.var["gene_symbol_upper"].astype(str).to_numpy()
    hits = np.flatnonzero(symbols == gene.upper())
    if len(hits) == 0:
        return np.zeros(obj.n_obs, dtype=np.float32)
    values = obj.X[:, int(hits[0])]
    if sp.issparse(values):
        values = values.toarray().ravel()
    return np.asarray(values, dtype=np.float32).ravel()


def annotate_and_normalize(obj: ad.AnnData, sample: dict, accession: str, cfg: dict) -> ad.AnnData:
    obj.obs_names = [f"{sample['sample_id']}:{x}" for x in obj.obs_names.astype(str)]
    obj.layers["counts"] = obj.X.copy()
    obj.obs["sample_id"] = sample["sample_id"]
    obj.obs["source_sample"] = sample.get("source_sample", sample["sample_id"])
    obj.obs["patient_id"] = sample["patient_id"]
    obj.obs["donor_id"] = sample["donor_id"]
    obj.obs["specimen_id"] = sample.get("specimen_id", sample["sample_id"])
    obj.obs["tissue"] = "lung"
    obj.obs["condition"] = "LAM"
    obj.obs["assay"] = sample["assay"]
    obj.obs["source_accession"] = accession
    obj.obs["evidence_type"] = sample.get("evidence_type", "independent_donor")
    obj.obs["independence_group"] = sample.get("independence_group", f"{accession}_{sample['patient_id']}")
    obj.obs["identity_status"] = sample.get("identity_status", "mapped_from_staged_metadata")
    obj.obs["treatment"] = "unknown"

    obj.var["mt"] = obj.var["gene_symbol_upper"].str.match(str(cfg.get("technical_patterns", {}).get("mitochondrial", "^MT-"))).to_numpy()
    obj.var["ribo"] = obj.var["gene_symbol_upper"].str.match(str(cfg.get("technical_patterns", {}).get("ribosomal", "^(RPS|RPL)"))).to_numpy()
    sc.pp.calculate_qc_metrics(obj, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

    marker_matrix = []
    for gene in cfg.get("marker_genes", MARKER_GENES):
        values = column_values(obj, gene)
        obj.obs[f"marker_expr_{gene}"] = values
        marker_matrix.append(values > 0)
    marker_matrix = np.column_stack(marker_matrix) if marker_matrix else np.zeros((obj.n_obs, 0), dtype=bool)
    obj.obs["known_marker_genes_detected"] = marker_matrix.sum(axis=1).astype(int)
    obj.obs["known_marker_combo_ge2"] = obj.obs["known_marker_genes_detected"] >= 2
    obj.obs["known_marker_score"] = marker_matrix.mean(axis=1).astype(float) if marker_matrix.shape[1] else 0.0
    obj.obs["lamcore_candidate_author_style"] = False
    obj.obs["lamcore_candidate_formal"] = False
    obj.obs["lamcore_candidate"] = False
    obj.obs["doublet_score"] = np.nan
    obj.obs["doublet_predicted"] = False

    sc.pp.normalize_total(obj, target_sum=10000.0, inplace=True)
    sc.pp.log1p(obj)
    obj.uns["external_provenance"] = {"accession": accession, "sample": sample, "processed_matrix_only": True}
    return obj


def convert_dataset(accession: str, spec: dict) -> tuple[ad.AnnData, list[dict]]:
    objects = []
    source_records = []
    for sample in spec["samples"]:
        path = ROOT / sample["path"]
        if not path.exists():
            raise FileNotFoundError(path)
        obj = read_matrix(path)
        obj = annotate_and_normalize(obj, sample, accession, spec)
        objects.append(obj)
        source_records.append({"accession": accession, "sample_id": sample["sample_id"], "source_sample": sample.get("source_sample"), "input": str(path.relative_to(ROOT)), "n_cells": int(obj.n_obs), "n_genes": int(obj.n_vars), "assay": sample["assay"], "patient_id": sample["patient_id"], "donor_id": sample["donor_id"]})
    combined = ad.concat(objects, join="outer", merge="first", index_unique=None)
    combined.layers["counts"] = (combined.layers["counts"].tocsr() if sp.issparse(combined.layers["counts"]) else sp.csr_matrix(combined.layers["counts"])).astype(np.int32)
    combined.var_names_make_unique()
    # Outer concatenation can introduce missing annotation values for genes
    # absent from the first object. Normalize these columns to HDF5-safe types.
    for col in ["mt", "ribo"]:
        if col in combined.var:
            combined.var[col] = combined.var[col].fillna(False).astype(bool)
    for col in ["gene_symbol", "gene_symbol_upper"]:
        if col in combined.var:
            fallback = pd.Series(combined.var_names.astype(str), index=combined.var.index)
            combined.var[col] = combined.var[col].fillna(fallback).astype(str)
    combined.uns["external_dataset"] = accession
    # AnnData/HDF5 cannot serialize a list of heterogeneous dictionaries in
    # ``uns`` directly. Keep the full provenance losslessly as JSON text.
    combined.uns["conversion"] = {"script": "scripts/prepare_discovery_external.py", "converted_at": datetime.now(timezone.utc).isoformat(), "source_records_json": json.dumps(source_records, ensure_ascii=False), "processed_matrix_only": True}
    return combined, source_records


def build_modality_manifest(root: Path) -> dict:
    return {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "independence_rule": "PatientID counts donors; same-donor RNA/ATAC/spatial adds orthogonal evidence only.",
        "datasets": {
            "GSE190260": {"rna_output": "data/processed/external/GSE190260.h5ad", "rna_samples": ["LAM1110", "LAM1158", "LAM1163", "LAM1164-1", "LAM1164-2", "LAM1164-3"], "patient_ids": ["LAM1110", "LAM1158", "LAM1163", "LAM1164"], "orthogonal_modalities": []},
            "GSE217108": {"rna_output": "data/processed/external/GSE217108.h5ad", "rna_samples": ["LAM32", "LAM44"], "patient_ids": ["LAM32", "LAM44"], "orthogonal_modalities": [{"assay": "snATAC", "patients": ["LAM32", "LAM44"]}]},
            "GSE302356": {"rna_output": "data/processed/external/GSE302356.h5ad", "rna_samples": ["LAM3", "LAM4", "LAM10", "LAM13"], "patient_ids": ["LAM32", "LAM18", "LAM3", "LAM50"], "orthogonal_modalities": [{"assay": "snATAC", "samples": ["LAM14", "LAM15", "LAM16"]}, {"assay": "spatial", "samples": ["LAM18", "LAM19", "LAM20"]}], "patient_id_pending": ["LAM16"]},
        },
        "note": "ATAC, spatial, and Xenium raw/staged files passed archive checks but remain separate from RNA AnnData until modality-specific analysis.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/external_conversion.yaml")
    parser.add_argument("--datasets", nargs="*", default=None, choices=["GSE190260", "GSE217108", "GSE302356"])
    args = parser.parse_args()
    cfg = read_yaml(ROOT / args.config)
    output_dir = ROOT / cfg["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = args.datasets or list(cfg["datasets"])
    all_records = []
    outputs = []
    for accession in selected:
        combined, records = convert_dataset(accession, cfg["datasets"][accession])
        output = output_dir / cfg["datasets"][accession]["output"]
        combined.write_h5ad(output, compression="gzip")
        outputs.append({"accession": accession, "output": str(output.relative_to(ROOT)), "shape": list(combined.shape), "patient_ids": sorted(combined.obs["patient_id"].astype(str).unique()), "n_cells": int(combined.n_obs), "n_genes": int(combined.n_vars)})
        all_records.extend(records)
        print(accession, combined.shape, output)

    manifest = build_modality_manifest(ROOT)
    manifest["converted_outputs"] = outputs
    (ROOT / cfg["modalities_manifest"]).write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    (output_dir / "external_conversion_manifest.json").write_text(json.dumps({"outputs": outputs, "source_records": all_records}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
