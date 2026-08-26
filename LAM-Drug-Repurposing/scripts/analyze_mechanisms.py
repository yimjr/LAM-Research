"""Prepare the STAT3 and SRPK2 mechanism-support datasets.

GSE84476 is a transcript-level kallisto export and therefore requires an
explicit transcript-to-gene map before gene-level claims are made. GSE104335
contains an Affymetrix HTA 2.0 archive with processed gene-level CHP values.
"""

from __future__ import annotations

import argparse
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from common import ROOT, write_json


def gse84476_design(path: Path) -> pd.DataFrame:
    # The GEO file is whitespace-delimited despite its .txt name.
    header = pd.read_csv(path, sep=r"\s+", engine="python", compression="infer", nrows=0).columns.tolist()
    sample_names = sorted({match.group(1) for column in header for match in [re.search(r"(?:tpm|est_counts)\.([^\.]+)$", column)] if match})
    return pd.DataFrame([
        {"sample": sample, "cell_context": sample.split("_")[0], "perturbation": sample.split("_")[1]}
        for sample in sample_names
    ])


def aggregate_if_mapped(path: Path, mapping_path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep=r"\s+", engine="python", compression="infer")
    mapping = pd.read_csv(mapping_path)
    required = {"target_id", "gene_symbol"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"Transcript map must contain {sorted(required)}")
    table = table.merge(mapping[list(required)], on="target_id", how="inner")
    tpm_cols = [column for column in table.columns if column.startswith("tpm.")]
    return table.groupby("gene_symbol", sort=False)[tpm_cols].sum()


def gse104335_sample_design(raw_path: Path) -> pd.DataFrame:
    """Recover the gene-expression/splice pairing and treatment groups."""
    with tarfile.open(raw_path) as archive:
        names = [member.name for member in archive.getmembers()]
    rows = []
    groups = [
        ("shGFP_vehicle", "shGFP", "vehicle", 1),
        ("shGFP_vehicle", "shGFP", "vehicle", 2),
        ("shGFP_vehicle", "shGFP", "vehicle", 3),
        ("shGFP_rapamycin", "shGFP", "rapamycin", 1),
        ("shGFP_rapamycin", "shGFP", "rapamycin", 2),
        ("shGFP_rapamycin", "shGFP", "rapamycin", 3),
        ("shSRPK2_vehicle", "shSRPK2", "vehicle", 1),
        ("shSRPK2_vehicle", "shSRPK2", "vehicle", 2),
        ("shSRPK2_vehicle", "shSRPK2", "vehicle", 3),
    ]
    gene_chps = sorted(name for name in names if name.endswith("sst-rma-gene-full.chp.gz"))
    splice_chps = sorted(name for name in names if name.endswith("sst-rma-alt-splice-sst-dabg.chp.gz"))
    for idx, (group, perturbation, treatment, replicate) in enumerate(groups):
        gene_name = gene_chps[idx] if idx < len(gene_chps) else ""
        splice_name = splice_chps[idx] if idx < len(splice_chps) else ""
        gsm = gene_name.split("_")[0] if gene_name else f"GSM279555{idx + 3}"
        rows.append({
            "sample_id": gsm,
            "group": group,
            "perturbation": perturbation,
            "treatment": treatment,
            "replicate": replicate,
            "layer": "gene_expression",
            "chp_file": gene_name,
        })
        rows.append({
            "sample_id": f"GSM27955{idx + 62}",
            "paired_gene_expression_sample_id": gsm,
            "group": group,
            "perturbation": perturbation,
            "treatment": treatment,
            "replicate": replicate,
            "layer": "splice_isoform",
            "chp_file": splice_name,
        })
    return pd.DataFrame(rows)


def build_gse104335_cross_dataset_summary(out_dir: Path) -> bool:
    """Join SRPK2/rapamycin mechanism contrasts to the human discovery table."""
    mechanism_path = out_dir / "GSE104335_gene_level_contrasts.csv"
    discovery_path = ROOT / "results/tables/GSE179044_factorial_contrasts.csv"
    if not mechanism_path.exists() or not discovery_path.exists():
        return False
    mechanism = pd.read_csv(mechanism_path)
    wide = mechanism.pivot(index="gene_symbol", columns="contrast", values=["logFC", "P.Value", "adj.P.Val"])
    wide.columns = [f"gse104335_{metric}_{contrast}" for metric, contrast in wide.columns]
    wide = wide.reset_index().rename(columns={"gene_symbol": "gene"})
    discovery = pd.read_csv(discovery_path)
    keep = [
        "gene",
        "tsc2_loss_hydrogel", "tsc2_loss_hydrogel_moderated_q",
        "residual_hydrogel", "residual_hydrogel_moderated_q",
        "hydrogel_specific_residual", "hydrogel_specific_residual_moderated_q",
        "escape_hydrogel", "escape_hydrogel_moderated_q",
        "environment_dependent_escape", "environment_dependent_escape_moderated_q",
    ]
    summary = wide.merge(discovery[keep], on="gene", how="outer")
    summary.to_csv(out_dir / "GSE104335_cross_dataset_summary.csv", index=False)
    rap_logfc = "gse104335_logFC_shGFP_rapamycin_minus_shGFP_vehicle"
    rap_q = "gse104335_adj.P.Val_shGFP_rapamycin_minus_shGFP_vehicle"
    hydrogel_logfc = "hydrogel_specific_residual"
    hydrogel_q = "hydrogel_specific_residual_moderated_q"
    overlap = summary.loc[summary[rap_q].lt(0.05) & summary[hydrogel_q].lt(0.05), [
        "gene", rap_logfc, rap_q, hydrogel_logfc, hydrogel_q,
        "gse104335_logFC_shSRPK2_vehicle_minus_shGFP_vehicle",
        "gse104335_adj.P.Val_shSRPK2_vehicle_minus_shGFP_vehicle",
    ]].copy()
    overlap["same_direction"] = overlap[rap_logfc] * overlap[hydrogel_logfc] > 0
    overlap = overlap.sort_values(["same_direction", rap_q], ascending=[False, True])
    overlap.to_csv(out_dir / "GSE104335_hydrogel_specific_overlap.csv", index=False)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcript-map", default="data/processed/mechanisms/ensembl_transcript_to_gene.csv")
    args = parser.parse_args()
    out_dir = ROOT / "results" / "mechanisms"
    out_dir.mkdir(parents=True, exist_ok=True)
    gse84476 = ROOT / "data/raw/GSE84476/GSE84476_VKMJ-STAT3_kallisto_abundance.txt.gz"
    design = gse84476_design(gse84476) if gse84476.exists() else pd.DataFrame()
    design.to_csv(out_dir / "GSE84476_sample_design.csv", index=False)
    mapped = ROOT / args.transcript_map
    status = {"GSE84476": "awaiting_explicit_transcript_to_gene_map", "GSE104335": "awaiting_processed_CHP_parse"}
    if mapped.exists() and gse84476.exists():
        gene_tpm = aggregate_if_mapped(gse84476, mapped)
        gene_tpm.to_csv(out_dir / "GSE84476_gene_level_tpm.csv")
        log_tpm = np.log2(gene_tpm.astype(float).clip(lower=0) + 1.0)
        contrasts = pd.DataFrame(index=log_tpm.index)
        for context in ("102cell", "103cell"):
            ctrl = log_tpm[f"tpm.{context}_siCtrl"]
            stat3 = log_tpm[f"tpm.{context}_siSTAT3"]
            rapa = log_tpm[f"tpm.{context}_Rap"]
            contrasts[f"{context}_siSTAT3_minus_siCtrl"] = stat3 - ctrl
            contrasts[f"{context}_rapamycin_minus_siCtrl"] = rapa - ctrl
            contrasts[f"{context}_siSTAT3_minus_rapamycin"] = stat3 - rapa
        contrasts.index.name = "gene"
        contrasts.to_csv(out_dir / "GSE84476_gene_level_log2_tpm_contrasts.csv")
        fold_change_path = ROOT / "data/raw/GSE84476/GSE84476_VKMJ-STAT3_fold_changes.txt.gz"
        if fold_change_path.exists():
            fold_changes = pd.read_csv(fold_change_path, sep="\t", compression="infer")
            fold_changes["gene"] = fold_changes["gene"].fillna("").astype(str)
            fold_summary = fold_changes[fold_changes["gene"].ne("")].groupby("gene", sort=False)["fold_change"].agg(["mean", "median", "max", "count"])
            fold_summary.to_csv(out_dir / "GSE84476_transcript_fold_change_by_gene.csv")
        status["GSE84476"] = "gene_level_table_created"
    raw104335 = ROOT / "data/raw/GSE104335/GSE104335_RAW.tar"
    if raw104335.exists() and tarfile.is_tarfile(raw104335):
        with tarfile.open(raw104335) as archive:
            members = pd.DataFrame({"member": [member.name for member in archive.getmembers()]})
        members.to_csv(out_dir / "GSE104335_archive_manifest.csv", index=False)
        gse104335_sample_design(raw104335).to_csv(out_dir / "GSE104335_sample_design.csv", index=False)
        if (out_dir / "GSE104335_gene_level_contrasts.csv").exists():
            status["GSE104335"] = "gene_level_CHP_contrasts_created"
            build_gse104335_cross_dataset_summary(out_dir)
        else:
            status["GSE104335"] = "archive_manifest_and_sample_design_created_gene_level_pending"
    write_json(ROOT / "manifests/mechanism_analysis.json", {
        "status": status,
        "GSE84476_role": "STAT3/TSC2/rapamycin mechanism support; transcript-level source requires mapping",
        "GSE104335_role": "SRPK2/rapamycin mechanism comparison; processed HTA 2.0 CHP values mapped to gene symbols and analyzed with limma",
        "caveat": "Mechanism-support datasets are not substitutes for the factorial discovery or external 2x2 validation.",
    })
    print({"status": status, "n_gse84476_samples": int(len(design))})


if __name__ == "__main__":
    main()
