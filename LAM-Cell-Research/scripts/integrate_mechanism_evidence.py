#!/usr/bin/env python3
"""Build a cross-modal evidence table for the unified LAM mechanism hypothesis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

PROTEIN_PRIORITY = {"MMP8", "PMEL", "VEGFD", "CCL14", "CTSK", "MMP2", "MMP9", "CTSS", "ELANE", "S100A4", "FAP"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/mechanism_integration")
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    retention = pd.read_csv(ROOT / "results/perturbation/GSE179044_gene_retention_effects.csv")
    retained = pd.read_csv(ROOT / "results/perturbation/GSE179044_cross_environment_retained_genes.csv")
    source = pd.read_csv(ROOT / "results/spatial/GSE302356/single_cell_protease_source_attribution.csv")
    adaptation = pd.read_csv(ROOT / "results/adaptation/four_group_lung_acquired_interactions.csv")
    known = pd.read_csv(ROOT / "results/program_discovery/benchmark/known_program_extended_benchmark_summary.csv")
    spatial_observed: dict[str, set[str]] = {}
    for modality in ["visium", "visium_hd", "xenium"]:
        path = ROOT / f"results/spatial/GSE302356/{modality}_source_attribution_and_balance.csv"
        if path.exists():
            columns = pd.read_csv(path, nrows=0).columns
            spatial_observed[modality] = {str(c).removeprefix("gene_").upper() for c in columns if str(c).startswith("gene_")}
        else:
            spatial_observed[modality] = set()
    # Gene membership is explicit and traceable rather than inferred from a
    # shared score name.
    program_membership = {
        "lung_adaptation": {"VEGFD", "FIGF", "CCL21", "PDPN", "LYVE1", "FLT4", "CTSK", "MMP2", "MMP9", "CXCL14"},
        "ECM_interaction": {"COL1A1", "COL1A2", "COL3A1", "FN1", "VCAN", "MMP2", "MMP9", "TIMP1", "TIMP2"},
        "survival_stress": {"BCL2", "BCL2L1", "LGALS3", "HIF1A", "VEGFA", "BNIP3", "DDIT4", "NDRG1", "TXNIP"},
        "immune_evasion": {"LGALS3", "CTSS", "TYROBP", "CCL2", "TGFB1", "CD274", "SERPINE1"},
    }
    genes = set(retained.loc[retained["cross_environment_retention_candidate"], "gene"].astype(str).str.upper())
    genes.update(source["gene"].astype(str).str.upper())
    genes.update(PROTEIN_PRIORITY)
    rows = []
    for gene in sorted(genes):
        gsource = source[source["gene"].astype(str).str.upper().eq(gene)]
        top = gsource.sort_values("transcript_fraction", ascending=False).iloc[0] if len(gsource) else None
        gret = retained[retained["gene"].astype(str).str.upper().eq(gene)]
        eligible = gret[gret["cross_environment_retention_candidate"]] if len(gret) else gret
        classes = sorted(set(gret["categories"].dropna().astype(str))) if len(gret) else []
        adapt_programs = [name for name, member_genes in program_membership.items() if gene in member_genes]
        rows.append({
            "gene": gene,
            "protein_priority_candidate": gene in PROTEIN_PRIORITY,
            "retention_eligible_across_environment": bool(len(eligible)),
            "retention_categories": ";".join(classes),
            "single_cell_top_source_state": str(top["source_state"]) if top is not None else "not_measured",
            "single_cell_top_source_transcript_fraction": float(top["transcript_fraction"]) if top is not None and np.isfinite(top["transcript_fraction"]) else np.nan,
            "spatial_gene_observed_visium": gene in spatial_observed["visium"],
            "spatial_gene_observed_visium_hd": gene in spatial_observed["visium_hd"],
            "spatial_gene_observed_xenium": gene in spatial_observed["xenium"],
            "adaptation_program_membership": ";".join(adapt_programs),
            "protein_matrix_status": "pending_public_matrix" if gene in PROTEIN_PRIORITY else "not_prioritized",
            "patient_matching": "not_assumed",
            "evidence_level": "candidate_cross_modal_link" if gene in PROTEIN_PRIORITY and (len(eligible) or len(gsource)) else "exploratory",
        })
    table = pd.DataFrame(rows)
    table.to_csv(out / "mechanism_evidence_crosswalk.csv", index=False)
    adaptation.to_csv(out / "lung_adaptation_interaction_snapshot.csv", index=False)
    known.to_csv(out / "known_program_benchmark_snapshot.csv", index=False)
    manifest = {
        "inputs": [
            "results/perturbation/GSE179044_cross_environment_retained_genes.csv",
            "results/spatial/GSE302356/single_cell_protease_source_attribution.csv",
            "results/adaptation/four_group_lung_acquired_interactions.csv",
            "results/program_discovery/benchmark/known_program_extended_benchmark_summary.csv",
        ],
        "protein_sources": ["MSV000099051", "LAM_SomaScan_2025"],
        "protein_status": "Protein matrices remain pending; this table does not claim a plasma/EV abundance result.",
        "integration_rule": "Cross-cohort evidence integration only; no patient-level RNA-protein correlation without matched PatientID.",
        "priority_genes": sorted(PROTEIN_PRIORITY),
        "interpretation": "A cross-modal candidate link is not a confirmed mechanism and requires independent validation.",
    }
    (out / "mechanism_evidence_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps({"genes": len(table), "priority_candidates": int(table["protein_priority_candidate"].sum()), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
