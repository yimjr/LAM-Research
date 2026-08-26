"""Map experimental programs to patient-level human cell states.

The primary interface is a processed h5ad with cell metadata and a gene-symbol
column. Scores are computed on patient-by-state pseudobulk profiles, then
converted to within-gene rank scores across profiles. This avoids treating
cell-level observations as independent patients.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.stats import rankdata

from common import ROOT, write_json


def load_sets() -> tuple[dict[str, dict[str, list[str]]], dict[str, list[str]]]:
    signatures = pd.read_csv(ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv")
    signed_sets: dict[str, dict[str, list[str]]] = {}
    for contrast, group in signatures.groupby("contrast"):
        signed_sets[contrast] = {
            "up": group.loc[group.direction.eq("up"), "gene"].astype(str).str.upper().unique().tolist(),
            "down": group.loc[group.direction.eq("down"), "gene"].astype(str).str.upper().unique().tolist(),
        }
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text())
    modules = {name: [str(gene).upper() for gene in genes] for name, genes in config["module_sets"].items()}
    return signed_sets, modules


def pseudobulk(adata: ad.AnnData, gene_symbols: list[str], state_column: str, patient_column: str, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    symbols = adata.var["gene_symbol_upper"].astype(str).str.upper() if "gene_symbol_upper" in adata.var else adata.var_names.astype(str).str.upper()
    symbol_to_indices: dict[str, list[int]] = {}
    for idx, symbol in enumerate(symbols):
        if symbol and symbol != "NAN":
            symbol_to_indices.setdefault(symbol, []).append(idx)
    selected_symbols = [symbol for symbol in dict.fromkeys(gene_symbols) if symbol in symbol_to_indices]
    selected_indices = [idx for symbol in selected_symbols for idx in symbol_to_indices[symbol]]
    if not selected_indices:
        raise ValueError("None of the requested genes are present in the h5ad gene universe")
    obs = adata.obs[[state_column, patient_column]].copy()
    obs[state_column] = obs[state_column].astype(str)
    obs[patient_column] = obs[patient_column].astype(str)
    obs["group"] = obs[patient_column] + "|" + obs[state_column]
    counts = obs.groupby("group", sort=True).size().rename("n_cells")
    keep = counts[counts >= min_cells].index
    obs = obs.loc[obs.group.isin(keep)]
    # Slice the expression matrix after filtering cells; backed AnnData keeps
    # this bounded by the selected gene set rather than loading all features.
    x = adata[obs.index, selected_indices].X
    if sp.issparse(x):
        x = x.tocsr()
    else:
        x = np.asarray(x)
    group_codes, groups = pd.factorize(obs["group"], sort=True)
    indicator = sp.csr_matrix((np.ones(len(obs)), (group_codes, np.arange(len(obs)))), shape=(len(groups), len(obs)))
    means = indicator @ x
    if sp.issparse(means):
        means = means.toarray()
    # Collapse duplicate platform features to the mean symbol value.
    collapsed = np.zeros((len(groups), len(selected_symbols)), dtype=float)
    for col, symbol in enumerate(selected_symbols):
        feature_cols = [j for j, idx in enumerate(selected_indices) if symbols.iloc[idx] == symbol]
        collapsed[:, col] = means[:, feature_cols].mean(axis=1)
    group_meta = obs.drop_duplicates("group").set_index("group").loc[groups, [patient_column, state_column]]
    group_meta["n_cells"] = counts.loc[groups].to_numpy()
    expression = pd.DataFrame(collapsed, index=groups, columns=selected_symbols)
    group_meta.index.name = "group"
    return expression, group_meta.reset_index()


def rank_scores(expression: pd.DataFrame, group_meta: pd.DataFrame, signed_sets: dict[str, dict[str, list[str]]], modules: dict[str, list[str]]) -> pd.DataFrame:
    ranks = pd.DataFrame(
        rankdata(expression.to_numpy(), axis=0, method="average"),
        index=expression.index,
        columns=expression.columns,
    )
    ranks = ranks.divide(len(expression), axis=0)
    rows = []
    for group in expression.index:
        meta = group_meta.loc[group_meta["group"].eq(group)].iloc[0]
        base = {"group": group, "patient_id": meta.iloc[1], "state": meta.iloc[2], "n_cells": int(meta.n_cells)}
        for contrast, sets in signed_sets.items():
            up = [gene for gene in sets["up"] if gene in ranks.columns]
            down = [gene for gene in sets["down"] if gene in ranks.columns]
            base[f"{contrast}_up_rank"] = float(ranks.loc[group, up].mean()) if up else np.nan
            base[f"{contrast}_down_rank"] = float(ranks.loc[group, down].mean()) if down else np.nan
            base[f"{contrast}_signed_enrichment"] = (
                base[f"{contrast}_up_rank"] - base[f"{contrast}_down_rank"]
                if up and down else np.nan
            )
        for module, genes in modules.items():
            available = [gene for gene in genes if gene in ranks.columns]
            base[f"module_{module}_rank"] = float(ranks.loc[group, available].mean()) if available else np.nan
        rows.append(base)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/GSE135851/GSE135851_lam_states_snapshot.h5ad")
    parser.add_argument("--state-column", default="lamcore_label")
    parser.add_argument("--patient-column", default="donor_id")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--dataset", default="GSE135851")
    args = parser.parse_args()
    input_path = ROOT / args.input
    adata = ad.read_h5ad(input_path, backed="r")
    if args.state_column not in adata.obs or args.patient_column not in adata.obs:
        raise ValueError(f"Missing required metadata columns: {args.state_column}, {args.patient_column}")
    signed_sets, modules = load_sets()
    all_genes = sorted({gene for sets in signed_sets.values() for genes in sets.values() for gene in genes} | {gene for genes in modules.values() for gene in genes})
    expression, group_meta = pseudobulk(adata, all_genes, args.state_column, args.patient_column, args.min_cells)
    scores = rank_scores(expression, group_meta, signed_sets, modules)
    out_dir = ROOT / "results" / "human_mapping"
    out_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(out_dir / f"{args.dataset}_patient_state_scores.csv", index=False)
    summary = scores.groupby("state", dropna=False).mean(numeric_only=True).reset_index()
    summary.to_csv(out_dir / f"{args.dataset}_state_summary.csv", index=False)
    write_json(ROOT / "manifests" / f"{args.dataset}_human_mapping.json", {
        "dataset": args.dataset,
        "input": str(args.input),
        "state_column": args.state_column,
        "patient_column": args.patient_column,
        "n_patient_state_profiles": int(len(scores)),
        "n_patients": int(scores.patient_id.nunique()),
        "n_states": int(scores.state.nunique()),
        "scoring": "patient-by-state pseudobulk followed by within-gene rank scoring; cells are not treated as independent replicates",
        "limitation": "GSE135851 snapshot contains candidate/other labels, not the formal LAMCORE1/2/3 and LAF-seed/LAF-niche states; formal state mapping requires GSE302356 processed state tables.",
    })
    print({"dataset": args.dataset, "profiles": int(len(scores)), "patients": int(scores.patient_id.nunique()), "states": sorted(scores.state.unique().tolist())})


if __name__ == "__main__":
    main()
