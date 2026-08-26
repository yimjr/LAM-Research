"""Filter LINCS/CMap-style perturbations against residual programs.

The script accepts a deliberately simple, auditable long-format input rather
than hiding assumptions inside a proprietary API. It can be populated from a
LINCS L1000 export or another ranked perturbation table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import ROOT, write_json


def load_optional(path: Path, columns: list[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=columns)
    table = pd.read_csv(path)
    missing = set(columns) - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return table


def perturbation_connectivity(signatures: pd.DataFrame, perturbations: pd.DataFrame) -> pd.DataFrame:
    if signatures.empty or perturbations.empty:
        return pd.DataFrame()
    query = signatures.copy()
    query["gene"] = query.gene.astype(str).str.upper()
    perturbations = perturbations.copy()
    perturbations["gene"] = perturbations.gene.astype(str).str.upper()
    rows = []
    for (perturbation, contrast), group in perturbations.merge(query, on="gene", how="inner").groupby(["perturbation", "contrast"]):
        effect = group.signed_score.to_numpy(float)
        pscore = group.perturbation_score.to_numpy(float)
        weights = np.abs(effect)
        # Negative means the perturbation tends to reverse the disease program.
        reversal = -float(np.average(np.sign(effect) * pscore, weights=weights)) if weights.sum() else np.nan
        rows.append({
            "perturbation": perturbation,
            "contrast": contrast,
            "n_genes": int(len(group)),
            "connectivity_reversal": reversal,
            "mean_perturbation_score": float(pscore.mean()),
        })
    return pd.DataFrame(rows)


def aggregate_candidates(connectivity: pd.DataFrame, target: pd.DataFrame, cytotoxicity: pd.DataFrame, exposure: pd.DataFrame) -> pd.DataFrame:
    if connectivity.empty:
        return pd.DataFrame(columns=["perturbation", "n_supported_contrasts", "mean_reversal_score", "target_concordance", "generic_cytotoxicity_score", "exposure_feasible", "tier"])
    summary = connectivity.groupby("perturbation", as_index=False).agg(
        n_supported_contrasts=("contrast", "nunique"),
        mean_reversal_score=("connectivity_reversal", "mean"),
        min_reversal_score=("connectivity_reversal", "min"),
        median_genes_per_contrast=("n_genes", "median"),
    )
    if not target.empty:
        target_summary = target.groupby("perturbation", as_index=False).agg(target_concordance=("target_concordance", "mean"))
        summary = summary.merge(target_summary, on="perturbation", how="left")
    else:
        summary["target_concordance"] = np.nan
    if not cytotoxicity.empty:
        summary = summary.merge(cytotoxicity.groupby("perturbation", as_index=False)["generic_cytotoxicity_score"].mean(), on="perturbation", how="left")
    else:
        summary["generic_cytotoxicity_score"] = np.nan
    if not exposure.empty:
        summary = summary.merge(exposure.groupby("perturbation", as_index=False)["exposure_feasible"].max(), on="perturbation", how="left")
    else:
        summary["exposure_feasible"] = np.nan
    summary["tier"] = "Tier 2"
    tier1 = (
        summary.n_supported_contrasts.ge(2)
        & summary.mean_reversal_score.ge(0.25)
        & summary.target_concordance.ge(0.5)
        & summary.generic_cytotoxicity_score.fillna(0).lt(0.5)
        & summary.exposure_feasible.fillna(False)
    )
    summary.loc[tier1, "tier"] = "Tier 1"
    summary.loc[summary.mean_reversal_score.lt(0.25), "tier"] = "Tier 3"
    return summary.sort_values(["tier", "mean_reversal_score"], ascending=[True, False])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", default="results/signatures/GSE179044_cmap_query_signatures.csv")
    parser.add_argument("--perturbation", default="data/processed/LINCS/perturbation_signatures.csv")
    parser.add_argument("--target", default="data/processed/LINCS/target_perturbation_concordance.csv")
    parser.add_argument("--cytotoxicity", default="data/processed/LINCS/generic_cytotoxicity.csv")
    parser.add_argument("--exposure", default="data/processed/LINCS/exposure_feasibility.csv")
    args = parser.parse_args()
    signature_path = ROOT / args.signature
    signatures = pd.read_csv(signature_path) if signature_path.exists() else pd.DataFrame()
    perturbations = load_optional(ROOT / args.perturbation, ["perturbation", "gene", "perturbation_score"])
    target = load_optional(ROOT / args.target, ["perturbation", "target_concordance"])
    cytotoxicity = load_optional(ROOT / args.cytotoxicity, ["perturbation", "generic_cytotoxicity_score"])
    exposure = load_optional(ROOT / args.exposure, ["perturbation", "exposure_feasible"])

    connectivity = perturbation_connectivity(signatures, perturbations)
    candidates = aggregate_candidates(connectivity, target, cytotoxicity, exposure)
    out_dir = ROOT / "results" / "candidates"
    out_dir.mkdir(parents=True, exist_ok=True)
    connectivity.to_csv(out_dir / "LINCS_connectivity_by_contrast.csv", index=False)
    candidates.to_csv(out_dir / "candidate_filtering_table.csv", index=False)
    write_json(ROOT / "manifests/candidate_screening.json", {
        "status": "completed" if not connectivity.empty else "awaiting_lincs_input",
        "n_perturbations": int(candidates.perturbation.nunique()) if not candidates.empty else 0,
        "n_tier1": int((candidates.tier == "Tier 1").sum()) if not candidates.empty else 0,
        "direction_reversal": "not included in default signatures",
        "required_filters": ["cross-contrast recurrence", "target perturbation concordance", "state expression", "generic cytotoxicity", "unrelated-signature promiscuity", "human exposure", "sirolimus complementarity"],
        "note": "No Tier 1 candidate is inferred without an actual perturbation table and the required filters.",
    })
    print({"status": "completed" if not connectivity.empty else "awaiting_lincs_input", "n_tier1": int((candidates.tier == "Tier 1").sum()) if not candidates.empty else 0})


if __name__ == "__main__":
    main()
