"""Patient-level state associations and antigen/presentation co-existence."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import PROJECT_ROOT, ensure_output_path, load_signatures, project_relative, write_json


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    pair = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(pair) < 3 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, len(pair)
    rho, p_value = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
    return float(rho), float(p_value), len(pair)


def patient_state_associations() -> pd.DataFrame:
    path = PROJECT_ROOT / "results" / "patient_summaries" / "patient_module_summary.csv"
    table = pd.read_csv(path)
    table = table[table["identity_pool"].eq("high_confidence")]
    rows = []
    state_modules = ["stem_like", "inflammatory"]
    visibility_modules = ["antigen_associated", "presentation_machinery", "ifn_response", "immune_evasion", "nk_ligands"]
    for dataset, group in table.groupby("dataset", observed=True):
        pivot = group.pivot_table(index="patient_id", columns="module", values="mean_score", aggfunc="mean")
        for state in state_modules:
            for visibility in visibility_modules:
                if state not in pivot or visibility not in pivot:
                    continue
                rho, p_value, n = safe_spearman(pivot[state], pivot[visibility])
                rows.append({
                    "dataset": dataset,
                    "state_module": state,
                    "visibility_module": visibility,
                    "n_patients": n,
                    "spearman_rho": rho,
                    "p_value_descriptive": p_value,
                    "interpretation": "patient-level descriptive association; not causal evidence",
                })
    return pd.DataFrame(rows)


def antigen_presentation_coexistence() -> pd.DataFrame:
    signatures = load_signatures()
    genes = [gene.upper() for gene in signatures["antigen_associated"]["genes"]]
    rows = []
    for path in sorted((PROJECT_ROOT / "results" / "cell_scores").glob("*_cell_visibility_scores.csv")):
        table = pd.read_csv(path)
        table = table[table["identity_pool"].eq("high_confidence")]
        if table.empty:
            continue
        for gene in genes:
            state_column = f"state__{gene}"
            if state_column not in table:
                continue
            states = table[state_column].astype(str)
            assayed = states.ne("not_assayed")
            antigen_detected = states.isin(["detected_low", "detected_high", "detected_unclassified"])
            for target in ["presentation_machinery", "ifn_response", "immune_evasion", "nk_ligands"]:
                score_column = f"module_{target}"
                if score_column not in table:
                    continue
                values = pd.to_numeric(table[score_column], errors="coerce")
                rows.append({
                    "dataset": table["dataset"].iloc[0],
                    "gene": gene,
                    "target_module": target,
                    "n_assayed_cells": int(assayed.sum()),
                    "n_antigen_detected_cells": int((antigen_detected & assayed).sum()),
                    "mean_target_when_antigen_detected": float(values[antigen_detected & assayed].mean()) if (antigen_detected & assayed).any() else np.nan,
                    "mean_target_when_antigen_not_detected": float(values[(~antigen_detected) & assayed].mean()) if ((~antigen_detected) & assayed).any() else np.nan,
                    "difference_detected_minus_not_detected": float(values[antigen_detected & assayed].mean() - values[(~antigen_detected) & assayed].mean()) if (antigen_detected & assayed).any() and ((~antigen_detected) & assayed).any() else np.nan,
                    "interpretation": "co-existence summary; expression is not proof of peptide presentation or immune recognition",
                })
    return pd.DataFrame(rows)


def main() -> None:
    output_dir = PROJECT_ROOT / "results" / "immune_context"
    states = patient_state_associations()
    coexistence = antigen_presentation_coexistence()
    states.to_csv(ensure_output_path(output_dir / "patient_state_visibility_associations.csv"), index=False)
    coexistence.to_csv(ensure_output_path(output_dir / "antigen_presentation_coexistence.csv"), index=False)
    write_json(PROJECT_ROOT / "manifests" / "state_association_manifest.json", {
        "patient_unit": "patient_id within dataset",
        "minimum_patients_for_descriptive_association": 3,
        "outputs": [project_relative(output_dir / "patient_state_visibility_associations.csv"), project_relative(output_dir / "antigen_presentation_coexistence.csv")],
    })
    print(f"patient associations: {len(states)} rows; coexistence: {len(coexistence)} rows")


if __name__ == "__main__":
    main()
