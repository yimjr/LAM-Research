"""Create donor/patient summaries from per-cell visibility scores."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from common import PROJECT_ROOT, ensure_output_path, load_signatures, project_relative, write_json


def load_scores() -> list[pd.DataFrame]:
    tables = []
    for path in sorted((PROJECT_ROOT / "results" / "cell_scores").glob("*_cell_visibility_scores.csv")):
        tables.append(pd.read_csv(path, index_col="cell_id"))
    if not tables:
        raise FileNotFoundError("Run score_visibility_modules.py before summarizing.")
    return tables


def module_summary(table: pd.DataFrame) -> pd.DataFrame:
    modules = sorted({column.removeprefix("module_") for column in table.columns if column.startswith("module_") and "__" not in column})
    rows = []
    group_columns = ["dataset", "patient_id", "identity_pool"]
    for keys, group in table.groupby(group_columns, dropna=False, observed=True):
        dataset, patient, pool = keys
        for module in modules:
            score = pd.to_numeric(group[f"module_{module}"], errors="coerce")
            available = pd.to_numeric(group.get(f"module_{module}__n_available"), errors="coerce")
            detected = pd.to_numeric(group.get(f"module_{module}__n_detected"), errors="coerce")
            status = group.get(f"module_{module}__status", pd.Series("unknown", index=group.index)).astype(str)
            rows.append({
                "dataset": dataset,
                "patient_id": patient,
                "identity_pool": pool,
                "n_cells": len(group),
                "module": module,
                "mean_score": score.mean(),
                "median_score": score.median(),
                "n_genes_available": available.iloc[0] if len(available) else np.nan,
                "mean_genes_detected": detected.mean(),
                "module_status": ";".join(sorted(status.dropna().unique())),
            })
    return pd.DataFrame(rows)


def antigen_summary(table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    group_columns = ["dataset", "patient_id", "identity_pool"]
    for keys, group in table.groupby(group_columns, dropna=False, observed=True):
        dataset, patient, pool = keys
        for gene in genes:
            state_column = f"state__{gene.upper()}"
            if state_column not in group:
                rows.append({"dataset": dataset, "patient_id": patient, "identity_pool": pool, "gene": gene.upper(), "n_cells": len(group), "n_assayed_cells": 0, "not_assayed_fraction": 1.0})
                continue
            states = group[state_column].astype(str)
            n = len(states)
            assayed = states.ne("not_assayed")
            detected = states.isin(["detected_low", "detected_high", "detected_unclassified"])
            rows.append({
                "dataset": dataset,
                "patient_id": patient,
                "identity_pool": pool,
                "gene": gene.upper(),
                "n_cells": n,
                "n_assayed_cells": int(assayed.sum()),
                "not_assayed_fraction": float((~assayed).mean()),
                "not_detected_fraction_among_assayed": float((states[assayed] == "not_detected").mean()) if assayed.any() else np.nan,
                "detected_fraction_among_assayed": float(detected[assayed].mean()) if assayed.any() else np.nan,
                "detected_low_fraction_among_assayed": float((states[assayed] == "detected_low").mean()) if assayed.any() else np.nan,
                "detected_high_fraction_among_assayed": float((states[assayed] == "detected_high").mean()) if assayed.any() else np.nan,
                "detected_unclassified_fraction_among_assayed": float((states[assayed] == "detected_unclassified").mean()) if assayed.any() else np.nan,
            })
    return pd.DataFrame(rows)


def main() -> None:
    signatures = load_signatures()
    tables = load_scores()
    modules = pd.concat([module_summary(table) for table in tables], ignore_index=True)
    antigens = pd.concat([antigen_summary(table, signatures["antigen_associated"]["genes"]) for table in tables], ignore_index=True)
    output_dir = PROJECT_ROOT / "results" / "patient_summaries"
    modules.to_csv(ensure_output_path(output_dir / "patient_module_summary.csv"), index=False)
    antigens.to_csv(ensure_output_path(output_dir / "patient_antigen_summary.csv"), index=False)

    dataset = modules.groupby(["dataset", "identity_pool"], dropna=False, observed=True).agg(
        patients=("patient_id", "nunique"),
        cells=("n_cells", "max"),
    ).reset_index()
    dataset.to_csv(ensure_output_path(output_dir / "dataset_pool_summary.csv"), index=False)
    write_json(PROJECT_ROOT / "manifests" / "patient_summary_manifest.json", {
        "module_summary": project_relative(output_dir / "patient_module_summary.csv"),
        "antigen_summary": project_relative(output_dir / "patient_antigen_summary.csv"),
        "denominator_rule": "not_assayed cells are excluded from detected/low/high fractions",
    })
    print(dataset.to_string(index=False))


if __name__ == "__main__":
    main()
