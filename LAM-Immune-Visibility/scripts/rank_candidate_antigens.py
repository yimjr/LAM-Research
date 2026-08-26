"""Rank antigen-associated expression candidates with explicit missingness."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import PROJECT_ROOT, ensure_output_path, load_signatures, load_source_manifest, project_relative, resolve_source, write_json


def state_flags(states: pd.Series) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    states = states.astype(str)
    assayed = states.ne("not_assayed")
    detected = states.isin(["detected_low", "detected_high", "detected_unclassified"])
    low = states.eq("detected_low")
    high = states.eq("detected_high")
    return assayed, detected, low, high


def independence_key(dataset: str, patient: str) -> str:
    # The existing donor registry explicitly maps LAM32 across GSE217108 and
    # GSE302356. Core GSE135851 donor labels are kept accession-scoped because
    # identical-looking labels such as LAM3 are not sufficient evidence of a
    # cross-study patient match.
    return str(patient) if dataset in {"GSE217108", "GSE302356"} else f"{dataset}:{patient}"


def main() -> None:
    signatures = load_signatures()
    genes = [str(g).upper() for g in signatures["antigen_associated"]["genes"]]
    tables = {path.stem.replace("_cell_visibility_scores", ""): pd.read_csv(path) for path in sorted((PROJECT_ROOT / "results" / "cell_scores").glob("*_cell_visibility_scores.csv"))}
    normal = tables.get("GSE122960_normal_lung")
    lam_tables = {name: table for name, table in tables.items() if name != "GSE122960_normal_lung"}
    rows = []
    patient_rows = []
    for gene in genes:
        dataset_fractions = []
        lam_detected_cells = []
        lam_assayed_cells = []
        lam_low_cells = []
        lam_high_cells = []
        lam_total_high_conf_cells = []
        presentation_detected_values = []
        presentation_not_detected_values = []
        evasion_detected_values = []
        evasion_not_detected_values = []
        patient_detected_flags: dict[str, bool] = {}
        for dataset, table in lam_tables.items():
            table = table[table["pool_high_confidence"].astype(bool)]
            state_col = f"state__{gene}"
            if state_col not in table:
                continue
            for (patient,), group in table.groupby(["patient_id"], observed=True):
                assayed, detected, low, high = state_flags(group[state_col])
                fraction = float(detected[assayed].mean()) if assayed.any() else np.nan
                patient_key = str(patient)
                independent_key = independence_key(dataset, patient_key)
                patient_detected_flags[independent_key] = patient_detected_flags.get(independent_key, False) or (bool(detected[assayed].any()) if assayed.any() else False)
                patient_rows.append({
                    "dataset": dataset,
                    "patient_id": patient,
                    "independence_key": independent_key,
                    "gene": gene,
                    "n_cells": len(group),
                    "n_assayed_cells": int(assayed.sum()),
                    "not_assayed_fraction": float((~assayed).mean()),
                    "detected_fraction_among_assayed": fraction,
                    "detected_low_fraction_among_assayed": float(low[assayed].mean()) if assayed.any() else np.nan,
                    "detected_high_fraction_among_assayed": float(high[assayed].mean()) if assayed.any() else np.nan,
                })
            all_assayed, all_detected, all_low, all_high = state_flags(table[state_col])
            lam_total_high_conf_cells.append(len(table))
            lam_assayed_cells.append(int(all_assayed.sum()))
            lam_detected_cells.append(int((all_detected & all_assayed).sum()))
            lam_low_cells.append(int((all_low & all_assayed).sum()))
            lam_high_cells.append(int((all_high & all_assayed).sum()))
            dataset_fractions.append(f"{dataset}:{float(all_detected[all_assayed].mean()) if all_assayed.any() else np.nan:.4g}")
            for target in ["presentation_machinery", "immune_evasion"]:
                score_col = f"module_{target}"
                if score_col not in table:
                    continue
                score = pd.to_numeric(table[score_col], errors="coerce")
                presentation_detected_values.extend(score[all_detected & all_assayed].dropna().tolist() if target == "presentation_machinery" else [])
                presentation_not_detected_values.extend(score[(~all_detected) & all_assayed].dropna().tolist() if target == "presentation_machinery" else [])
                evasion_detected_values.extend(score[all_detected & all_assayed].dropna().tolist() if target == "immune_evasion" else [])
                evasion_not_detected_values.extend(score[(~all_detected) & all_assayed].dropna().tolist() if target == "immune_evasion" else [])

        normal_fraction = np.nan
        normal_assayed = 0
        normal_detected = 0
        if normal is not None and f"state__{gene}" in normal:
            assayed, detected, _, _ = state_flags(normal[f"state__{gene}"])
            normal_assayed = int(assayed.sum())
            normal_detected = int((detected & assayed).sum())
            normal_fraction = float(detected[assayed].mean()) if assayed.any() else np.nan
        lam_fraction = (sum(lam_detected_cells) / sum(lam_assayed_cells)) if sum(lam_assayed_cells) else np.nan
        rows.append({
            "gene": gene,
            "n_lam_datasets_with_high_confidence_cells": len(dataset_fractions),
            "n_lam_patients": len(patient_detected_flags),
            "n_lam_patients_with_any_detected": int(sum(patient_detected_flags.values())),
            "patient_consistency_any_detected": float(np.mean(list(patient_detected_flags.values()))) if patient_detected_flags else np.nan,
            "lam_detected_fraction_among_assayed": lam_fraction,
            "lam_detected_low_fraction_among_assayed": (sum(lam_low_cells) / sum(lam_assayed_cells)) if sum(lam_assayed_cells) else np.nan,
            "lam_detected_high_fraction_among_assayed": (sum(lam_high_cells) / sum(lam_assayed_cells)) if sum(lam_assayed_cells) else np.nan,
            "lam_not_assayed_cells": int(sum(lam_total_high_conf_cells) - sum(lam_assayed_cells)) if lam_assayed_cells else np.nan,
            "normal_lung_assayed_cells": normal_assayed,
            "normal_lung_detected_cells": normal_detected,
            "normal_lung_detected_fraction_among_assayed": normal_fraction,
            "lam_minus_normal_lung_detected_fraction": lam_fraction - normal_fraction if np.isfinite(lam_fraction) and np.isfinite(normal_fraction) else np.nan,
            "mean_presentation_when_antigen_detected": float(np.mean(presentation_detected_values)) if presentation_detected_values else np.nan,
            "mean_presentation_when_antigen_not_detected": float(np.mean(presentation_not_detected_values)) if presentation_not_detected_values else np.nan,
            "mean_evasion_when_antigen_detected": float(np.mean(evasion_detected_values)) if evasion_detected_values else np.nan,
            "mean_evasion_when_antigen_not_detected": float(np.mean(evasion_not_detected_values)) if evasion_not_detected_values else np.nan,
            "dataset_detected_fraction": ";".join(dataset_fractions),
            "evidence_class": "antigen-associated expression candidate; not confirmed immunogenicity",
        })
    ranking = pd.DataFrame(rows)
    if ranking.empty:
        raise RuntimeError("No antigen candidates were available in scored tables.")
    retention_path = resolve_source(load_source_manifest()["sirolimus"]["retention_cross_environment"], load_source_manifest()["source_root"])
    if retention_path.exists():
        retention = pd.read_csv(retention_path)
        retention["gene"] = retention["gene"].astype(str).str.upper()
        retention = retention[["gene", "environments", "categories", "cross_environment_retention_candidate"]].drop_duplicates("gene")
        ranking = ranking.merge(retention, on="gene", how="left")
    output_dir = PROJECT_ROOT / "results" / "candidate_antigens"
    ranking.to_csv(ensure_output_path(output_dir / "candidate_antigen_ranking.csv"), index=False)
    pd.DataFrame(patient_rows).to_csv(ensure_output_path(output_dir / "candidate_antigen_by_patient.csv"), index=False)
    write_json(PROJECT_ROOT / "manifests" / "candidate_antigen_manifest.json", {
        "positive_control": "PMEL/gp100 is retained as an expression and prior-evidence positive control, not as proof of current HLA peptide presentation.",
        "normal_reference": "GSE122960 normal lung only; no whole-body normal tissue safety conclusion.",
        "denominator_rule": "not_assayed cells are excluded from detected/low/high fractions",
        "independence_rule": "LAM32 is deduplicated across GSE217108/GSE302356; core GSE135851 donor labels remain accession-scoped unless the donor registry confirms a match.",
        "outputs": [project_relative(output_dir / "candidate_antigen_ranking.csv"), project_relative(output_dir / "candidate_antigen_by_patient.csv")],
    })
    print(ranking[["gene", "patient_consistency_any_detected", "lam_detected_fraction_among_assayed", "normal_lung_detected_fraction_among_assayed"]].to_string(index=False))


if __name__ == "__main__":
    main()
