"""Evaluate locally measured pharmacological references for target modules.

These compounds are not assumed to be selective.  The output explicitly
labels their multi-target limitations and asks the useful first question:
does a pharmacological perturbation reproduce the stable gene module?
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from analyze_lincs_gene_programs import (
    DATASETS,
    load_gene_mapping,
    load_panel,
    load_sig_info,
    summarize_release,
)
from common import CANDIDATE_PROGRAMS, CANDIDATE_VALIDATION, ROOT


REFERENCE_COMPOUNDS = {
    "ibrutinib": {
        "axis": "BTK",
        "selectivity_note": "BTK-directed pharmacology but known off-target kinases; not a clean BTK-only perturbation",
    },
    "crizotinib": {
        "axis": "NTRK/ALK/MET",
        "selectivity_note": "multi-target kinase inhibitor; NTRK3 relevance is indirect, not a clean NTRK3-only perturbation",
    },
    "cabozantinib": {
        "axis": "RET/MET/VEGFR",
        "selectivity_note": "multi-target RET/MET/VEGFR pharmacology; not a clean RET-only perturbation",
    },
    "vandetanib": {
        "axis": "RET/EGFR/VEGFR",
        "selectivity_note": "multi-target RET/EGFR/VEGFR pharmacology; not a clean RET-only perturbation",
    },
}


def load_reference_entities() -> pd.DataFrame:
    rows = []
    for name, info in REFERENCE_COMPOUNDS.items():
        rows.append({
            "entity_id": "reference::" + name,
            "entity_type": "compound_reference",
            "pert_iname": name,
            "scope": "selective_reference",
            "axis": info["axis"],
            "selectivity_note": info["selectivity_note"],
        })
    return pd.DataFrame(rows)


def load_stable_modules() -> dict[str, dict]:
    matches = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_cross_release_module_matches.csv")
    matches = matches.loc[matches["same_module_first_pass"].eq(True)].reset_index(drop=True)
    modules = {}
    for idx, row in matches.iterrows():
        module_id = f"{row['scope']}__stable_{idx + 1}"
        modules[module_id] = {
            "module_id": module_id,
            "scope": row["scope"],
            "genes": sorted(set(str(row["common_genes"]).split(";")) - {"", "nan"}),
            "jaccard": float(row["jaccard"]),
        }
    return modules


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-size", type=int, default=256)
    args = parser.parse_args()
    panel = load_panel(ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv")
    entities = load_reference_entities()
    modules = load_stable_modules()
    all_summary = []
    audits = []
    for dataset in DATASETS:
        mapped, _, _, _ = load_gene_mapping(dataset, panel)
        selected, _ = load_sig_info(dataset, entities.rename(columns={"pert_iname": "pert_iname", "scope": "scope"}), {})
        summary, audit = summarize_release(dataset, mapped, selected, entities, args.chunk_size)
        all_summary.append(summary)
        audits.append(audit.assign(dataset=dataset))
    summary = pd.concat(all_summary, ignore_index=True)
    audit = pd.concat(audits, ignore_index=True)
    modules_rows = []
    for dataset, group in summary.groupby("dataset", sort=True):
        for reference in entities.itertuples():
            entity_id = reference.entity_id
            profile = group.loc[group["entity_id"].eq(entity_id)].set_index("gene")
            n_signatures = int(audit.loc[(audit.dataset == dataset) & (audit.entity_id == entity_id), "n_selected_signatures"].sum()) if not audit.empty else 0
            for module_id, module in modules.items():
                genes = [gene for gene in module["genes"] if gene in profile.index]
                if n_signatures == 0 or not genes:
                    modules_rows.append({
                        "dataset": dataset,
                        "reference_compound": reference.pert_iname,
                        "axis": reference.axis,
                        "selectivity_note": reference.selectivity_note,
                        "module_id": module_id,
                        "module_scope": module["scope"],
                        "n_module_genes_compared": 0,
                        "n_selected_signatures": n_signatures,
                        "data_status": "not_available",
                        "reference_reversal_fraction": np.nan,
                        "reference_mimic_fraction": np.nan,
                        "reference_module_direction_fraction": np.nan,
                        "module_jaccard_cross_release": module["jaccard"],
                        "interpretation": "reference not measured or no analyzable module genes; not a negative biological result",
                    })
                    continue
                weights = profile.loc[genes, "disease_weight"].to_numpy(float)
                effects = profile.loc[genes, "drug_effect_median"].to_numpy(float)
                signed = weights * effects
                reversal_fraction = float(np.mean(signed < 0))
                mimic_fraction = float(np.mean(signed > 0))
                wanted = reversal_fraction if module["scope"] == "reversal_only" else mimic_fraction
                modules_rows.append({
                    "dataset": dataset,
                    "reference_compound": reference.pert_iname,
                    "axis": reference.axis,
                    "selectivity_note": reference.selectivity_note,
                    "module_id": module_id,
                    "module_scope": module["scope"],
                    "n_module_genes_compared": len(genes),
                    "n_selected_signatures": n_signatures,
                    "data_status": "not_available" if n_signatures == 0 else ("available_but_weak" if n_signatures < 3 else "available"),
                    "reference_reversal_fraction": reversal_fraction,
                    "reference_mimic_fraction": mimic_fraction,
                    "reference_module_direction_fraction": wanted,
                    "module_jaccard_cross_release": module["jaccard"],
                    "interpretation": "pharmacological module reproduction reference; not independent single-target proof",
                })
    availability = entities[["entity_id", "pert_iname", "axis", "selectivity_note"]].merge(
        audit.groupby(["dataset", "entity_id"], as_index=False).agg(
            n_selected_signatures=("n_selected_signatures", "sum"),
            n_selected_cells=("n_selected_cells", "sum"),
            perturbation_types=("perturbation_types", lambda x: ";".join(sorted(set(";".join(x.astype(str)).split(";"))))),
        ),
        on="entity_id",
        how="outer",
    )
    availability["measured"] = availability["n_selected_signatures"].fillna(0).gt(0)
    out = CANDIDATE_VALIDATION
    out.mkdir(parents=True, exist_ok=True)
    availability.to_csv(out / "LINCS_selective_reference_availability.csv", index=False)
    pd.DataFrame(modules_rows).to_csv(out / "LINCS_selective_reference_module_reproduction.csv", index=False)
    summary.to_csv(out / "LINCS_selective_reference_gene_summary.csv.gz", index=False, compression="gzip")
    (ROOT / "manifests" / "LINCS_selective_reference_manifest.json").write_text(json.dumps({
        "references": REFERENCE_COMPOUNDS,
        "interpretation": "Local LINCS pharmacological references are multi-target and do not replace selective drug or independent CRISPR/shRNA validation.",
        "priority_axes": ["RET", "BTK", "NTRK3", "MKNK1/2"],
    }, indent=2))
    print(availability.to_string(index=False))


if __name__ == "__main__":
    main()
