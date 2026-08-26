"""Generic-stress deconvolution, module enrichment, and human mapping.

The script works on the release-specific drug x gene summaries created by
``analyze_lincs_gene_programs.py``.  Generic stress programs are treated as
confounder axes, not as proof of mechanism.  GO/Reactome/MSigDB enrichment is
interpretive and uses the actual 204-gene common analyzable background.

GSE135851 and GSE302356 are mapped with their available metadata.  The former
has preliminary candidate/other labels; the latter has operational paper-
derived marker scores but no formal cell-state annotation in the staged raw
archives.  The output preserves those limitations.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.stats import fisher_exact
from sklearn.linear_model import Ridge
from statsmodels.stats.multitest import multipletests

from common import CANDIDATE_DECONVOLUTION, CANDIDATE_PROGRAMS, CANDIDATE_AUDIT, ROOT


STRESS_PROGRAMS = {
    "generic_apoptosis": {
        "up": ["BAX", "BBC3", "PMAIP1", "CASP3", "CASP7", "CASP8", "CASP9", "FAS", "TNFRSF10B", "APAF1", "DIABLO", "BCL2L11"],
        "down": ["BCL2", "BCL2L1", "MCL1", "XIAP", "CFLAR", "BIRC2", "BIRC3"],
    },
    "upr_er_stress": {
        "up": ["HSPA5", "DDIT3", "ATF4", "XBP1", "ATF6", "EIF2AK3", "ERN1", "DNAJB9", "HERPUD1", "MANF", "EDEM1", "SEL1L", "PPP1R15A", "TRIB3", "HYOU1", "CRELD2"],
        "down": [],
    },
    "heat_shock": {
        "up": ["HSPA1A", "HSPA1B", "HSP90AA1", "HSPH1", "DNAJB1", "DNAJB4", "HSPB1", "HSF1", "HSPA8", "HSPA9", "HSP90AB1"],
        "down": [],
    },
    "cell_cycle_arrest": {
        "up": ["CDKN1A", "CDKN1B", "GADD45A", "GADD45B", "BTG1", "BTG2", "SESN1", "SESN2", "DDIT4", "RPRM", "KLF4"],
        "down": ["MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "TOP2A", "CCNB1", "CCNE1", "CDK1", "CDC20", "PLK1", "BUB1"],
    },
    "proteostasis_translation_stress": {
        "up": ["SQSTM1", "NFE2L2", "ATF4", "DDIT3", "HSPA9", "GDF15", "LAMP3", "NDUFA4L2", "TXN", "PSME2"],
        "down": ["RPLP0", "RPLP1", "RPS6", "EIF4E", "EEF1A1", "EIF3A", "EIF4G1"],
    },
}


GMT_FILES = {
    "GO_Biological_Process_2023": "data/processed/gene_sets/GO_Biological_Process_2023.gmt",
    "Reactome_2022": "data/processed/gene_sets/Reactome_2022.gmt",
    "MSigDB_Hallmark_2020": "data/processed/gene_sets/MSigDB_Hallmark_2020.gmt",
}

SAMPLE_META = {
    "LAM3": {"modality": "scRNA-seq"},
    "LAM4": {"modality": "scRNA-seq"},
    "LAM18": {"modality": "Visium HD"},
    "LAM20": {"modality": "Visium"},
}


def read_10x_h5_local(path: Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Read a standard 10x filtered feature-barcode HDF5 without scanpy."""
    with h5py.File(path, "r") as h5:
        group = h5["matrix"] if "matrix" in h5 else h5["raw/matrix"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        shape = tuple(int(x) for x in group["shape"][:])
        matrix_gene_cell = sp.csc_matrix((data, indices, indptr), shape=shape)
        matrix_cell_gene = matrix_gene_cell.T.tocsr()
        features = group["features"]
        feature_key = "name" if "name" in features else "id"
        genes = [value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value) for value in features[feature_key][:]]
        barcodes = [value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value) for value in group["barcodes"][:]]
    return matrix_cell_gene, genes, barcodes


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    center = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - center)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale < 1e-8:
        scale = float(np.nanstd(values))
    if not np.isfinite(scale) or scale < 1e-8:
        scale = 1.0
    return (values - center) / scale


def load_common_genes() -> set[str]:
    # The comparison table also contains gene-level rows outside the primary
    # 300-gene query panel.  The enrichment/deconvolution background must be
    # the genes that are both in that panel and analyzable in both releases.
    audit = pd.read_csv(
        CANDIDATE_AUDIT / "LINCS_gene_panel_mapping_audit.csv",
        usecols=["gene", "dataset", "gene_available"],
    )
    available = audit.loc[audit["gene_available"].astype(bool)]
    counts = available.groupby("gene")["dataset"].nunique()
    return set(counts.index[counts.eq(2)].astype(str).str.upper())


def load_stable_modules(common_genes: set[str]) -> dict[str, dict]:
    matches = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_cross_release_module_matches.csv")
    matches = matches.loc[matches["same_module_first_pass"].eq(True)].copy()
    clusters = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_drug_gene_clusters.csv")
    modules = {}
    for idx, row in matches.reset_index(drop=True).iterrows():
        scope = row["scope"]
        module_id = f"{scope}__stable_{idx + 1}"
        genes = {gene for gene in str(row["common_genes"]).split(";") if gene and gene != "nan"} & common_genes
        left = clusters.loc[
            clusters["dataset"].eq("GSE92742")
            & clusters["scope"].eq(scope)
            & clusters["object_type"].eq("gene")
            & clusters["cluster_id"].eq(row["module_GSE92742"]),
            "object_id",
        ]
        right = clusters.loc[
            clusters["dataset"].eq("GSE70138")
            & clusters["scope"].eq(scope)
            & clusters["object_type"].eq("gene")
            & clusters["cluster_id"].eq(row["module_GSE70138"]),
            "object_id",
        ]
        modules[module_id] = {
            "module_id": module_id,
            "scope": scope,
            "genes": sorted(genes),
            "genes_GSE92742": sorted(set(left.astype(str)) & common_genes),
            "genes_GSE70138": sorted(set(right.astype(str)) & common_genes),
            "jaccard": float(row["jaccard"]),
            "overlap_coefficient": float(row["overlap_coefficient"]),
        }
    return modules


def template_matrix(common_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    vectors = {}
    audit = []
    for name, sets in STRESS_PROGRAMS.items():
        vector = pd.Series(0.0, index=common_genes)
        for gene in sets["up"]:
            if gene in vector.index:
                vector.loc[gene] += 1.0
        for gene in sets["down"]:
            if gene in vector.index:
                vector.loc[gene] -= 1.0
        vectors[name] = vector
        audit.append({
            "program": name,
            "n_up_in_common_panel": int(sum(gene in common_genes for gene in sets["up"])),
            "n_down_in_common_panel": int(sum(gene in common_genes for gene in sets["down"])),
            "n_nonzero_template_genes": int((vector != 0).sum()),
            "definition": "operational local generic stress reference; not a formal MSigDB signature",
        })
    matrix = pd.DataFrame(vectors, index=common_genes)
    return matrix, pd.DataFrame(audit)


def deconvolve(summary: pd.DataFrame, common_genes: set[str], modules: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = sorted(common_genes)
    templates, template_audit = template_matrix(genes)
    valid_programs = [column for column in templates.columns if np.linalg.norm(templates[column].to_numpy()) > 0]
    x_template = templates[valid_programs].to_numpy(float)
    x_template = x_template / np.maximum(np.linalg.norm(x_template, axis=0, keepdims=True), 1e-8)
    compounds = summary.loc[
        summary["entity_type"].eq("compound") & summary["gene"].isin(common_genes)
    ].copy()
    effect_rows = []
    burden_rows = []
    module_rows = []
    coefficient_rows = []
    for dataset, group in compounds.groupby("dataset", sort=True):
        pivot = group.pivot_table(index="entity_id", columns="gene", values="drug_effect_median", aggfunc="first").reindex(columns=genes)
        for entity_id, vector_series in pivot.iterrows():
            vector = vector_series.to_numpy(float)
            valid = np.isfinite(vector)
            if valid.sum() < 20:
                continue
            z = robust_z(vector)
            model = Ridge(alpha=1.0, fit_intercept=True)
            model.fit(x_template[valid], z[valid])
            predicted = model.predict(x_template)
            residual = z - predicted
            total_ss = float(np.sum((z[valid] - np.mean(z[valid])) ** 2))
            residual_ss = float(np.sum((residual[valid] - np.mean(residual[valid])) ** 2))
            r2 = 1.0 - residual_ss / total_ss if total_ss > 0 else np.nan
            meta = group.loc[group["entity_id"].eq(entity_id)].iloc[0]
            burden_rows.append({
                "dataset": dataset,
                "entity_id": entity_id,
                "pert_iname": meta["pert_iname"],
                "scope": meta["scope"],
                "generic_component_r2": r2,
                "generic_component_rms": float(np.sqrt(np.mean(predicted[valid] ** 2))),
                "generic_stress_burden": float(np.mean(np.abs(predicted[valid]))),
                "residual_rms": float(np.sqrt(np.mean(residual[valid] ** 2))),
                "n_genes_used": int(valid.sum()),
                "stress_reference": "operational local apoptosis/UPR/heat-shock/cell-cycle/proteostasis templates",
            })
            for program, coefficient in zip(valid_programs, model.coef_):
                coefficient_rows.append({
                    "dataset": dataset,
                    "entity_id": entity_id,
                    "pert_iname": meta["pert_iname"],
                    "scope": meta["scope"],
                    "program": program,
                    "coefficient": float(coefficient),
                })
            for gene, raw_effect, effect_z, generic_component, adjusted_z, weight in zip(genes, vector, z, predicted, residual, group.loc[group["entity_id"].eq(entity_id)].set_index("gene").reindex(genes)["disease_weight"].to_numpy(float)):
                before_sign = weight * raw_effect if np.isfinite(raw_effect) else np.nan
                after_sign = weight * adjusted_z if np.isfinite(adjusted_z) else np.nan
                effect_rows.append({
                    "dataset": dataset,
                    "entity_id": entity_id,
                    "pert_iname": meta["pert_iname"],
                    "scope": meta["scope"],
                    "gene": gene,
                    "disease_weight": weight,
                    "drug_effect_median": raw_effect,
                    "effect_robust_z": effect_z,
                    "generic_component": generic_component,
                    "stress_adjusted_effect_z": adjusted_z,
                    "direction_before": "reversal" if before_sign < 0 else ("mimic" if before_sign > 0 else "neutral"),
                    "direction_after": "reversal" if after_sign < 0 else ("mimic" if after_sign > 0 else "neutral"),
                })
            for module_id, module in modules.items():
                module_genes = [gene for gene in module["genes"] if gene in genes]
                if not module_genes:
                    continue
                indices = [genes.index(gene) for gene in module_genes]
                before_signs = vector[indices] * np.asarray([group.loc[(group.entity_id == entity_id) & (group.gene == gene), "disease_weight"].iloc[0] for gene in module_genes])
                after_signs = residual[indices] * np.asarray([group.loc[(group.entity_id == entity_id) & (group.gene == gene), "disease_weight"].iloc[0] for gene in module_genes])
                wanted = before_signs < 0 if module["scope"] == "reversal_only" else before_signs > 0
                wanted_after = after_signs < 0 if module["scope"] == "reversal_only" else after_signs > 0
                before_fraction = float(np.mean(wanted))
                after_fraction = float(np.mean(wanted_after))
                module_rows.append({
                    "dataset": dataset,
                    "entity_id": entity_id,
                    "pert_iname": meta["pert_iname"],
                    "scope": meta["scope"],
                    "module_id": module_id,
                    "module_scope": module["scope"],
                    "n_module_genes": len(module_genes),
                    "before_direction_fraction": before_fraction,
                    "after_stress_adjusted_direction_fraction": after_fraction,
                    "direction_fraction_retained": after_fraction / before_fraction if before_fraction > 0 else np.nan,
                    "n_before_direction": int(wanted.sum()),
                    "n_after_direction": int(wanted_after.sum()),
                    "module_jaccard_cross_release": module["jaccard"],
                })
    return pd.DataFrame(effect_rows), pd.DataFrame(burden_rows), pd.DataFrame(module_rows), pd.concat([template_audit, pd.DataFrame(coefficient_rows)], axis=0, ignore_index=True)


def parse_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = {gene.upper() for gene in parts[2:] if gene}
            if genes:
                gene_sets[name] = genes
    return gene_sets


def enrich_modules(modules: dict[str, dict], common_genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    audit = []
    background = set(common_genes)
    for library, relative_path in GMT_FILES.items():
        path = ROOT / relative_path
        if not path.exists():
            audit.append({"library": library, "status": "not_available", "path": relative_path})
            continue
        gene_sets = parse_gmt(path)
        audit.append({"library": library, "status": "loaded", "path": relative_path, "n_terms": len(gene_sets), "background_size": len(background)})
        for module_id, module in modules.items():
            foreground = set(module["genes"]) & background
            for term, term_genes in gene_sets.items():
                term_background = term_genes & background
                overlap = foreground & term_background
                if len(term_background) < 3 or not overlap:
                    continue
                a = len(overlap)
                b = len(foreground - term_background)
                c = len(term_background - foreground)
                d = len(background - foreground - term_background)
                _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
                rows.append({
                    "module_id": module_id,
                    "module_scope": module["scope"],
                    "module_size": len(foreground),
                    "library": library,
                    "term": term,
                    "term_size_in_background": len(term_background),
                    "overlap_genes": len(overlap),
                    "leading_genes": ";".join(sorted(overlap)),
                    "p_value": float(p_value),
                    "background_size": len(background),
                    "background_definition": "204 common analyzable genes entering cross-release module comparison",
                })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = np.nan
        for (module_id, library), indices in result.groupby(["module_id", "library"]).groups.items():
            result.loc[indices, "fdr"] = multipletests(result.loc[indices, "p_value"], method="fdr_bh")[1]
        result["fdr_significant"] = result["fdr"] < 0.05
        result = result.sort_values(["module_id", "library", "fdr", "p_value", "overlap_genes"], ascending=[True, True, True, True, False])
    return result, pd.DataFrame(audit)


def group_module_scores(expression: pd.DataFrame, meta: pd.DataFrame, modules: dict[str, dict], dataset: str) -> pd.DataFrame:
    ranks = expression.rank(axis=0, pct=True)
    rows = []
    for group in expression.index:
        group_meta = meta.loc[meta["group"].eq(group)].iloc[0]
        base = {
            "dataset": dataset,
            "group": group,
            "patient_id": group_meta.get("patient_id", ""),
            "state_label": group_meta.get("state_label", ""),
            "n_cells": int(group_meta.get("n_cells", 0)),
        }
        for module_id, module in modules.items():
            genes = [gene for gene in module["genes"] if gene in ranks.columns]
            base[f"{module_id}_rank_score"] = float(ranks.loc[group, genes].mean()) if genes else np.nan
            base[f"{module_id}_n_genes"] = len(genes)
        rows.append(base)
    return pd.DataFrame(rows)


def map_gse135851(modules: dict[str, dict]) -> pd.DataFrame:
    adata = ad.read_h5ad(ROOT / "data/processed/GSE135851/GSE135851_lam_states_snapshot.h5ad", backed="r")
    symbols = pd.Series(adata.var_names.astype(str).str.upper())
    lookup = {}
    for idx, symbol in symbols.items():
        lookup.setdefault(symbol, []).append(idx)
    genes = sorted({gene for module in modules.values() for gene in module["genes"]})
    selected_symbols = [gene for gene in genes if gene in lookup]
    indices = [idx for gene in selected_symbols for idx in lookup[gene]]
    obs = adata.obs[["donor_id", "lamcore_label", "analysis_pass"]].copy()
    obs["analysis_pass"] = obs["analysis_pass"].astype(bool)
    obs = obs.loc[obs["analysis_pass"]].copy()
    obs["group"] = obs["donor_id"].astype(str) + "|" + obs["lamcore_label"].astype(str)
    counts = obs.groupby("group").size().rename("n_cells")
    keep = counts[counts >= 20].index
    obs = obs.loc[obs["group"].isin(keep)]
    x = adata[obs.index, indices].X
    if sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=float)
    group_codes, groups = pd.factorize(obs["group"], sort=True)
    means = np.zeros((len(groups), len(selected_symbols)))
    for group_idx in range(len(groups)):
        group_values = x[group_codes == group_idx]
        for symbol_idx, symbol in enumerate(selected_symbols):
            feature_cols = [j for j, idx in enumerate(indices) if symbols.iloc[idx] == symbol]
            means[group_idx, symbol_idx] = group_values[:, feature_cols].mean()
    expression = pd.DataFrame(means, index=groups, columns=selected_symbols)
    meta = obs.drop_duplicates("group").set_index("group").loc[groups, ["donor_id", "lamcore_label"]].reset_index()
    if "group" not in meta.columns:
        meta = meta.rename(columns={meta.columns[0]: "group"})
    meta = meta.rename(columns={"donor_id": "patient_id", "lamcore_label": "state_label"})
    meta["n_cells"] = counts.loc[groups].to_numpy()
    return group_module_scores(expression, meta, modules, "GSE135851")


def process_raw_gse302356_sample(sample_id: str, modules: dict[str, dict]) -> dict:
    path = ROOT / "data/raw/GSE302356/unpacked" / sample_id / f"{sample_id}_filtered_feature_bc_matrix.h5"
    x, var_names, obs_names = read_10x_h5_local(path)
    x = x.astype(np.float32)
    total = np.asarray(x.sum(axis=1)).ravel()
    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    valid = total > 0
    if SAMPLE_META[sample_id]["modality"].startswith("Visium"):
        positions = ROOT / "data/raw/GSE302356/unpacked" / sample_id / f"{sample_id}_tissue_positions.parquet"
        if not positions.exists():
            positions = positions.with_suffix(".csv")
        table = pd.read_parquet(positions) if positions.suffix == ".parquet" else pd.read_csv(positions)
        tissue = set(table.loc[table["in_tissue"].astype(int).eq(1), "barcode"].astype(str))
        valid &= np.array([str(barcode) in tissue for barcode in obs_names])
    x = x[valid]
    total = total[valid]
    x = sp.diags(1e4 / total) @ x
    x = x.tocsr()
    x.data = np.log1p(x.data)
    lookup = {gene.upper(): idx for idx, gene in enumerate(var_names)}
    rows = []
    for module_id, module in modules.items():
        indices = [lookup[gene] for gene in module["genes"] if gene in lookup]
        score = np.asarray(x[:, indices].mean(axis=1)).ravel() if indices else np.full(x.shape[0], np.nan)
        rows.append({
            "sample_id": sample_id,
            "modality": SAMPLE_META[sample_id]["modality"],
            "module_id": module_id,
            "module_score_mean_log1p": float(np.nanmean(score)),
            "module_score_median_log1p": float(np.nanmedian(score)),
            "module_score_p90_log1p": float(np.nanquantile(score, 0.9)),
            "n_genes_available": len(indices),
            "n_cells_or_spots": int(x.shape[0]),
        })
    return rows


def map_gse302356(modules: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for sample_id in SAMPLE_META:
        rows.extend(process_raw_gse302356_sample(sample_id, modules))
    result = pd.DataFrame(rows)
    state_scores = pd.read_csv(ROOT / "results/human_mapping/GSE302356_raw_sample_scores.csv")
    state_cols = [
        column
        for column in state_scores.columns
        if (column.startswith("state_") and not column.startswith("state_signature"))
        or column.startswith("fraction_z_state_")
    ]
    state_subset = state_scores[["sample_id", "modality", *state_cols]].copy()
    return result.merge(state_subset, on=["sample_id", "modality"], how="left")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert robust_z(np.array([1.0, 2.0, 3.0])).shape == (3,)
        assert "generic_apoptosis" in STRESS_PROGRAMS
        print("self_test: PASS")
        return

    out_dir = CANDIDATE_DECONVOLUTION
    out_dir.mkdir(parents=True, exist_ok=True)
    common_genes = load_common_genes()
    modules = load_stable_modules(common_genes)
    summary = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_drug_gene_response_summary.csv.gz")
    adjusted, burden, module_retention, coefficient_audit = deconvolve(summary, common_genes, modules)
    enrichment, enrichment_audit = enrich_modules(modules, common_genes)
    human_gse135851 = map_gse135851(modules)
    human_gse302356 = map_gse302356(modules)

    adjusted.to_csv(out_dir / "LINCS_stress_adjusted_gene_effects.csv.gz", index=False, compression="gzip")
    burden.to_csv(out_dir / "LINCS_generic_stress_burden.csv", index=False)
    module_retention.to_csv(out_dir / "LINCS_module_stress_adjusted_retention.csv", index=False)
    coefficient_audit.to_csv(out_dir / "LINCS_generic_stress_program_coefficients_and_audit.csv", index=False)
    enrichment.to_csv(out_dir / "LINCS_stable_module_GO_Reactome_MSigDB_enrichment.csv.gz", index=False, compression="gzip")
    enrichment_audit.to_csv(out_dir / "LINCS_enrichment_library_audit.csv", index=False)
    human_gse135851.to_csv(ROOT / "results/human_mapping/GSE135851_stable_module_scores.csv", index=False)
    human_gse302356.to_csv(ROOT / "results/human_mapping/GSE302356_stable_module_sample_scores.csv", index=False)

    manifest = {
        "common_gene_background": len(common_genes),
        "stable_modules": {module_id: {"scope": module["scope"], "n_genes": len(module["genes"]), "jaccard": module["jaccard"]} for module_id, module in modules.items()},
        "generic_stress_reference": "operational local signed templates; deconvolution by ridge regression on robust-z drug effects",
        "generic_programs": list(STRESS_PROGRAMS),
        "enrichment_background": "204 genes shared by GSE92742/GSE70138 and used for cross-release modules",
        "enrichment_libraries": list(GMT_FILES),
        "human_GSE135851": "preliminary candidate/other labels; patient-state pseudobulk; not formal LAMCORE/LAF taxonomy",
        "human_GSE302356": "sample/modality-level mapping using raw matrices and paper-derived marker-score columns; no formal state labels in staged archives",
        "interpretation": "generic stress burden and adjusted module retention are hypothesis-generating; no cytotoxicity cutoff is used as a hard candidate filter in this pass",
    }
    (ROOT / "manifests" / "LINCS_deconvolution_human_mapping_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
