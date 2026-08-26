"""Analyze GSE179044's 2x2x2 design and produce signed residual tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
import yaml

from common import ROOT, classify_ratio, collapse_duplicate_genes, signed_ratio, write_json


def load_gse179044(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    sample_cols = [column for column in raw.columns if str(column).startswith("Pietrobon_S")]
    if len(sample_cols) != 16:
        raise ValueError(f"Expected 16 Pietrobon sample columns, found {len(sample_cols)}")
    gene_column = "gene_name" if "gene_name" in raw.columns else "gene_id"
    table = raw.set_index(gene_column)[sample_cols]
    table = collapse_duplicate_genes(table)
    library_size = table.sum(axis=0).replace(0, np.nan)
    return np.log2(table.divide(library_size, axis=1) * 1_000_000 + 1.0)


def sample_table() -> pd.DataFrame:
    rows = []
    for rep in (1, 2):
        rows += [
            {"sample_id": f"Pietrobon_S{1 + (rep - 1) * 8:02d}", "genotype": "WT", "treatment": "vehicle", "environment": "hydrogel", "replicate": rep},
            {"sample_id": f"Pietrobon_S{2 + (rep - 1) * 8:02d}", "genotype": "KO", "treatment": "vehicle", "environment": "hydrogel", "replicate": rep},
            {"sample_id": f"Pietrobon_S{3 + (rep - 1) * 8:02d}", "genotype": "WT", "treatment": "rapamycin", "environment": "hydrogel", "replicate": rep},
            {"sample_id": f"Pietrobon_S{4 + (rep - 1) * 8:02d}", "genotype": "KO", "treatment": "rapamycin", "environment": "hydrogel", "replicate": rep},
            {"sample_id": f"Pietrobon_S{5 + (rep - 1) * 8:02d}", "genotype": "WT", "treatment": "vehicle", "environment": "plastic", "replicate": rep},
            {"sample_id": f"Pietrobon_S{6 + (rep - 1) * 8:02d}", "genotype": "KO", "treatment": "vehicle", "environment": "plastic", "replicate": rep},
            {"sample_id": f"Pietrobon_S{7 + (rep - 1) * 8:02d}", "genotype": "WT", "treatment": "rapamycin", "environment": "plastic", "replicate": rep},
            {"sample_id": f"Pietrobon_S{8 + (rep - 1) * 8:02d}", "genotype": "KO", "treatment": "rapamycin", "environment": "plastic", "replicate": rep},
        ]
    return pd.DataFrame(rows)


def group_mean(table: pd.DataFrame, samples: pd.DataFrame, **filters: str) -> pd.Series:
    ids = samples.loc[np.logical_and.reduce([samples[key].eq(value).to_numpy() for key, value in filters.items()]), "sample_id"]
    return table.loc[:, list(ids)].mean(axis=1)


def simple_interactions(table: pd.DataFrame, samples: pd.DataFrame) -> dict[str, pd.Series]:
    result: dict[str, pd.Series] = {}
    for environment in ("plastic", "hydrogel"):
        d_vehicle = group_mean(table, samples, environment=environment, treatment="vehicle", genotype="KO") - group_mean(table, samples, environment=environment, treatment="vehicle", genotype="WT")
        d_rapa = group_mean(table, samples, environment=environment, treatment="rapamycin", genotype="KO") - group_mean(table, samples, environment=environment, treatment="rapamycin", genotype="WT")
        ko_response = group_mean(table, samples, environment=environment, treatment="rapamycin", genotype="KO") - group_mean(table, samples, environment=environment, treatment="vehicle", genotype="KO")
        wt_response = group_mean(table, samples, environment=environment, treatment="rapamycin", genotype="WT") - group_mean(table, samples, environment=environment, treatment="vehicle", genotype="WT")
        result[f"tsc2_loss_{environment}"] = d_vehicle
        result[f"residual_{environment}"] = d_rapa
        result[f"escape_{environment}"] = ko_response - wt_response
    result["hydrogel_specific_residual"] = result["residual_hydrogel"] - result["residual_plastic"]
    result["environment_dependent_escape"] = result["escape_hydrogel"] - result["escape_plastic"]
    return result


def cell_vector(genotype: int, treatment: int, environment: int) -> np.ndarray:
    """Treatment-coded vector for one cell of the saturated 2x2x2 model."""
    g, r, e = genotype, treatment, environment
    return np.array([1.0, g, r, e, g * r, g * e, r * e, g * r * e])


def moderated_statistics(table: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Fit the full factorial model and moderate gene-wise residual variances.

    The two biological replicates per cell give only eight residual degrees of
    freedom. We therefore retain the ordinary least-squares effect estimates,
    but use an empirical-Bayes-style pooled prior for the residual variance to
    stabilize standard errors, t statistics and FDR values. Ratios remain
    descriptive and are never used as p-values.
    """
    design_rows = []
    for row in samples.itertuples(index=False):
        design_rows.append(cell_vector(
            int(row.genotype == "KO"),
            int(row.treatment == "rapamycin"),
            int(row.environment == "hydrogel"),
        ))
    x = np.vstack(design_rows)
    y = table.loc[:, samples.sample_id].to_numpy(dtype=float).T
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = (xtx_inv @ x.T @ y).T
    residual = y - x @ beta.T
    df_resid = x.shape[0] - x.shape[1]
    s2 = np.sum(residual ** 2, axis=0) / df_resid
    prior_df = 4.0
    prior_s2 = float(np.nanmedian(s2))
    moderated_s2 = (prior_df * prior_s2 + df_resid * s2) / (prior_df + df_resid)

    vectors: dict[str, np.ndarray] = {}
    for environment_name, e in (("plastic", 0), ("hydrogel", 1)):
        d0 = cell_vector(1, 0, e) - cell_vector(0, 0, e)
        d1 = cell_vector(1, 1, e) - cell_vector(0, 1, e)
        escape = d1 - d0
        vectors[f"tsc2_loss_{environment_name}"] = d0
        vectors[f"residual_{environment_name}"] = d1
        vectors[f"escape_{environment_name}"] = escape
    vectors["hydrogel_specific_residual"] = vectors["residual_hydrogel"] - vectors["residual_plastic"]
    vectors["environment_dependent_escape"] = vectors["escape_hydrogel"] - vectors["escape_plastic"]

    index = table.index
    result = pd.DataFrame(index=index)
    contrast_cov = xtx_inv
    for name, vector in vectors.items():
        effect = beta @ vector
        se = np.sqrt(np.maximum(0.0, moderated_s2 * (vector @ contrast_cov @ vector)))
        t_value = np.divide(effect, se, out=np.full_like(effect, np.nan), where=se > 0)
        p_value = 2.0 * student_t.sf(np.abs(t_value), df=df_resid + prior_df)
        q_value = multipletests(np.nan_to_num(p_value, nan=1.0), method="fdr_bh")[1]
        result[name] = effect
        result[f"{name}_moderated_se"] = se
        result[f"{name}_moderated_t"] = t_value
        result[f"{name}_moderated_p"] = p_value
        result[f"{name}_moderated_q"] = q_value
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/raw/GSE179044/GSE179044_raw_counts_matrix.csv.gz")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    input_path = ROOT / args.input
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    table = load_gse179044(input_path)
    samples = sample_table()
    if set(samples.sample_id) != set(table.columns):
        raise ValueError(f"Sample identifiers do not match GSE179044 matrix: {list(table.columns)}")
    table = table.loc[:, samples.sample_id]
    contrasts = simple_interactions(table, samples)
    output = moderated_statistics(table, samples)
    # Keep the transparent group-mean contrasts as a consistency check. In a
    # balanced saturated design they equal the model estimates up to floating
    # point precision.
    mean_contrasts = pd.DataFrame(contrasts)
    for column in mean_contrasts.columns:
        output[f"{column}_group_mean"] = mean_contrasts[column]
    baseline = output["tsc2_loss_hydrogel"]
    rapa = output["residual_hydrogel"]
    output["signed_residual_ratio_hydrogel"] = signed_ratio(rapa, baseline, float(config["min_baseline_effect"]))
    output["absolute_residual_ratio_hydrogel"] = output["signed_residual_ratio_hydrogel"].abs()
    output["residual_class_hydrogel"] = classify_ratio(output["signed_residual_ratio_hydrogel"])
    baseline_p = output["tsc2_loss_plastic"]
    rapa_p = output["residual_plastic"]
    output["signed_residual_ratio_plastic"] = signed_ratio(rapa_p, baseline_p, float(config["min_baseline_effect"]))
    output["absolute_residual_ratio_plastic"] = output["signed_residual_ratio_plastic"].abs()
    output["residual_class_plastic"] = classify_ratio(output["signed_residual_ratio_plastic"])
    output.index.name = "gene"
    result_dir = ROOT / "results" / "tables"
    result_dir.mkdir(parents=True, exist_ok=True)
    output.to_csv(result_dir / "GSE179044_factorial_contrasts.csv")
    summary = {
        "dataset": "GSE179044",
        "shape": list(table.shape),
        "n_genes": int(table.shape[0]),
        "n_samples": int(table.shape[1]),
        "contrast_columns": list(output.columns),
        "notes": [
            "GSE179044 contrasts are orthogonal comparisons within one 2x2x2 experiment, not independent replication.",
            "Hydrogel-specific residual is rapamycin-condition GxE; environment-dependent escape is GxRxE.",
            "GxR is genotype-dependent rapamycin response and is not automatically escape.",
            "Moderated statistics stabilize uncertainty using a pooled residual-variance prior; ratios are descriptive and gated by min_baseline_effect.",
        ],
    }
    write_json(ROOT / "manifests" / "GSE179044_analysis.json", summary)
    print(summary)


if __name__ == "__main__":
    main()
