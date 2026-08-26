"""Build non-redundant query signatures and downstream functional summaries.

Experimental contrasts become candidate CMap/LINCS query objects. Functional
modules are summarized separately as interpretation layers, so the same genes
are not counted once as a signature and again as an independent signature.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from common import ROOT, write_json


QUERY_CONTRASTS = [
    "tsc2_loss_plastic",
    "tsc2_loss_hydrogel",
    "residual_plastic",
    "residual_hydrogel",
    "escape_plastic",
    "escape_hydrogel",
    "hydrogel_specific_residual",
    "environment_dependent_escape",
]


def top_signed_genes(table: pd.DataFrame, contrast: str, fdr: float, min_effect: float, top_n: int) -> pd.DataFrame:
    q_col = f"{contrast}_moderated_q"
    eligible = table.loc[
        table[q_col].le(fdr) & table[contrast].abs().ge(min_effect),
        [contrast, q_col],
    ].dropna()
    up = eligible.loc[eligible[contrast] > 0].sort_values([q_col, contrast], ascending=[True, False]).head(top_n)
    down = eligible.loc[eligible[contrast] < 0].sort_values([q_col, contrast], ascending=[True, True]).head(top_n)
    pieces = []
    for direction, subset in (("up", up), ("down", down)):
        if subset.empty:
            continue
        part = subset.reset_index().rename(columns={"index": "gene", contrast: "signed_score", q_col: "moderated_q"})
        part["direction"] = direction
        part["contrast"] = contrast
        part["default_cmap_query"] = True
        pieces.append(part[["contrast", "gene", "direction", "signed_score", "moderated_q", "default_cmap_query"]])
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def module_summary(table: pd.DataFrame, module_sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for contrast in QUERY_CONTRASTS:
        for module, genes in module_sets.items():
            available = [gene for gene in genes if gene in table.index]
            values = table.loc[available, contrast] if available else pd.Series(dtype=float)
            q_col = f"{contrast}_moderated_q"
            q_values = table.loc[available, q_col] if available else pd.Series(dtype=float)
            rows.append({
                "contrast": contrast,
                "module": module,
                "n_genes": len(available),
                "mean_effect": float(values.mean()) if len(values) else np.nan,
                "median_effect": float(values.median()) if len(values) else np.nan,
                "n_fdr_supported": int((q_values <= 0.1).sum()) if len(q_values) else 0,
                "fraction_fdr_supported": float((q_values <= 0.1).mean()) if len(q_values) else np.nan,
            })
    return pd.DataFrame(rows)


def direction_reversal_diagnostics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for environment in ("plastic", "hydrogel"):
        ratio = table[f"signed_residual_ratio_{environment}"]
        baseline = table[f"tsc2_loss_{environment}"]
        residual = table[f"residual_{environment}"]
        rows.append({
            "environment": environment,
            "n_ratio_eligible": int(ratio.notna().sum()),
            "n_direction_reversal": int((ratio < 0).sum()),
            "n_reversal_with_baseline_q_le_0.1": int(((ratio < 0) & (table[f"tsc2_loss_{environment}_moderated_q"] <= 0.1)).sum()),
            "n_reversal_with_residual_q_le_0.1": int(((ratio < 0) & (table[f"residual_{environment}_moderated_q"] <= 0.1)).sum()),
            "median_abs_baseline": float(baseline.abs().median()),
            "median_abs_residual": float(residual.abs().median()),
        })
    return pd.DataFrame(rows)


def cross_model_overlap(factorial: pd.DataFrame, external: pd.DataFrame, top_n: int = 300) -> pd.DataFrame:
    """Compare top signed effects after a conservative case-normalized symbol join.

    This is supportive only: human-cell and mouse-MEF experiments are not
    treated as interchangeable or as a formal meta-analysis.
    """
    f = factorial.copy()
    e = external.copy()
    f["join_gene"] = f.index.astype(str).str.upper()
    e["join_gene"] = e.index.astype(str).str.upper()
    rows = []
    for contrast in ("tsc2_loss_hydrogel", "residual_hydrogel", "escape_hydrogel"):
        f_top = f.loc[
            f[f"{contrast}_moderated_q"].le(0.1),
            ["join_gene", contrast],
        ].dropna().sort_values(contrast, key=lambda s: s.abs(), ascending=False).head(top_n)
        if contrast == "tsc2_loss_hydrogel":
            e_col = "tsc2_loss"
        elif contrast == "residual_hydrogel":
            e_col = "residual_after_rapamycin"
        else:
            e_col = "genotype_dependent_rapamycin_response"
        e_top = e[["join_gene", e_col]].dropna().sort_values(e_col, key=lambda s: s.abs(), ascending=False).head(top_n)
        merged = f_top.merge(e_top, on="join_gene", suffixes=("_human", "_mouse"))
        if merged.empty:
            continue
        rows.append({
            "factorial_contrast": contrast,
            "external_contrast": e_col,
            "n_common_top_features": int(len(merged)),
            "sign_concordance": float((np.sign(merged[contrast]) == np.sign(merged[e_col])).mean()),
            "spearman_like_rank_correlation": float(merged[contrast].rank().corr(merged[e_col].rank())),
        })
    return pd.DataFrame(rows)


def main() -> None:
    config = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    table = pd.read_csv(ROOT / "results/tables/GSE179044_factorial_contrasts.csv", index_col="gene")
    module_sets = {name: [str(gene) for gene in genes] for name, genes in config["module_sets"].items()}
    signature_dir = ROOT / "results" / "signatures"
    signature_dir.mkdir(parents=True, exist_ok=True)

    signatures = [
        top_signed_genes(table, contrast, float(config["fdr_cutoff"]), float(config["min_baseline_effect"]), 150)
        for contrast in QUERY_CONTRASTS
    ]
    signatures = [part for part in signatures if not part.empty]
    signature_table = pd.concat(signatures, ignore_index=True) if signatures else pd.DataFrame()
    signature_table.to_csv(signature_dir / "GSE179044_cmap_query_signatures.csv", index=False)

    module_summary(table, module_sets).to_csv(ROOT / "results/tables/GSE179044_functional_module_summary.csv", index=False)
    direction_reversal_diagnostics(table).to_csv(ROOT / "results/tables/GSE179044_direction_reversal_diagnostics.csv", index=False)

    external_path = ROOT / "results/tables/GSE27982_external_response.csv"
    overlap = cross_model_overlap(table, pd.read_csv(external_path, index_col="probe_or_gene")) if external_path.exists() else pd.DataFrame()
    overlap.to_csv(ROOT / "results/tables/GSE179044_GSE27982_supportive_overlap.csv", index=False)

    write_json(ROOT / "manifests/signature_build.json", {
        "query_contrasts": QUERY_CONTRASTS,
        "top_n_per_direction": 150,
        "fdr_cutoff": float(config["fdr_cutoff"]),
        "min_effect_gate": float(config["min_baseline_effect"]),
        "direction_reversal": "diagnostic only; excluded from default CMap query set",
        "functional_modules": "downstream summaries, not independent signatures",
        "external_overlap": "supportive case-normalized human/mouse comparison, not formal replication/meta-analysis",
    })
    print({"n_signature_rows": int(len(signature_table)), "n_module_rows": len(QUERY_CONTRASTS) * len(module_sets), "n_overlap_rows": int(len(overlap))})


if __name__ == "__main__":
    main()
