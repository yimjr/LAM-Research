"""Analyze TSC2-loss translational efficiency and its residual programs.

GSE277844 provides paired total and polysomal RNA-seq count columns for
isogenic TSC2-WT and TSC2-null neural progenitor cells, with DMSO, RMC-6272
and eFT-508 treatment.  This script uses library-size-normalized log2 CPM and
defines translational efficiency (TE) as log2 CPM(polysome) - log2 CPM(total)
within each genotype/treatment/replicate pair.

The analysis is intentionally effect-size first.  The TSC2-loss TE effect is
compared with the existing GSE179044 residual tables, then the two
translation-targeting perturbations are evaluated by whether they reduce the
KO-vs-WT TE distance.  Ratio classes are descriptive and gated by a minimum
baseline effect; they are not used as p-values.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from scipy.stats import fisher_exact
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from common import ROOT, sha256_file, write_json


DEFAULT_INPUT = ROOT / "data/raw/GSE277844/GSE277844_raw_counts.txt.gz"
FACTORIAL_INPUT = ROOT / "results/tables/GSE179044_factorial_contrasts.csv"
OUTPUT_DIR = ROOT / "data/processed/translation_analysis"
REPORT_DIR = ROOT / "reports/07_translation_analysis"
MANIFEST_DIR = ROOT / "manifests"

SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277844/suppl/"
    "GSE277844_raw_counts.txt.gz"
)

SAMPLE_PATTERN = re.compile(
    r"^(?P<genotype>WT|TSC2null)_(?P<treatment>DMSO|eFT508|RMC6272)_"
    r"(?P<fraction>Poly|Total)_rep(?P<replicate>\d+)$"
)

DRUGS = {
    "eFT508": "MNK1/2 inhibitor eFT-508",
    "RMC6272": "mTORC1 inhibitor RMC-6272",
}


def parse_sample_metadata(columns: list[str]) -> pd.DataFrame:
    rows = []
    for sample_id in columns:
        match = SAMPLE_PATTERN.match(str(sample_id))
        if not match:
            raise ValueError(f"Unrecognized GSE277844 sample column: {sample_id}")
        row = match.groupdict()
        row["sample_id"] = sample_id
        row["replicate"] = int(row["replicate"])
        row["pair_id"] = (
            f"{row['genotype']}_{row['treatment']}_rep{row['replicate']}"
        )
        rows.append(row)
    metadata = pd.DataFrame(rows)
    expected = metadata.groupby(["genotype", "treatment", "replicate"])["fraction"].agg(set)
    if not expected.map(lambda values: values == {"Poly", "Total"}).all():
        raise ValueError("Every GSE277844 genotype/treatment/replicate must have Poly and Total counts")
    return metadata


def load_counts(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    raw = pd.read_csv(path, sep="\t", compression="infer")
    if raw.shape[1] < 3:
        raise ValueError(f"Unexpected GSE277844 count matrix shape: {raw.shape}")
    gene_column = "geneID" if "geneID" in raw.columns else raw.columns[0]
    sample_columns = [column for column in raw.columns if column != gene_column]
    metadata = parse_sample_metadata(sample_columns)
    counts = raw.set_index(gene_column)[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # These are raw counts, so duplicate symbols (if present) are summed rather
    # than averaged.  Symbols are used as the cross-dataset comparison key.
    counts.index = counts.index.astype(str).str.strip().str.upper()
    counts = counts.groupby(level=0, sort=False).sum()
    counts = counts.loc[(counts >= 10).sum(axis=1) >= 4].copy()
    library_sizes = {}
    normalized = counts.copy()
    for fraction in ("Poly", "Total"):
        columns = metadata.loc[metadata["fraction"].eq(fraction), "sample_id"].tolist()
        sizes = counts[columns].sum(axis=0)
        library_sizes[fraction] = {str(key): float(value) for key, value in sizes.items()}
        normalized[columns] = counts[columns].divide(sizes.replace(0, np.nan), axis=1) * 1_000_000
    log_cpm = np.log2(normalized + 0.5)
    poly_columns = metadata.loc[metadata["fraction"].eq("Poly"), "sample_id"].tolist()
    total_columns = metadata.loc[metadata["fraction"].eq("Total"), "sample_id"].tolist()
    poly_log_cpm = log_cpm[poly_columns].copy()
    total_log_cpm = log_cpm[total_columns].copy()
    poly_meta = metadata.loc[metadata["fraction"].eq("Poly")].set_index("pair_id")
    total_meta = metadata.loc[metadata["fraction"].eq("Total")].set_index("pair_id")
    pair_ids = sorted(set(poly_meta.index) & set(total_meta.index))
    te = pd.DataFrame(index=log_cpm.index)
    for pair_id in pair_ids:
        poly_sample = poly_meta.loc[pair_id, "sample_id"]
        total_sample = total_meta.loc[pair_id, "sample_id"]
        te[pair_id] = log_cpm[poly_sample] - log_cpm[total_sample]
    audit = {
        "n_raw_genes": int(raw.shape[0]),
        "n_genes_after_count_filter": int(counts.shape[0]),
        "n_samples": int(len(sample_columns)),
        "n_te_pairs": int(len(pair_ids)),
        "count_filter": "counts >= 10 in at least 4 of 32 sample columns",
        "normalization": "library-size normalized CPM separately within Poly and Total fractions",
        "te_definition": "log2(CPM Poly + 0.5) - log2(CPM Total + 0.5) within matched genotype/treatment/replicate",
        "translation_effect_model": "conditional model log2 CPM Poly ~ centered log2 CPM Total + genotype, analogous to the translation component of anota2seq; not a claim to reproduce official anota2seq output",
        "library_sizes": library_sizes,
    }
    return counts, te, metadata, poly_log_cpm, total_log_cpm, audit


def condition_values(te: pd.DataFrame, metadata: pd.DataFrame, genotype: str, treatment: str) -> pd.DataFrame:
    pairs = metadata.loc[
        metadata["genotype"].eq(genotype)
        & metadata["treatment"].eq(treatment)
        & metadata["fraction"].eq("Poly"),
        "pair_id",
    ].drop_duplicates().tolist()
    return te.reindex(columns=pairs)


def fraction_values(matrix: pd.DataFrame, metadata: pd.DataFrame, genotype: str, treatment: str, fraction: str) -> pd.DataFrame:
    selected = metadata.loc[
        metadata["genotype"].eq(genotype)
        & metadata["treatment"].eq(treatment)
        & metadata["fraction"].eq(fraction),
        ["sample_id", "replicate"],
    ].sort_values("replicate")
    return matrix.reindex(columns=selected["sample_id"].tolist())


def condition_pair_map(metadata: pd.DataFrame, genotype: str, treatment: str) -> dict[int, str]:
    selected = metadata.loc[
        metadata["genotype"].eq(genotype)
        & metadata["treatment"].eq(treatment)
        & metadata["fraction"].eq("Poly"),
        ["replicate", "pair_id"],
    ].drop_duplicates()
    return {int(row.replicate): str(row.pair_id) for row in selected.itertuples(index=False)}


def safe_ttest(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if len(left) < 2 or len(right) < 2:
        return np.nan, np.nan
    test = ttest_ind(left, right, equal_var=False, nan_policy="omit")
    return float(test.statistic), float(test.pvalue)


def q_values(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    q = multipletests(p.fillna(1.0).to_numpy(float), method="fdr_bh")[1]
    return pd.Series(q, index=values.index)


def conditional_translation_effects(
    poly_log_cpm: pd.DataFrame,
    total_log_cpm: pd.DataFrame,
    metadata: pd.DataFrame,
    treatment: str,
) -> pd.DataFrame:
    """Estimate genotype-dependent translation conditional on total mRNA.

    This is a transparent Python analogue of the relevant anota2seq logic:
    the genotype coefficient in a centered ANCOVA captures polysome changes
    not explained by total mRNA.  It is deliberately labelled anota2seq-like
    because it does not implement the package's full regulation-mode and
    quality-control machinery.
    """
    poly_wt = fraction_values(poly_log_cpm, metadata, "WT", treatment, "Poly")
    poly_ko = fraction_values(poly_log_cpm, metadata, "TSC2null", treatment, "Poly")
    total_wt = fraction_values(total_log_cpm, metadata, "WT", treatment, "Total")
    total_ko = fraction_values(total_log_cpm, metadata, "TSC2null", treatment, "Total")
    rows = []
    for gene in poly_log_cpm.index:
        y = np.concatenate([poly_wt.loc[gene].to_numpy(float), poly_ko.loc[gene].to_numpy(float)])
        x = np.concatenate([total_wt.loc[gene].to_numpy(float), total_ko.loc[gene].to_numpy(float)])
        genotype = np.concatenate([np.zeros(len(poly_wt.columns)), np.ones(len(poly_ko.columns))])
        valid = np.isfinite(y) & np.isfinite(x)
        if valid.sum() < 4 or np.unique(genotype[valid]).size < 2:
            rows.append({"gene": gene, "translation_effect_TSC2null_minus_WT": np.nan, "translation_separation_t": np.nan, "translation_separation_p": np.nan})
            continue
        x_centered = x[valid] - np.mean(x[valid])
        design = np.column_stack([np.ones(valid.sum()), x_centered, genotype[valid]])
        covariance = np.linalg.pinv(design.T @ design)
        beta = covariance @ design.T @ y[valid]
        residual = y[valid] - design @ beta
        df_resid = max(int(valid.sum() - design.shape[1]), 1)
        residual_variance = float(np.sum(residual ** 2) / df_resid)
        se = float(np.sqrt(max(0.0, residual_variance * covariance[2, 2])))
        statistic = beta[2] / se if se > 0 else np.nan
        p_value = float(2.0 * student_t.sf(abs(statistic), df=df_resid)) if np.isfinite(statistic) else np.nan
        rows.append({
            "gene": gene,
            "translation_effect_TSC2null_minus_WT": float(beta[2]),
            "translation_separation_t": float(statistic) if np.isfinite(statistic) else np.nan,
            "translation_separation_p": p_value,
            "translation_total_effect_TSC2null_minus_WT": float(np.mean(total_ko.loc[gene]) - np.mean(total_wt.loc[gene])),
            "translation_polysome_effect_TSC2null_minus_WT": float(np.mean(poly_ko.loc[gene]) - np.mean(poly_wt.loc[gene])),
            "n_TSC2null": int(np.sum(valid & (genotype == 1))),
            "n_WT": int(np.sum(valid & (genotype == 0)),),
        })
    return pd.DataFrame(rows).set_index("gene")


def translation_effects(
    poly_log_cpm: pd.DataFrame,
    total_log_cpm: pd.DataFrame,
    metadata: pd.DataFrame,
    min_effect: float,
    fdr_cutoff: float,
) -> pd.DataFrame:
    result = conditional_translation_effects(poly_log_cpm, total_log_cpm, metadata, "DMSO")
    result["translation_separation_q"] = q_values(result["translation_separation_p"])
    effect = result["translation_effect_TSC2null_minus_WT"]
    selected = result["translation_separation_q"].le(fdr_cutoff) & effect.abs().ge(min_effect) & effect.ne(0)
    result["translation_direction"] = "not_selected"
    result.loc[selected & effect.gt(0), "translation_direction"] = "translation_up"
    result.loc[selected & effect.lt(0), "translation_direction"] = "translation_down"
    result["translation_selected"] = selected
    return result.reset_index()


def classify_rescue_ratio(ratio: float) -> str:
    if not np.isfinite(ratio):
        return "unclassified"
    if ratio < 0:
        return "direction_reversal"
    if ratio < 0.2:
        return "near_complete_rescue"
    if ratio < 0.8:
        return "partial_rescue_residual"
    if ratio <= 1.2:
        return "persistent_residual"
    return "worsened_residual"


def drug_translation_effects(
    poly_log_cpm: pd.DataFrame,
    total_log_cpm: pd.DataFrame,
    te: pd.DataFrame,
    metadata: pd.DataFrame,
    translation: pd.DataFrame,
    min_residual_effect: float,
    fdr_cutoff: float,
) -> pd.DataFrame:
    rows = []
    baseline_ko = condition_values(te, metadata, "TSC2null", "DMSO")
    baseline_wt = condition_values(te, metadata, "WT", "DMSO")
    baseline_conditional = translation.set_index("gene")
    baseline_ko_pairs = condition_pair_map(metadata, "TSC2null", "DMSO")
    baseline_wt_pairs = condition_pair_map(metadata, "WT", "DMSO")
    for drug in DRUGS:
        ko_after = condition_values(te, metadata, "TSC2null", drug)
        wt_after = condition_values(te, metadata, "WT", drug)
        after_conditional = conditional_translation_effects(poly_log_cpm, total_log_cpm, metadata, drug)
        after_conditional["translation_separation_q"] = q_values(after_conditional["translation_separation_p"])
        after_ko_pairs = condition_pair_map(metadata, "TSC2null", drug)
        after_wt_pairs = condition_pair_map(metadata, "WT", drug)
        for gene in te.index:
            ko_base = baseline_ko.loc[gene].to_numpy(float)
            wt_base = baseline_wt.loc[gene].to_numpy(float)
            ko_post = ko_after.loc[gene].to_numpy(float)
            wt_post = wt_after.loc[gene].to_numpy(float)
            baseline_effect = float(baseline_conditional.loc[gene, "translation_effect_TSC2null_minus_WT"])
            after_effect = float(after_conditional.loc[gene, "translation_effect_TSC2null_minus_WT"])
            ko_response = float(np.nanmean(ko_post) - np.nanmean(ko_base))
            wt_response = float(np.nanmean(wt_post) - np.nanmean(wt_base))
            # Replicate-matched response vectors provide the most direct
            # genotype-by-drug comparison available in this small dataset.
            common_reps_ko = sorted(set(baseline_ko_pairs) & set(after_ko_pairs))
            ko_response_vector = np.asarray(
                [ko_after.loc[gene, after_ko_pairs[rep]] - baseline_ko.loc[gene, baseline_ko_pairs[rep]] for rep in common_reps_ko],
                dtype=float,
            )
            common_reps_wt = sorted(set(baseline_wt_pairs) & set(after_wt_pairs))
            wt_response_vector = np.asarray(
                [wt_after.loc[gene, after_wt_pairs[rep]] - baseline_wt.loc[gene, baseline_wt_pairs[rep]] for rep in common_reps_wt],
                dtype=float,
            )
            interaction_t, interaction_p = safe_ttest(ko_response_vector, wt_response_vector)
            ratio = after_effect / baseline_effect if abs(baseline_effect) >= min_residual_effect else np.nan
            rows.append({
                "gene": gene,
                "drug": drug,
                "drug_description": DRUGS[drug],
                "baseline_translation_effect": baseline_effect,
                "baseline_translation_t": float(baseline_conditional.loc[gene, "translation_separation_t"]),
                "baseline_translation_p": float(baseline_conditional.loc[gene, "translation_separation_p"]),
                "baseline_translation_q": float(baseline_conditional.loc[gene, "translation_separation_q"]),
                "after_translation_effect": after_effect,
                "after_translation_t": float(after_conditional.loc[gene, "translation_separation_t"]),
                "after_translation_p": float(after_conditional.loc[gene, "translation_separation_p"]),
                "after_translation_q": float(after_conditional.loc[gene, "translation_separation_q"]),
                "KO_translation_response": ko_response,
                "WT_translation_response": wt_response,
                "genotype_drug_interaction": ko_response - wt_response,
                "genotype_drug_interaction_t": interaction_t,
                "genotype_drug_interaction_p": interaction_p,
                "signed_residual_ratio": ratio,
                "absolute_residual_ratio": abs(ratio) if np.isfinite(ratio) else np.nan,
                "distance_to_WT_reduced": abs(after_effect) < abs(baseline_effect) if np.isfinite(after_effect) else False,
                "rescue_class": classify_rescue_ratio(ratio),
                "n_KO_after": int(np.isfinite(ko_post).sum()),
                "n_WT_after": int(np.isfinite(wt_post).sum()),
            })
    result = pd.DataFrame(rows)
    for drug in DRUGS:
        mask = result["drug"].eq(drug)
        result.loc[mask, "genotype_drug_interaction_q"] = q_values(result.loc[mask, "genotype_drug_interaction_p"])
    result = result.merge(
        translation[["gene", "translation_direction", "translation_selected"]],
        on="gene",
        how="left",
        validate="many_to_one",
    )
    result["baseline_effect_eligible"] = result["baseline_translation_effect"].abs().ge(min_residual_effect)
    result["after_significant_and_effect_sized"] = (
        result["after_translation_q"].le(fdr_cutoff) & result["after_translation_effect"].abs().ge(min_residual_effect)
    )
    return result


def load_residual_sets(path: Path, min_effect: float, fdr_cutoff: float) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    factorial = pd.read_csv(path)
    factorial["gene"] = factorial["gene"].astype(str).str.upper().str.strip()
    factorial = factorial.drop_duplicates("gene").set_index("gene")
    sets: dict[str, set[str]] = {}
    for environment in ("plastic", "hydrogel"):
        class_column = f"residual_class_{environment}"
        residual_column = f"residual_{environment}"
        q_column = f"{residual_column}_moderated_q"
        sets[f"persistent_residual_{environment}"] = set(
            factorial.index[factorial[class_column].eq("persistent_residual")]
        )
        sets[f"persistent_residual_{environment}_q10"] = set(
            factorial.index[
                factorial[class_column].eq("persistent_residual")
                & factorial[residual_column].abs().ge(min_effect)
                & factorial[q_column].le(fdr_cutoff)
            ]
        )
    sets["hydrogel_residual_q10"] = set(
        factorial.index[
            factorial["residual_hydrogel"].abs().ge(min_effect)
            & factorial["residual_hydrogel_moderated_q"].le(fdr_cutoff)
        ]
    )
    sets["hydrogel_specific_residual_q10"] = set(
        factorial.index[
            factorial["hydrogel_specific_residual"].abs().ge(min_effect)
            & factorial["hydrogel_specific_residual_moderated_q"].le(fdr_cutoff)
        ]
    )
    return factorial.reset_index(), sets


def make_overlap_tables(translation: pd.DataFrame, factorial: pd.DataFrame, residual_sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = translation.loc[translation["translation_selected"]].copy()
    selected["gene"] = selected["gene"].astype(str).str.upper().str.strip()
    fact = factorial.set_index("gene")
    overlap_rows = []
    for _, row in selected.iterrows():
        gene = row["gene"]
        for category, genes in residual_sets.items():
            if gene not in genes:
                continue
            environment = "hydrogel" if "hydrogel" in category else "plastic"
            residual_column = f"residual_{environment}"
            q_column = f"{residual_column}_moderated_q"
            residual_effect = float(fact.loc[gene, residual_column]) if gene in fact.index else np.nan
            relation = "unknown"
            if np.isfinite(residual_effect) and residual_effect != 0:
                relation = "same_sign" if np.sign(residual_effect) == np.sign(row["translation_effect_TSC2null_minus_WT"]) else "opposite_sign"
            overlap_rows.append({
                "gene": gene,
                "translation_direction": row["translation_direction"],
                "translation_effect": row["translation_effect_TSC2null_minus_WT"],
                "translation_q": row["translation_separation_q"],
                "residual_category": category,
                "residual_effect": residual_effect,
                "residual_q": float(fact.loc[gene, q_column]) if gene in fact.index else np.nan,
                "residual_class": fact.loc[gene, f"residual_class_{environment}"] if gene in fact.index else "not_available",
                "translation_residual_direction_relation": relation,
            })
    long = pd.DataFrame(overlap_rows)
    summary_rows = []
    for direction in ("translation_up", "translation_down"):
        genes = set(selected.loc[selected["translation_direction"].eq(direction), "gene"])
        for category, residual_genes in residual_sets.items():
            intersection = genes & residual_genes
            summary_rows.append({
                "translation_direction": direction,
                "residual_category": category,
                "n_translation_genes": len(genes),
                "n_residual_genes": len(residual_genes),
                "n_overlap": len(intersection),
                "overlap_fraction_of_translation_genes": len(intersection) / len(genes) if genes else np.nan,
                "overlap_genes": ";".join(sorted(intersection)),
            })
    return long, pd.DataFrame(summary_rows)


def make_drug_summary(drug_effects: pd.DataFrame, overlap_membership: pd.DataFrame) -> pd.DataFrame:
    """Summarize treatment rescue separately for each residual category."""
    if overlap_membership.empty:
        return pd.DataFrame()
    subset = drug_effects.merge(
        overlap_membership[["gene", "translation_direction", "residual_category"]].drop_duplicates(),
        on=["gene", "translation_direction"],
        how="inner",
        validate="many_to_many",
    )
    rows = []
    for (drug, category), group in subset.groupby(["drug", "residual_category"], sort=True):
        eligible = group["baseline_effect_eligible"]
        rows.append({
            "drug": drug,
            "drug_description": group["drug_description"].iloc[0],
            "residual_category": category,
            "n_translation_residual_overlap_genes": int(len(group)),
            "n_baseline_effect_eligible": int(eligible.sum()),
            "n_distance_to_WT_reduced": int(group.loc[eligible, "distance_to_WT_reduced"].sum()),
            "fraction_distance_to_WT_reduced": float(group.loc[eligible, "distance_to_WT_reduced"].mean()) if eligible.any() else np.nan,
            "median_abs_baseline_effect": float(group.loc[eligible, "baseline_translation_effect"].abs().median()) if eligible.any() else np.nan,
            "median_abs_after_effect": float(group.loc[eligible, "after_translation_effect"].abs().median()) if eligible.any() else np.nan,
            "median_signed_residual_ratio": float(group.loc[eligible, "signed_residual_ratio"].median()) if eligible.any() else np.nan,
            "n_near_complete_rescue": int((group["rescue_class"] == "near_complete_rescue").sum()),
            "n_partial_rescue_residual": int((group["rescue_class"] == "partial_rescue_residual").sum()),
            "n_persistent_residual": int((group["rescue_class"] == "persistent_residual").sum()),
            "n_worsened_residual": int((group["rescue_class"] == "worsened_residual").sum()),
            "n_direction_reversal": int((group["rescue_class"] == "direction_reversal").sum()),
            "n_after_significant_and_effect_sized": int(group["after_significant_and_effect_sized"].sum()),
        })
    return pd.DataFrame(rows)


def make_background_comparison(
    drug_effects: pd.DataFrame,
    overlap_membership: pd.DataFrame,
) -> pd.DataFrame:
    """Compare overlap recovery with the full selected translation background.

    The primary comparison uses the same baseline-effect eligibility gate for
    both groups.  The all-translation background includes the overlap genes,
    while the non-overlap comparison is retained as a less circular sensitivity
    analysis.
    """
    selected_all = drug_effects.loc[drug_effects["translation_selected"]].copy()
    selected = selected_all.loc[selected_all["baseline_effect_eligible"]].copy()
    categories = sorted(overlap_membership["residual_category"].dropna().unique())
    rows = []

    def recovery_stats(group: pd.DataFrame) -> tuple[int, int, float]:
        n = int(len(group))
        recovered = int(group["distance_to_WT_reduced"].sum())
        rate = float(recovered / n) if n else np.nan
        return n, recovered, rate

    for drug, drug_group in selected.groupby("drug", sort=True):
        all_selected_group = selected_all.loc[selected_all["drug"].eq(drug)]
        all_selected_n, all_selected_recovered, all_selected_rate = recovery_stats(all_selected_group)
        all_n, all_recovered, all_rate = recovery_stats(drug_group)
        for category in categories:
            category_genes = set(
                overlap_membership.loc[
                    overlap_membership["residual_category"].eq(category), "gene"
                ]
            )
            overlap = drug_group.loc[drug_group["gene"].isin(category_genes)]
            nonoverlap = drug_group.loc[~drug_group["gene"].isin(category_genes)]
            overlap_all_selected = all_selected_group.loc[all_selected_group["gene"].isin(category_genes)]
            nonoverlap_all_selected = all_selected_group.loc[~all_selected_group["gene"].isin(category_genes)]
            overlap_all_n, overlap_all_recovered, overlap_all_rate = recovery_stats(overlap_all_selected)
            nonoverlap_all_n, nonoverlap_all_recovered, nonoverlap_all_rate = recovery_stats(nonoverlap_all_selected)
            overlap_n, overlap_recovered, overlap_rate = recovery_stats(overlap)
            nonoverlap_n, nonoverlap_recovered, nonoverlap_rate = recovery_stats(nonoverlap)
            fisher_table = [
                [overlap_recovered, overlap_n - overlap_recovered],
                [nonoverlap_recovered, nonoverlap_n - nonoverlap_recovered],
            ]
            fisher_odds_ratio, fisher_p = fisher_exact(fisher_table, alternative="two-sided")
            rows.append({
                "drug": drug,
                "drug_description": drug_group["drug_description"].iloc[0],
                "residual_category": category,
                "background_definition": "translation_selected and abs baseline translation effect >= 0.5",
                "n_all_selected_translation_genes": all_selected_n,
                "n_all_selected_recovered": all_selected_recovered,
                "all_selected_recovery_rate": all_selected_rate,
                "n_overlap_genes_all_selected": overlap_all_n,
                "n_overlap_recovered_all_selected": overlap_all_recovered,
                "overlap_recovery_rate_all_selected": overlap_all_rate,
                "overlap_minus_all_selected_rate": overlap_all_rate - all_selected_rate if np.isfinite(overlap_all_rate) else np.nan,
                "n_nonoverlap_genes_all_selected": nonoverlap_all_n,
                "n_nonoverlap_recovered_all_selected": nonoverlap_all_recovered,
                "nonoverlap_recovery_rate_all_selected": nonoverlap_all_rate,
                "n_all_selected_translation_genes_effect_eligible": all_n,
                "n_all_recovered": all_recovered,
                "all_recovery_rate": all_rate,
                "n_overlap_genes_effect_eligible": overlap_n,
                "n_overlap_recovered": overlap_recovered,
                "overlap_recovery_rate": overlap_rate,
                "overlap_minus_all_rate": overlap_rate - all_rate if np.isfinite(overlap_rate) else np.nan,
                "n_nonoverlap_genes_effect_eligible": nonoverlap_n,
                "n_nonoverlap_recovered": nonoverlap_recovered,
                "nonoverlap_recovery_rate": nonoverlap_rate,
                "overlap_minus_nonoverlap_rate": overlap_rate - nonoverlap_rate if np.isfinite(overlap_rate) and np.isfinite(nonoverlap_rate) else np.nan,
                "fisher_odds_ratio_overlap_vs_nonoverlap": float(fisher_odds_ratio) if np.isfinite(fisher_odds_ratio) else np.nan,
                "fisher_p_overlap_vs_nonoverlap": float(fisher_p),
            })
    return pd.DataFrame(rows)


def make_drug_gene_concordance(overlap_drug_effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Report gene-level agreement between the two translation-targeting drugs."""
    key_columns = ["gene", "translation_direction", "residual_category"]
    base = overlap_drug_effects[key_columns].drop_duplicates().copy()
    for drug in DRUGS:
        columns = [
            "baseline_effect_eligible",
            "distance_to_WT_reduced",
            "signed_residual_ratio",
            "absolute_residual_ratio",
            "after_translation_effect",
            "after_significant_and_effect_sized",
            "rescue_class",
        ]
        subset = overlap_drug_effects.loc[
            overlap_drug_effects["drug"].eq(drug), key_columns + columns
        ].drop_duplicates(key_columns)
        subset = subset.rename(columns={column: f"{drug}_{column}" for column in columns})
        base = base.merge(subset, on=key_columns, how="left", validate="one_to_one")

    rmc_eligible = base["RMC6272_baseline_effect_eligible"].fillna(False).astype(bool)
    eft_eligible = base["eFT508_baseline_effect_eligible"].fillna(False).astype(bool)
    rmc_reduced = base["RMC6272_distance_to_WT_reduced"].fillna(False).astype(bool)
    eft_reduced = base["eFT508_distance_to_WT_reduced"].fillna(False).astype(bool)
    both_eligible = rmc_eligible & eft_eligible
    base["both_drugs_distance_reduced"] = both_eligible & rmc_reduced & eft_reduced
    base["drug_support_pattern"] = "not_baseline_effect_eligible"
    base.loc[both_eligible & rmc_reduced & eft_reduced, "drug_support_pattern"] = "both_reduced"
    base.loc[both_eligible & rmc_reduced & ~eft_reduced, "drug_support_pattern"] = "RMC6272_only_reduced"
    base.loc[both_eligible & ~rmc_reduced & eft_reduced, "drug_support_pattern"] = "eFT508_only_reduced"
    base.loc[both_eligible & ~rmc_reduced & ~eft_reduced, "drug_support_pattern"] = "neither_reduced"

    summary_rows = []
    for category, group in base.groupby("residual_category", sort=True):
        eligible = group["drug_support_pattern"].ne("not_baseline_effect_eligible")
        counts = group.loc[eligible, "drug_support_pattern"].value_counts()
        n_eligible = int(eligible.sum())
        both = int(counts.get("both_reduced", 0))
        rmc_only = int(counts.get("RMC6272_only_reduced", 0))
        eft_only = int(counts.get("eFT508_only_reduced", 0))
        neither = int(counts.get("neither_reduced", 0))
        union = both + rmc_only + eft_only
        both_genes = sorted(group.loc[group["drug_support_pattern"].eq("both_reduced"), "gene"])
        rmc_only_genes = sorted(group.loc[group["drug_support_pattern"].eq("RMC6272_only_reduced"), "gene"])
        eft_only_genes = sorted(group.loc[group["drug_support_pattern"].eq("eFT508_only_reduced"), "gene"])
        neither_genes = sorted(group.loc[group["drug_support_pattern"].eq("neither_reduced"), "gene"])
        summary_rows.append({
            "residual_category": category,
            "n_overlap_genes": int(len(group)),
            "n_both_drug_effect_eligible": n_eligible,
            "n_both_reduced": both,
            "fraction_both_reduced": both / n_eligible if n_eligible else np.nan,
            "n_RMC6272_only_reduced": rmc_only,
            "n_eFT508_only_reduced": eft_only,
            "n_neither_reduced": neither,
            "fraction_any_drug_reduced": (both + rmc_only + eft_only) / n_eligible if n_eligible else np.nan,
            "reduced_set_jaccard": both / union if union else np.nan,
            "both_reduced_genes": ";".join(both_genes),
            "RMC6272_only_reduced_genes": ";".join(rmc_only_genes),
            "eFT508_only_reduced_genes": ";".join(eft_only_genes),
            "neither_reduced_genes": ";".join(neither_genes),
        })
    return base, pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(DEFAULT_INPUT.relative_to(ROOT)))
    parser.add_argument("--factorial", default=str(FACTORIAL_INPUT.relative_to(ROOT)))
    parser.add_argument("--min-effect", type=float, default=0.5, help="effect-size gate for residual ratios and treatment rescue summaries")
    parser.add_argument("--min-translation-effect", type=float, default=0.0, help="effect-size gate for selecting translation-up/down genes")
    parser.add_argument("--fdr-cutoff", type=float, default=0.15, help="FDR cutoff, matching the published anota2seq discovery convention")
    args = parser.parse_args()

    input_path = ROOT / args.input
    factorial_path = ROOT / args.factorial
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if not factorial_path.exists():
        raise FileNotFoundError(factorial_path)

    counts, te, metadata, poly_log_cpm, total_log_cpm, audit = load_counts(input_path)
    translation = translation_effects(poly_log_cpm, total_log_cpm, metadata, args.min_translation_effect, args.fdr_cutoff)
    drug_effects = drug_translation_effects(poly_log_cpm, total_log_cpm, te, metadata, translation, args.min_effect, args.fdr_cutoff)
    factorial, residual_sets = load_residual_sets(factorial_path, args.min_effect, args.fdr_cutoff)
    overlap_long, overlap_summary = make_overlap_tables(translation, factorial, residual_sets)
    overlap_genes = set(overlap_long["gene"]) if not overlap_long.empty else set()
    overlap_membership = overlap_long[["gene", "translation_direction", "residual_category", "residual_effect", "residual_q", "residual_class"]].drop_duplicates() if not overlap_long.empty else pd.DataFrame(columns=["gene", "translation_direction", "residual_category", "residual_effect", "residual_q", "residual_class"])
    overlap_drug_effects = drug_effects.merge(overlap_membership, on=["gene", "translation_direction"], how="inner", validate="many_to_many")
    drug_summary = make_drug_summary(drug_effects, overlap_membership)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    paths = {
        "sample_metadata": OUTPUT_DIR / "GSE277844_sample_metadata.csv",
        "translation_te_matrix": OUTPUT_DIR / "GSE277844_translation_efficiency_log2cpm.csv.gz",
        "translation_effects": OUTPUT_DIR / "GSE277844_tsc2_loss_translation_effects.csv",
        "residual_overlap_genes": OUTPUT_DIR / "GSE277844_translation_residual_overlap_genes.csv",
        "residual_overlap_summary": OUTPUT_DIR / "GSE277844_translation_residual_overlap_summary.csv",
        "drug_translation_effects": OUTPUT_DIR / "GSE277844_translation_targeting_drug_effects.csv.gz",
        "overlap_drug_effects": OUTPUT_DIR / "GSE277844_translation_residual_overlap_drug_effects.csv",
        "drug_summary": OUTPUT_DIR / "GSE277844_translation_residual_overlap_drug_summary.csv",
        "background_comparison": OUTPUT_DIR / "GSE277844_translation_residual_recovery_background_comparison.csv",
        "drug_gene_concordance": OUTPUT_DIR / "GSE277844_translation_residual_drug_gene_concordance.csv",
        "drug_gene_concordance_summary": OUTPUT_DIR / "GSE277844_translation_residual_drug_concordance_summary.csv",
        "report": REPORT_DIR / "GSE277844_translation_residual_analysis.md",
    }
    metadata.to_csv(paths["sample_metadata"], index=False)
    te.to_csv(paths["translation_te_matrix"], compression="gzip")
    translation.to_csv(paths["translation_effects"], index=False)
    overlap_long.to_csv(paths["residual_overlap_genes"], index=False)
    overlap_summary.to_csv(paths["residual_overlap_summary"], index=False)
    drug_effects.to_csv(paths["drug_translation_effects"], index=False, compression="gzip")
    overlap_drug_effects.to_csv(paths["overlap_drug_effects"], index=False)
    drug_summary.to_csv(paths["drug_summary"], index=False)
    background_comparison = make_background_comparison(drug_effects, overlap_membership)
    drug_gene_concordance, drug_gene_concordance_summary = make_drug_gene_concordance(overlap_drug_effects)
    background_comparison.to_csv(paths["background_comparison"], index=False)
    drug_gene_concordance.to_csv(paths["drug_gene_concordance"], index=False)
    drug_gene_concordance_summary.to_csv(paths["drug_gene_concordance_summary"], index=False)

    selected_up = set(translation.loc[translation["translation_direction"].eq("translation_up"), "gene"])
    selected_down = set(translation.loc[translation["translation_direction"].eq("translation_down"), "gene"])
    report_summary = {
        "n_selected_translation_up": len(selected_up),
        "n_selected_translation_down": len(selected_down),
        "n_translation_residual_overlap_genes": len(overlap_genes),
        "residual_set_sizes": {key: len(value) for key, value in residual_sets.items()},
        "drug_summary": drug_summary.to_dict("records"),
        "background_comparison": background_comparison.to_dict("records"),
        "drug_gene_concordance_summary": drug_gene_concordance_summary.to_dict("records"),
    }
    overlap_display = overlap_summary[[
        "translation_direction", "residual_category", "n_translation_genes", "n_residual_genes", "n_overlap", "overlap_fraction_of_translation_genes"
    ]]
    background_display = background_comparison[[
        "drug", "residual_category",
        "n_all_selected_translation_genes", "n_overlap_genes_all_selected", "n_overlap_recovered_all_selected",
        "all_selected_recovery_rate", "overlap_recovery_rate_all_selected", "overlap_minus_all_selected_rate",
        "n_all_selected_translation_genes_effect_eligible", "n_overlap_genes_effect_eligible",
        "n_overlap_recovered", "all_recovery_rate", "overlap_recovery_rate", "overlap_minus_all_rate",
        "overlap_minus_nonoverlap_rate", "fisher_p_overlap_vs_nonoverlap",
    ]]
    report_lines = [
        "# GSE277844 translation-loss and residual-program analysis",
        "",
        "## Research question",
        "",
        "先定义 TSC2 loss 是否改变 polysome-associated mRNA 相对于 total mRNA 的翻译效率，再与 GSE179044 的 persistent/hydrogel residual 比较，最后检查 mTORC1/MNK1/2 translation-targeting treatment 是否把重叠基因的 KO-vs-WT 翻译效率差异拉回 WT 附近。",
        "",
        "## Data and method",
        "",
        f"- input: `{input_path.relative_to(ROOT)}`; GEO supplementary source: {SOURCE_URL}",
        f"- raw genes: {audit['n_raw_genes']}; genes after count filter: {audit['n_genes_after_count_filter']}; TE pairs: {audit['n_te_pairs']}",
        "- TE = log2(CPM polysome + 0.5) − log2(CPM total + 0.5), paired by genotype × treatment × replicate.",
        "- TSC2-loss translation genes are called from the conditional polysome-vs-total model at the published-style FDR cutoff; the translation effect size is retained for ranking and interpretation rather than used as the primary significance test.",
        "- residual overlap uses the existing moderated GSE179044 contrasts. `persistent_residual_plastic` and `persistent_residual_hydrogel` are kept as separate categories; `_q10` additionally requires absolute residual effect ≥ 0.5 and moderated FDR ≤ 0.10. `hydrogel_residual_q10` is the post-rapamycin hydrogel residual with the same effect/FDR gate. No intersection between ordinary and hydrogel residual is required.",
        "- drug rescue is assessed by signed residual ratio = post-treatment KO-vs-WT TE effect / DMSO KO-vs-WT TE effect. Ratios are descriptive, especially for RMC-6272 where each group has two replicates.",
        "",
        "## Step 1: TSC2-loss translation program",
        "",
        f"- selected translation-up genes: {len(selected_up)}",
        f"- selected translation-down genes: {len(selected_down)}",
        "",
        "## Step 2: overlap with GSE179044 residual programs",
        "",
        overlap_display.to_string(index=False),
        "",
        f"The union of selected translation genes overlapping at least one residual category contains {len(overlap_genes)} genes for recordkeeping. Treatment summaries are reported separately by residual category, so this union is not used to require ordinary and hydrogel residual to overlap first.",
        "",
        "## Step 3: translation-targeting treatment",
        "",
        "The following summaries concern genes in each translation × residual overlap category separately. `distance_to_WT_reduced` asks whether the KO-vs-WT TE distance became smaller; it does not by itself prove a mechanistic rescue.",
        "",
        drug_summary.to_string(index=False),
        "",
        "## Background comparison",
        "",
        "The direct comparison uses all selected GSE277844 translation-abnormal genes; the effect-eligible comparison applies the same |baseline translation effect| ≥ 0.5 gate used for rescue ratios. The all-gene comparison includes overlap genes, while the non-overlap comparison and Fisher exact p-value are sensitivity diagnostics, not gene-independent proof.",
        "",
        background_display.to_string(index=False),
        "",
        "## Gene-level agreement between translation-targeting drugs",
        "",
        "`both_reduced` means both drugs reduced the KO-vs-WT translation distance for that gene. `RMC6272_only_reduced` and `eFT508_only_reduced` identify drug-specific support; genes failing the baseline effect gate are not assigned a recovery-support pattern.",
        "",
        drug_gene_concordance_summary.to_string(index=False),
        "",
        "## Interpretation limits",
        "",
        "GSE277844 is a human NPC model and is biologically distinct from the LAM cell model in GSE179044. The comparison is therefore a cross-model program test, not a direct LAM replication. Translation efficiency is estimated from normalized bulk count fractions, and treatment-specific rescue ratios are unstable when the baseline TE effect is small. RMC-6272/eFT-508 results should be followed by gene/module-level confirmation and, where possible, independent translational or genetic perturbation data.",
        "",
        "## Outputs",
        "",
        "- `GSE277844_tsc2_loss_translation_effects.csv`: all filtered genes with TE effect, p/q values and up/down selection.",
        "- `GSE277844_translation_residual_overlap_genes.csv`: gene-level overlap with persistent/hydrogel residual categories and direction relation.",
        "- `GSE277844_translation_residual_overlap_drug_effects.csv`: treatment effects for each gene × residual category overlap, including rescue ratios and classes.",
        "- `GSE277844_translation_residual_overlap_drug_summary.csv`: treatment × residual-category summary of recovery toward WT.",
        "- `GSE277844_translation_residual_recovery_background_comparison.csv`: overlap recovery compared with all selected translation-abnormal genes and the non-overlap sensitivity background.",
        "- `GSE277844_translation_residual_drug_gene_concordance.csv`: gene-level RMC-6272/eFT-508 recovery pattern for every residual-category overlap.",
        "- `GSE277844_translation_residual_drug_concordance_summary.csv`: per-category counts and gene lists for both-drug, RMC-6272-only and eFT-508-only recovery.",
    ]
    paths["report"].write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "GSE277844 translation efficiency → GSE179044 residual overlap → translation-targeting rescue",
        "input": str(input_path.relative_to(ROOT)),
        "factorial_input": str(factorial_path.relative_to(ROOT)),
        "source_url": SOURCE_URL,
        "input_sha256": sha256_file(input_path),
        "parameters": {
            "min_effect": args.min_effect,
            "min_translation_effect": args.min_translation_effect,
            "fdr_cutoff": args.fdr_cutoff,
            "translation_definition": audit["te_definition"],
            "drug_rescue_definition": "distance_to_WT_reduced = abs(post-treatment KO-WT TE) < abs(DMSO KO-WT TE)",
            "background_recovery_definition": "compare residual-overlap genes with all selected translation-abnormal genes; primary gated comparison uses abs baseline translation effect >= 0.5, with all-selected sensitivity comparison",
            "drug_gene_concordance_definition": "both_reduced, RMC6272_only_reduced, eFT508_only_reduced, neither_reduced, or not_baseline_effect_eligible based on distance_to_WT_reduced",
            "background_test": "two-sided Fisher exact test for overlap versus non-overlap among baseline-effect-eligible genes; descriptive only",
        },
        "audit": audit,
        "summary": report_summary,
        "residual_set_definitions": {
            "persistent_residual_plastic": "GSE179044 residual_class_plastic == persistent_residual",
            "persistent_residual_hydrogel": "GSE179044 residual_class_hydrogel == persistent_residual",
            "persistent_residual_*_q10": "persistent class plus abs residual effect >= 0.5 and moderated q <= 0.10",
            "hydrogel_residual_q10": "abs residual_hydrogel >= 0.5 and moderated q <= 0.10",
            "hydrogel_specific_residual_q10": "abs hydrogel_specific_residual >= 0.5 and moderated q <= 0.10",
        },
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in paths.items()},
    }
    write_json(MANIFEST_DIR / "GSE277844_translation_residual_analysis.json", manifest)
    print(json.dumps(report_summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
