"""Annotate the fixed hydrogel translation-residual core.

The analysis object is deliberately locked to genes satisfying all of:
``residual_category == hydrogel_residual_q10``, both drug baseline-effect
eligibility flags, and ``both_drugs_distance_reduced == True``.  This keeps
the functional interpretation separate from the broader residual categories.

The 13-gene set is annotated against local GO Biological Process, GO Cellular
Component, Reactome and MSigDB Hallmark GMT files.  Term-level overlap and
over-representation are retained, but repeated themes across genes are the
main interpretation target because the gene set is small.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from common import ROOT, sha256_file, write_json


DEFAULT_CONCORDANCE = ROOT / "data/processed/translation_analysis/GSE277844_translation_residual_drug_gene_concordance.csv"
DEFAULT_TRANSLATION = ROOT / "data/processed/translation_analysis/GSE277844_tsc2_loss_translation_effects.csv"
DEFAULT_GENE_SET_DIR = ROOT / "data/processed/gene_sets"
OUTPUT_DIR = ROOT / "data/processed/translation_analysis"
REPORT_DIR = ROOT / "reports/07_translation_analysis"
MANIFEST_DIR = ROOT / "manifests"

EXPECTED_CORE_ORDER = [
    "CDC42EP3", "RND3", "SERPINE2", "GPR27", "WWTR1", "FBN2", "REEP2",
    "ZNF354C", "PNMA2", "CACFD1", "NFATC4", "SPIN4", "GPC4",
]
EXPECTED_CORE = set(EXPECTED_CORE_ORDER)

LIBRARY_FILES = {
    "GO_Biological_Process": "GO_Biological_Process_2023.gmt",
    "GO_Cellular_Component": "GO_Cellular_Component_2023.gmt",
    "Reactome": "Reactome_2022.gmt",
    "MSigDB_Hallmark": "MSigDB_Hallmark_2020.gmt",
}

THEME_PATTERNS = {
    "extracellular_matrix": [
        r"extracellular matrix", r"collagen", r"elastic fibre", r"microfibril",
    ],
    "actin_cytoskeleton": [r"actin", r"cytoskeleton", r"pseudopod"],
    "cell_adhesion": [r"adhesion", r"cell-substrate junction", r"junction"],
    "Rho_GTPase": [r"rho gtpase", r"cdc42 gtpase", r"rhoq gtpase"],
    "focal_adhesion": [r"focal adhesion", r"cell-substrate junction"],
    "Hippo_YAP_TAZ": [r"hippo", r"yap", r"taz", r"transcriptional coactivator"],
    "mechanotransduction": [r"mechanotransduction", r"mechanical", r"cell-substrate junction"],
    "TGF_beta": [r"tgf[- ]?beta", r"transforming growth factor", r"smad"],
    "migration": [r"migration", r"motility", r"pseudopod", r"cell projection"],
}


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def parse_term(term: str) -> tuple[str, str]:
    match = re.search(r"\s*\((GO:\d+|R-HSA-\d+)\)\s*$", term)
    if not match:
        return term.strip(), ""
    return term[: match.start()].strip(), match.group(1)


def parse_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term = fields[0].strip()
            genes = {gene.strip().upper() for gene in fields[2:] if gene.strip()}
            if term and genes:
                gene_sets.setdefault(term, set()).update(genes)
    return gene_sets


def themes_for_term(term_name: str) -> list[str]:
    text = term_name.lower()
    return sorted(
        theme for theme, patterns in THEME_PATTERNS.items()
        if any(re.search(pattern, text) for pattern in patterns)
    )


def load_core(concordance_path: Path, translation_path: Path) -> pd.DataFrame:
    concordance = pd.read_csv(concordance_path)
    concordance["gene"] = concordance["gene"].astype(str).str.strip().str.upper()
    mask = (
        concordance["residual_category"].eq("hydrogel_residual_q10")
        & parse_bool(concordance["RMC6272_baseline_effect_eligible"])
        & parse_bool(concordance["eFT508_baseline_effect_eligible"])
        & parse_bool(concordance["both_drugs_distance_reduced"])
    )
    core = concordance.loc[mask].drop_duplicates("gene").copy()
    observed = set(core["gene"])
    if observed != EXPECTED_CORE:
        raise ValueError(
            "Fixed hydrogel translation-residual core changed: "
            f"expected {sorted(EXPECTED_CORE)}, observed {sorted(observed)}"
        )
    translation = pd.read_csv(translation_path)
    translation["gene"] = translation["gene"].astype(str).str.strip().str.upper()
    keep_translation = [
        "gene", "translation_effect_TSC2null_minus_WT",
        "translation_separation_q",
    ]
    core = core.merge(
        translation[keep_translation], on="gene", how="left", validate="one_to_one"
    )
    core["baseline_effect_eligible"] = True
    core["selection_criteria"] = (
        "residual_category=hydrogel_residual_q10; "
        "baseline_effect_eligible=True for RMC6272 and eFT508; "
        "both_drugs_distance_reduced=True"
    )
    core["fixed_core_order"] = core["gene"].map({gene: i for i, gene in enumerate(EXPECTED_CORE_ORDER)})
    return core.sort_values("fixed_core_order").reset_index(drop=True)


def enrich_libraries(
    core_genes: set[str],
    background_genes: set[str],
    libraries: dict[str, dict[str, set[str]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    enrichment_rows = []
    annotation_rows = []
    n_core = len(core_genes)
    n_background = len(background_genes)
    for library, gene_sets in libraries.items():
        for raw_term, genes in gene_sets.items():
            term_name, term_id = parse_term(raw_term)
            term_background = genes & background_genes
            hits = core_genes & term_background
            if not hits:
                continue
            themes = themes_for_term(term_name)
            overlap = len(hits)
            p_value = float(hypergeom.sf(overlap - 1, n_background, len(term_background), n_core))
            enrichment_rows.append({
                "library": library,
                "term": raw_term,
                "term_name": term_name,
                "term_id": term_id,
                "term_size_in_background": len(term_background),
                "core_size": n_core,
                "background_size": n_background,
                "overlap_count": overlap,
                "overlap_fraction_of_core": overlap / n_core,
                "overlap_genes": ";".join(sorted(hits)),
                "focus_themes": ";".join(themes),
                "p_value": p_value,
            })
            for gene in sorted(hits):
                annotation_rows.append({
                    "gene": gene,
                    "library": library,
                    "term": raw_term,
                    "term_name": term_name,
                    "term_id": term_id,
                    "term_size_in_background": len(term_background),
                    "core_overlap_count": overlap,
                    "core_overlap_genes": ";".join(sorted(hits)),
                    "focus_themes": ";".join(themes),
                })
    enrichment = pd.DataFrame(enrichment_rows)
    if not enrichment.empty:
        enrichment["fdr"] = np.nan
        for library, indices in enrichment.groupby("library").groups.items():
            enrichment.loc[indices, "fdr"] = multipletests(
                enrichment.loc[indices, "p_value"].to_numpy(float), method="fdr_bh"
            )[1]
        enrichment = enrichment.sort_values(["library", "fdr", "p_value", "term_name"])
    annotations = pd.DataFrame(annotation_rows)
    if not annotations.empty:
        annotations = annotations.sort_values(["gene", "library", "term_name"])
    return enrichment, annotations


def make_gene_summary(core: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    summary = core[[
        "gene", "translation_direction", "translation_effect_TSC2null_minus_WT",
        "translation_separation_q", "residual_category", "RMC6272_distance_to_WT_reduced",
        "eFT508_distance_to_WT_reduced", "both_drugs_distance_reduced",
    ]].copy()
    summary = summary.rename(columns={
        "translation_effect_TSC2null_minus_WT": "translation_effect",
        "translation_separation_q": "translation_q",
    })
    for library in LIBRARY_FILES:
        if annotations.empty:
            group = pd.DataFrame(columns=["gene", "term_name"])
        else:
            group = annotations.loc[annotations["library"].eq(library)]
        term_counts = group.groupby("gene").size().rename(f"{library}_n_terms")
        summary = summary.merge(term_counts, on="gene", how="left")
    summary = summary.fillna({column: 0 for column in summary.columns if column.endswith("_n_terms")})
    theme_rows = []
    for gene in summary["gene"]:
        group = annotations.loc[annotations["gene"].eq(gene)] if not annotations.empty else pd.DataFrame()
        themes = sorted({theme for value in group["focus_themes"].fillna("") for theme in value.split(";") if theme}) if not group.empty else []
        theme_rows.append({
            "gene": gene,
            "focus_themes": ";".join(themes),
            "focus_term_count": int(group["focus_themes"].astype(str).ne("").sum()) if not group.empty else 0,
        })
    summary = summary.merge(pd.DataFrame(theme_rows), on="gene", validate="one_to_one")
    return summary


def make_theme_tables(core_genes: set[str], annotations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    theme_by_gene_rows = []
    for gene in sorted(core_genes, key=lambda value: EXPECTED_CORE_ORDER.index(value)):
        gene_annotations = annotations.loc[annotations["gene"].eq(gene)] if not annotations.empty else pd.DataFrame()
        for theme in THEME_PATTERNS:
            hits = gene_annotations.loc[gene_annotations["focus_themes"].fillna("").str.contains(rf"(?:^|;){re.escape(theme)}(?:;|$)", regex=True)] if not gene_annotations.empty else pd.DataFrame()
            theme_by_gene_rows.append({
                "gene": gene,
                "theme": theme,
                "observed": not hits.empty,
                "n_supporting_terms": int(len(hits)),
                "supporting_libraries": ";".join(sorted(hits["library"].unique())) if not hits.empty else "",
                "supporting_terms": ";".join(hits["term_name"].drop_duplicates().tolist()) if not hits.empty else "",
            })
    by_gene = pd.DataFrame(theme_by_gene_rows)
    summary_rows = []
    for theme, group in by_gene.groupby("theme", sort=False):
        observed = group.loc[group["observed"]]
        summary_rows.append({
            "theme": theme,
            "n_core_genes": int(observed["gene"].nunique()),
            "genes": ";".join(observed["gene"].drop_duplicates().tolist()),
            "n_supporting_gene_term_links": int(observed["n_supporting_terms"].sum()),
            "n_supporting_libraries": int(len({library for value in observed["supporting_libraries"] for library in value.split(";") if library})),
            "supporting_libraries": ";".join(sorted({library for value in observed["supporting_libraries"] for library in value.split(";") if library})),
            "supporting_terms": ";".join(observed["supporting_terms"].drop_duplicates().tolist()),
            "repeated_across_genes": int(observed["gene"].nunique()) >= 2,
        })
    return by_gene, pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--concordance", default=str(DEFAULT_CONCORDANCE.relative_to(ROOT)))
    parser.add_argument("--translation", default=str(DEFAULT_TRANSLATION.relative_to(ROOT)))
    parser.add_argument("--gene-set-dir", default=str(DEFAULT_GENE_SET_DIR.relative_to(ROOT)))
    args = parser.parse_args()

    concordance_path = ROOT / args.concordance
    translation_path = ROOT / args.translation
    gene_set_dir = ROOT / args.gene_set_dir
    for path in (concordance_path, translation_path):
        if not path.exists():
            raise FileNotFoundError(path)
    missing_libraries = [name for name, filename in LIBRARY_FILES.items() if not (gene_set_dir / filename).exists()]
    if missing_libraries:
        raise FileNotFoundError(f"Missing local gene-set libraries: {missing_libraries}")

    core = load_core(concordance_path, translation_path)
    translation = pd.read_csv(translation_path)
    translation["gene"] = translation["gene"].astype(str).str.strip().str.upper()
    background = set(translation.loc[translation["translation_effect_TSC2null_minus_WT"].notna(), "gene"])
    libraries = {name: parse_gmt(gene_set_dir / filename) for name, filename in LIBRARY_FILES.items()}
    enrichment, annotations = enrich_libraries(set(core["gene"]), background, libraries)
    gene_summary = make_gene_summary(core, annotations)
    theme_by_gene, theme_summary = make_theme_tables(set(core["gene"]), annotations)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    paths = {
        "core_genes": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_genes.csv",
        "functional_annotations": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_functional_annotations.csv.gz",
        "gene_function_summary": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_gene_function_summary.csv",
        "theme_by_gene": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_functional_theme_by_gene.csv",
        "theme_summary": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_functional_theme_summary.csv",
        "functional_enrichment": OUTPUT_DIR / "GSE277844_hydrogel_translation_core_functional_enrichment.csv.gz",
        "report": REPORT_DIR / "GSE277844_hydrogel_translation_core_functional_annotation.md",
    }
    core.to_csv(paths["core_genes"], index=False)
    annotations.to_csv(paths["functional_annotations"], index=False, compression="gzip")
    gene_summary.to_csv(paths["gene_function_summary"], index=False)
    theme_by_gene.to_csv(paths["theme_by_gene"], index=False)
    theme_summary.to_csv(paths["theme_summary"], index=False)
    enrichment.to_csv(paths["functional_enrichment"], index=False, compression="gzip")

    annotation_display = gene_summary[[
        "gene", "translation_direction", "translation_effect", "translation_q",
        "GO_Biological_Process_n_terms", "GO_Cellular_Component_n_terms",
        "Reactome_n_terms", "MSigDB_Hallmark_n_terms", "focus_themes",
    ]]
    repeated_themes = theme_summary.loc[theme_summary["repeated_across_genes"]].copy()
    repeated_display = repeated_themes[[
        "theme", "n_core_genes", "genes", "n_supporting_libraries", "supporting_libraries",
    ]] if not repeated_themes.empty else theme_summary.head(0)
    enrichment_display = enrichment.sort_values(["library", "fdr", "p_value"]).groupby("library", sort=False).head(10) if not enrichment.empty else enrichment
    enrichment_columns = ["library", "term_name", "overlap_count", "overlap_genes", "focus_themes", "p_value", "fdr"]

    report_lines = [
        "# GSE277844 hydrogel translation-residual core functional annotation",
        "",
        "## Fixed analysis object",
        "",
        "The analysis object is fixed by `residual_category=hydrogel_residual_q10`, baseline-effect eligibility for both RMC-6272 and eFT-508, and `both_drugs_distance_reduced=True`. It contains exactly 13 genes; ordinary/plastic residual genes are not used in this functional annotation step.",
        "",
        "; ".join(EXPECTED_CORE_ORDER),
        "",
        "## Data and method",
        "",
        f"- concordance input: `{concordance_path.relative_to(ROOT)}`",
        f"- translation input/background: `{translation_path.relative_to(ROOT)}`; background = {len(background)} genes with finite conditional translation effect",
        "- libraries: GO Biological Process, GO Cellular Component, Reactome and MSigDB Hallmark from local GMT files",
        "- term-level enrichment uses a hypergeometric over-representation test with the analyzable GSE277844 background and Benjamini–Hochberg FDR within each library.",
        "- the primary interpretation is gene-level annotation and repeated themes across the 13 genes; non-significant FDR does not erase a coherent small-module signal.",
        "",
        "## Per-gene functional annotation summary",
        "",
        annotation_display.to_string(index=False),
        "",
        "## Repeated functional themes",
        "",
        repeated_display.to_string(index=False),
        "",
        "Themes are treated as repeated when at least two of the 13 genes have one or more supporting terms. Keyword-based theme labels are operational summaries; the full term-level annotations are retained in the long output.",
        "",
        "## Term enrichment (descriptive)",
        "",
        enrichment_display[enrichment_columns].to_string(index=False) if not enrichment_display.empty else "No term-level overlap was found.",
        "",
        "## Interpretation",
        "",
        "The repeated-theme table should be read as a functional convergence check, not as proof that every annotated term is active in the same cell. GO/Reactome/Hallmark terms are overlapping, so counts across libraries are not independent evidence. The 13 genes were selected using residual and treatment-response criteria, and GSE277844 remains a cross-model human NPC dataset rather than a direct LAM experiment.",
        "",
        "## Outputs",
        "",
        "- `GSE277844_hydrogel_translation_core_genes.csv`: the locked 13-gene analysis object with treatment-response fields.",
        "- `GSE277844_hydrogel_translation_core_functional_annotations.csv.gz`: gene × library × term annotations.",
        "- `GSE277844_hydrogel_translation_core_gene_function_summary.csv`: one row per fixed gene with library term counts and focus themes.",
        "- `GSE277844_hydrogel_translation_core_functional_theme_by_gene.csv`: gene-level theme support matrix.",
        "- `GSE277844_hydrogel_translation_core_functional_theme_summary.csv`: repeated-theme summary across genes.",
        "- `GSE277844_hydrogel_translation_core_functional_enrichment.csv.gz`: descriptive term enrichment with FDR.",
    ]
    paths["report"].write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "fixed 13-gene hydrogel residual × translation core functional annotation",
        "selection": {
            "residual_category": "hydrogel_residual_q10",
            "baseline_effect_eligible": True,
            "both_drugs_distance_reduced": True,
            "expected_genes": EXPECTED_CORE_ORDER,
        },
        "inputs": {
            "concordance": str(concordance_path.relative_to(ROOT)),
            "translation": str(translation_path.relative_to(ROOT)),
            "concordance_sha256": sha256_file(concordance_path),
            "translation_sha256": sha256_file(translation_path),
        },
        "gene_set_libraries": {
            name: {
                "path": str((gene_set_dir / filename).relative_to(ROOT)),
                "sha256": sha256_file(gene_set_dir / filename),
            }
            for name, filename in LIBRARY_FILES.items()
        },
        "background": "all genes in GSE277844_tsc2_loss_translation_effects.csv with finite translation_effect_TSC2null_minus_WT",
        "methods": {
            "term_enrichment": "hypergeometric over-representation test; Benjamini-Hochberg FDR within library",
            "focus_themes": THEME_PATTERNS,
            "interpretation": "gene-level and repeated-theme evidence prioritized because n=13; FDR is descriptive and overlapping libraries are not independent",
        },
        "summary": {
            "n_core_genes": int(len(core)),
            "n_annotations": int(len(annotations)),
            "n_enriched_terms_with_overlap": int(len(enrichment)),
            "repeated_themes": repeated_themes["theme"].tolist(),
        },
        "outputs": {key: str(path.relative_to(ROOT)) for key, path in paths.items()},
    }
    write_json(MANIFEST_DIR / "GSE277844_hydrogel_translation_core_functional_annotation.json", manifest)
    print(json.dumps(manifest["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
