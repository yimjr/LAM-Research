"""Rank-based module enrichment utilities for future GSE135851/GSE302356 inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from common import ROOT


def rank_module_score(expression: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Return a simple rank-based module score, robust to platform scale."""
    available = [gene for gene in genes if gene in expression.index]
    if len(available) < 2:
        return pd.Series(np.nan, index=expression.columns)
    ranks = expression.rank(axis=0, pct=True)
    return ranks.loc[available].mean(axis=0)


def score_modules(expression: pd.DataFrame, module_sets: dict[str, list[str]]) -> pd.DataFrame:
    return pd.DataFrame({name: rank_module_score(expression, genes) for name, genes in module_sets.items()})


def run_from_tables(expression_path: str, modules_path: str, output_path: str) -> None:
    expression = pd.read_csv(ROOT / expression_path, index_col=0)
    modules = pd.read_csv(ROOT / modules_path)
    module_sets = {name: group["gene"].astype(str).tolist() for name, group in modules.groupby("module")}
    scores = score_modules(expression, module_sets)
    scores.to_csv(ROOT / output_path)
