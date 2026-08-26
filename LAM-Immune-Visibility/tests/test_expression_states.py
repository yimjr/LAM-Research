import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from common import ExpressionData, expression_intervals, gene_state_matrix, score_modules  # noqa: E402


def make_expression() -> ExpressionData:
    raw = np.asarray([
        [0, 0],
        [1, 1],
        [2, 2],
        [4, 4],
        [8, 8],
        [16, 16],
    ], dtype=float)
    normalized = np.log1p(raw)
    obs = pd.DataFrame({"assay": ["scRNA"] * 6}, index=[f"c{i}" for i in range(6)])
    return ExpressionData(
        obs=obs,
        target_genes=["A", "B", "C"],
        available_genes=["A", "B"],
        raw_counts=raw,
        normalized=normalized,
        gene_to_col={"A": 0, "B": 1},
        raw_counts_available=True,
        source_path="fixture",
    )


class ExpressionStateTests(unittest.TestCase):
    def test_nonzero_count_is_detected_and_missing_gene_is_not_assayed(self):
        expr = make_expression()
        intervals = expression_intervals(expr, positive_quantile=0.25, min_positive_observations=3)
        states = gene_state_matrix(expr, intervals)
        self.assertEqual(states[0, 0], "not_detected")
        self.assertTrue(states[1, 0].startswith("detected_"))
        self.assertEqual(states[0, 2], "not_assayed")
        self.assertIn("detected_low", set(states[:, 0]))


    def test_low_expression_is_detected_and_classified_by_expression_interval(self):
        expr = make_expression()
        intervals = expression_intervals(expr, positive_quantile=0.25, min_positive_observations=3)
        states = gene_state_matrix(expr, intervals)
        self.assertIn("detected_low", set(states[:, 0]))
        self.assertIn("detected_high", set(states[:, 0]))
        self.assertIn("not_detected", set(states[:, 0]))


    def test_module_missingness_is_explicit(self):
        expr = make_expression()
        scores, availability = score_modules(expr, {"test": {"genes": ["A", "C"]}}, min_available_genes=1)
        self.assertEqual(scores["module_test__n_available"].iloc[0], 1)
        self.assertEqual(scores["module_test__status"].iloc[0], "partial_coverage")
        self.assertEqual(availability.loc[0, "missing_genes"], "C")
