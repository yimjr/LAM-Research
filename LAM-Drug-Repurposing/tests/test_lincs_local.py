"""Hard algorithm acceptance tests for the local LINCS/CMap scorer."""

from __future__ import annotations

import math
import unittest
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_lincs_local import (  # noqa: E402
    _normalise_perturbation_key,
    weighted_ks_from_positions,
    weighted_pearson,
    wtcs,
)


class TestLocalLINCSAlgorithms(unittest.TestCase):
    def test_weighted_ks_hand_calculation(self):
        self.assertTrue(math.isclose(weighted_ks_from_positions(np.array([0]), np.array([1.0]), 4), 1.0))
        self.assertTrue(math.isclose(weighted_ks_from_positions(np.array([3]), np.array([1.0]), 4), -1.0))


    def test_wtcs_same_sign_is_zero_and_swap_reverses_direction(self):
        self.assertEqual(wtcs(0.5, 0.5), 0.0)
        self.assertEqual(wtcs(-0.5, -0.5), 0.0)
        self.assertEqual(wtcs(0.5, -0.5), -wtcs(-0.5, 0.5))


    def test_weighted_correlation_uses_nonnegative_weights(self):
        x = np.array([1.0, 2.0, -1.0])
        y = np.array([-1.0, -2.0, 1.0])
        self.assertLess(weighted_pearson(x, y, np.abs(x)), -0.99)


    def test_perturbation_class_separation(self):
        self.assertEqual(_normalise_perturbation_key("trt_cp", "sirolimus")[1], "compound")
        self.assertEqual(_normalise_perturbation_key("trt_sh", "MTOR")[1], "genetic")


if __name__ == "__main__":
    unittest.main()
