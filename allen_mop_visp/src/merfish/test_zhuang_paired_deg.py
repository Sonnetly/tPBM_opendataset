"""Checks for the animal-paired regional DEG helpers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zhuang_paired_deg import bh, paired_table  # noqa: E402


def test_bh_matches_manual():
    q = bh([0.01, 0.04, 0.03, 0.5])
    assert np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.5], atol=1e-6)


def test_paired_table_signs_and_counts():
    rows = []
    for animal, shift in zip("ABCD", [0.20, 0.25, 0.22, 0.18]):
        rows.append({"animal": animal, "subclass_clean": "X", "cell_class": "Glutamatergic", "gene": "up",
                     "MOp_minus_VISp": shift})
        rows.append({"animal": animal, "subclass_clean": "X", "cell_class": "Glutamatergic", "gene": "mixed",
                     "MOp_minus_VISp": shift * (-1) ** "ABCD".index(animal)})
    out = paired_table(pd.DataFrame(rows), "MOp_minus_VISp").set_index("gene")
    assert out.loc["up", "all_same_sign"] and out.loc["up", "n_animals"] == 4
    assert np.isclose(out.loc["up", "mean_difference"], 0.2125)
    assert not out.loc["mixed", "all_same_sign"]
    assert out.loc["up", "p"] < out.loc["mixed", "p"]
