"""Synthetic checks for the supertype standardization in zhuang_supertype_check.py."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zhuang_supertype_check import supertype_decomposition  # noqa: E402


def _make(n_mop, n_visp, level_by_supertype, shift_by_region=None):
    shift_by_region = shift_by_region or {"MOp": 0.0, "VISp": 0.0}
    labels, regions, supertypes, values = [], [], [], []
    for region, counts in (("MOp", n_mop), ("VISp", n_visp)):
        for supertype, n in counts.items():
            for i in range(n):
                labels.append(f"{region}_{supertype}_{i}")
                regions.append(region)
                supertypes.append(supertype)
                values.append(level_by_supertype[supertype] + shift_by_region[region])
    cells = pd.DataFrame({"cell_label": labels, "animal": "A", "region": regions,
                          "subclass_clean": "L4/5 IT CTX Glut", "cell_class": "Glutamatergic",
                          "supertype": supertypes})
    expression = pd.DataFrame({"G": values}, index=pd.Index(labels))
    return cells, expression


def test_pure_composition_difference_is_removed():
    # Same per-supertype expression in both regions, different proportions.
    cells, expression = _make({"a": 300, "b": 100}, {"a": 100, "b": 300}, {"a": 1.0, "b": 0.0})
    out, _ = supertype_decomposition(cells, expression, min_cells=20, min_subclass_cells=50)
    row = out.iloc[0]
    assert np.isclose(row.raw_restricted, 0.5)
    assert np.isclose(row.standardized, 0.0)
    assert np.isclose(row.composition_part, 0.5)
    assert np.isclose(row.coverage, 1.0)


def test_pure_expression_difference_is_kept():
    cells, expression = _make({"a": 200, "b": 200}, {"a": 200, "b": 200}, {"a": 1.0, "b": 0.0},
                              {"MOp": -0.1, "VISp": 0.0})
    out, _ = supertype_decomposition(cells, expression, min_cells=20, min_subclass_cells=50)
    row = out.iloc[0]
    assert np.isclose(row.standardized, -0.1)
    assert np.isclose(row.composition_part, 0.0, atol=1e-7)


def test_unshared_supertype_is_excluded_and_reported_in_coverage():
    cells, expression = _make({"a": 200, "b": 100}, {"a": 200, "b": 5}, {"a": 1.0, "b": 0.0})
    out, shares = supertype_decomposition(cells, expression, min_cells=20, min_subclass_cells=50)
    row = out.iloc[0]
    assert row.supertypes_kept == 1 and row.supertypes_total == 2
    assert np.isclose(row.coverage, 400 / 505)
    assert set(shares.supertype) == {"a", "b"}
