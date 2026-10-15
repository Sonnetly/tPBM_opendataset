#!/usr/bin/env python3
"""Synthetic geometry tests for src/merfish/depth.py. No Allen data is touched."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from depth import (assign_shell, cortex_depth_limit, cortical_depth, depth_cumulative,
                   depth_layer_profile, depth_shell_difference, depth_shell_ratios)

LAYER_DEPTH_MM = {"1": 0.0, "2/3": 0.2, "5": 0.6, "6a": 1.0, "6b": 1.4}


def flat_cortex(region, hemisphere, thickness_scale, n_per_layer=40, section="S1", seed=0):
    """A flat sheet of cortex. Pia is the y = 0 line, depth is exactly y.

    thickness_scale compresses every layer, which is how MOp and VISp differ.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for layer, depth in LAYER_DEPTH_MM.items():
        for i in range(n_per_layer):
            rows.append({
                "cell_label": f"{region}-{layer}-{i}",
                "brain_section_label": section, "region": region, "hemisphere": hemisphere,
                "anatomical_layer": layer,
                "y_ccf": depth * thickness_scale,
                "z_ccf": 1.0 + 0.01 * i,
                "major_group": "GABAergic neuron" if i % 5 == 0 else "Glutamatergic neuron",
                "raw_subclass": "Pvalb" if i % 5 == 0 else "L5 IT CTX",
                "gaba_family": "Pvalb" if i % 5 == 0 else "Other GABA",
                "is_neuron": True, "is_gaba": i % 5 == 0, "qc_annotated": True,
            })
    # A few glia so the per-100-neurons columns are exercised.
    for i in range(10):
        rows.append({
            "cell_label": f"{region}-astro-{i}", "brain_section_label": section, "region": region,
            "hemisphere": hemisphere, "anatomical_layer": "2/3",
            "y_ccf": 0.2 * thickness_scale, "z_ccf": 1.0 + 0.01 * i,
            "major_group": "Astrocyte", "raw_subclass": "Astro-TE NN", "gaba_family": "Other GABA",
            "is_neuron": False, "is_gaba": False, "qc_annotated": True,
        })
    _ = rng
    return pd.DataFrame(rows)


class DepthTests(unittest.TestCase):
    def setUp(self):
        thick = flat_cortex("MOp", "left", thickness_scale=1.0, seed=1)
        thin = flat_cortex("VISp", "left", thickness_scale=0.5, seed=2)
        self.cells = pd.concat([thick, thin], ignore_index=True)
        self.depth = cortical_depth(self.cells, min_reference_cells=5)

    def test_depth_matches_known_geometry(self):
        """On a flat sheet, depth must equal the y offset in micrometres."""
        for layer, mm in LAYER_DEPTH_MM.items():
            mop = self.depth[self.depth.region.eq("MOp") & self.depth.anatomical_layer.eq(layer)]
            self.assertTrue(np.allclose(mop.cortical_depth_um, mm * 1000, atol=11),
                            f"MOp {layer}: {mop.cortical_depth_um.describe()}")
            visp = self.depth[self.depth.region.eq("VISp") & self.depth.anatomical_layer.eq(layer)]
            self.assertTrue(np.allclose(visp.cortical_depth_um, mm * 500, atol=11))

    def test_layer_labels_are_not_depth_matched(self):
        """The failure mode this module exists for: same label, different physical depth."""
        profile = depth_layer_profile(self.depth).set_index(["region", "anatomical_layer"])
        self.assertGreater(profile.loc[("MOp", "5"), "depth_um_median"],
                           profile.loc[("VISp", "5"), "depth_um_median"] * 1.8)

    def test_input_frame_is_not_mutated(self):
        self.assertNotIn("cortical_depth_um", self.cells.columns)

    def test_insufficient_pia_reference_gives_nan_not_a_guess(self):
        sparse = self.cells[~(self.cells.region.eq("MOp") & self.cells.anatomical_layer.eq("1"))]
        out = cortical_depth(sparse, min_reference_cells=5)
        self.assertTrue(out[out.region.eq("MOp")].cortical_depth_um.isna().all())
        self.assertFalse(out[out.region.eq("VISp")].cortical_depth_um.isna().any())

    def test_shell_edges_are_half_open_and_ordered(self):
        shells = assign_shell(pd.Series([0.0, 99.9, 100.0, 1199.0, 1200.0, np.nan]))
        self.assertEqual(list(shells), ["0000-0100", "0000-0100", "0100-0200",
                                        "1000-1200", "1200+", "unknown"])

    def test_cortex_limit_flags_subcortical_shells(self):
        limits = cortex_depth_limit(self.depth)
        self.assertLess(limits["VISp"], limits["MOp"])
        shells = depth_shell_ratios(self.depth, None, min_shell_cells=1)
        deep_visp = shells[shells.region.eq("VISp") & shells.depth_shell.eq("1000-1200")]
        if len(deep_visp):
            self.assertFalse(bool(deep_visp.within_cortex_depth.iloc[0]))

    def test_cumulative_is_monotonic_and_reaches_100(self):
        table = depth_cumulative(self.depth)
        for _, frame in table.groupby("region"):
            values = frame.sort_values("depth_um").cumulative_percentage.to_numpy()
            self.assertTrue(np.all(np.diff(values) >= -1e-9))
            self.assertAlmostEqual(values[-1], 100.0, places=6)

    def test_empty_depth_does_not_raise(self):
        """Tiny or malformed inputs must give typed empty tables, not a KeyError."""
        blank = cortical_depth(self.cells.head(3), min_reference_cells=1000)
        self.assertTrue(depth_layer_profile(blank).empty)
        shells = depth_shell_ratios(blank, None)
        self.assertTrue(shells.empty)
        self.assertTrue(depth_shell_difference(shells).empty)

    def test_per_100_neurons_removes_the_shared_denominator(self):
        shells = depth_shell_ratios(self.depth, None, min_shell_cells=1)
        row = shells[shells.region.eq("MOp") & shells.depth_shell.eq("0200-0300")].iloc[0]
        neurons = row.shell_neurons
        self.assertAlmostEqual(row["Astrocyte per 100 neurons"], 100 * 10 / neurons, places=6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
