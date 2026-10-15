#!/usr/bin/env python3
"""Synthetic tests for the Zhuang analysis logic. No network and no Allen data.

The point of these cases is the claim the pipeline is built to support: a region-level
technical offset shifts every cell class the same way, so it cannot produce a high
cell class discordance score, while a genuine class-specific pattern can.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_zhuang import (cell_class, clean_subclass, cross_animal_agreement, discordance_score,
                        within_subclass_difference)

GABA = ["Pvalb Gaba", "Sst Gaba", "Vip Gaba", "Lamp5 Gaba"]
GLUT = ["L2/3 IT CTX Glut", "L4/5 IT CTX Glut", "L5 IT CTX Glut", "L6 CT CTX Glut"]


def synthetic(animals=("A1", "A2"), n=60, batch_shift=0.0, seed=0):
    """Three genes.

    dissociated   inhibitory subclasses up in MOp, excitatory subclasses down.
    uniform       every subclass up in MOp by the same amount, a region-level offset.
    flat          no regional difference at all.
    """
    rng = np.random.default_rng(seed)
    meta, values = [], []
    for animal in animals:
        for subclass in GABA + GLUT:
            is_gaba = subclass in GABA
            for region in ("MOp", "VISp"):
                bump = 1.0 if region == "MOp" else 0.0
                for i in range(n):
                    meta.append({"cell_label": f"{animal}-{subclass}-{region}-{i}", "animal": animal,
                                 "region": region, "subclass": subclass,
                                 "class": "01 GABA" if is_gaba else "02 Glut",
                                 "subclass_clean": clean_subclass(subclass),
                                 "cell_class": "GABAergic" if is_gaba else "Glutamatergic"})
                    values.append({
                        "dissociated": (0.4 if is_gaba else -0.4) * bump + rng.normal(0, 0.01),
                        "uniform": 0.4 * bump + batch_shift * bump + rng.normal(0, 0.01),
                        "flat": rng.normal(0, 0.01),
                    })
    cells = pd.DataFrame(meta)
    expression = pd.DataFrame(values, index=pd.Index(cells.cell_label, name="cell_label"))
    return cells, expression


class ZhuangLogicTests(unittest.TestCase):
    def setUp(self):
        self.cells, self.expression = synthetic()
        self.difference, self.dropped = within_subclass_difference(
            self.cells, self.expression, min_cells=10)
        self.scores = discordance_score(self.difference)

    def test_class_labels(self):
        self.assertEqual(cell_class("01 CTX-CGE GABA", "053 Lamp5 Gaba"), "GABAergic")
        self.assertEqual(cell_class("02 IT-ET Glut", "005 L2/3 IT CTX Glut"), "Glutamatergic")
        self.assertEqual(clean_subclass("053 Lamp5 Gaba"), "Lamp5 Gaba")

    def test_only_shared_subclasses_are_compared(self):
        """A subclass present in one region only must be dropped, never imputed."""
        trimmed = self.cells[~(self.cells.subclass.eq("Vip Gaba") & self.cells.region.eq("VISp"))]
        expression = self.expression.loc[trimmed.cell_label]
        difference, dropped = within_subclass_difference(trimmed, expression, min_cells=10)
        self.assertNotIn("Vip Gaba", set(difference.subclass_clean))
        self.assertIn("Vip Gaba", set(dropped.subclass_clean))

    def test_small_strata_are_dropped_with_a_reason(self):
        difference, dropped = within_subclass_difference(self.cells, self.expression, min_cells=1000)
        self.assertTrue(difference.empty)
        self.assertTrue((dropped.reason.str.contains("fewer than")).all())

    def test_dissociated_gene_scores_high_uniform_gene_scores_zero(self):
        pooled = self.scores.groupby("gene").discordance.mean()
        self.assertAlmostEqual(pooled["dissociated"], 1.0, places=6)
        self.assertAlmostEqual(pooled["uniform"], 0.0, places=6)

    def test_region_level_offset_cannot_create_discordance(self):
        """The confound argument, as an executable check.

        A batch effect that shifts one whole region must leave the discordance score of a
        uniformly expressed gene at zero no matter how large the shift is.
        """
        for shift in (0.0, 2.0, 50.0):
            cells, expression = synthetic(batch_shift=shift, seed=3)
            difference, _ = within_subclass_difference(cells, expression, min_cells=10)
            scores = discordance_score(difference)
            value = scores[scores.gene.eq("uniform")].discordance.mean()
            self.assertAlmostEqual(value, 0.0, places=6, msg=f"shift={shift}")

    def test_ranks_are_computed_against_the_whole_panel(self):
        self.assertTrue((self.scores.panel_genes_ranked == 3).all())
        best = self.scores[self.scores.gene.eq("dissociated")]
        self.assertTrue((best.discordance_rank == 1).all())

    def test_cross_animal_agreement_counts_animals(self):
        agreement = cross_animal_agreement(self.difference)
        row = agreement[agreement.gene.eq("dissociated") & agreement.subclass_clean.eq("Pvalb Gaba")].iloc[0]
        self.assertEqual(row.animals, 2)
        self.assertTrue(row.all_animals_agree)
        self.assertEqual(row.animals_positive, 2)

    def test_flat_gene_does_not_agree_across_animals(self):
        agreement = cross_animal_agreement(self.difference)
        flat = agreement[agreement.gene.eq("flat")]
        self.assertLess(flat.all_animals_agree.mean(), 1.0)

    def test_difference_direction_matches_construction(self):
        row = self.difference[self.difference.gene.eq("dissociated")
                              & self.difference.cell_class.eq("GABAergic")]
        self.assertTrue((row.MOp_minus_VISp > 0).all())
        row = self.difference[self.difference.gene.eq("dissociated")
                              & self.difference.cell_class.eq("Glutamatergic")]
        self.assertTrue((row.MOp_minus_VISp < 0).all())


if __name__ == "__main__":
    unittest.main(verbosity=2)
