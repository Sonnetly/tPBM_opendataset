"""Synthetic-only correctness checks; these cells are never study results."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("run_merfish", Path(__file__).with_name("run_merfish.py"))
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)


class SyntheticMerfishTests(unittest.TestCase):
    def test_raw_ontology_mapping_is_not_marker_reclassification(self):
        self.assertEqual(pipeline.classify("31 OPC-Oligo", "326 OPC NN", "NA"), "OPC")
        self.assertEqual(pipeline.classify("31 OPC-Oligo", "327 Oligo NN", "NA"), "Oligodendrocyte")
        self.assertEqual(pipeline.classify("30 Astro-Epen", "321 Astroependymal NN", "NA"), "Other non-neuronal")
        self.assertEqual(pipeline.classify("07 CTX-MGE GABA", "051 Pvalb chandelier Gaba", ""), "GABAergic neuron")
        self.assertEqual(pipeline.classify("01 IT-ET Glut", "007 L2/3 IT CTX Glut", ""), "Glutamatergic neuron")
        self.assertEqual(pipeline.gaba_family("056 Sst Chodl Gaba"), "Sst")
        self.assertEqual(pipeline.layer_from_acronym("MOp", "MOp2/3"), "2/3")
        self.assertEqual(pipeline.layer_from_acronym("MOp", "MOs2/3"), "unknown")

    def test_denominators_and_zero_categories(self):
        cells = pd.DataFrame({"region": ["MOp"]*4 + ["VISp"]*2,
                              "major_group": ["Glutamatergic neuron", "GABAergic neuron", "GABAergic neuron", "Astrocyte", "Glutamatergic neuron", "Astrocyte"],
                              "raw_subclass": ["Glut", "Sst", "Vip", "Astro", "Glut", "Astro"],
                              "is_neuron": [True, True, True, False, True, False],
                              "is_gaba": [False, True, True, False, False, False],
                              "gaba_family": ["Other GABA", "Sst", "Vip", "Other GABA", "Other GABA", "Other GABA"]})
        whole = pipeline.proportions(cells, ["region"], "raw_subclass")
        sst = whole[whole.raw_subclass.eq("Sst")].set_index("region")
        self.assertEqual(sst.loc["MOp", "fraction"], .25)
        self.assertEqual(sst.loc["VISp", "fraction"], 0)
        neurons = pipeline.proportions(cells, ["region"], "raw_subclass", "neurons")
        self.assertAlmostEqual(neurons.query('region == "MOp" and raw_subclass == "Sst"').iloc[0].fraction, 1/3)
        gaba = pipeline.proportions(cells, ["region"], "gaba_family", "GABAergic_neurons")
        self.assertEqual(gaba.query('gaba_family == "Sst"').iloc[0].fraction, .5)
        np.testing.assert_allclose(whole.groupby("region").fraction.sum(), [1, 1])

    def make_fixture(self, folder):
        parcel = pd.DataFrame({"parcellation_index": [1, 2, 3], "structure": ["MOp", "VISp", "MOs"], "substructure": ["MOp2/3", "VISp4", "MOs2/3"]})
        reconstructed = pd.DataFrame({"cell_label": ["a", "b", "c", "d", "outside"], "parcellation_index": [1, 1, 2, 2, 3], "x": [1.]*5, "y": [2.]*5, "z": [3.]*5})
        ccf = pd.DataFrame({"cell_label": ["a", "b", "c", "d", "outside"], "parcellation_index": [1, 1, 2, 2, 3], "x": [1.]*5, "y": [2.]*5, "z": [4., 6., 4., 6., 4.]})
        annotation = pd.DataFrame({"cell_label": ["a", "b", "c", "d", "outside"], "brain_section_label": ["section1", "section1", "section2", "section2", "section3"],
                                   "cluster_alias": [1, 2, 1, 2, 1], "average_correlation_score": [.2, .9, .8, .8, .9],
                                   "donor_label": ["C57BL6J-638850"]*5, "donor_genotype": ["wt/wt"]*5, "donor_sex": ["M"]*5,
                                   "class": ["01 IT-ET Glut", "07 CTX-MGE GABA", "01 IT-ET Glut", "07 CTX-MGE GABA", "01 IT-ET Glut"],
                                   "subclass": ["007 L2/3 IT CTX Glut", "053 Sst Gaba", "006 L4/5 IT CTX Glut", "052 Pvalb Gaba", "007 L2/3 IT CTX Glut"],
                                   "neurotransmitter": ["Glut", "GABA", "Glut", "GABA", "Glut"], "supertype": ["synthetic"]*5, "cluster": ["synthetic"]*5})
        paths = {}
        for name, frame in {"parcellation": parcel, "reconstructed": reconstructed, "ccf": ccf, "annotation": annotation}.items():
            paths[name] = folder / f"synthetic_{name}.csv"
            frame.to_csv(paths[name], index=False)
        return paths

    def test_exact_region_join_layer_hemisphere_and_figures(self):
        with tempfile.TemporaryDirectory(prefix="merfish_synthetic_") as temp:
            folder = Path(temp)
            cells, audit = pipeline.prepare_cells(self.make_fixture(folder), chunksize=2)
            self.assertEqual(set(cells.cell_label), {"a", "b", "c", "d"})
            self.assertEqual(set(cells.anatomical_layer), {"2/3", "4"})
            self.assertEqual(set(cells.hemisphere), {"left", "right"})
            self.assertEqual(audit["selected_QC_annotated_cells"], 4)
            # A low correlation score is retained: no invented threshold.
            self.assertTrue(cells.set_index("cell_label").loc["a", "qc_annotated"])
            pipeline.write_composition(cells, folder/"synthetic_results", folder/"synthetic_figures")
            self.assertTrue((folder/"synthetic_figures"/"regional_major_composition.png").exists())
            self.assertTrue((folder/"synthetic_results"/"composition_layer_hemisphere_section_subclass_neurons.csv").exists())

    def test_measured_raw_penk_normalization_and_missing_not_zero(self):
        with tempfile.TemporaryDirectory(prefix="merfish_synthetic_") as temp:
            folder = Path(temp)
            cells, _ = pipeline.prepare_cells(self.make_fixture(folder), chunksize=2)
            raw = folder/"synthetic_raw.h5ad"
            with h5py.File(raw, "w") as f:
                dt = h5py.string_dtype("utf-8")
                obs = f.create_group("obs"); obs.attrs["_index"] = "_index"
                obs.create_dataset("_index", data=np.array(["c", "a", "d", "b"], dtype=object), dtype=dt)
                var = f.create_group("var"); var.attrs["_index"] = "_index"
                var.create_dataset("_index", data=np.array(["gene1", "gene2"], dtype=object), dtype=dt)
                var.create_dataset("gene_symbol", data=np.array(["Penk", "Other"], dtype=object), dtype=dt)
                f.create_dataset("X", data=np.array([[2., 2.], [0., 2.], [0., 0.], [1., 3.]]))
            status = pipeline.measured_penk(cells, raw, folder, expected_panel_size=2)
            self.assertEqual(status["independent_animals"], 1)
            expression = pd.read_csv(folder/"measured_Penk_cells.csv.gz").set_index("cell_label")
            self.assertEqual(expression.loc["a", "Penk_raw_count"], 0)
            self.assertAlmostEqual(expression.loc["c", "Penk_log1p_panel_CP10K"], np.log1p(5000))
            self.assertTrue(np.isnan(expression.loc["d", "Penk_log1p_panel_CP10K"]))


if __name__ == "__main__":
    unittest.main()
