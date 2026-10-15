#!/usr/bin/env python3
"""Measured B6J MERFISH MOp/VISp composition, single-animal description only."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from depth import (DEPTH_CAUTION, cortical_depth, depth_cumulative, depth_layer_profile,
                   depth_shell_difference, depth_shell_ratios)

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "input" / "merfish"
RESULTS = ROOT / "results" / "merfish"
FIGURES = ROOT / "figures" / "merfish"
RELEASE = "20250531"
BASE = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
PANEL = ["Slc29a2", "Slc29a1", "Entpd1", "Nt5e", "Adk", "Ada", "Adora1",
         "Adora2a", "Grm2", "Ryr2", "Trpc5", "Atp2b3", "Cacna1c", "Grin2a",
         "Mrps6", "Penk", "Kcnj6", "Shisa7", "Slitrk3", "Dlgap1", "Bcs1l", "Ppargc1a"]
REQUIREMENTS = {
    "annotation": ("MERFISH-C57BL6J-638850", "metadata", "cell_metadata_with_cluster_annotation"),
    "reconstructed": ("MERFISH-C57BL6J-638850-CCF", "metadata", "reconstructed_coordinates"),
    "ccf": ("MERFISH-C57BL6J-638850-CCF", "metadata", "ccf_coordinates"),
    "parcellation": ("Allen-CCF-2020", "metadata", "parcellation_to_parcellation_term_membership_acronym"),
    "gene": ("MERFISH-C57BL6J-638850", "metadata", "gene"),
}
RAW_REQUIREMENT = ("MERFISH-C57BL6J-638850", "expression_matrices", "C57BL6J-638850")
GROUPS = ["Glutamatergic neuron", "GABAergic neuron", "Other neuron", "Astrocyte",
          "Oligodendrocyte", "OPC", "Microglia", "Endothelial", "Pericyte/SMC", "VLMC", "Other non-neuronal"]
LAYERS = ["1", "2/3", "4", "5", "6a", "6b", "unknown"]
CAUTION = "One C57BL/6J male WT donor (638850); sections, hemispheres and cells are not independent animals. No animal-level p values or density estimates."


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def checksum(path: Path, algorithm="md5") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_spec(manifest, triplet, kind=None):
    entry = manifest["file_listing"][triplet[0]][triplet[1]][triplet[2]]
    if kind:
        entry = entry[kind]
    options = entry["files"]
    if len(options) != 1:
        raise ValueError(f"Ambiguous formats for {triplet}: {list(options)}")
    spec = next(iter(options.values())).copy()
    relative = Path(spec["relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Unsafe manifest path")
    if not spec["url"].startswith(BASE + "/"):
        raise ValueError("Manifest URL is not the official Allen bucket")
    spec["local_path"] = str(INPUT / "cache" / relative)
    return spec


def download_verified(spec):
    """Atomic HTTPS download with the unchanged Allen-published MD5 and size."""
    path = Path(spec["local_path"])
    if path.exists() and path.stat().st_size == spec["size"] and checksum(path) == spec["file_hash"]:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    print(f"Download {spec['url']} ({spec['size']:,} bytes)", flush=True)
    try:
        with urllib.request.urlopen(spec["url"], timeout=60) as response, temporary.open("wb") as out:
            shutil.copyfileobj(response, out, length=4 * 1024 * 1024)
        if temporary.stat().st_size != spec["size"]:
            raise ValueError(f"File size mismatch: {path.name}")
        if checksum(temporary) != spec["file_hash"]:
            raise ValueError(f"Allen MD5 mismatch: {path.name}")
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def classify(raw_class, raw_subclass, raw_neurotransmitter):
    """Display-only mapping of Allen 2023 labels. Raw labels are retained."""
    cls, sub, nt = str(raw_class), str(raw_subclass), str(raw_neurotransmitter)
    sub_clean = re.sub(r"^\d+\s+", "", sub)
    non_neuronal = {"OPC NN": "OPC", "Oligo NN": "Oligodendrocyte", "Microglia NN": "Microglia",
                    "Endo NN": "Endothelial", "Peri NN": "Pericyte/SMC", "SMC NN": "Pericyte/SMC", "VLMC NN": "VLMC"}
    if sub_clean in non_neuronal:
        return non_neuronal[sub_clean]
    if sub_clean.startswith("Astro-"):
        return "Astrocyte"
    match = re.match(r"^(\d+)\s", cls)
    is_neuronal_class = bool(match and 1 <= int(match[1]) <= 29)
    # Use the released class ontology for neuronal status. A missing NT is not a glial call.
    if is_neuronal_class:
        if cls.endswith(" Glut") or ("Glut" in nt and "Gaba" not in nt and "GABA" not in nt):
            return "Glutamatergic neuron"
        if cls.endswith(" GABA") or "Gaba" in nt or "GABA" in nt:
            return "GABAergic neuron"
        return "Other neuron"
    return "Other non-neuronal"


def gaba_family(raw_subclass):
    clean = re.sub(r"^\d+\s+", "", str(raw_subclass))
    for family in ["Pvalb", "Sst", "Vip", "Lamp5", "Sncg"]:
        if clean.startswith(family + " "):
            return family
    return "Other GABA"


def layer_from_acronym(region, acronym):
    suffix = str(acronym)[len(region):] if str(acronym).startswith(region) else ""
    return suffix if suffix in LAYERS else "unknown"


def prepare_cells(paths, chunksize=200_000):
    """Filter by CCF before reading wide all-brain annotation into memory."""
    parcel = pd.read_csv(paths["parcellation"], keep_default_na=False)
    parcel = parcel.rename(columns={c: f"parcellation_{c}" for c in parcel if c != "parcellation_index"})
    required = {"parcellation_index", "parcellation_structure", "parcellation_substructure"}
    if not required.issubset(parcel):
        raise ValueError(f"Parcellation schema missing {required - set(parcel)}")
    if parcel["parcellation_index"].duplicated().any():
        raise ValueError("Duplicate parcellation index")
    parcel = parcel[parcel.parcellation_structure.isin(["MOp", "VISp"])].copy()
    indices = set(parcel.parcellation_index)
    selected = []
    n_reconstructed = 0
    for frame in pd.read_csv(paths["reconstructed"], chunksize=chunksize, dtype={"cell_label": str}):
        n_reconstructed += len(frame)
        hit = frame[frame.parcellation_index.isin(indices)].copy()
        selected.append(hit)
    cells = pd.concat(selected, ignore_index=True)
    if cells.empty:
        raise ValueError("No exact MOp/VISp parcellation cells")
    if cells.cell_label.duplicated().any():
        raise ValueError("Duplicate reconstructed cell labels")
    cells = cells.rename(columns={c: f"{c}_reconstructed" for c in ["x", "y", "z"]})
    cells = cells.merge(parcel, on="parcellation_index", how="left", validate="many_to_one")
    selected_ids = set(cells.cell_label)
    frames = []
    for frame in pd.read_csv(paths["ccf"], chunksize=chunksize, dtype={"cell_label": str}):
        frames.append(frame.loc[frame.cell_label.isin(selected_ids), ["cell_label", "x", "y", "z"]])
    ccf = pd.concat(frames, ignore_index=True).rename(columns={c: f"{c}_ccf" for c in ["x", "y", "z"]})
    cells = cells.merge(ccf, on="cell_label", how="left", validate="one_to_one")
    annotation_columns = ["cell_label", "brain_section_label", "cluster_alias", "average_correlation_score",
                          "donor_label", "donor_genotype", "donor_sex", "neurotransmitter", "class",
                          "subclass", "supertype", "cluster", "x", "y", "z"]
    header = pd.read_csv(paths["annotation"], nrows=0)
    required = {"cell_label", "brain_section_label", "cluster_alias", "donor_label", "donor_genotype", "donor_sex", "class", "subclass", "neurotransmitter"}
    if not required.issubset(header):
        raise ValueError(f"Annotation schema missing {required - set(header)}")
    frames = []
    total_annotation = 0
    for frame in pd.read_csv(paths["annotation"], chunksize=chunksize, usecols=[c for c in annotation_columns if c in header],
                             dtype={"cell_label": str}, keep_default_na=False):
        total_annotation += len(frame)
        frames.append(frame[frame.cell_label.isin(selected_ids)])
    annotation = pd.concat(frames, ignore_index=True)
    annotation = annotation.rename(columns={c: f"raw_{c}" for c in ["neurotransmitter", "class", "subclass", "supertype", "cluster"]})
    annotation = annotation.rename(columns={c: f"{c}_section" for c in ["x", "y", "z"]})
    cells = cells.merge(annotation, on="cell_label", how="left", validate="one_to_one", indicator="annotation_join")
    cells["region"] = cells.parcellation_structure
    cells["anatomical_layer"] = [layer_from_acronym(a, b) for a, b in zip(cells.region, cells.parcellation_substructure)]
    # Layer-derived grouping, retained for continuity. It is NOT depth matched between
    # regions: MOp cortex is about twice as thick as VISp, so the same label spans a
    # different physical depth range in each. Use depth_shell for anything dose related.
    cells["depth_group"] = cells.anatomical_layer.map({"1": "superficial", "2/3": "superficial", "4": "middle", "5": "deep", "6a": "deep", "6b": "deep"}).fillna("unknown")
    # Allen CCF ASL, z axis left-to-right, 1140 voxels × 10 μm, midline z=5.7 mm.
    cells["hemisphere"] = np.select([cells.z_ccf < 5.7, cells.z_ccf > 5.7, cells.z_ccf == 5.7], ["left", "right", "midline"], default="unknown")
    # Physical depth below the pial surface, in micrometres. See src/merfish/depth.py.
    cells = cortical_depth(cells)
    for column in ["raw_class", "raw_subclass", "raw_neurotransmitter"]:
        cells[column] = cells[column].fillna("")
    cells["major_group"] = [classify(a, b, c) for a, b, c in zip(cells.raw_class, cells.raw_subclass, cells.raw_neurotransmitter)]
    cells["is_neuron"] = cells.major_group.str.endswith("neuron")
    cells["is_gaba"] = cells.major_group.eq("GABAergic neuron")
    cells["gaba_family"] = cells.raw_subclass.map(gaba_family)
    # No new correlation-score threshold: use released classified metadata; record exclusions.
    cells["qc_annotated"] = cells.annotation_join.eq("both") & cells.raw_class.ne("") & cells.raw_subclass.ne("") & cells.cluster_alias.notna() & cells.cluster_alias.astype(str).ne("")
    n = {"released_annotation_cells": total_annotation, "reconstructed_coordinate_cells": n_reconstructed,
         "selected_MOp_VISp_coordinate_cells": len(cells), "missing_annotation": int(cells.annotation_join.ne("both").sum()),
         "missing_CCF_coordinates": int(cells.x_ccf.isna().sum()), "selected_QC_annotated_cells": int(cells.qc_annotated.sum())}
    retained = cells[cells.qc_annotated].copy()
    donors = retained[["donor_label", "donor_genotype", "donor_sex"]].drop_duplicates()
    if len(donors) != 1 or donors.iloc[0].donor_label != "C57BL6J-638850":
        raise ValueError(f"Expected exactly the specified B6J donor; observed {donors.to_dict('records')}")
    if donors.iloc[0].donor_genotype != "wt/wt" or donors.iloc[0].donor_sex != "M":
        raise ValueError("Unexpected genotype/sex for B6J638850 dataset")
    if retained.groupby("region").size().reindex(["MOp", "VISp"], fill_value=0).min() == 0:
        raise ValueError("Both exact MOp and VISp need retained annotated cells")
    return cells, n


def proportions(cells, dimensions, category, denominator="all_QC_annotated"):
    if denominator == "neurons":
        cells = cells[cells.is_neuron]
    elif denominator == "GABAergic_neurons":
        cells = cells[cells.is_gaba]
    elif denominator != "all_QC_annotated":
        raise ValueError(denominator)
    totals = cells.groupby(dimensions, observed=True, dropna=False).size().rename("denominator_cells")
    counts = cells.groupby(dimensions + [category], observed=True, dropna=False).size().rename("cell_count")
    # Include true zero categories observed elsewhere, rather than silently dropping them.
    categories = sorted(cells[category].dropna().unique())
    expanded = totals.reset_index()[dimensions].merge(pd.DataFrame({category: categories}), how="cross")
    result = expanded.merge(counts.reset_index(), on=dimensions + [category], how="left")
    result["cell_count"] = result.cell_count.fillna(0).astype(int)
    result = result.merge(totals.reset_index(), on=dimensions, validate="many_to_one")
    result["fraction"] = result.cell_count / result.denominator_cells
    result["percentage"] = 100 * result.fraction
    result["denominator"] = denominator
    result["donor_label"] = "C57BL6J-638850"
    result["independent_animals"] = 1
    return result


def write_composition(cells, outdir, figdir):
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    retained = cells[cells.qc_annotated].copy()
    cells.to_csv(outdir / "selected_cells_raw_labels_and_coordinates.csv.gz", index=False)
    mapping = retained[["raw_class", "raw_subclass", "raw_neurotransmitter", "major_group", "is_neuron", "is_gaba", "gaba_family"]].drop_duplicates().sort_values(["raw_class", "raw_subclass"])
    mapping["taxonomy_version"] = "WMB-taxonomy 20231215"
    mapping.to_csv(outdir / "raw_label_display_group_mapping.csv", index=False)
    samples = retained.groupby(["donor_label", "donor_genotype", "donor_sex", "region", "anatomical_layer", "hemisphere", "brain_section_label"], dropna=False).size().reset_index(name="QC_annotated_cells")
    samples.to_csv(outdir / "sample_coverage.csv", index=False)
    dimensions = {"region": ["region"], "layer": ["region", "anatomical_layer"], "hemisphere": ["region", "hemisphere"],
                  "section": ["region", "brain_section_label"], "layer_hemisphere": ["region", "anatomical_layer", "hemisphere"],
                  "layer_section": ["region", "anatomical_layer", "brain_section_label"],
                  "layer_hemisphere_section": ["region", "anatomical_layer", "hemisphere", "brain_section_label"],
                  "depth_group": ["region", "depth_group"],
                  "depth_shell": ["region", "depth_shell"]}
    for name, dim in dimensions.items():
        for category, denominator, suffix in [("major_group", "all_QC_annotated", "major"), ("raw_subclass", "all_QC_annotated", "subclass_all_cells"),
                                              ("raw_subclass", "neurons", "subclass_neurons"), ("gaba_family", "GABAergic_neurons", "gaba_within_gaba")]:
            proportions(retained, dim, category, denominator).to_csv(outdir / f"composition_{name}_{suffix}.csv", index=False)
    gaba_all = proportions(retained, ["region"], "raw_subclass")
    gaba_labels = set(retained.loc[retained.is_gaba, "raw_subclass"])
    gaba_all[gaba_all.raw_subclass.isin(gaba_labels)].to_csv(outdir / "gaba_subclasses_fraction_of_all_QC_cells.csv", index=False)
    for category, denominator, label in [("major_group", "all_QC_annotated", "major"), ("raw_subclass", "neurons", "neuron_subclass"), ("gaba_family", "GABAergic_neurons", "gaba_family")]:
        table = proportions(retained, ["region"], category, denominator)
        wide = table.pivot(index=category, columns="region", values="percentage").reindex(columns=["MOp", "VISp"])
        wide["MOp_minus_VISp_percentage_points"] = wide.MOp - wide.VISp
        wide["denominator"] = denominator
        wide.to_csv(outdir / f"regional_descriptive_difference_{label}.csv")
    write_depth_tables(retained, outdir)
    plot_composition(retained, figdir)


def write_depth_tables(retained, outdir):
    """Physical depth tables. Percentages of all QC cells share one denominator, so a
    shift in the neuron share moves every other class mechanically; the per-100-neurons
    table removes that closure effect and is the one to read for dose arguments."""
    depth_layer_profile(retained).to_csv(outdir / "depth_layer_profile.csv", index=False)
    depth_cumulative(retained).to_csv(outdir / "depth_cumulative_distribution.csv", index=False)
    shells = depth_shell_ratios(retained, proportions)
    shells.to_csv(outdir / "composition_depth_shell_per100_neurons.csv", index=False)
    depth_shell_difference(shells).to_csv(outdir / "depth_shell_regional_difference.csv", index=False)
    (outdir / "depth_caution.txt").write_text(DEPTH_CAUTION + "\n")


def plot_composition(cells, figdir):
    cache = ROOT / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    table = proportions(cells, ["region"], "major_group")
    wide = table.pivot(index="region", columns="major_group", values="percentage").reindex(index=["MOp", "VISp"], columns=GROUPS).fillna(0)
    fig, ax = plt.subplots(figsize=(10, 5))
    wide.plot.bar(stacked=True, ax=ax, colormap="tab20", rot=0)
    totals = cells.groupby("region").size()
    ax.set_xticklabels([f"{r}\nn={totals[r]:,} sampled cells" for r in wide.index])
    ax.set(ylabel="% of released QC annotated cells", xlabel="", ylim=(0, 100), title="B6J638850: regional sampled-cell composition (one animal)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout(); fig.savefig(figdir / "regional_major_composition.png", dpi=180); plt.close(fig)
    table = proportions(cells, ["region", "anatomical_layer"], "major_group")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    for ax, region in zip(axes, ["MOp", "VISp"]):
        data = table[table.region.eq(region)].pivot(index="major_group", columns="anatomical_layer", values="percentage").reindex(index=GROUPS, fill_value=0).reindex(columns=LAYERS)
        im = ax.imshow(np.ma.masked_invalid(data.to_numpy(dtype=float)), vmin=0, vmax=100, cmap="viridis", aspect="auto")
        totals = cells[cells.region.eq(region)].groupby("anatomical_layer").size()
        ax.set_xticks(range(len(LAYERS)), [f"L{x}\nn={totals.get(x, 0):,}" for x in LAYERS], rotation=40, ha="right")
        ax.set_yticks(range(len(GROUPS)), GROUPS if region == "MOp" else [""]*len(GROUPS))
        ax.set_title(region); ax.set_xlabel("CCF anatomical layer (blank = not sampled)")
    fig.colorbar(im, ax=axes, label="% of QC annotated cells within layer", shrink=.8)
    fig.suptitle("Single-animal layer composition; transcriptomic subclasses are separate")
    fig.savefig(figdir / "layer_major_composition.png", dpi=180); plt.close(fig)
    table = proportions(cells, ["region"], "raw_subclass", "neurons")
    wide = table.pivot(index="raw_subclass", columns="region", values="percentage").reindex(columns=["MOp", "VISp"]).fillna(0)
    wide = wide.loc[wide.max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(10, max(5, len(wide)*.27)))
    y = np.arange(len(wide))
    ax.plot(wide.MOp, y, "o", label="MOp"); ax.plot(wide.VISp, y, "s", label="VISp")
    ax.set_yticks(y, wide.index, fontsize=8); ax.invert_yaxis()
    ax.set(xlabel="% of annotated neurons in region", title="Released raw neuronal subclasses; one B6J animal")
    ax.legend(); fig.tight_layout(); fig.savefig(figdir / "regional_neuron_subclasses.png", dpi=180); plt.close(fig)


def read_h5_strings(group, key):
    obj = group[key]
    if hasattr(obj, "shape"):
        return np.asarray(obj.asstr()[:], dtype=str)
    codes = obj["codes"][:]
    categories = np.asarray(obj["categories"].asstr()[:], dtype=str)
    return np.where(codes >= 0, categories[np.maximum(codes, 0)], "")


def measured_penk(cells, raw_path, outdir, expected_panel_size=500):
    """Raw detection and explicitly defined measured-panel CP10K log1p."""
    import h5py
    from scipy.sparse import csc_matrix
    selected = cells[cells.qc_annotated].copy().set_index("cell_label")
    if not selected.index.is_unique:
        raise ValueError("Duplicate selected expression observations")
    with h5py.File(raw_path, "r") as f:
        obskey = f["obs"].attrs.get("_index", "_index")
        varkey = f["var"].attrs.get("_index", "_index")
        all_obs = read_h5_strings(f["obs"], obskey)
        var_ids = read_h5_strings(f["var"], varkey)
        genes = read_h5_strings(f["var"], "gene_symbol") if "gene_symbol" in f["var"] else var_ids
        if len(genes) != expected_panel_size or len(set(genes)) != len(genes):
            raise ValueError(f"Expected {expected_panel_size} unique measured panel genes, found {len(genes)} features")
        if np.count_nonzero(genes == "Penk") != 1:
            raise ValueError("Raw matrix must have exactly one measured Penk symbol")
        gene_index = int(np.flatnonzero(genes == "Penk")[0])
        positions = pd.Series(np.arange(len(all_obs)), index=all_obs)
        if positions.index.duplicated().any():
            raise ValueError("Duplicate raw matrix observations")
        missing = selected.index.difference(positions.index)
        if len(missing):
            raise ValueError(f"{len(missing)} selected cells absent from raw expression")
        wanted = positions.loc[selected.index].to_numpy(dtype=np.int64)
        sorted_order = np.argsort(wanted)
        sorted_rows = wanted[sorted_order]
        values = np.empty(len(wanted), dtype=float); totals = np.empty(len(wanted), dtype=float)
        x = f["X"]
        if hasattr(x, "shape"):
            for offset in range(0, len(sorted_rows), 1024):
                rows = sorted_rows[offset:offset+1024]
                batch = x[rows, :]
                target = sorted_order[offset:offset+len(rows)]
                values[target] = batch[:, gene_index]
                totals[target] = batch.sum(axis=1)
        else:
            encoding = x.attrs.get("encoding-type", "")
            if isinstance(encoding, bytes): encoding = encoding.decode()
            shape = tuple(x.attrs["shape"])
            if encoding == "csr_matrix":
                pointers = x["indptr"][:]
                for offset in range(0, len(sorted_rows), 1024):
                    rows = sorted_rows[offset:offset+1024]
                    for local, row in enumerate(rows):
                        start, end = int(pointers[row]), int(pointers[row+1])
                        dat = x["data"][start:end]; col = x["indices"][start:end]
                        target = sorted_order[offset+local]
                        totals[target] = dat.sum()
                        values[target] = dat[col == gene_index].sum()
            elif encoding == "csc_matrix":
                matrix = csc_matrix((x["data"][:], x["indices"][:], x["indptr"][:]), shape=shape)
                subset = matrix[wanted, :]
                values[:] = subset[:, gene_index].toarray().ravel()
                totals[:] = np.asarray(subset.sum(axis=1)).ravel()
            else:
                raise ValueError(f"Unsupported H5AD X encoding: {encoding}")
    if not np.isfinite(values).all() or (values < 0).any() or not np.isfinite(totals).all():
        raise ValueError("Nonfinite/negative measured raw expression")
    selected["Penk_raw_count"] = values
    selected["measured_panel_total_count"] = totals
    selected["Penk_detected"] = values > 0
    selected["Penk_log1p_panel_CP10K"] = np.where(totals > 0, np.log1p(values / np.maximum(totals, 1) * 10000), np.nan)
    selected.to_csv(outdir / "measured_Penk_cells.csv.gz")
    for name, dimensions in {"region": ["region"], "major_group": ["region", "major_group"], "subclass": ["region", "raw_subclass"], "layer": ["region", "anatomical_layer"], "subclass_layer": ["region", "raw_subclass", "anatomical_layer"]}.items():
        aggregate = selected.groupby(dimensions).agg(cells=("Penk_detected", "size"), detected_cells=("Penk_detected", "sum"), detection_fraction=("Penk_detected", "mean"), mean_raw_count=("Penk_raw_count", "mean"), mean_log1p_panel_CP10K=("Penk_log1p_panel_CP10K", "mean"), zero_panel_cells=("measured_panel_total_count", lambda x: int(x.eq(0).sum())))
        aggregate.to_csv(outdir / f"measured_Penk_{name}.csv")
    return {"status": "completed", "measurement": "raw MERFISH counts", "normalization": f"log1p(Penk / summed measured {expected_panel_size}-gene panel counts * 10000)", "measured_panel_features": expected_panel_size, "detection": "raw Penk count > 0", "independent_animals": 1}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Download required official files after network access is configured")
    parser.add_argument("--with-penk", action="store_true", help="Also download/analyze measured raw 7.63 GB H5AD")
    parser.add_argument("--manifest", type=Path, default=INPUT / "official_sdk_manifest_snapshot.json")
    args = parser.parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True)
    report = {"run_time_UTC": utc_now(), "status": "pending", "dataset": "MERFISH-C57BL6J-638850", "caution": CAUTION,
              "actual_analysis_computed": False, "manifest_origin": "Official AllenInstitute/abc_atlas_access GitHub test-resource snapshot; not yet fetched from S3", "pending": [], "files": {}}
    try:
        if not args.manifest.exists():
            raise FileNotFoundError(f"Manifest not found: {args.manifest}")
        manifest = json.loads(args.manifest.read_text())
        report["manifest_version"] = manifest["version"]
        report["manifest_sha256"] = checksum(args.manifest, "sha256")
        specs = {name: file_spec(manifest, triplet) for name, triplet in REQUIREMENTS.items()}
        if args.with_penk:
            specs["raw"] = file_spec(manifest, RAW_REQUIREMENT, kind="raw")
        report["download_bytes_if_all_absent"] = sum(s["size"] for s in specs.values())
        if args.download:
            report["manifest_endpoint"] = f"{BASE}/releases/{manifest['version']}/manifest.json"
            # Require the official live manifest to agree before downloading multi-GB artifacts.
            with urllib.request.urlopen(report["manifest_endpoint"], timeout=30) as response:
                live = json.load(response)
            live_path = INPUT / f"live_manifest_{live['version']}.json"
            live_path.write_text(json.dumps(live, indent=2))
            specs = {name: file_spec(live, triplet) for name, triplet in REQUIREMENTS.items()}
            if args.with_penk: specs["raw"] = file_spec(live, RAW_REQUIREMENT, kind="raw")
            report["manifest_origin"] = str(live_path.relative_to(ROOT))
            report["manifest_sha256"] = checksum(live_path, "sha256")
            required_disk = sum(s["size"] for s in specs.values() if not Path(s["local_path"]).exists()) + 512 * 1024**2
            if shutil.disk_usage(INPUT).free < required_disk:
                raise OSError(f"Need at least {required_disk:,} free bytes for requested verified files")
            for spec in specs.values(): download_verified(spec)
        for name, spec in specs.items():
            path = Path(spec["local_path"])
            present = path.exists()
            verified = present and path.stat().st_size == spec["size"] and checksum(path) == spec["file_hash"]
            report["files"][name] = {**spec, "local_path": str(path.relative_to(ROOT)), "present": present, "verified_published_MD5": verified}
            if not verified: report["pending"].append(f"Missing or unverified input: {name}")
        if report["pending"]:
            return 2
        paths = {name: Path(spec["local_path"]) for name, spec in specs.items()}
        genes = pd.read_csv(paths["gene"], keep_default_na=False)
        if "gene_symbol" not in genes: raise ValueError("Gene metadata missing gene_symbol")
        pd.DataFrame({"gene": PANEL, "measured_in_B6J_panel": [g in set(genes.gene_symbol) for g in PANEL], "not_measured_is_not_zero": True}).to_csv(RESULTS / "candidate_panel_measured_coverage.csv", index=False)
        cells, audit = prepare_cells(paths)
        write_composition(cells, RESULTS, FIGURES)
        report["cell_audit"] = audit
        report["Penk"] = measured_penk(cells, paths["raw"], RESULTS) if args.with_penk else {"status": "pending_optional_raw_matrix", "not_measured_genes": "All other 21 candidates are unmeasured in B6J500 panel"}
        report["status"] = "composition_completed" if not args.with_penk else "composition_and_measured_Penk_completed"
        report["actual_analysis_computed"] = True
        return 0
    except (OSError, ValueError, KeyError, urllib.error.URLError) as error:
        report["pending"].append(f"{type(error).__name__}: {error}")
        print(report["pending"][-1], file=sys.stderr)
        return 2
    finally:
        (RESULTS / "readiness_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
        print(json.dumps({"status": report["status"], "actual_analysis_computed": report["actual_analysis_computed"], "pending": report["pending"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
