#!/usr/bin/env python3
"""Zhuang-ABCA MERFISH MOp/VISp within-subclass expression, with a panel-wide null.

Why this dataset
----------------
MERFISH-C57BL6J-638850 is one animal and its 500-gene panel measures only 7 of the
62 PBM candidates. It does not measure Adora1, which is the gene where the CTX-HPF
10x analysis found the clearest pattern. The Zhuang-ABCA MERFISH release covers four
animals on a roughly 1,122-gene panel aligned to the same WMB taxonomy and the same
CCF, and that panel does contain Adora1 and Adora2a. It is the only public option
that relaxes the single-animal limit and the measurement gap at the same time.

What this script computes
-------------------------
1. Panel coverage. Which of the 62 candidates are measured, stated before any result.
2. Within-subclass MOp minus VISp difference for EVERY measured gene, per animal.
   Only cell subclasses present in both regions are compared, so a difference in
   subclass composition cannot masquerade as a difference in expression.
3. A cell-class discordance score per gene. The 10x observation was that Adora1 runs
   one way in every GABAergic subclass and the other way in every glutamatergic
   subclass. A region-level batch or detection-efficiency difference is a scalar and
   cannot flip sign between two cell classes measured in the same sections, so this
   statistic is partly protected against the confound that blocks a plain regional
   comparison. The score is computed for all measured genes, which makes the 62
   candidates a prespecified subset inside a panel-wide null rather than a hand-picked
   list. A gene's rank in that null is reportable without an animal-level p value.
4. Between-animal agreement. Four animals allow sign concordance to be counted
   directly, which the single-animal dataset never could.

Inference scope
---------------
Descriptive. Four animals is enough to report per-animal agreement and is not enough
to claim a population effect size. No PBM treatment is present in this atlas, so
nothing here measures a PBM response. Panel absence means unmeasured, never zero.

Usage
-----
    python src/merfish/run_zhuang.py --list                 # show manifest entries only
    python src/merfish/run_zhuang.py --download             # fetch and analyse
    python src/merfish/run_zhuang.py --download --animals Zhuang-ABCA-1 Zhuang-ABCA-2
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "input" / "zhuang"
RESULTS = ROOT / "results" / "zhuang"
FIGURES = ROOT / "figures" / "zhuang"
BASE = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
RELEASE = "20250531"
# Coronal animals first: MOp and VISp are both sampled in coronal series.
DEFAULT_ANIMALS = ("Zhuang-ABCA-1", "Zhuang-ABCA-2")
CANDIDATE_CSV = ROOT / "results" / "pbm" / "exploration_genes62.csv"
# The MERFISH run already stores the official release manifest snapshot and it covers
# every directory in the release, Zhuang included. Planning works offline from it; a
# --download run still refetches the live manifest and verifies every published MD5.
MANIFEST_SNAPSHOT = ROOT / "input" / "merfish" / "official_sdk_manifest_snapshot.json"
GABA_SUBCLASS_HINTS = ("Pvalb", "Sst", "Vip", "Lamp5", "Sncg", "Sst Chodl")
MIN_CELLS_PER_STRATUM = 50

CAUTION = (
    "Zhuang-ABCA MERFISH, normal untreated brains. Descriptive within-subclass regional "
    "comparison only. No PBM condition is present, so nothing here estimates a PBM response. "
    "Panel absence means unmeasured, not zero. Per-animal agreement is reported instead of an "
    "animal-level p value."
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def checksum(path: Path, algorithm="md5") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_manifest(release=RELEASE, snapshot=MANIFEST_SNAPSHOT, live=True):
    """Live manifest when the network allows it, otherwise the stored snapshot.

    Returns (manifest, path, origin). Planning and entry resolution are identical either
    way. Only a --download run needs the live copy, because the published MD5 that every
    download is checked against must come from the release itself.
    """
    if live:
        url = f"{BASE}/releases/{release}/manifest.json"
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                manifest = json.load(response)
            INPUT.mkdir(parents=True, exist_ok=True)
            path = INPUT / f"live_manifest_{manifest['version']}.json"
            path.write_text(json.dumps(manifest, indent=2))
            return manifest, path, "live S3 release manifest"
        except (urllib.error.URLError, OSError, ValueError) as error:
            if not snapshot.exists():
                raise
            print(f"Live manifest unavailable ({type(error).__name__}); using stored snapshot",
                  file=sys.stderr)
    if not snapshot.exists():
        raise FileNotFoundError(f"No manifest snapshot at {snapshot}")
    return json.loads(snapshot.read_text()), snapshot, "stored snapshot, not refetched"


def describe_manifest(manifest, animals=None):
    """What each Zhuang directory offers, with byte sizes.

    Entry names have changed across releases, so the script resolves them from the
    manifest rather than hard-coding a path that may silently point at the wrong file.
    """
    listing = manifest["file_listing"]
    out = {}
    for directory in sorted(k for k in listing if k.startswith("Zhuang")):
        sections = {}
        for section, entries in listing[directory].items():
            for name, entry in entries.items():
                kinds = [k for k in ("raw", "log2") if k in entry] or [None]
                for kind in kinds:
                    spec = next(iter((entry[kind] if kind else entry)["files"].values()))
                    label = f"{name}/{kind}" if kind else name
                    sections[f"{section}/{label}"] = {"bytes": spec["size"],
                                                      "format": Path(spec["relative_path"]).suffix}
        out[directory] = sections
    return out


def download_plan(manifest, animals):
    """Bytes this run would fetch, per animal, before anything is downloaded."""
    plan, total = {}, 0
    for animal in animals:
        specs = {
            "annotation": resolve(manifest, animal, "metadata", "cell_metadata_with_cluster_annotation"),
            "ccf": resolve(manifest, f"{animal}-CCF", "metadata", "ccf_coordinates",
                           "cell_metadata_with_parcellation_annotation"),
            "expression": resolve(manifest, animal, "expression_matrices", animal, kind="log2"),
        }
        plan[animal] = {name: {"entry": spec["manifest_entry"], "bytes": spec["size"],
                               "cached": Path(spec["local_path"]).exists()}
                        for name, spec in specs.items()}
        total += sum(spec["size"] for spec in specs.values())
    parcellation = resolve(manifest, "Allen-CCF-2020", "metadata",
                           "parcellation_to_parcellation_term_membership_acronym")
    plan["Allen-CCF-2020"] = {"parcellation": {"entry": parcellation["manifest_entry"],
                                               "bytes": parcellation["size"],
                                               "cached": Path(parcellation["local_path"]).exists()}}
    total += parcellation["size"]
    plan["total_bytes_if_nothing_cached"] = total
    plan["total_GB_if_nothing_cached"] = round(total / 1024 ** 3, 2)
    return plan


def resolve(manifest, directory, section, *name_candidates, kind=None):
    """Pick the first manifest entry that exists, and say clearly what was available."""
    entries = manifest["file_listing"].get(directory, {}).get(section, {})
    for name in name_candidates:
        if name in entries:
            entry = entries[name]
            if kind:
                if kind not in entry:
                    continue
                entry = entry[kind]
            options = entry["files"]
            spec = next(iter(options.values())).copy()
            relative = Path(spec["relative_path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Unsafe manifest path")
            if not spec["url"].startswith(BASE + "/"):
                raise ValueError("Manifest URL is not the official Allen bucket")
            spec["local_path"] = str(INPUT / "cache" / relative)
            spec["manifest_entry"] = f"{directory}/{section}/{name}" + (f"/{kind}" if kind else "")
            return spec
    raise KeyError(
        f"None of {name_candidates} under {directory}/{section}. Available: {sorted(entries)}")


def download_verified(spec):
    path = Path(spec["local_path"])
    if path.exists() and path.stat().st_size == spec["size"] and checksum(path) == spec["file_hash"]:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    print(f"Download {spec['url']} ({spec['size']:,} bytes)", flush=True)
    try:
        with urllib.request.urlopen(spec["url"], timeout=120) as response, temporary.open("wb") as out:
            shutil.copyfileobj(response, out, length=4 * 1024 * 1024)
        if temporary.stat().st_size != spec["size"]:
            raise ValueError(f"Size mismatch: {path.name}")
        if checksum(temporary) != spec["file_hash"]:
            raise ValueError(f"Allen MD5 mismatch: {path.name}")
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def candidate_genes():
    frame = pd.read_csv(CANDIDATE_CSV)
    column = next(c for c in frame.columns if c.lstrip("﻿").lower() == "gene")
    return [str(g) for g in frame[column]]


def clean_subclass(value):
    return str(value).split(" ", 1)[1] if str(value)[:1].isdigit() else str(value)


def cell_class(raw_class, raw_subclass):
    """Coarse excitatory / inhibitory call from the released WMB labels."""
    cls, sub = str(raw_class), clean_subclass(raw_subclass)
    if cls.endswith(" GABA") or sub.startswith(GABA_SUBCLASS_HINTS):
        return "GABAergic"
    if cls.endswith(" Glut"):
        return "Glutamatergic"
    return "Other"


def select_cells(manifest, animal, chunksize=200_000):
    """MOp and VISp cells of one Zhuang animal, joined to taxonomy and CCF parcellation."""
    parcellation = resolve(manifest, "Allen-CCF-2020", "metadata",
                           "parcellation_to_parcellation_term_membership_acronym")
    annotation = resolve(manifest, animal, "metadata", "cell_metadata_with_cluster_annotation")
    ccf = resolve(manifest, f"{animal}-CCF", "metadata", "ccf_coordinates", "cell_metadata_with_parcellation_annotation")
    paths = {name: download_verified(spec) for name, spec in
             {"parcellation": parcellation, "annotation": annotation, "ccf": ccf}.items()}

    parcel = pd.read_csv(paths["parcellation"], keep_default_na=False)
    parcel = parcel.rename(columns={c: f"parcellation_{c}" for c in parcel if c != "parcellation_index"})
    parcel = parcel[parcel.parcellation_structure.isin(["MOp", "VISp"])]
    indices = set(parcel.parcellation_index)

    frames = []
    for frame in pd.read_csv(paths["ccf"], chunksize=chunksize, dtype={"cell_label": str}, low_memory=False):
        if "parcellation_index" not in frame.columns:
            raise ValueError(f"{paths['ccf'].name} has no parcellation_index; columns {list(frame.columns)[:12]}")
        frames.append(frame[frame.parcellation_index.isin(indices)])
    cells = pd.concat(frames, ignore_index=True)
    if cells.empty:
        raise ValueError(f"{animal}: no MOp/VISp cells in the CCF parcellation")
    cells = cells.merge(parcel, on="parcellation_index", how="left", validate="many_to_one")

    wanted = set(cells.cell_label)
    frames = []
    for frame in pd.read_csv(paths["annotation"], chunksize=chunksize, dtype={"cell_label": str},
                             keep_default_na=False, low_memory=False):
        frames.append(frame[frame.cell_label.isin(wanted)])
    taxonomy = pd.concat(frames, ignore_index=True)
    overlap = [c for c in taxonomy.columns if c in cells.columns and c != "cell_label"]
    cells = cells.merge(taxonomy.drop(columns=overlap), on="cell_label", how="inner", validate="one_to_one")

    cells["animal"] = animal
    cells["region"] = cells.parcellation_structure
    cells["subclass_clean"] = cells["subclass"].map(clean_subclass)
    cells["cell_class"] = [cell_class(a, b) for a, b in zip(cells["class"], cells["subclass"])]
    cells["qc_annotated"] = cells["class"].ne("") & cells["subclass"].ne("")
    return cells[cells.qc_annotated].copy()


def load_expression(manifest, animal, cell_labels, genes=None):
    """Return a cells-by-genes frame of the official log2 expression for the given cells."""
    import anndata

    spec = resolve(manifest, animal, "expression_matrices", animal, kind="log2")
    path = download_verified(spec)
    adata = anndata.read_h5ad(path, backed="r")
    symbols = (adata.var["gene_symbol"].astype(str).to_numpy()
               if "gene_symbol" in adata.var else adata.var_names.to_numpy().astype(str))
    if genes is not None:
        keep = np.flatnonzero(np.isin(symbols, list(genes)))
    else:
        keep = np.arange(len(symbols))
    index = pd.Index(adata.obs_names.astype(str))
    present = index.intersection(pd.Index(list(cell_labels)))
    rows = index.get_indexer(present)
    order = np.argsort(rows)
    block = adata[rows[order], keep].to_memory().X
    block = np.asarray(block.todense() if hasattr(block, "todense") else block, dtype=np.float32)
    frame = pd.DataFrame(block, index=present[order], columns=symbols[keep])
    adata.file.close()
    return frame, symbols


def within_subclass_difference(cells, expression, min_cells=MIN_CELLS_PER_STRATUM):
    """Mean MOp minus VISp inside each subclass present in BOTH regions.

    Restricting to shared subclasses is what separates an expression difference from a
    composition difference. A subclass that exists in only one region is dropped and
    counted, never imputed.
    """
    meta = cells.set_index("cell_label").loc[expression.index, ["animal", "region", "subclass_clean", "cell_class"]]
    grouped = expression.groupby([meta.animal, meta.subclass_clean, meta.cell_class, meta.region], observed=True)
    means = grouped.mean()
    sizes = grouped.size().rename("cells")
    rows = []
    dropped = []
    for (animal, subclass, klass), frame in means.groupby(level=[0, 1, 2], observed=True):
        regions = frame.index.get_level_values("region")
        counts = {r: int(sizes.loc[(animal, subclass, klass, r)]) for r in regions}
        if not {"MOp", "VISp"}.issubset(counts) or min(counts.values()) < min_cells:
            dropped.append({"animal": animal, "subclass_clean": subclass, "cell_class": klass,
                            "n_MOp": counts.get("MOp", 0), "n_VISp": counts.get("VISp", 0),
                            "reason": "subclass absent in one region" if not {"MOp", "VISp"}.issubset(counts)
                            else f"fewer than {min_cells} cells"})
            continue
        mop = frame.xs("MOp", level="region").iloc[0]
        visp = frame.xs("VISp", level="region").iloc[0]
        for gene in expression.columns:
            rows.append({"animal": animal, "subclass_clean": subclass, "cell_class": klass, "gene": gene,
                         "n_MOp": counts["MOp"], "n_VISp": counts["VISp"],
                         "mean_MOp": mop[gene], "mean_VISp": visp[gene],
                         "MOp_minus_VISp": mop[gene] - visp[gene]})
    table = pd.DataFrame(rows)
    table["units"] = "difference of mean official log2 MERFISH expression; not log2 fold change"
    return table, pd.DataFrame(dropped)


def discordance_score(difference):
    """Per gene and animal: how cleanly GABAergic and glutamatergic subclasses disagree.

    For each gene the sign of the regional difference is counted separately in the two
    cell classes. `discordance` is 1 when every inhibitory subclass goes one way and every
    excitatory subclass the other, and 0 when both classes point the same way. Because a
    region-level technical offset shifts both classes together, a high score cannot be
    produced by that offset alone.
    """
    rows = []
    neuronal = difference[difference.cell_class.isin(["GABAergic", "Glutamatergic"])]
    for (animal, gene), frame in neuronal.groupby(["animal", "gene"], observed=True):
        gaba = frame[frame.cell_class.eq("GABAergic")].MOp_minus_VISp
        glut = frame[frame.cell_class.eq("Glutamatergic")].MOp_minus_VISp
        if len(gaba) < 2 or len(glut) < 2:
            continue
        gaba_sign = np.sign(gaba).mean()
        glut_sign = np.sign(glut).mean()
        rows.append({
            "animal": animal, "gene": gene,
            "n_gaba_subclasses": len(gaba), "n_glut_subclasses": len(glut),
            "gaba_mean_difference": gaba.mean(), "glut_mean_difference": glut.mean(),
            "gaba_sign_consistency": gaba_sign, "glut_sign_consistency": glut_sign,
            "discordance": abs(gaba_sign - glut_sign) / 2.0,
            "separation": abs(gaba.mean() - glut.mean()),
        })
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["discordance_rank"] = table.groupby("animal").discordance.rank(ascending=False, method="min")
    table["separation_rank"] = table.groupby("animal").separation.rank(ascending=False, method="min")
    table["panel_genes_ranked"] = table.groupby("animal").gene.transform("size")
    table["discordance_percentile"] = 100 * (1 - (table.discordance_rank - 1) / table.panel_genes_ranked)
    table["inference"] = ("Rank within the whole measured panel. Reported so a candidate gene is read "
                          "against a panel-wide null rather than on its own.")
    return table


def cross_animal_agreement(difference):
    """Sign agreement across animals for each gene and subclass. Four animals, counted."""
    rows = []
    for (subclass, klass, gene), frame in difference.groupby(
            ["subclass_clean", "cell_class", "gene"], observed=True):
        signs = np.sign(frame.MOp_minus_VISp)
        rows.append({"subclass_clean": subclass, "cell_class": klass, "gene": gene,
                     "animals": len(frame), "animals_positive": int((signs > 0).sum()),
                     "animals_negative": int((signs < 0).sum()),
                     "all_animals_agree": bool(len(frame) > 1 and abs(signs.sum()) == len(frame)),
                     "mean_difference": frame.MOp_minus_VISp.mean(),
                     "min_difference": frame.MOp_minus_VISp.min(),
                     "max_difference": frame.MOp_minus_VISp.max()})
    return pd.DataFrame(rows)


def plot_results(scores, difference, candidates, figdir):
    import os
    cache = ROOT / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir.mkdir(parents=True, exist_ok=True)

    # Where the candidates sit in the panel-wide null.
    pooled = scores.groupby("gene").agg(discordance=("discordance", "mean"),
                                        separation=("separation", "mean")).reset_index()
    pooled["is_candidate"] = pooled.gene.isin(candidates)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    other = pooled[~pooled.is_candidate]
    ax.scatter(other.separation, other.discordance, s=6, color="lightgrey", label="other panel genes")
    hit = pooled[pooled.is_candidate]
    ax.scatter(hit.separation, hit.discordance, s=26, color="#1f5fa8", label="62 PBM candidates")
    top = pooled.sort_values(["discordance", "separation"], ascending=False).head(12)
    for _, row in top.iterrows():
        ax.annotate(row.gene, (row.separation, row.discordance), fontsize=7,
                    textcoords="offset points", xytext=(4, 3))
    ax.set(xlabel="separation between GABAergic and glutamatergic mean difference (log2 units)",
           ylabel="cell class discordance (1 = classes point opposite ways)",
           title="MOp minus VISp within shared subclasses, all measured panel genes")
    ax.legend(); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(figdir / "zhuang_panel_discordance_null.png", dpi=180); plt.close(fig)

    # Per-subclass detail for the candidates that clear the panel.
    focus = [g for g in ["Adora1", "Adora2a"] if g in set(difference.gene)]
    if focus:
        fig, axes = plt.subplots(1, len(focus), figsize=(5.5 * len(focus), 5), squeeze=False)
        for ax, gene in zip(axes[0], focus):
            sub = difference[difference.gene.eq(gene) & difference.cell_class.ne("Other")]
            order = sub.groupby("subclass_clean").MOp_minus_VISp.mean().sort_values().index
            for i, subclass in enumerate(order):
                frame = sub[sub.subclass_clean.eq(subclass)]
                color = "#d2691e" if frame.cell_class.iloc[0] == "GABAergic" else "#1f5fa8"
                ax.scatter(frame.MOp_minus_VISp, [i] * len(frame), color=color, s=22)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_yticks(range(len(order)), order, fontsize=8)
            ax.set_title(gene); ax.set_xlabel("MOp minus VISp (log2 units), one point per animal")
            ax.grid(axis="x", alpha=0.25)
        fig.suptitle("Orange = GABAergic subclass, blue = glutamatergic subclass")
        fig.tight_layout(); fig.savefig(figdir / "zhuang_candidate_by_subclass.png", dpi=180); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--animals", nargs="+", default=list(DEFAULT_ANIMALS))
    parser.add_argument("--release", default=RELEASE)
    parser.add_argument("--download", action="store_true", help="Fetch verified files and analyse")
    parser.add_argument("--list", action="store_true",
                        help="Print Zhuang manifest entries and the download plan, then stop")
    parser.add_argument("--offline", action="store_true",
                        help="Skip the live manifest and plan from the stored snapshot")
    parser.add_argument("--all-panel-genes", action="store_true", default=True,
                        help="Compute the panel-wide null, not just the 62 candidates")
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_STRATUM)
    args = parser.parse_args(argv)

    RESULTS.mkdir(parents=True, exist_ok=True)
    report = {"run_time_UTC": utc_now(), "status": "pending", "caution": CAUTION,
              "animals_requested": args.animals, "release": args.release, "pending": []}
    try:
        manifest, manifest_path, origin = fetch_manifest(args.release, live=not args.offline)
        report["manifest_version"] = manifest["version"]
        report["manifest_origin"] = origin
        report["manifest_sha256"] = checksum(manifest_path, "sha256")
        report["zhuang_manifest_entries"] = describe_manifest(manifest)
        report["download_plan"] = download_plan(manifest, args.animals)
        candidates = candidate_genes()
        report["candidate_genes"] = len(candidates)
        if args.list:
            report["status"] = "manifest_listed"
            print(json.dumps(report["download_plan"], indent=2))
            return 0

        if not args.download:
            report["pending"].append("Pass --download to fetch the verified Allen files")
            return 2
        if origin != "live S3 release manifest":
            raise ValueError("Refusing to download against a stored snapshot. The published MD5 "
                             "must come from the live release manifest; retry when S3 is reachable.")
        required = report["download_plan"]["total_bytes_if_nothing_cached"] + 8 * 1024 ** 3
        INPUT.mkdir(parents=True, exist_ok=True)
        free = shutil.disk_usage(INPUT).free
        report["free_bytes"] = free
        if free < required:
            raise OSError(f"Need about {required:,} free bytes including working room, have {free:,}")

        per_animal = []
        coverage_written = False
        for animal in args.animals:
            cells = select_cells(manifest, animal)
            genes = None if args.all_panel_genes else candidates
            expression, symbols = load_expression(manifest, animal, cells.cell_label, genes)
            if not coverage_written:
                measured = set(symbols)
                pd.DataFrame({"gene": candidates,
                              "measured_in_zhuang_panel": [g in measured for g in candidates],
                              "not_measured_is_not_zero": True}).to_csv(
                    RESULTS / "candidate_panel_measured_coverage.csv", index=False)
                report["panel_size"] = len(symbols)
                report["candidates_measured"] = sorted(g for g in candidates if g in measured)
                report["candidates_not_measured"] = sorted(g for g in candidates if g not in measured)
                coverage_written = True
            difference, dropped = within_subclass_difference(cells, expression, args.min_cells)
            dropped.to_csv(RESULTS / f"{animal}_dropped_strata.csv", index=False)
            per_animal.append(difference)
            report.setdefault("cells", {})[animal] = {
                "MOp": int((cells.region == "MOp").sum()), "VISp": int((cells.region == "VISp").sum()),
                "shared_strata": int(difference[["subclass_clean"]].drop_duplicates().shape[0]),
                "dropped_strata": int(len(dropped))}

        difference = pd.concat(per_animal, ignore_index=True)
        difference.to_csv(RESULTS / "within_subclass_region_difference.csv.gz", index=False)
        scores = discordance_score(difference)
        scores.to_csv(RESULTS / "panel_discordance_scores.csv.gz", index=False)
        agreement = cross_animal_agreement(difference)
        agreement.to_csv(RESULTS / "cross_animal_sign_agreement.csv.gz", index=False)
        scores[scores.gene.isin(candidates)].to_csv(RESULTS / "candidate62_discordance_ranks.csv", index=False)
        plot_results(scores, difference, set(candidates), FIGURES)
        report["status"] = "zhuang_within_subclass_analysis_completed"
        return 0
    except (OSError, ValueError, KeyError, ImportError, urllib.error.URLError) as error:
        report["pending"].append(f"{type(error).__name__}: {error}")
        print(report["pending"][-1], file=sys.stderr)
        return 2
    finally:
        (RESULTS / "zhuang_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, default=str))
        print(json.dumps({"status": report["status"], "pending": report["pending"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
