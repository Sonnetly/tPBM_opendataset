#!/usr/bin/env python3
"""Supertype-level composition check for the Zhuang-ABCA MOp minus VISp comparison.

Why
---
run_zhuang.py compares MOp and VISp inside each WMB subclass. A subclass is still a
mixture of supertypes, and MOp and VISp can hold different supertype proportions (for
example the L4/5 IT subclass, since MOp has little classical L4). A within-subclass
difference can therefore be partly a within-subclass composition difference.

What this computes, per animal, subclass and gene
-------------------------------------------------
- Keep supertypes present in both regions with at least --min-cells cells in each.
- raw_restricted: plain MOp mean minus VISp mean over the cells of the kept supertypes.
- standardized: weighted mean over kept supertypes of (MOp mean minus VISp mean), with
  weights equal to the pooled supertype cell count. This removes the supertype
  composition difference between the regions.
- composition_part = raw_restricted minus standardized.
- coverage = share of the subclass cells (both regions) that fall in kept supertypes.

It runs on every panel gene so candidate genes can be read against the panel-wide
distribution of composition_part. Descriptive only, no PBM condition is present.

Usage
-----
    python src/merfish/zhuang_supertype_check.py --animals Zhuang-ABCA-1 Zhuang-ABCA-2
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_zhuang as rz  # noqa: E402

OUT = rz.RESULTS / "supertype"
MIN_SUPERTYPE_CELLS = 20


def supertype_decomposition(cells, expression, min_cells=MIN_SUPERTYPE_CELLS,
                            min_subclass_cells=rz.MIN_CELLS_PER_STRATUM):
    """Return per subclass x gene raw, restricted, standardized and composition parts."""
    meta = cells.set_index("cell_label").loc[
        expression.index, ["animal", "region", "subclass_clean", "cell_class", "supertype"]]
    rows = []
    shares = []
    for (animal, subclass, klass), idx in meta.groupby(
            ["animal", "subclass_clean", "cell_class"], observed=True).groups.items():
        sub = meta.loc[idx]
        n_region = sub.region.value_counts()
        if not {"MOp", "VISp"}.issubset(n_region.index) or n_region.min() < min_subclass_cells:
            continue
        counts = sub.groupby(["supertype", "region"]).size().unstack(fill_value=0)
        for region in ("MOp", "VISp"):
            if region not in counts:
                counts[region] = 0
        for supertype, row in counts.iterrows():
            shares.append({"animal": animal, "subclass_clean": subclass, "cell_class": klass,
                           "supertype": supertype, "n_MOp": int(row.MOp), "n_VISp": int(row.VISp),
                           "share_MOp": row.MOp / n_region["MOp"],
                           "share_VISp": row.VISp / n_region["VISp"]})
        kept = counts[(counts.MOp >= min_cells) & (counts.VISp >= min_cells)]
        if kept.empty:
            continue
        expr = expression.loc[idx]
        raw = expr[sub.region.eq("MOp").to_numpy()].mean() - expr[sub.region.eq("VISp").to_numpy()].mean()
        in_kept = sub.supertype.isin(kept.index).to_numpy()
        ksub, kexpr = sub[in_kept], expr[in_kept]
        raw_restricted = (kexpr[ksub.region.eq("MOp").to_numpy()].mean()
                          - kexpr[ksub.region.eq("VISp").to_numpy()].mean())
        means = kexpr.groupby([ksub.supertype.to_numpy(), ksub.region.to_numpy()]).mean()
        diffs = means.xs("MOp", level=1) - means.xs("VISp", level=1)
        weights = (kept.MOp + kept.VISp).reindex(diffs.index).astype(float)
        standardized = diffs.mul(weights, axis=0).sum() / weights.sum()
        coverage = in_kept.mean()
        frame = pd.DataFrame({"raw_subclass": raw, "raw_restricted": raw_restricted,
                              "standardized": standardized})
        frame["composition_part"] = frame.raw_restricted - frame.standardized
        frame.index.name = "gene"
        frame = frame.reset_index()
        frame.insert(0, "animal", animal)
        frame.insert(1, "subclass_clean", subclass)
        frame.insert(2, "cell_class", klass)
        frame["supertypes_kept"] = len(kept)
        frame["supertypes_total"] = len(counts)
        frame["coverage"] = coverage
        rows.append(frame)
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()), pd.DataFrame(shares)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--animals", nargs="+", default=list(rz.DEFAULT_ANIMALS))
    parser.add_argument("--release", default=rz.RELEASE)
    parser.add_argument("--min-cells", type=int, default=MIN_SUPERTYPE_CELLS)
    args = parser.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    manifest, _, origin = rz.fetch_manifest(args.release, live=True)
    candidates = rz.candidate_genes()
    parts, share_tables, report = [], [], {
        "run_time_UTC": rz.utc_now(), "animals": args.animals, "release": args.release,
        "manifest_origin": origin, "min_cells_per_supertype_per_region": args.min_cells,
        "caution": rz.CAUTION + " Supertype standardization removes within-subclass composition "
                   "differences only for supertypes sampled in both regions; coverage is reported."}
    for animal in args.animals:
        cells = rz.select_cells(manifest, animal)
        if "supertype" not in cells:
            raise KeyError(f"{animal}: no supertype column in the cluster annotation")
        expression, _ = rz.load_expression(manifest, animal, cells.cell_label, None)
        decomposition, shares = supertype_decomposition(cells, expression, args.min_cells)
        parts.append(decomposition)
        share_tables.append(shares)
        report.setdefault("cells", {})[animal] = {
            "MOp": int((cells.region == "MOp").sum()), "VISp": int((cells.region == "VISp").sum()),
            "subclasses_decomposed": int(decomposition.subclass_clean.nunique()) if len(decomposition) else 0}
        del expression
        print(f"{animal}: decomposed {report['cells'][animal]['subclasses_decomposed']} subclasses", flush=True)

    allparts = pd.concat(parts, ignore_index=True)
    allparts.to_csv(OUT / "supertype_standardized_difference.csv.gz", index=False)
    pd.concat(share_tables, ignore_index=True).to_csv(OUT / "supertype_shares_by_region.csv.gz", index=False)

    neuronal = allparts[allparts.cell_class.ne("Other")].copy()
    neuronal["abs_composition"] = neuronal.composition_part.abs()
    neuronal["abs_raw"] = neuronal.raw_restricted.abs()
    # Panel-wide reference for how large composition parts usually are, per subclass.
    reference = (neuronal.groupby(["animal", "subclass_clean"]).abs_composition
                 .quantile([0.5, 0.95]).unstack().rename(columns={0.5: "panel_median_abs_composition",
                                                                  0.95: "panel_p95_abs_composition"})
                 .reset_index())
    focus = neuronal[neuronal.gene.isin(candidates)].merge(reference, on=["animal", "subclass_clean"])
    focus["composition_above_panel_p95"] = focus.abs_composition > focus.panel_p95_abs_composition
    focus["sign_kept_after_standardization"] = (
        np.sign(focus.raw_restricted) == np.sign(focus.standardized)) | focus.raw_restricted.abs().lt(0.02)
    focus = focus.drop(columns=["abs_composition", "abs_raw"])
    focus.to_csv(OUT / "candidate_supertype_decomposition.csv", index=False)

    agreement = (neuronal.pivot_table(index=["subclass_clean", "gene"], columns="animal",
                                      values="standardized"))
    report["standardized_difference_between_animal_pearson"] = (
        agreement.corr().round(3).to_dict() if agreement.shape[1] > 1 else None)
    report["panel_median_abs_composition_part_neuronal"] = float(neuronal.abs_composition.median())
    report["panel_median_abs_raw_restricted_neuronal"] = float(neuronal.abs_raw.median())
    report["status"] = "zhuang_supertype_check_completed"
    (OUT / "supertype_check_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps({"status": report["status"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
