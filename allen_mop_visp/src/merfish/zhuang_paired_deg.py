#!/usr/bin/env python3
"""Animal-paired MOp versus VISp differential expression inside each neuronal subclass.

Input is the per-animal within-subclass table written by run_zhuang.py
(results/zhuang/within_subclass_region_difference.csv.gz), and, as a robustness
column, the supertype-standardized table from zhuang_supertype_check.py.

Why paired by animal
--------------------
In Zhuang-ABCA MERFISH both regions sit in the same sections of the same animal, so
the animal is the replicate and region is a within-animal factor. Each animal gives one
pseudobulk difference (mean MOp log2 minus mean VISp log2) per subclass and gene. A
one-sample t test on those differences across animals is the paired test. With 3 or 4
animals the test has little power, so the sign agreement and the effect size are
reported next to the p value and the BH q value, and a gene is called only when all
three criteria agree.

Call rule (prespecified)
------------------------
- animals with the subclass sampled >= 3
- all animals share the sign of the difference
- |mean difference| >= 0.10 log2 units
- BH q < 0.10 within subclass across the panel

Descriptive inference about these four B6N animals. No PBM condition is present.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "zhuang"
OUT = RES / "paired_deg"
MIN_ANIMALS = 3
MIN_EFFECT = 0.10
MAX_Q = 0.10


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * ok.sum() / (np.arange(ok.sum()) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    out = np.empty(ok.sum())
    out[order] = ranked
    q[ok] = out
    return q


def paired_table(frame, value, min_animals=MIN_ANIMALS):
    """frame has animal, subclass_clean, cell_class, gene and a difference column."""
    wide = frame.pivot_table(index=["cell_class", "subclass_clean", "gene"], columns="animal", values=value)
    values = wide.to_numpy()
    n = np.isfinite(values).sum(axis=1)
    mean = np.nanmean(values, axis=1)
    sd = np.nanstd(values, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = mean / (sd / np.sqrt(n))
    p = np.where(n >= 2, 2 * stats.t.sf(np.abs(t), df=np.maximum(n - 1, 1)), np.nan)
    pos = (values > 0).sum(axis=1)
    neg = (values < 0).sum(axis=1)
    out = pd.DataFrame({"n_animals": n, "mean_difference": mean, "sd_difference": sd, "t": t, "p": p,
                        "animals_positive": pos, "animals_negative": neg,
                        "min_difference": np.nanmin(values, axis=1),
                        "max_difference": np.nanmax(values, axis=1)}, index=wide.index).reset_index()
    for animal in wide.columns:
        out[f"diff_{animal}"] = wide[animal].to_numpy()
    out = out[out.n_animals >= min_animals].copy()
    out["q_BH_within_subclass"] = out.groupby("subclass_clean").p.transform(lambda s: bh(s.to_numpy()))
    out["all_same_sign"] = (out.animals_positive == out.n_animals) | (out.animals_negative == out.n_animals)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-effect", type=float, default=MIN_EFFECT)
    parser.add_argument("--max-q", type=float, default=MAX_Q)
    args = parser.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(RES / "within_subclass_region_difference.csv.gz")
    raw = raw[raw.cell_class.ne("Other")]
    table = paired_table(raw, "MOp_minus_VISp")

    std_path = RES / "supertype" / "supertype_standardized_difference.csv.gz"
    if std_path.exists():
        std = pd.read_csv(std_path)
        std = std[std.cell_class.ne("Other")]
        st = paired_table(std, "standardized")[["subclass_clean", "gene", "n_animals", "mean_difference",
                                                 "all_same_sign", "q_BH_within_subclass"]]
        st = st.rename(columns={"n_animals": "std_n_animals", "mean_difference": "std_mean_difference",
                                "all_same_sign": "std_all_same_sign", "q_BH_within_subclass": "std_q"})
        table = table.merge(st, on=["subclass_clean", "gene"], how="left")
    table["called"] = (table.all_same_sign & table.mean_difference.abs().ge(args.min_effect)
                       & table.q_BH_within_subclass.lt(args.max_q))
    if "std_mean_difference" in table:
        table["survives_supertype_standardization"] = (
            table.called & table.std_all_same_sign.fillna(False)
            & (np.sign(table.std_mean_difference) == np.sign(table.mean_difference))
            & table.std_mean_difference.abs().ge(args.min_effect / 2))
    table["direction"] = np.where(table.mean_difference > 0, "MOp_higher", "VISp_higher")
    table = table.sort_values(["subclass_clean", "p"])
    table.to_csv(OUT / "paired_region_deg_all.csv.gz", index=False)
    called = table[table.called]
    called.to_csv(OUT / "paired_region_deg_called.csv", index=False)

    summary = (called.groupby(["cell_class", "subclass_clean", "direction"]).size()
               .unstack(fill_value=0).reset_index())
    if "survives_supertype_standardization" in called:
        robust = called[called.survives_supertype_standardization]
        summary = summary.merge(robust.groupby("subclass_clean").size().rename("survive_supertype").reset_index(),
                                on="subclass_clean", how="left").fillna({"survive_supertype": 0})
    summary.to_csv(OUT / "paired_region_deg_summary.csv", index=False)
    report = {"rule": {"min_animals": MIN_ANIMALS, "all_same_sign": True, "min_abs_mean_difference": args.min_effect,
                       "max_q_BH_within_subclass": args.max_q},
              "units": "difference of mean official log2 MERFISH expression; not log2 fold change",
              "tested_rows": int(len(table)), "called_rows": int(len(called)),
              "called_genes": int(called.gene.nunique()),
              "caution": "Four B6N animals, paired by animal. Panel of 1,122 genes chosen for cell typing, "
                         "so absent genes are unmeasured. No PBM condition is present."}
    (OUT / "paired_region_deg_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
