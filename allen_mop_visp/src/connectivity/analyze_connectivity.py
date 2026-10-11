#!/usr/bin/env python3
"""Describe anatomically labeled axons; do not infer synapses or sleep effects."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from acquire_connectivity import ROOT, INPUT, RESULTS, SOURCE_IDS, source_region, read_response

FIGURES = ROOT / "figures/connectivity"
HEMISPHERES = {1: "left", 2: "right", 3: "bilateral_aggregate"}
METRICS = ["normalized_projection_volume", "projection_volume", "projection_density", "projection_energy", "projection_intensity", "volume", "max_voxel_density", "sum_pixels", "sum_projection_pixels"]
CAUTION = "Anterograde viral tracer segmented axonal signal: not proof of synapses, postsynaptic identity, monosynaptic functional coupling, or NREM promotion. No physiological response is measured."


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False))


def load_inputs():
    structures = {r["id"]: r for r in read_response(INPUT / "ontology_structure_graph1.json")["msg"]}
    experiments = read_response(INPUT / "api_connectivity_experiments.json")["msg"]
    manifest = json.loads((INPUT / "acquisition_manifest.json").read_text())
    for response in manifest["base_responses"]:
        path = ROOT / response["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != response["sha256"]:
            raise ValueError(f"Base official API response changed: {path.name}")
    batches = {"injection": [], "target": [], "manual": []}
    for response in manifest["batch_responses"]:
        path = ROOT / response["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != response["sha256"]:
            raise ValueError(f"Raw API response changed: {path.name}")
        if path.name.startswith("injection_unionizes_"): category = "injection"
        elif path.name.startswith("target_unionizes_"): category = "target"
        elif path.name.startswith("manual_metadata_"): category = "manual"
        else: raise ValueError(f"Unexpected batch: {path.name}")
        batches[category].extend(read_response(path)["msg"])
    definitions = json.loads((INPUT / "nrem_target_definitions.json").read_text())
    if isinstance(definitions, dict): definitions = definitions["targets"]
    for target in definitions:
        sid = target.get("structure_id")
        if sid is not None and structures[sid]["acronym"] != target["acronym"]:
            raise ValueError(f"Target definition ontology mismatch: {target}")
    return structures, experiments, batches, definitions


def metadata_table(structures, experiments, batches):
    manual = {m["id"]: m for m in batches["manual"]}
    if len(manual) != len(batches["manual"]): raise ValueError("Duplicate manual metadata")
    injection = pd.DataFrame(batches["injection"])
    if injection.duplicated(["section_data_set_id", "structure_id", "hemisphere_id", "is_injection"]).any():
        raise ValueError("Duplicate injection unionize keys")
    lookup = injection.set_index(["section_data_set_id", "structure_id", "hemisphere_id"])
    rows, selection_audit = [], []
    for raw in experiments:
        source = source_region(raw, structures)
        selected = raw["product_id"] == 5 and source is not None
        declared = [int(x) for x in raw.get("injection_structures", "").split("/") if x]
        selection_audit.append({"experiment_id": raw["data_set_id"], "product_id": raw["product_id"], "primary_structure_id": raw["structure_id"],
                                "primary_acronym": raw["structure_abbrev"], "source_region_from_primary_ontology": source,
                                "primary_MOp_VISp_product5": selected, "transgenic_line_name": (raw.get("transgenic_line") or {}).get("name"),
                                "strain": raw.get("strain"), "declared_injection_structures": raw.get("injection_structures"),
                                "contains_MOp_VISp_declared_injection": any(source_region({"structure_id": sid}, structures) for sid in declared)})
        if not selected: continue
        eid = raw["data_set_id"]
        if eid not in manual: raise ValueError(f"Missing manual quality metadata {eid}")
        m = manual[eid]
        specimen = m["specimen"]; donor = specimen.get("donor") or {}; injections = specimen.get("injections") or []
        lineage = raw.get("transgenic_line") or {}
        raw_cre = bool(lineage or raw.get("transgenic_line_id"))
        donor_cre = donor.get("transgenic_mouse_id") is not None
        source_ids = {sid for sid, s in structures.items() if f"/{SOURCE_IDS[source]}/" in s["structure_id_path"]}
        outside = [sid for sid in declared if sid not in source_ids]
        values = {}
        for sid, label in [(997, "whole_brain"), (SOURCE_IDS[source], "source_parent")]:
            for hemi in [1, 2, 3]:
                key = (eid, sid, hemi)
                values[f"{label}_injection_signal_volume_h{hemi}"] = float(lookup.loc[key, "projection_volume"]) if key in lookup.index else np.nan
        total = values["whole_brain_injection_signal_volume_h3"]
        source_total = values["source_parent_injection_signal_volume_h3"]
        metadata_total = raw.get("injection_volume")
        if metadata_total is None or not np.isfinite(metadata_total) or metadata_total <= 0:
            raise ValueError(f"Invalid metadata injection signal volume: {eid}")
        left, right = values["whole_brain_injection_signal_volume_h1"], values["whole_brain_injection_signal_volume_h2"]
        # Missing injection hemisphere records are missing, rather than silently imputed zero.
        if np.isfinite(left) and np.isfinite(right):
            dominant = 1 if left > right else (2 if right > left else None)
        elif np.isfinite(left) and left > 0: dominant = 1
        elif np.isfinite(right) and right > 0: dominant = 2
        else: dominant = None
        fraction = source_total / total if np.isfinite(total) and total > 0 and np.isfinite(source_total) else np.nan
        row = {"experiment_id": eid, "source_region": source, "primary_structure_id": raw["structure_id"], "primary_structure_name": raw["structure_name"],
               "product_id": raw["product_id"], "product_name": "MouseConnectivityProjection", "specimen_id": m["specimen_id"], "specimen_name": raw["specimen_name"],
               "donor_id": donor.get("id"), "donor_strain": donor.get("strain"), "api_strain": raw.get("strain"), "donor_sex": donor.get("sex"),
               "transgenic_line_name": lineage.get("name"), "transgenic_line_id": raw.get("transgenic_line_id"),
               "transgenic_mouse_present": donor_cre, "is_Cre_line": raw_cre,
               "cohort": "Cre_line" if raw_cre or donor_cre else ("WT_C57BL6J" if raw.get("strain") == "C57BL/6J" and donor.get("strain") == "C57BL/6J" else "non_Cre_other_or_unknown_strain"),
               "failed": m["failed"], "atlas_quality_pass": m["failed"] is False, "qc_date": m.get("qc_date"), "reference_space_id": m.get("reference_space_id"),
               "injection_x_um": raw.get("injection_x"), "injection_y_um": raw.get("injection_y"), "injection_z_um": raw.get("injection_z"),
               "injection_signal_volume_metadata_mm3": metadata_total,
               "injection_structures_raw": raw.get("injection_structures"), "off_source_declared_structure_ids": "/".join(map(str, outside)),
               "off_source_declared_acronyms": "/".join(structures.get(sid, {}).get("acronym", str(sid)) for sid in outside),
               "source_only_annotated_injection_structures": len(outside) == 0 and len(declared) > 0,
               "injection_source_parent_fraction": fraction, "off_source_fraction": 1 - fraction if np.isfinite(fraction) else np.nan,
               "injection_dominant_hemisphere_id": dominant, "injection_dominant_hemisphere": HEMISPHERES.get(dominant, "undetermined"),
               "manual_injection_records": len(injections), "manual_injection_quality_codes": json.dumps(sorted({str(i.get('injection_quality')) for i in injections})),
               "manual_injection_age_days": json.dumps(sorted({i.get('age', {}).get('days') for i in injections if i.get('age', {}).get('days') is not None})),
               "injection_age_P56_only": bool(injections and all(i.get('age', {}).get('days') == 56 for i in injections)),
               "manual_days_post_injection": json.dumps(sorted({i.get('days_post_injection') for i in injections if i.get('days_post_injection') is not None})),
               "manual_injection_methods": json.dumps(sorted({str(i.get('injection_method')) for i in injections})),
               "manual_injection_materials": json.dumps(sorted({str(i.get('injection_materials')) for i in injections})),
               "manual_registration_points": json.dumps(sorted({str(i.get('registration_point')) for i in injections})),
               "injection_volume_root_matches_metadata": bool(np.isfinite(total) and np.isclose(total, metadata_total, rtol=2e-6, atol=1e-12)),
               "experiment_viewer_url": f"https://connectivity.brain-map.org/projection/experiment/{eid}", **values}
        rows.append(row)
    metadata = pd.DataFrame(rows)
    pd.DataFrame(selection_audit).to_csv(RESULTS / "all_experiments_selection_audit.csv", index=False)
    injection.to_csv(RESULTS / "injection_unionizes_all_structures.csv.gz", index=False)
    metadata.to_csv(RESULTS / "selected_experiment_metadata_and_injection_audit.csv", index=False)
    return metadata


def target_table(metadata, raw, definitions, structures):
    unionizes = pd.DataFrame(raw)
    if unionizes.duplicated(["section_data_set_id", "structure_id", "hemisphere_id", "is_injection"]).any(): raise ValueError("Duplicate target unionize keys")
    if unionizes.is_injection.any(): raise ValueError("Injection record in target projection query")
    unionizes.to_csv(RESULTS / "target_unionizes_raw_metrics.csv.gz", index=False)
    lookup = unionizes.set_index(["section_data_set_id", "structure_id", "hemisphere_id"])
    columns = ["experiment_id", "source_region", "cohort", "transgenic_line_name", "donor_id", "atlas_quality_pass", "source_only_annotated_injection_structures", "injection_source_parent_fraction", "injection_age_P56_only", "injection_signal_volume_metadata_mm3", "injection_dominant_hemisphere_id"]
    rows = []
    for experiment in metadata.to_dict("records"):
        for target in definitions:
            sid = target.get("structure_id")
            for hemi in [1, 2, 3]:
                key = (experiment["experiment_id"], sid, hemi)
                observed = sid is not None and key in lookup.index
                row = {k: experiment[k] for k in columns}
                row.update({"target_key": target["target_key"], "target_acronym": target["acronym"], "target_structure_id": sid,
                            "target_atlas_name": target.get("atlas_name"), "target_tier": target["tier"], "target_interpretation": target["interpretation"],
                            "mapping_status": target["mapping_status"], "postsynaptic_celltype_resolved": False,
                            "hemisphere_id": hemi, "hemisphere": HEMISPHERES[hemi],
                            "hemisphere_relative_to_dominant_injection": "both" if hemi == 3 else ("ipsilateral" if hemi == experiment["injection_dominant_hemisphere_id"] else ("contralateral" if pd.notna(experiment["injection_dominant_hemisphere_id"]) else "undetermined")),
                            "record_observed": observed, "measurement_status": "unmapped_anatomical_target" if sid is None else ("observed" if observed else "missing_API_record")})
                for metric in METRICS: row[metric] = float(lookup.loc[key, metric]) if observed else np.nan
                row["segmented_signal_positive"] = bool(row["projection_volume"] > 0) if observed else None
                row["observed_zero_signal"] = bool(row["projection_volume"] == 0) if observed else None
                row["computed_normalized_projection_volume"] = row["projection_volume"] / experiment["injection_signal_volume_metadata_mm3"] if observed else np.nan
                row["normalization_arithmetic_matches"] = bool(np.isclose(row["normalized_projection_volume"], row["computed_normalized_projection_volume"], rtol=2e-6, atol=1e-14)) if observed else None
                rows.append(row)
    output = pd.DataFrame(rows)
    output.to_csv(RESULTS / "projection_by_experiment_target_hemisphere.csv", index=False)
    overlaps = []
    mapped = [d for d in definitions if d.get("structure_id") is not None]
    for a in mapped:
        for b in mapped:
            if a is not b and f"/{a['structure_id']}/" in structures[b["structure_id"]]["structure_id_path"]:
                overlaps.append({"ancestor": a["acronym"], "descendant": b["acronym"]})
    pd.DataFrame(overlaps, columns=["ancestor", "descendant"]).to_csv(RESULTS / "target_hierarchy_overlap_audit.csv", index=False)
    return output


def summarize(output, metadata):
    bilateral = output[output.hemisphere_id.eq(3) & output.atlas_quality_pass].copy()
    def summary(frame, dimensions):
        return frame.groupby(dimensions, dropna=False).agg(n_experiments=("experiment_id", "size"),
            n_records_observed=("record_observed", "sum"), n_positive_segmented_signal=("segmented_signal_positive", "sum"),
            n_observed_zero=("observed_zero_signal", "sum"), median_NPV=("normalized_projection_volume", "median"),
            q25_NPV=("normalized_projection_volume", lambda s: s.quantile(.25)), q75_NPV=("normalized_projection_volume", lambda s: s.quantile(.75)),
            min_NPV=("normalized_projection_volume", "min"), max_NPV=("normalized_projection_volume", "max"),
            median_projection_volume_mm3=("projection_volume", "median"), median_projection_density=("projection_density", "median"),
            median_injection_source_parent_fraction=("injection_source_parent_fraction", "median")).reset_index()
    dimensions = ["source_region", "target_acronym", "target_tier"]
    wt = bilateral[bilateral.cohort.eq("WT_C57BL6J")]
    all_wt = summary(wt, dimensions); all_wt["injection_scope"] = "all_primary_source_injections"
    pure_wt = summary(wt[wt.source_only_annotated_injection_structures], dimensions); pure_wt["injection_scope"] = "source_only_annotated_injection_structures"
    age_wt = summary(wt[wt.injection_age_P56_only], dimensions); age_wt["injection_scope"] = "exact_P56_age"
    purity_wt = summary(wt[wt.injection_source_parent_fraction >= .8], dimensions); purity_wt["injection_scope"] = "source_parent_signal_fraction_at_least_0.8_exploratory"
    both_wt = summary(wt[wt.injection_age_P56_only & (wt.injection_source_parent_fraction >= .8)], dimensions); both_wt["injection_scope"] = "exact_P56_and_source_parent_signal_fraction_at_least_0.8_exploratory"
    combined = pd.concat([all_wt, pure_wt, age_wt, purity_wt, both_wt], ignore_index=True)
    combined.to_csv(RESULTS / "WT_target_projection_descriptive_summary.csv", index=False)
    summary(bilateral[bilateral.cohort.eq("Cre_line")], dimensions + ["transgenic_line_name"]).to_csv(RESULTS / "Cre_line_separate_exploratory_summary.csv", index=False)
    comparisons = []
    for scope, group in combined.groupby("injection_scope"):
        wide = group.pivot(index=["target_acronym", "target_tier"], columns="source_region", values="median_NPV")
        for (target, tier), row in wide.iterrows():
            mop, visp = row.get("MOp", np.nan), row.get("VISp", np.nan)
            comparisons.append({"target_acronym": target, "target_tier": tier, "injection_scope": scope,
                                "MOp_median_NPV": mop, "VISp_median_NPV": visp,
                                "MOp_minus_VISp_median_NPV": mop - visp,
                                "MOp_to_VISp_median_NPV_ratio": mop / visp if pd.notna(visp) and visp > 0 else np.nan,
                                "animal_or_region_effect_inference": "not estimated; heterogeneous injection sampling, spread and tracer conditions"})
    pd.DataFrame(comparisons).to_csv(RESULTS / "WT_exploratory_regional_descriptive_difference.csv", index=False)
    metadata.groupby(["source_region", "cohort", "source_only_annotated_injection_structures", "atlas_quality_pass"], dropna=False).agg(experiments=("experiment_id", "size"), donors=("donor_id", "nunique"),
                min_injection_volume_mm3=("injection_signal_volume_metadata_mm3", "min"), max_injection_volume_mm3=("injection_signal_volume_metadata_mm3", "max"),
                median_source_parent_fraction=("injection_source_parent_fraction", "median"), min_source_parent_fraction=("injection_source_parent_fraction", "min")).reset_index().to_csv(RESULTS / "source_cohort_coverage_and_spread.csv", index=False)
    return combined


def figures(output, metadata, definitions):
    cache = ROOT / ".cache"; cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib")); os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    FIGURES.mkdir(parents=True, exist_ok=True)
    target_order = [t["acronym"] for t in definitions if t.get("structure_id") is not None]
    wt = output[output.cohort.eq("WT_C57BL6J") & output.hemisphere_id.eq(3) & output.atlas_quality_pass]
    nonzero = wt.loc[wt.normalized_projection_volume > 0, "normalized_projection_volume"]
    vmin, vmax = nonzero.min(), nonzero.max()
    cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("white")
    fig, axes = plt.subplots(1, 2, figsize=(14, 12), constrained_layout=True, gridspec_kw={"width_ratios": [1, 2]})
    for ax, source in zip(axes, ["MOp", "VISp"]):
        subset = wt[wt.source_region.eq(source)]
        order = metadata[metadata.cohort.eq("WT_C57BL6J") & metadata.source_region.eq(source)].sort_values(["source_only_annotated_injection_structures", "injection_source_parent_fraction"], ascending=False)
        values = subset.pivot(index="experiment_id", columns="target_acronym", values="normalized_projection_volume").reindex(index=order.experiment_id, columns=target_order)
        masked = np.ma.masked_where(values.to_numpy() <= 0, values.to_numpy())
        image = ax.imshow(masked, aspect="auto", cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
        for row, eid in enumerate(order.experiment_id):
            for column, target in enumerate(target_order):
                if values.loc[eid, target] == 0: ax.text(column, row, "0", ha="center", va="center", fontsize=6)
        labels = [f"{r.experiment_id}{'*' if r.source_only_annotated_injection_structures else ''}  f={r.injection_source_parent_fraction:.3f}" for r in order.itertuples()]
        ax.set_yticks(range(len(labels)), labels, fontsize=7)
        ax.set_xticks(range(len(target_order)), target_order, rotation=40)
        ax.set_title(f"{source} WT; n={len(order)} tracer experiments")
    fig.colorbar(image, ax=axes, label="Bilateral normalized projection volume (log scale)", shrink=.75)
    fig.suptitle("Axonal segmented signal into atlas regions; functional NREM effects unmeasured\n* source-only annotated injection labels; f = injection signal fraction in source parent; 0 = observed zero")
    fig.savefig(FIGURES / "WT_individual_target_projection_heatmap.png", dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, source in zip(axes, ["MOp", "VISp"]):
        subset = metadata[metadata.cohort.eq("WT_C57BL6J") & metadata.source_region.eq(source)].sort_values("injection_source_parent_fraction", ascending=False)
        colors = np.where(subset.source_only_annotated_injection_structures, "#225ea8", "#fc8d59")
        ax.bar(np.arange(len(subset)), subset.injection_source_parent_fraction, color=colors)
        ax.set_xticks(np.arange(len(subset)), subset.experiment_id, rotation=90, fontsize=7)
        ax.set(title=source, ylabel="Fraction of total injection signal in source parent", ylim=(0,1.02))
    fig.suptitle("Injection spread audit: blue = source-only declared labels; orange = regional spread labels")
    fig.savefig(FIGURES / "WT_injection_source_purity.png", dpi=180); plt.close(fig)


def run():
    RESULTS.mkdir(parents=True, exist_ok=True)
    structures, experiments, batches, definitions = load_inputs()
    metadata = metadata_table(structures, experiments, batches)
    output = target_table(metadata, batches["target"], definitions, structures)
    summary = summarize(output, metadata)
    figures(output, metadata, definitions)
    observed = output[output.record_observed]
    mapped_count = sum(t.get("structure_id") is not None for t in definitions)
    expected = len(metadata) * mapped_count * 3
    if len(observed) != expected: raise ValueError(f"Incomplete target coverage: {len(observed)} of {expected}")
    if not observed.normalization_arithmetic_matches.all(): raise ValueError("Target NPV arithmetic mismatch")
    if not metadata.injection_volume_root_matches_metadata.all(): raise ValueError("Injection volume root mismatch")
    for cohort in ["WT_C57BL6J", "Cre_line"]:
        subset = metadata[metadata.cohort.eq(cohort)]
        if subset.donor_id.duplicated().any():
            subset[subset.donor_id.duplicated(False)].to_csv(RESULTS / f"{cohort}_duplicate_donor_audit.csv", index=False)
    status = {"completed_at_UTC": datetime.now(timezone.utc).isoformat(), "status": "actual_projection_analysis_completed",
              "actual_data_analyzed": True, "selected_experiments": len(metadata), "WT_experiments_by_source": metadata[metadata.cohort.eq("WT_C57BL6J")].source_region.value_counts().to_dict(),
              "Cre_experiments_by_source": metadata[metadata.cohort.eq("Cre_line")].source_region.value_counts().to_dict(),
              "WT_source_only_annotation_by_source": metadata[metadata.cohort.eq("WT_C57BL6J") & metadata.source_only_annotated_injection_structures].source_region.value_counts().to_dict(),
              "failed_flag_counts": {str(k): int(v) for k,v in metadata.failed.value_counts().items()}, "mapped_targets": mapped_count,
              "unmapped_targets": [d["acronym"] for d in definitions if d.get("structure_id") is None], "expected_mapped_target_hemisphere_records": expected,
              "observed_mapped_target_hemisphere_records": len(observed), "missing_mapped_target_hemisphere_records": int((output.measurement_status=="missing_API_record").sum()),
              "normalization_arithmetic_verified": True, "injection_root_total_arithmetic_verified": True, "no_inference_p_values": True,
              "sensitivity_definitions": ["source-only annotated labels", "exact P56 injection age", "source-parent injection signal fraction >= 0.8", "P56 and fraction >= 0.8"],
              "purity_cutoff_interpretation": "0.8 is an exploratory injection-spread sensitivity condition, not a signal detection cutoff, official quality score, or functional-evidence threshold",
              "summary_unit": "tracer experiment; WT and Cre lines separated", "caution": CAUTION,
              "source_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [INPUT / "nrem_target_definitions.json", Path(__file__)]}}
    json_write(RESULTS / "status.json", status)
    print(json.dumps(status, ensure_ascii=False))
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    run()
