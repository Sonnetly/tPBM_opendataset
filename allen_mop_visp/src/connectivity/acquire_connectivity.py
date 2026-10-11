#!/usr/bin/env python3
"""Acquire public Allen projection metadata with HTTPS and raw response provenance."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "input/connectivity"
RESULTS = ROOT / "results/connectivity"
BASE = "https://api.brain-map.org/api/v2/data/query.json"
SOURCE_IDS = {"MOp": 985, "VISp": 385}
# Exact ontology targets, including separately labeled BF regional-context controls.
TARGET_IDS = [689, 452, 797, 342, 262, 298, 596]


def now():
    return datetime.now(timezone.utc).isoformat()


def make_url(criteria):
    return BASE + "?" + urllib.parse.urlencode({"criteria": criteria})


def read_response(path):
    data = json.loads(path.read_text())
    if not data.get("success") or not isinstance(data.get("msg"), list):
        raise ValueError(f"Unsuccessful Allen RMA response: {path.name}")
    if data.get("num_rows") != len(data["msg"]) or data.get("total_rows") != len(data["msg"]):
        raise ValueError(f"Incomplete Allen RMA response: {path.name}")
    return data


def fetch(criteria, filename, download=False):
    path = INPUT / filename
    provenance_path = path.with_suffix(".provenance.json")
    url = make_url(criteria)
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"Missing official input {path}; use --download")
        request = urllib.request.Request(url, headers={"User-Agent": "tPBM-opendataset-connectivity-analysis/1.0"})
        with urllib.request.urlopen(request, timeout=90) as response:
            if not response.url.startswith("https://"):
                raise ValueError("Unexpected non-HTTPS redirect")
            content = response.read()
            observation = {"url": url, "final_url": response.url, "status_code": response.status,
                           "retrieved_at_UTC": now(), "content_type": response.headers.get("Content-Type"),
                           "response_date": response.headers.get("Date")}
        temporary = path.with_suffix(path.suffix + ".partial")
        temporary.write_bytes(content)
        read_response(temporary)
        temporary.replace(path)
    else:
        if provenance_path.exists():
            observation = json.loads(provenance_path.read_text())
            current_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            if current_hash != observation.get("sha256"):
                raise ValueError(f"Cached official response SHA256 changed; preserve and investigate: {path.name}")
            if observation.get("url") != url:
                raise ValueError(f"Cached response criteria do not match requested official URL: {path.name}")
            observation["last_loaded_existing_response_at_UTC"] = now()
            observation["cached_sha256_matches_initial_record"] = True
        else:
            observation = {"url": url, "loaded_existing_raw_response_at_UTC": now(),
                           "initial_acquisition": "Unknown initial retrieval: existing response had no provenance record; current bytes preserved"}
    data = read_response(path)
    content = path.read_bytes()
    observation.update({"path": str(path.relative_to(ROOT)), "response_bytes": len(content),
                        "sha256": hashlib.sha256(content).hexdigest(), "rows": len(data["msg"]),
                        "total_rows": data["total_rows"], "TLS_verification": "preserved",
                        "publisher_checksum": "not supplied for dynamic API response; local SHA256 records exact response"})
    provenance_path.write_text(json.dumps(observation, indent=2))
    return data["msg"], observation


def source_region(experiment, structures):
    path = structures.get(experiment["structure_id"], {}).get("structure_id_path", "")
    hits = [name for name, sid in SOURCE_IDS.items() if f"/{sid}/" in path]
    if len(hits) > 1:
        raise ValueError("Ambiguous source ontology membership")
    return hits[0] if hits else None


def acquire(download=False):
    INPUT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    structure_criteria = "model::Structure,rma::criteria,[graph_id$eq1],rma::options[num_rows$eqall]"
    experiment_criteria = "model::ApiConnectivity,rma::options[num_rows$eqall]"
    structure_rows, structure_provenance = fetch(structure_criteria, "ontology_structure_graph1.json", download)
    experiments, experiment_provenance = fetch(experiment_criteria, "api_connectivity_experiments.json", download)
    structures = {s["id"]: s for s in structure_rows}
    # MouseConnectivityProjection product 5 only. Product 31 is retinal projection.
    selected = [e for e in experiments if e["product_id"] == 5 and source_region(e, structures)]
    ids = sorted(e["data_set_id"] for e in selected)
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError("Empty or duplicate selected experiment identifiers")
    requests = []
    for number, start in enumerate(range(0, len(ids), 80)):
        chunk = ids[start:start+80]
        id_filter = "[section_data_set_id$in" + ",".join(map(str, chunk)) + "]"
        options = ",rma::options[num_rows$eqall]"
        requests.append((f"injection_unionizes_{number:02d}.json", "model::ProjectionStructureUnionize,rma::criteria," + id_filter + "[is_injection$eqtrue]" + options))
        target_filter = "[structure_id$in" + ",".join(map(str, TARGET_IDS)) + "]"
        requests.append((f"target_unionizes_{number:02d}.json", "model::ProjectionStructureUnionize,rma::criteria," + id_filter + "[is_injection$eqfalse]" + target_filter + options))
        details_filter = "[id$in" + ",".join(map(str, chunk)) + "]"
        include = ",rma::include,specimen(donor(transgenic_mouse(transgenic_lines)),injections(structure,age)),equalization,products"
        requests.append((f"manual_metadata_{number:02d}.json", "model::SectionDataSet,rma::criteria," + details_filter + include + options))
    # Independent read-only RMA requests; different raw response files.
    def perform(request):
        filename, criteria = request
        print(f"Collect {filename}", flush=True)
        rows, provenance = fetch(criteria, filename, download)
        return {"filename": filename, "rows": rows, "provenance": provenance}
    with ThreadPoolExecutor(max_workers=3) as pool:
        collected = list(pool.map(perform, requests))
    provenance = {"completed_at_UTC": now(), "status": "download_complete", "ontology_graph": 1,
                  "source_selection": "Product 5; primary injection structure is MOp/VISp or a graph-1 descendant",
                  "source_region_ids": SOURCE_IDS, "requested_target_ids": TARGET_IDS,
                  "api_experiments_returned": len(experiments), "selected_primary_source_experiments": len(ids),
                  "independent_unit": "tracer experiment; donor identifiers retained for duplicate-animal audit",
                  "base_responses": [structure_provenance, experiment_provenance],
                  "batch_responses": [c["provenance"] for c in collected]}
    (INPUT / "acquisition_manifest.json").write_text(json.dumps(provenance, indent=2))
    return selected, structures, collected


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    selected, structures, collected = acquire(args.download)
    print(json.dumps({"selected_experiments": len(selected), "response_batches": len(collected)}))
