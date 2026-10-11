#!/usr/bin/env python3
"""Download a public GSE185862 archive file with TLS and size verification."""
import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/"
DEFAULT_ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", choices=["GSE185862_expression_matrix_10x.hdf5", "GSE185862_expression_matrix_SSv4.hdf5", "GSE185862_metadata_10x.csv.gz", "GSE185862_metadata_ssv4.csv.gz", "GSE185862_sample_group_mapping.csv.gz", "GSE185862_CTX_HPF_2021_10x_readme.txt", "GSE185862_CTX_HPF_2021_SS_readme.txt", "GSE185862_CTX-HPF-73538sample_metadata.xlsx"])
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.project_root.resolve()
    path = root / "input/rnaseq/raw" / args.file
    path.parent.mkdir(parents=True, exist_ok=True)
    url = BASE + args.file
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as response:
        expected = int(response.headers["Content-Length"])
    if path.exists():
        if path.stat().st_size != expected:
            raise RuntimeError("Existing final file has wrong length; preserve and investigate it before replacement")
        digest = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
        print(json.dumps({"status": "existing_complete_size_verified", "bytes": expected, "sha256": digest, "file": str(path)}), flush=True)
        return
    temp = path.with_name(path.name + ".part")
    if temp.exists():
        raise RuntimeError("A partial file exists; preserve and investigate/resume it explicitly")
    digest = hashlib.sha256()
    size = 0
    start = time.monotonic()
    with urllib.request.urlopen(url, timeout=60) as response, temp.open("wb") as target:
        headers = dict(response.headers)
        while block := response.read(8 * 1024 * 1024):
            target.write(block)
            digest.update(block)
            size += len(block)
            if size // (256 * 1024 * 1024) != (size - len(block)) // (256 * 1024 * 1024):
                print(json.dumps({"file": args.file, "downloaded_bytes": size, "expected_bytes": expected, "seconds": round(time.monotonic() - start, 1)}), flush=True)
    if size != expected:
        raise RuntimeError(f"Truncated download: {size} bytes, expected {expected}; partial file retained")
    temp.rename(path)
    record = {"status": "download_complete", "source_url": url, "local_path": str(path.relative_to(root)), "bytes": size, "sha256": digest.hexdigest(), "response_headers": headers, "TLS_verification": "preserved", "source": "official GEO GSE185862 archive", "publisher_checksum_available": False}
    manifest = root / "input/rnaseq" / (args.file + ".download.json")
    manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
