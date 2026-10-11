#!/usr/bin/env python3
"""Discover exact links on Allen's official CTX-HPF download page.

Use standard TLS verification. Candidate links are a discovery inventory, not
proof that an expression artifact has a particular scale or cell population.
"""
import argparse
import csv
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

DEFAULT_ROOT = Path(__file__).resolve().parents[2]
PAGES = {
    "10x": "https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/mouse-whole-cortex-and-hippocampus-10x",
    "smartseq": "https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/mouse-whole-cortex-and-hippocampus-smart-seq",
}


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.items = []
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and attrs.get("href"):
            self.items.append(attrs["href"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--page", choices=PAGES, default="10x")
    p.add_argument("--url", help="Official Allen URL override after verified portal relocation")
    args = p.parse_args()
    root = args.project_root.resolve()
    out = root / "input/rnaseq"
    logs = root / "logs/rnaseq"
    out.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    url = args.url or PAGES[args.page]
    host = urllib.parse.urlsplit(url).hostname
    if not (host and (host == "brain-map.org" or host.endswith(".brain-map.org"))):
        p.error("Discovery accepts official brain-map.org pages; verify other sources separately")
    try:
        with urllib.request.urlopen(url, timeout=45) as response:
            data = response.read()
            final_url = response.url
        (out / f"official_{args.page}_download_page.html").write_bytes(data)
        parser = Links()
        parser.feed(data.decode("utf-8", errors="replace"))
        # Portal scripts can contain download URLs outside ordinary anchors.
        links = parser.items + re.findall(r'https?://[^\s<>"\\]+', data.decode("utf-8", errors="replace"))
        records = []
        for href in sorted(set(urllib.parse.urljoin(final_url, x) for x in links)):
            if re.search(r"metadata|matrix|\.csv|\.h5|\.loom|\.zip|\.tar|\.gz|download|GSE[0-9]+", href, re.I):
                records.append({"url": href, "host": urllib.parse.urlsplit(href).hostname, "source_page": final_url, "status": "discovered_not_validated"})
        with (out / f"official_{args.page}_download_links.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["url", "host", "source_page", "status"])
            writer.writeheader(); writer.writerows(records)
        status = {"status": "official_page_downloaded", "url": final_url, "sha256": hashlib.sha256(data).hexdigest(), "candidate_links": len(records)}
    except (urllib.error.URLError, OSError) as error:
        status = {"status": "blocked_official_page", "url": url, "error": str(error), "tls_verification": "preserved", "required_domain": host}
    (logs / f"download_discovery_{args.page}.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    main()
