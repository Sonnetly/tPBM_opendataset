#!/usr/bin/env python3
"""Acquire or reuse exact public responses and run the projection description."""
import argparse
from datetime import datetime, timezone
import json
import sys
from acquire_connectivity import acquire, RESULTS
from analyze_connectivity import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Fetch absent official RMA JSON responses with HTTPS")
    args = parser.parse_args()
    try:
        acquire(download=args.download)
        run()
        return 0
    except Exception as error:
        RESULTS.mkdir(parents=True, exist_ok=True)
        status_path = RESULTS / "status.json"
        previous = json.loads(status_path.read_text()) if status_path.exists() else None
        status = {"run_at_UTC": datetime.now(timezone.utc).isoformat(), "status": "failed_current_run",
                  "actual_data_analyzed_this_run": False, "error": f"{type(error).__name__}: {error}",
                  "previous_status": previous,
                  "note": "Existing result files may belong to the prior successful run; this failed run does not validate them. Raw cached responses were preserved."}
        status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False))
        print(status["error"], file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
