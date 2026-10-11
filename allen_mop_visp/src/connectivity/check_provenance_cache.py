"""Meaningful synthetic check: altered cached API responses must not be trusted."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import acquire_connectivity as acquisition


class CacheIntegrityTests(unittest.TestCase):
    def test_initial_metadata_preserved_and_mutated_response_rejected(self):
        with tempfile.TemporaryDirectory(prefix="connectivity_cache_integrity_") as directory:
            root = Path(directory); inputs = root / "input"; inputs.mkdir()
            raw = inputs / "synthetic.json"
            response = {"success": True, "num_rows": 1, "total_rows": 1, "msg": [{"synthetic": 1}]}
            raw.write_text(json.dumps(response))
            with patch.object(acquisition, "ROOT", root), patch.object(acquisition, "INPUT", inputs):
                _, first = acquisition.fetch("model::Synthetic", raw.name)
                provenance = raw.with_suffix(".provenance.json")
                first["retrieved_at_UTC"] = "2026-10-11T00:00:00+00:00"
                first["status_code"] = 200
                provenance.write_text(json.dumps(first))
                _, reused = acquisition.fetch("model::Synthetic", raw.name)
                self.assertEqual(reused["retrieved_at_UTC"], first["retrieved_at_UTC"])
                self.assertEqual(reused["status_code"], 200)
                self.assertTrue(reused["cached_sha256_matches_initial_record"])
                response["msg"][0]["synthetic"] = 2
                raw.write_text(json.dumps(response))
                with self.assertRaisesRegex(ValueError, "SHA256 changed"):
                    acquisition.fetch("model::Synthetic", raw.name)
                self.assertEqual(json.loads(provenance.read_text())["sha256"], first["sha256"])


if __name__ == "__main__":
    unittest.main()
