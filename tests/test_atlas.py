"""Tiny synthetic atlas construction through the public command line."""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd


class AtlasTests(unittest.TestCase):
    def test_small_reference_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            n = 16
            with gzip.open(root / "index.tsv.gz", "wt") as stream:
                stream.write("".join("synthetic\t%d\n" % i for i in range(n)))
            rows = []
            for component, count in [("component_a", 0), ("component_b", 10)]:
                path = root / (component + ".beta")
                np.tile([count, 10], (n, 1)).astype("uint8").tofile(path)
                rows.append({"cell_type": component, "beta_file": str(path)})
            pd.DataFrame(rows).to_csv(root / "reference.tsv", sep="\t", index=False)
            pd.DataFrame({"beta_file": [rows[0]["beta_file"]]}).to_csv(
                root / "coverage.tsv", sep="\t", index=False)
            env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
            command = [sys.executable, "-m", "cell_deconvolution.atlas",
                "--reference-manifest", str(root / "reference.tsv"),
                "--coverage-manifest", str(root / "coverage.tsv"),
                "--cpg-index", str(root / "index.tsv.gz"),
                "--output", str(root / "atlas.tsv")]
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            frame = pd.read_csv(root / "atlas.tsv", sep="\t")
            self.assertEqual(len(frame), 1)
            self.assertEqual(frame.loc[0, "component_a"], 0)
            self.assertEqual(frame.loc[0, "component_b"], 1)
            self.assertNotIn(str(root), result.stdout)


if __name__ == "__main__":
    unittest.main()
