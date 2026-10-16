"""Numeric checks using only artificial data."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cell_deconvolution import V1Model
from cell_deconvolution.v1 import extract_block_amf_from_bedgraph, _nnls_sqrt_cov, _r2_sqrt_cov
from cell_deconvolution.selection import compute_specificity_scores, select_blocks


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        spec = importlib.util.spec_from_file_location("synthetic", ROOT / "examples" / "generate_synthetic.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.generate(self.directory)
        self.model = V1Model(self.directory / "atlas.tsv")

    def test_known_mixture(self):
        result = self.model.deconv_bedgraph(self.directory / "mixture.bedGraph")
        self.assertAlmostEqual(result["component_a"], 0.25)
        self.assertAlmostEqual(result["component_b"], 0.75)
        self.assertAlmostEqual(result["raw_sum"], 1.0)
        self.assertAlmostEqual(result["R2"], 1.0)
        self.assertEqual(result["n_blocks"], 60)

    def test_coverage_and_missing_values(self):
        betas = self.model.A @ np.array([0.25, 0.75])
        cov = np.full(60, 3.0)
        cov[:10] = 2
        self.assertEqual(self.model.deconv(betas, cov)["n_blocks"], 50)
        betas[10] = np.nan
        self.assertIsNone(self.model.deconv(betas, cov))

    def test_zero_signal(self):
        result = self.model.deconv(np.zeros(60), np.ones(60) * 10)
        self.assertEqual(result["raw_sum"], 0)
        self.assertEqual(result["R2"], 0)
        self.assertEqual(result["component_a"], 0)

    def test_half_open_intervals(self):
        blocks = pd.DataFrame({"chr": ["synthetic"] * 3, "start": [0, 5, 20], "end": [5, 10, 30]})
        path = self.directory / "boundaries.bedGraph"
        path.write_text("track type=bedGraph\nsynthetic\t0\t1\t0\t1\t3\n"
                        "synthetic\t4\t5\t0\t3\t1\nsynthetic\t5\t6\t0\t9\t1\n"
                        "synthetic\t10\t11\t0\t50\t50\n")
        beta, cov = extract_block_amf_from_bedgraph(path, blocks)
        np.testing.assert_allclose(beta[:2], [0.5, 0.9])
        np.testing.assert_array_equal(cov, [8, 10, 0])
        self.assertTrue(np.isnan(beta[2]))

    def test_varying_coverage(self):
        matrix = np.vstack([np.eye(2), [[0.5, 0.2], [0.1, 0.8]]])
        cov = np.array([3., 8., 20., 60.])
        observed = matrix @ np.array([0.3, 0.6])
        coef = _nnls_sqrt_cov(matrix, cov, observed)
        np.testing.assert_allclose(coef, [0.3, 0.6])
        self.assertAlmostEqual(_r2_sqrt_cov(matrix, cov, observed, coef), 1.)

    def test_reference_order_guard(self):
        with self.assertRaises(ValueError):
            V1Model(self.directory / "atlas.tsv", {"celltypes": ["component_b", "component_a"]})

    def test_marker_selection(self):
        matrix = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0.8, 0.1, 0.1]])
        _, owner, score = compute_specificity_scores(matrix, ["a", "b", "c"])
        np.testing.assert_array_equal(select_blocks(owner, score, 3, 1), [0, 1, 2])

    def test_cli(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
        result = subprocess.run([sys.executable, "-m", "cell_deconvolution",
            "--atlas", str(self.directory / "atlas.tsv"),
            "--bedgraph", str(self.directory / "mixture.bedGraph")],
            env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertAlmostEqual(json.loads(result.stdout)["component_a"], 0.25)
        self.assertNotIn(str(self.directory), result.stdout)


if __name__ == "__main__":
    unittest.main()
