"""Regression tests for the Python replacements of the R output tools."""

import csv
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_script(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pvalue = load_script("csv_pvalue", "IntaRNA_CSV_p-value.py")
plot = load_script("plot_regions", "IntaRNA_plotRegions.py")


class OutputToolsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)

    def write(self, name, text):
        path = self.path / name
        path.write_text(text, encoding="utf-8")
        return path

    def run_cli(self, script, *arguments, success=True):
        result = subprocess.run([sys.executable, str(ROOT / script), *map(str, arguments)],
                                capture_output=True, text=True,
                                env={**os.environ, "MPLCONFIGDIR": str(self.path / "mpl"),
                                     "DISPLAY": ""})
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        return result

    def energy_table(self, column="E"):
        # Deterministic GEV quantiles: independent inverse-CDF formula, xi=0.2.
        probabilities = (np.arange(1000) + 0.5) / 1000
        energies = -(10 + 2 * ((-np.log(probabilities)) ** -0.2 - 1) / 0.2)
        source = self.path / "energies.csv"
        with source.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter=";")
            writer.writerow(["id", column, "annotation"])
            for index, energy in enumerate(energies):
                writer.writerow(["{:04}".format(index), energy, 'RNA; "α"'])
        return source, energies

    def test_gev_tail_sign_gumbel_and_endpoints(self):
        energies = np.array([-3.0, -8.0, -12.0, -40.0])
        for xi in (-0.2, 0.2):
            # Original R/evir formula, using explicit xi rather than SciPy's c.
            y = 1 + xi * (-energies[:3] - 10) / 2
            expected = -np.expm1(-y ** (-1 / xi))
            np.testing.assert_allclose(pvalue.gev_pvalues(energies[:3], (-xi, 10, 2)), expected)
        expected = -np.expm1(-np.exp((energies + 10) / 2))
        np.testing.assert_allclose(pvalue.gev_pvalues(energies, (0, 10, 2)), expected)
        np.testing.assert_equal(pvalue.gev_pvalues([-100, 100], (0.2, 10, 2)), [0, 1])
        # Tiny upper tails must not be lost to 1-cdf cancellation.
        self.assertGreater(pvalue.gev_pvalues([-60], (0, 0, 1))[0], 0)

    def test_fit_recovers_known_distribution(self):
        _, energies = self.energy_table()
        params = pvalue.fit_gev(energies)
        np.testing.assert_allclose(params, [-0.2, 10, 2], atol=0.015, rtol=0)
        expected = 1 - (np.arange(1000) + 0.5) / 1000
        np.testing.assert_allclose(pvalue.gev_pvalues(energies, params), expected, atol=0.002)

    def test_bh_unsorted_ties_and_bounds(self):
        np.testing.assert_allclose(pvalue.adjust_bh([0.04, 0.01, 0.03, 0.01, 1, 0]),
                                   [0.048, 0.02, 0.045, 0.02, 1, 0])

    def test_legacy_r_reference(self):
        # R 4.4.3, original R/IntaRNA_CSV_p-value.R at 2c4d0ec, 1000 inverse-CDF
        # quantiles per shape, location=10, scale=2. Compare rows 0/99/499/899/999.
        # Optimizers differ; BH can amplify the seven-decimal tail differences.
        references = {
            -0.2: ([0.9994943, 0.9006078, 0.5007587, 0.1001753, 0.0004651],
                   [0.9994943, 0.9994715, 0.9994715, 0.9918347, 0.4651]),
            0.0: ([0.9995020, 0.9006536, 0.5005593, 0.1003176, 0.0004918],
                  [0.9995020, 0.9995020, 0.9991204, 0.9932436, 0.4918]),
            0.2: ([0.9995087, 0.9006690, 0.5003434, 0.1003388, 0.0004986],
                  [0.9995087, 0.9995087, 0.9986894, 0.9934535, 0.4986]),
        }
        quantiles = (np.arange(1000) + 0.5) / 1000
        for xi, (expected_p, expected_fdr) in references.items():
            with self.subTest(xi=xi):
                energies = (-(10 + 2 * ((-np.log(quantiles)) ** -xi - 1) / xi)
                            if xi else -(10 - 2 * np.log(-np.log(quantiles))))
                actual = np.round(pvalue.gev_pvalues(energies, pvalue.fit_gev(energies)), 7)
                fdr = np.round(pvalue.adjust_bh(actual), 7)
                indices = [0, 99, 499, 899, 999]
                np.testing.assert_allclose(actual[indices], expected_p, atol=1e-4, rtol=0)
                np.testing.assert_allclose(fdr[indices], expected_fdr, atol=0.002, rtol=0)

    def test_csv_default_in_place_preserves_fields(self):
        source, energies = self.energy_table()
        source.chmod(0o640)
        self.run_cli("IntaRNA_CSV_p-value.py", source)
        with source.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=";")
            rows = list(reader)
            self.assertEqual(reader.fieldnames, ["id", "E", "annotation", "p-value", "fdr"])
        self.assertEqual(len(rows), 1000)
        self.assertEqual(rows[0]["id"], "0000")
        self.assertEqual(rows[0]["annotation"], 'RNA; "α"')
        self.assertEqual(float(rows[0]["E"]), energies[0])
        self.assertLess(float(rows[-1]["p-value"]), float(rows[0]["p-value"]))
        self.assertTrue(all(0 <= float(row["p-value"]) <= float(row["fdr"]) <= 1 for row in rows))
        self.assertEqual(source.stat().st_mode & 0o777, 0o640)
        self.assertNotIn(b"\r", source.read_bytes())

    def test_csv_separate_output_custom_energy_and_rounding(self):
        source, _ = self.energy_table("normalized energy")
        before = source.read_bytes()
        output = self.path / "extended.csv"
        self.run_cli("IntaRNA_CSV_p-value.py", source, output, "normalized energy")
        self.assertEqual(source.read_bytes(), before)
        with output.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=";"))
        expected = pvalue.adjust_bh([float(row["p-value"]) for row in rows])
        self.assertEqual([row["fdr"] for row in rows], [format(value, ".7f") for value in expected])

    def test_invalid_csv_does_not_replace_input_or_output(self):
        invalid = ["", "id;score\na;-5\n", "E;E\n-5;-6\n", "E\n",
                   "E\n-5\n", "E\n-5\n-5\n-5\n", "E\n-5\n-6\nNA\n",
                   "E\n-5\n-6\nnan\n", "E\n-5\n-6\ninf\n", "E\n-5;extra\n",
                   "E;p-value\n-5;0.1\n"]
        for index, content in enumerate(invalid):
            with self.subTest(content=content):
                source = self.write("invalid{}.csv".format(index), content)
                output = self.write("output.csv", "keep me\n")
                for destination in (source, output):
                    before = destination.read_bytes()
                    with self.assertRaises(ValueError):
                        pvalue.extend_csv(source, destination)
                    self.assertEqual(destination.read_bytes(), before)

    def test_fit_failure_does_not_overwrite(self):
        source, _ = self.energy_table()
        before = source.read_bytes()
        with patch.object(pvalue, "fit_gev", side_effect=ValueError("fit failed")):
            with self.assertRaisesRegex(ValueError, "fit failed"):
                pvalue.extend_csv(source, source)
        self.assertEqual(source.read_bytes(), before)

    def test_publish_failure_cleans_temporary_output(self):
        source, _ = self.energy_table()
        before = source.read_bytes()
        with patch.object(pvalue.os, "replace", side_effect=OSError("cannot replace")):
            with self.assertRaises(OSError):
                pvalue.extend_csv(source, source)
        self.assertEqual(source.read_bytes(), before)
        self.assertEqual(list(self.path.iterdir()), [source])

    def test_cli_errors_and_help(self):
        for script in ("IntaRNA_CSV_p-value.py", "IntaRNA_plotRegions.py"):
            self.assertIn("usage:", self.run_cli(script, "--help").stdout)
            self.run_cli(script, success=False)
        bad = self.write("bad.csv", "E\nNaN\n")
        self.assertIn("nonfinite", self.run_cli("IntaRNA_CSV_p-value.py", bad, success=False).stderr)
        self.run_cli("IntaRNA_plotRegions.py", bad, "1", self.path / "out.png", success=False)
        self.run_cli("IntaRNA_plotRegions.py", bad, "1", self.path / "out.xyz", success=False)

    def test_region_selectors_and_default_row_limit(self):
        source = ROOT / "examples/regions.csv"
        self.assertEqual(plot.read_regions(source, plot.read_parameters("1"))[0], ("b0001", 266, 273))
        self.assertEqual(plot.read_regions(source, plot.read_parameters("2"))[0], ("query", 116, 123))
        source = self.write("many.csv", "id1;start1;end1\n" + "001;2;5\n" * 205)
        regions = plot.read_regions(source, plot.read_parameters("1"))
        self.assertEqual(len(regions), 200)
        self.assertEqual(regions[0][0], "001")

    def test_parameters_custom_columns_comments_defaults_and_limits(self):
        path = self.write("plot.params", "# comment\nid = name\nstart = left\nend = right\n"
                          "title = custom title\nxmin = -5\nxmax = 8\nxvline = 0\nrows = 2\n")
        parameters = plot.read_parameters(str(path))
        source = self.write("custom.csv", "name;left;right\nx;4;-2\nx;0;0\ny;10;15\n")
        regions = plot.read_regions(source, parameters)
        self.assertEqual(regions, [("x", 4, -2), ("x", 0, 0)])
        figure = plot.create_figure(regions, parameters)
        import matplotlib.pyplot as plt
        self.addCleanup(plt.close, figure)
        coverage, intervals = figure.axes
        self.assertEqual(coverage.get_title(), "custom title")
        self.assertEqual(intervals.get_xlim(), (-5, 8))
        self.assertEqual([tick.get_text() for tick in intervals.get_yticklabels()], ["x", "x"])
        self.assertEqual(list(intervals.lines[0].get_xdata()), [4, -2])
        self.assertEqual(list(intervals.lines[-1].get_xdata()), [0, 0])
        minimal = self.write("minimal.params", "id=name\nstart=left\nend=right\n")
        self.assertEqual(plot.read_parameters(str(minimal))["rows"], 200)

    def test_invalid_parameters_and_regions(self):
        for content in ("id=x\n", "id=x\nstart=s\nend=e\nrows=0\n",
                        "id=x\nstart=s\nend=e\nrows=2.5\n", "id=x\nid=y\n",
                        "id=x\nstart=s\nend=e\nxmin=nan\n", "id=x\nunknown=3\n"):
            with self.subTest(content=content):
                path = self.write("invalid.params", content)
                with self.assertRaises(ValueError):
                    plot.read_parameters(str(path))
        for content in ("", "id1;start1;end1\n", "id1;start1\na;3\n",
                        "id1;start1;end1\na;3;NA\n", "id1;start1;end1\na;3.5;4\n",
                        "id1;start1;end1\na;3;4;extra\n"):
            with self.subTest(content=content):
                with self.assertRaises(ValueError):
                    plot.read_regions(self.write("invalid.csv", content), plot.read_parameters("1"))

    def test_density_inclusive_reversed_clipped_and_count_scaled(self):
        regions = [("a", -1, 1), ("b", 1, 0)]
        x, actual = plot.coverage_density(regions, -1, 1)
        # Five bases including both endpoints, with repeated covered positions.
        bases = np.array([-1, 0, 1, 1, 0])
        bandwidth = 0.9 * min(np.std(bases, ddof=1), 1 / 1.34) * len(bases) ** -0.2
        expected = sum(np.exp(-0.5 * ((x - base) / bandwidth) ** 2)
                       / (math.sqrt(2 * math.pi) * bandwidth) for base in bases)
        np.testing.assert_allclose(actual, expected)
        # X limits remove out-of-range samples before density estimation.
        np.testing.assert_allclose(plot.coverage_density(regions, 0, 1)[1],
                                   plot.coverage_density([("a", 0, 1), ("b", 1, 0)], 0, 1)[1])
        self.assertTrue(np.all(plot.coverage_density(regions, 10, 20)[1] == 0))

    def test_single_base_region(self):
        parameters = plot.read_parameters("1")
        figure = plot.create_figure([("one", 0, 0)], parameters)
        import matplotlib.pyplot as plt
        self.addCleanup(plt.close, figure)
        self.assertEqual(figure.axes[0].get_xlim(), (-0.5, 0.5))
        self.assertTrue(np.all(np.isfinite(figure.axes[0].lines[0].get_ydata())))

    def test_every_output_format_without_display(self):
        signatures = {".png": b"\x89PNG", ".pdf": b"%PDF", ".svg": b"<?xml",
                      ".eps": b"%!PS", ".ps": b"%!PS", ".jpeg": b"\xff\xd8", ".tiff": b"II"}
        for extension, signature in signatures.items():
            with self.subTest(extension=extension):
                output = self.path / ("plot" + extension)
                self.run_cli("IntaRNA_plotRegions.py", ROOT / "examples/regions.csv",
                             ROOT / "examples/regions.params", output)
                data = output.read_bytes()
                self.assertGreater(len(data), 500)
                self.assertTrue(data.startswith(signature))


if __name__ == "__main__":
    unittest.main()
