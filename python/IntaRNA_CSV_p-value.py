#!/usr/bin/env python3
"""Append genome-wide GEV p-values and Benjamini-Hochberg FDRs to IntaRNA CSV."""

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import tempfile


def fit_gev(energies):
    """Fit negated energies; return SciPy's (shape=-xi, location, scale)."""
    import numpy as np
    from scipy.optimize import minimize
    from scipy.stats import genextreme

    # The background must contain enough distinct, finite observations to fit.
    values = np.asarray(energies, dtype=float)
    if not np.all(np.isfinite(values)) or np.unique(values).size < 3:
        raise ValueError("GEV fitting requires at least three distinct finite energies")
    data = -values
    scale = float(np.std(data, ddof=1) * math.sqrt(6) / math.pi)
    location = float(np.mean(data) - 0.57722 * scale)
    if not math.isfinite(scale) or scale <= 0 or not math.isfinite(location):
        raise ValueError("energy range is too large for GEV fitting")

    # Retain the R script's initial xi=0.1 and moment-based scale/location.
    # SciPy uses the opposite shape sign and handles the Gumbel limit (xi=0).
    def optimizer(function, initial, args=(), disp=0):
        result = minimize(function, initial, args=args, method="Nelder-Mead",
                          options={"maxiter": 10000, "maxfev": 20000,
                                   "xatol": 1e-8, "fatol": 1e-8})
        if not result.success:
            raise ValueError("GEV fit did not converge: " + result.message)
        return result.x

    params = genextreme.fit(data, -0.1, loc=location, scale=scale, optimizer=optimizer)
    if (not np.all(np.isfinite(params)) or params[2] <= 0
            or not np.all(np.isfinite(genextreme.logpdf(data, *params)))):
        raise ValueError("GEV fit is not finite or does not cover all energies")
    return params


def gev_pvalues(energies, params):
    """Return P(E <= observed E), using the upper tail of the negated energies."""
    import numpy as np
    from scipy.stats import genextreme

    # sf avoids cancellation for significant interactions and handles support ends.
    values = genextreme.sf(-np.asarray(energies, dtype=float), *params)
    if not np.all(np.isfinite(values)):
        raise ValueError("GEV fit produced nonfinite p-values")
    return values


def adjust_bh(pvalues):
    """Adjust p-values by Benjamini-Hochberg, retaining the original row order."""
    result = [0.0] * len(pvalues)
    previous = 1.0
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        previous = min(previous, pvalues[index] * len(order) / rank)
        result[index] = previous
    return result


def extend_csv(source, destination, energy_column="E"):
    """Validate and fit the entire table before atomically publishing the output."""
    with open(source, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter=";", strict=True)
        header = next(reader, [])
        if not header or len(set(header)) != len(header):
            raise ValueError("CSV must have a nonempty header with unique column names")
        if energy_column not in header:
            raise ValueError("energy column {!r} is missing".format(energy_column))
        if "p-value" in header or "fdr" in header:
            raise ValueError("CSV already contains a 'p-value' or 'fdr' column")
        column = header.index(energy_column)
        rows, energies = [], []
        for row in reader:
            if not row:  # Allow blank lines, as in read.csv().
                continue
            if len(row) != len(header):
                raise ValueError("CSV line {} has the wrong number of fields".format(reader.line_num))
            try:
                energy = float(row[column])
            except ValueError:
                raise ValueError("invalid energy on CSV line {}".format(reader.line_num)) from None
            if not math.isfinite(energy):
                raise ValueError("nonfinite energy on CSV line {}".format(reader.line_num))
            rows.append(row)
            energies.append(energy)

    # Preserve the historical seven-decimal rounding *before* BH adjustment.
    pvalues = [round(float(p), 7) for p in gev_pvalues(energies, fit_gev(energies))]
    fdrs = adjust_bh(pvalues)
    destination = Path(destination).resolve()
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8",
                                         dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.writer(handle, delimiter=";", lineterminator="\n")
            writer.writerow(header + ["p-value", "fdr"])
            for row, pvalue, fdr in zip(rows, pvalues, fdrs):
                writer.writerow(row + [format(pvalue, ".7f"), format(fdr, ".7f")])
        if destination.exists():
            shutil.copymode(destination, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        "Use a sufficiently large, diverse genome-wide prediction set as the background. "
        "If output is omitted, input is overwritten after a successful fit. "
        "Requires Python 3.10+, NumPy and SciPy (see python/requirements-output-tools.txt)."))
    parser.add_argument("input", help="semicolon-separated IntaRNA CSV")
    parser.add_argument("output", nargs="?", help="output CSV (default: overwrite input)")
    parser.add_argument("energy_column", nargs="?", default="E", help="energy column (default: E)")
    args = parser.parse_args(argv)
    try:
        extend_csv(args.input, args.output or args.input, args.energy_column)
    except ImportError as error:
        parser.exit(1, "Missing Python dependency: {}. Install numpy and scipy.\n".format(error))
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        parser.exit(1, "Error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
