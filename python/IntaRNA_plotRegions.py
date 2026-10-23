#!/usr/bin/env python3
"""Plot IntaRNA interaction regions and their smoothed nucleotide coverage."""

import argparse
import csv
import math
from pathlib import Path


from intarna_tools import csv_header, csv_rows


FORMATS = {".pdf", ".png", ".svg", ".eps", ".ps", ".jpeg", ".tiff"}


def read_parameters(selector):
    """Read a strand selector or the original key=value plotting parameters."""
    parameters = {"title": None, "xvline": None, "xmin": None, "xmax": None, "rows": 200}
    if selector in ("1", "2"):
        parameters.update({key: key + selector for key in ("id", "start", "end")})
        return parameters

    # Comments and optional values retain the original parameter-file syntax.
    seen = set()
    with open(selector, encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            key, separator, value = line.partition("=")
            key, value = key.strip(), value.strip()
            if not separator or key not in {"id", "start", "end", *parameters}:
                raise ValueError("unknown or malformed parameter on line {}".format(line_number))
            if key in seen:
                raise ValueError("duplicate parameter: " + key)
            seen.add(key)
            parameters[key] = value
    for key in ("id", "start", "end"):
        if not parameters.get(key):
            raise ValueError("parameter file must define " + key)
    for key in ("xmin", "xmax", "xvline"):
        value = parameters[key]
        if value is None or value in ("", "NA"):
            parameters[key] = None
        else:
            parameters[key] = float(value)
            if not math.isfinite(parameters[key]):
                raise ValueError(key + " must be finite")
    try:
        parameters["rows"] = int(parameters["rows"])
    except ValueError:
        raise ValueError("rows must be a positive integer") from None
    if parameters["rows"] < 1:
        raise ValueError("rows must be a positive integer")
    return parameters


def read_regions(source, parameters):
    """Read the selected rows in input order, keeping inclusive signed positions."""
    regions = []
    with open(source, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter=";", strict=True)
        header = csv_header(reader)
        for key in ("id", "start", "end"):
            if parameters[key] not in header:
                raise ValueError("missing CSV column: " + parameters[key])
        columns = [header.index(parameters[key]) for key in ("id", "start", "end")]
        for row in csv_rows(reader, header):
            try:
                start, end = (int(row[index]) for index in columns[1:])
            except ValueError:
                raise ValueError("positions must be integers on CSV line {}".format(reader.line_num)) from None
            regions.append((row[columns[0]], start, end))
            if len(regions) == parameters["rows"]:
                break
    if not regions:
        raise ValueError("CSV contains no interaction regions")
    return regions


def coverage_density(regions, xmin, xmax):
    """Gaussian count density of inclusive covered positions inside the x limits."""
    import numpy as np
    from scipy.stats import norm

    # ggplot's original scale limits filter positions before estimating density.
    positions = [np.arange(max(min(start, end), math.ceil(xmin)),
                           min(max(start, end), math.floor(xmax)) + 1, dtype=float)
                 for _, start, end in regions]
    positions = np.concatenate(positions)
    x = np.linspace(xmin, xmax, 512)
    if positions.size == 0:
        return x, np.zeros_like(x)

    # R's default bw.nrd0 bandwidth, including its zero-spread fallback.
    spread = float(np.std(positions, ddof=1)) if positions.size > 1 else 0.0
    iqr = float(np.subtract(*np.percentile(positions, [75, 25])))
    bandwidth = 0.9 * (min(spread, iqr / 1.34) or spread or abs(positions[0]) or 1.0)
    bandwidth *= positions.size ** -0.2
    # Evaluate in chunks so the number of plotted bases does not multiply memory by 512.
    density = np.zeros_like(x)
    for offset in range(0, len(positions), 4096):
        distances = (x[:, None] - positions[None, offset:offset + 4096]) / bandwidth
        density += norm.pdf(distances).sum(axis=1) / bandwidth
    return x, density


def create_figure(regions, parameters):
    """Create aligned coverage and region panels without requiring a display."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Preserve the input ordering and independent rows even for duplicate IDs.
    xmin = parameters["xmin"]
    xmax = parameters["xmax"]
    if xmin is None:
        xmin = min(min(start, end) for _, start, end in regions)
    if xmax is None:
        xmax = max(max(start, end) for _, start, end in regions)
    if xmin > xmax:
        raise ValueError("xmin must not exceed xmax")
    if xmin == xmax:
        xmin, xmax = xmin - 0.5, xmax + 0.5
    x, density = coverage_density(regions, xmin, xmax)
    region_height = max(2, len(regions) / 9)
    figure, (coverage, intervals) = plt.subplots(
        2, 1, figsize=(6, 2 + region_height), sharex=True,
        gridspec_kw={"height_ratios": [2, region_height]}, layout="constrained")
    coverage.plot(x, density, color="black", linewidth=1)
    coverage.set_ylabel("coverage")
    coverage.set_ylim(bottom=0)
    if parameters["title"]:
        coverage.set_title(parameters["title"], fontweight="bold")

    for index, (_, start, end) in enumerate(regions):
        intervals.plot([start, end], [index, index], color="dodgerblue",
                       linewidth=2, marker="o", markersize=3)
    intervals.set_yticks(range(len(regions)), labels=[region[0] for region in regions], fontsize=7)
    intervals.set_ylim(len(regions) - 0.5, -0.5)
    intervals.grid(axis="y", color="lightgray", linewidth=0.7)
    intervals.set_axisbelow(True)
    intervals.set_xlabel("position")
    intervals.axvline(xmin, color="black", linewidth=0.7)
    for axis in (coverage, intervals):
        axis.set_xlim(xmin, xmax)
        axis.yaxis.tick_right()
        axis.yaxis.set_label_position("right")
        axis.spines["top"].set_visible(False)
        axis.spines["left"].set_visible(False)
        line = parameters["xvline"]
        if line is not None and xmin <= line <= xmax:
            axis.axvline(line, color="black", linewidth=0.7)
    return figure


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        "Parameter files accept id, start, end, title, xvline, xmin, xmax and rows (default: 200). "
        "Requires Python 3.10+, NumPy, SciPy and Matplotlib; see python/requirements-output-tools.txt."))
    parser.add_argument("input", help="semicolon-separated IntaRNA CSV")
    parser.add_argument("selector", help="1 (target), 2 (query), or a key=value parameter file")
    parser.add_argument("output", help="plot filename: " + ", ".join(sorted(FORMATS)))
    args = parser.parse_args(argv)
    try:
        if Path(args.output).suffix.lower() not in FORMATS:
            raise ValueError("unsupported output extension; choose " + ", ".join(sorted(FORMATS)))
        parameters = read_parameters(args.selector)
        regions = read_regions(args.input, parameters)
        figure = create_figure(regions, parameters)
        import matplotlib.pyplot as plt
        try:
            figure.savefig(args.output, dpi=150)
        finally:
            plt.close(figure)
    except ImportError as error:
        parser.exit(1, "Missing Python dependency: {}. Install numpy, scipy and matplotlib.\n".format(error))
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        parser.exit(1, "Error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
