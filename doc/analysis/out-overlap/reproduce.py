#!/usr/bin/env python3
"""Capture issue 212 observations without changing IntaRNA.

Historical defects and their repairs are recorded, not blessed as expectations.
Regional examples accept either historical output or the revised input rejection.
Only the small independent oracle and the non-exhausting controls are asserted.
Requires Python 3 and an already built IntaRNA executable.
"""

import argparse
import csv
import io
import itertools
import json
import pathlib
import subprocess


COLS = "start1,end1,start2,end2,E,ED1,ED2"
COMMON = ["--threads=1", "--outMode=C", "--outCsvCols=" + COLS]
TOY = ["--energy=B", "--acc=N", "--outDeltaE=100"]
FORBIDDEN = {"N": "12", "T": "2", "Q": "1", "B": ""}


def overlaps(a, b, axis):
    return max(a["start" + axis], b["start" + axis]) <= min(
        a["end" + axis], b["end" + axis])


def conflicts(rows, mode):
    return [
        [i + 1, j + 1, axis]
        for (i, a), (j, b) in itertools.combinations(enumerate(rows), 2)
        for axis in FORBIDDEN[mode] if overlaps(a, b, axis)
    ]


def site(row):
    return tuple(row[key] for key in ("start1", "end1", "start2", "end2"))


def base_pair_oracle(target, query):
    """Enumerate every antiparallel structure for the four-base toy example.

    No accessibility, seed or lonely-pair restrictions apply. At this length
    every interior loop fits the default loop limit. Each pair contributes -1.
    Keep the best energy for each pair of inclusive sequence intervals.
    """
    best = {}
    pairs = [(i, j) for i, t in enumerate(target, 1)
             for j, q in enumerate(query, 1) if t + q in
             {"AU", "UA", "CG", "GC", "GU", "UG"}]

    def extend(structure):
        i, j = structure[-1]
        bounds = (structure[0][0], i, j, structure[0][1])
        best[bounds] = min(best.get(bounds, 0), -len(structure))
        for ni, nj in pairs:
            if ni > i and nj < j:
                extend(structure + [(ni, nj)])

    for pair in pairs:
        extend([pair])
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    observations = {}

    def run(name, options, mode="B", expected_status=0):
        command = COMMON + options + ["--outOverlap=" + mode]
        process = subprocess.run([binary] + command, capture_output=True,
                                 text=True, timeout=30)
        if expected_status is not None and (process.returncode == 0) != (expected_status == 0):
            raise RuntimeError(f"{name}: unexpected exit {process.returncode}\n"
                               + process.stdout + process.stderr)
        lines = [line for line in process.stdout.splitlines()
                 if line and not line.startswith("#")]
        rows = []
        if process.returncode == 0:
            reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter=";")
            if reader.fieldnames != COLS.split(","):
                raise RuntimeError(f"{name}: unexpected CSV header {reader.fieldnames}")
            for row in reader:
                rows.append({key: (int(value) if key.startswith(("start", "end"))
                                   else float(value)) for key, value in row.items()})
        observations[name] = {
            "arguments": command, "returncode": process.returncode,
            "stdout": process.stdout, "stderr": process.stderr,
            "rows": rows, "forbidden_overlaps": conflicts(rows, mode),
        }
        return rows

    # Every CLI predictor combination, with enough disjoint two-pair sites to
    # avoid the separately investigated ensemble exhaustion defect.
    configurations = []
    for model, modes, seeded in [
        ("S", "HM", False), ("S", "HMS", True), ("X", "HMRS", True),
        ("P", "HM", False), ("P", "HMS", True),
        ("B", "H", False), ("B", "H", True),
    ]:
        for mode in modes:
            configurations.append((model, mode, seeded))
    for model, mode, seeded in configurations:
        for overlap in "BNTQ":
            name = f"control_{model}_{mode}_{'seed' if seeded else 'noSeed'}_{overlap}"
            options = TOY + ["-t", "CCAACC", "-q", "GGAAGG", "-n", "2",
                             "--intLenMax=2", "--model=" + model, "--mode=" + mode,
                             "--helixMinBP=2", "--helixMaxBP=2",
                             "--seedBP=2" if seeded else "--noSeed"]
            rows = run(name, options, overlap)
            assert len(rows) == 2, (name, rows)
            assert all(row["E"] == -2 for row in rows), (name, rows)
            assert not conflicts(rows, overlap), (name, rows)

    for overlap in "BNTQ":
        run("target_regions_" + overlap, TOY + ["-t", "CCAACC", "-q", "GG",
            "--seedBP=2", "--tRegion=1-2,5-6", "-n", "10"], overlap,
            expected_status=None if overlap in "NT" else 0)
        run("query_regions_" + overlap, TOY + ["-t", "CC", "-q", "GGAAGG",
            "--seedBP=2", "--qRegion=1-2,5-6", "-n", "10"], overlap,
            expected_status=None if overlap in "NQ" else 0)
        run("regional_delta_" + overlap, ["--energy=B", "--acc=N", "-t", "CCAAC",
            "-q", "GG", "--noSeed", "--model=S", "--tRegion=1-2,5-5",
            "--outDeltaE=0", "-n", "10"], overlap,
            expected_status=None if overlap in "NT" else 0)
        run("window_rejection_" + overlap, TOY + ["-t", "CCAACCAACCAA", "-q", "GGAAGGAAGGAA",
            "--seedBP=2", "--intLenMax=2", "--windowWidth=10", "--windowOverlap=2",
            "-n", "2"], overlap, expected_status=0 if overlap == "B" else 1)

    run("automatic_regions", ["--energy=B", "--acc=C", "-t", "UAUCGGCC", "-q", "GG",
        "--seedBP=2", "--tRegionLenMax=4", "--outDeltaE=100", "-n", "10"], "N", expected_status=None)
    run("regions_vienna", ["--acc=N", "-t", "CCCACCC", "-q", "GGG", "--seedBP=2",
        "--tRegion=1-3,5-7", "--outDeltaE=100", "-n", "10"], "N", expected_status=None)
    run("explicit_per_region", TOY + ["-t", "CCAACC", "-q", "GG", "--seedBP=2",
        "--tRegion=1-2,5-6", "--outPerRegion", "-n", "1"], "N")

    seeded_missing = TOY + ["-t", "CCCAAU", "-q", "GGGUU", "--seedBP=2", "-n", "100"]
    for mode in "HM":
        for overlap in "BNTQ":
            run(f"missing_seed_{mode}_{overlap}", seeded_missing + ["--mode=" + mode], overlap)
        run("blocked_seed_" + mode, seeded_missing + ["--mode=" + mode,
            "--tAccConstr=b:1-3", "--qAccConstr=b:1-3"], "N")

    for overlap, target, query in [("N", "CCAU", "CGGU"), ("T", "GUA", "GUA"),
                                   ("Q", "GCGC", "AGUG")]:
        options = TOY + ["-t", target, "-q", query, "--model=S", "--mode=M",
                         "--noSeed", "--outNoLP=false", "-n", "1000"]
        all_rows = run("oracle_B_" + overlap, options)
        oracle = base_pair_oracle(target, query)
        assert {site(row): row["E"] for row in all_rows} == oracle
        assert len(all_rows) == len(oracle)
        selected = run("oracle_" + overlap, options, overlap)
        observations["oracle_" + overlap]["omitted_compatible_sites"] = [
            row for row in all_rows if all(not overlaps(row, old, axis)
                for old in selected for axis in FORBIDDEN[overlap])]

    for model, mode, seed in [("S", "H", ["--noSeed"]),
                              ("S", "M", ["--noSeed"]),
                              ("S", "H", ["--seedBP=2", "--seedMaxUP=2"]),
                              ("X", "H", ["--seedBP=2", "--seedMaxUP=2"])]:
        for overlap in "BT":
            name = f"gu_filter_{model}_{mode}_{seed[0]}_{overlap}"
            rows = run(name, TOY + ["-t", "UUGA", "-q", "CAUU", "--model=" + model,
                "--mode=" + mode, "--outNoGUend", "--outNoLP=false", "-n", "10"] + seed, overlap)
            observations[name]["gu_ended_rows"] = [i + 1 for i, row in enumerate(rows)
                if any("UUGA"[row[t] - 1] + "CAUU"[row[q] - 1] in ("GU", "UG")
                       for t, q in [("start1", "end2"), ("end1", "start2")])]

    for mode in "HM":
        for overlap in "BNTQ":
            run(f"ensemble_exhaustion_{mode}_{overlap}", TOY + ["-t", "CC", "-q", "GG",
                "--model=P", "--mode=" + mode, "--noSeed", "-n", "2"], overlap)
    for overlap in "BN":
        run("ensemble_accessibility_" + overlap, ["--energy=B", "--acc=C", "-t", "AGAGC",
            "-q", "GAUUC", "--model=P", "--mode=H", "--noSeed", "--outDeltaE=100",
            "-n", "100" if overlap == "B" else "2", "--outMaxE=100" if overlap == "B"
            else "--outMaxE=0"], overlap)
    run("blocking_can_bridge", TOY + ["-t", "CCGG", "-q", "CCGG", "--noSeed",
        "--model=S", "--mode=M", "--tAccConstr=b:2-3", "-n", "2"], "N")

    data = {"version": subprocess.check_output([binary, "--version"], text=True).strip(),
            "control_configurations": len(configurations),
            "control_runs": len(configurations) * 4,
            "independent_oracles": 3, "invocations": len(observations),
            "observations": observations}
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"Captured {len(observations)} runs; {len(configurations) * 4} controls and "
          f"3 independent site oracles passed. Observations: {args.output}")


if __name__ == "__main__":
    main()
