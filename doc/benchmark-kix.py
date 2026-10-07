#!/usr/bin/env python3
"""Small, single-thread comparison of IntaRNAsnap and default IntaRNA.

Requires Python 3, GNU time and a release IntaRNA binary. Run from any directory:
  python3 doc/benchmark-kix.py /path/to/IntaRNA --output results.json
"""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
import time


PAIRS = [("fhlA", "OxyS"), ("phoB", "GcvB"), ("ilvE", "GcvB.ST")]
COLUMNS = "start1,end1,start2,end2,E,hybridDB"


def read_fasta(path):
    lines = path.read_text().splitlines()
    assert sum(line.startswith(">") for line in lines) == 1, path
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return {"file": "handson/" + path.name, "header": lines[0][1:],
            "length_nt": len(sequence), "sequence": sequence,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def prediction(stdout):
    rows = list(csv.DictReader(io.StringIO(stdout), delimiter=";"))
    if not rows:
        return None
    assert len(rows) == 1, rows
    row = rows[0]
    result = {key: int(row[key]) for key in ("start1", "end1", "start2", "end2")}
    result.update(E_kcal_mol=float(row["E"]), hybridDB=row["hybridDB"])
    # All inputs use the default one-based, ascending coordinates.
    result["length_nt"] = max(result["end1"] - result["start1"] + 1,
                              result["end2"] - result["start2"] + 1)
    return result


def distribution(samples):
    return {"median": statistics.median(samples), "min": min(samples), "max": max(samples)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--time", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--build-description", default="unspecified")
    args = parser.parse_args()
    if args.repetitions < 1 or args.warmups < 0:
        parser.error("repetitions must be positive and warmups nonnegative")
    binary = args.binary.resolve()
    timer = args.time.resolve()
    fixtures = Path(__file__).resolve().parent / "handson"
    common = ["--threads=1", "--outMode=C", "--outCsvCols=" + COLUMNS,
              "--outNumber=1", "--default-log-file=/dev/null"]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    report = {
        "schema": 1,
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "version": subprocess.check_output([str(binary), "--version"], text=True).strip(),
        "build": args.build_description,
        "platform": platform.platform(), "machine": platform.machine(),
        "python": platform.python_version(),
        "timer": subprocess.check_output([str(timer), "--version"], text=True).splitlines()[0],
        "common_arguments": common,
        "repetitions": args.repetitions, "warmups_per_configuration": args.warmups,
        "timing": "perf_counter wall seconds around GNU time + binary; includes startup, folding, seed search and prediction",
        "memory": "GNU time %M: maximum resident set size of each child process, KiB",
        "noGU": "--outNoGUend=true; seedNoGU remains false for both programs",
        "baseline": "default IntaRNA (model X, mode H, outNoLP false)",
        "comparison": "IntaRNAsnap defaults (model X, mode K, outNoLP true, kineticScore A)",
        "length": "max(end1-start1+1, end2-start2+1), nt; default one-based coordinates",
        "deviations": "signed IntaRNAsnap minus default IntaRNA, within each GU setting; not a global-optimum error bound",
        "cases": [],
    }
    with tempfile.TemporaryDirectory(prefix="intarna-kix-benchmark-") as directory:
        rss_file = Path(directory) / "time.txt"
        for target, query in PAIRS:
            inputs = ["--target=" + str(fixtures / (target + ".fasta")),
                      "--query=" + str(fixtures / (query + ".fasta"))]
            case = {"name": target + "/" + query,
                    "target": read_fasta(fixtures / (target + ".fasta")),
                    "query": read_fasta(fixtures / (query + ".fasta")), "runs": [], "comparisons": []}
            configurations = []
            for no_gu in (False, True):
                for program in ("IntaRNA", "IntaRNAsnap"):
                    options = [] if program == "IntaRNA" else ["--personality=IntaRNAsnap"]
                    options.append("--outNoGUend=" + str(no_gu).lower())
                    configurations.append((program, no_gu, options))
                    case["runs"].append({"program": program, "noGU": no_gu,
                                         "arguments": options, "samples": []})
            # Warm each configuration; rotate execution order every round.
            for iteration in range(args.warmups + args.repetitions):
                for offset in range(len(configurations)):
                    index = (iteration + offset) % len(configurations)
                    program, no_gu, options = configurations[index]
                    run = case["runs"][index]
                    command = [str(binary), *common, *inputs, *options]
                    start = time.perf_counter()
                    process = subprocess.run([str(timer), "-f", "%M", "-o", str(rss_file),
                                              *command], env=env, text=True, capture_output=True, check=True)
                    elapsed = time.perf_counter() - start
                    result = prediction(process.stdout)
                    if "prediction" in run:
                        assert run["prediction"] == result, (case["name"], options, process.stdout)
                    run["prediction"] = result
                    if iteration >= args.warmups:
                        run["samples"].append({"wall_s": elapsed, "peak_rss_kib": int(rss_file.read_text())})
            for run in case["runs"]:
                run["wall_s"] = distribution([sample["wall_s"] for sample in run["samples"]])
                run["peak_rss_kib"] = distribution([sample["peak_rss_kib"] for sample in run["samples"]])
                # Outside timing, reevaluate each reported structure independently.
                if run["prediction"] is not None:
                    evaluated = subprocess.check_output([str(binary), *common, *inputs,
                                                         "--rri=" + run["prediction"]["hybridDB"]],
                                                        env=env, text=True)
                    assert prediction(evaluated) == run["prediction"], (case["name"], run)
                    run["structure_reevaluation_matches"] = True
            for index, no_gu in ((0, False), (2, True)):
                baseline, kix = case["runs"][index:index + 2]
                available = baseline["prediction"] is not None and kix["prediction"] is not None
                case["comparisons"].append({
                    "noGU": no_gu,
                    "delta_E_kcal_mol": round(kix["prediction"]["E_kcal_mol"] - baseline["prediction"]["E_kcal_mol"], 2) if available else None,
                    "delta_length_nt": kix["prediction"]["length_nt"] - baseline["prediction"]["length_nt"] if available else None,
                    "wall_ratio_kix_over_default": kix["wall_s"]["median"] / baseline["wall_s"]["median"],
                    "rss_ratio_kix_over_default": kix["peak_rss_kib"]["median"] / baseline["peak_rss_kib"]["median"],
                })
            report["cases"].append(case)
            print(case["name"] + ": " + json.dumps(case["comparisons"]), file=sys.stderr, flush=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
