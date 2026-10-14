#!/usr/bin/env python3
"""Compare actual kernel executables, then benchmark each in fresh processes.

Example: python3 run.py --dense /tmp/dense --hash /tmp/hash --output /tmp/results
Correctness runs write and byte-compare every callback, z and all pair masses.
Timing runs have only a callback counter, with no dump I/O in the core timer.
GNU time measures peak RSS including startup, mask/domain setup and the result.
Its small resident parent avoids inheriting the Python runner's RSS high-water
floor when the executable is launched directly through fork/exec.
All processes run sequentially, with alternating variant order per repetition.
"""

import argparse
import csv
import filecmp
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile


def cases(suite):
    """A full matrix plus larger sparse workloads and a multi-start subset."""
    if suite == "smoke":
        return [(32, d, length, "single", pairs, no_lp)
                for d, length, pairs, no_lp in
                [(100, 1, 0, 0), (1000, 7, 1, 0), (3750, 20, 1, 1), (10000, 7, 0, 1)]]
    matrix = [(n, d, length, "single", pairs, no_lp)
              for n, d, length, pairs, no_lp in itertools.product(
                  (64, 192), (100, 1000, 3750, 10000), (1, 7, 20), (0, 1), (0, 1))]
    matrix += [(64, d, 7, "multiple", pairs, no_lp)
               for d, pairs, no_lp in itertools.product((100, 1000, 3750, 10000), (0, 1), (0, 1))]
    matrix += [(1024, d, 20, "single", pairs, 0)
               for d, pairs in itertools.product((100, 1000), (0, 1))]
    return matrix


def invoke(executable, case, dump=None):
    n, density, length, starts, pairs, no_lp = case
    command = [str(executable), "--n", str(n), "--density", str(density),
               "--length", str(length), "--starts", starts, "--pairs", str(pairs),
               "--noLP", str(no_lp)]
    if dump is not None:
        command += ["--dump", str(dump)]
    rss_tag = "INTARNA_STORAGE_PEAK_RSS_KIB="
    wrapped = ["/usr/bin/time", "-f", rss_tag + "%M", "--"] + command
    result = subprocess.run(wrapped, check=True, text=True, capture_output=True,
                            env={**os.environ, "OMP_NUM_THREADS": "1"})
    rss_lines = [line for line in result.stderr.splitlines() if line.startswith(rss_tag)]
    if len(rss_lines) != 1:
        raise RuntimeError("expected exactly one GNU time peak-RSS measurement: " + result.stderr)
    measured = json.loads(result.stdout)
    # A directly forked child can retain Python's resident-size floor across
    # exec. GNU time forks the executable from its own much smaller process.
    measured["peak_rss_kib"] = int(rss_lines[0][len(rss_tag):])
    return measured


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense", required=True, type=Path)
    parser.add_argument("--hash", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", default=3, type=int)
    parser.add_argument("--suite", choices=("full", "smoke"), default="full")
    parser.add_argument("--cpu", type=int, help="pin this runner and its children to one CPU (Linux)")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.cpu is not None:
        if not hasattr(os, "sched_setaffinity"):
            parser.error("--cpu requires sched_setaffinity; use the platform's affinity tool instead")
        os.sched_setaffinity(0, {args.cpu})
    variants = {"dense": args.dense.resolve(), "hash": args.hash.resolve()}
    args.output.mkdir(parents=True, exist_ok=True)
    workloads = cases(args.suite)
    checks, samples = [], []
    with tempfile.TemporaryDirectory(prefix="seed-storage-") as temporary:
        dumps = {name: Path(temporary) / (name + ".txt") for name in variants}
        # Finish correctness comparison for every workload before collecting times.
        for index, case in enumerate(workloads):
            print(f"check {index + 1}/{len(workloads)}: {case}", flush=True)
            reference = None
            for name, executable in variants.items():
                result = invoke(executable, case, dumps[name])
                invariant = {key: value for key, value in result.items()
                             if key not in ("core_ms", "peak_rss_kib")}
                if reference is None:
                    reference = invariant
                elif invariant != reference:
                    raise RuntimeError(f"metadata/result mismatch in {case}: {invariant} != {reference}")
            if not filecmp.cmp(dumps["dense"], dumps["hash"], shallow=False):
                # Retain failing evidence outside the temporary directory.
                for name in variants:
                    dumps[name].replace(args.output / (name + "-mismatch.txt"))
                raise RuntimeError(f"exact callback/mass mismatch in {case}; dumps retained in {args.output}")
            checks.append({"case": index, **reference, "dump_sha256": sha256(dumps["dense"]),
                           "dump_bytes": dumps["dense"].stat().st_size})
        write_csv(args.output / "correctness.csv", checks)
        for index, case in enumerate(workloads):
            print(f"time {index + 1}/{len(workloads)}: {case}", flush=True)
            for executable in variants.values():
                invoke(executable, case)  # One unrecorded warmup per variant/workload.
            for repeat in range(args.repeats):
                order = list(variants) if (index + repeat) % 2 == 0 else list(reversed(variants))
                for name in order:
                    result = invoke(variants[name], case)
                    if any(result[key] != value for key, value in checks[index].items()
                           if key not in ("case", "dump_sha256", "dump_bytes")):
                        raise RuntimeError(f"result changed during timing: {name}, {case}")
                    samples.append({"case": index, "variant": name, "repeat": repeat, **result})
            # Preserve completed measurements if a later run fails.
            write_csv(args.output / "samples.csv", samples)
    summary = []
    for index, case in enumerate(workloads):
        row = {key: value for key, value in checks[index].items()
               if key not in ("dump_sha256", "dump_bytes")}
        for name in variants:
            selected = [r for r in samples if r["case"] == index and r["variant"] == name]
            row[name + "_median_ms"] = statistics.median(r["core_ms"] for r in selected)
            row[name + "_min_ms"] = min(r["core_ms"] for r in selected)
            row[name + "_max_ms"] = max(r["core_ms"] for r in selected)
            row[name + "_median_rss_kib"] = statistics.median(r["peak_rss_kib"] for r in selected)
        row["hash_time_over_dense"] = row["hash_median_ms"] / row["dense_median_ms"]
        row["hash_rss_over_dense"] = row["hash_median_rss_kib"] / row["dense_median_rss_kib"]
        summary.append(row)
    write_csv(args.output / "summary.csv", summary)
    (args.output / "run.json").write_text(json.dumps({
        "executables": {name: {"path": str(path), "sha256": sha256(path)}
                        for name, path in variants.items()},
        "suite": args.suite, "repeats": args.repeats, "workloads": len(workloads),
        "warmups_per_variant_and_workload": 1,
        "platform": platform.platform(), "python": platform.python_version(),
        "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "OMP_NUM_THREADS": "1", "correctness": "exact ordered callbacks, partition, all pair masses",
        "runtime": "steady_clock around compute only; callback counter, no dump I/O",
        "memory": "GNU time %M peak RSS including setup/result; small launcher avoids Python fork RSS floor"
    }, indent=2) + "\n")
    print(f"Compared {len(checks)} workloads exactly; wrote {len(samples)} samples to {args.output}")


if __name__ == "__main__":
    main()
