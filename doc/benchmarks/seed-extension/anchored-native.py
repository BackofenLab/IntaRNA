#!/usr/bin/env python3
"""Compare native OxyS/fhlA CLI runs; requires GNU time and keeps raw outputs."""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys

QUERY = "GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGC"
TARGET = "AGUUAGUCAAUGACCUUUUGCACCGCUUUGCGGUGCUUUCCUGGAACAACAAAAUGUCAUAUACACCGAUGAGUGAUCUCGGACAACAAGGGUUGUUCGACAUCACUCGGACA"
CASES = ("M-seeded", "M-unseeded", "H-seeded", "H-unseeded")
RELATIVE_TOLERANCE = 2e-12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read_pairs(path):
    """Keep labels and NA cells exact; compare every finite probability numerically."""
    with path.open(newline="") as stream:
        rows = list(csv.reader(stream, delimiter=";"))
    if not rows or len(rows) != len(TARGET) + 1:
        raise ValueError(f"unexpected probability matrix shape: {path}")
    if any(len(row) != len(QUERY) + 1 for row in rows):
        raise ValueError(f"unexpected probability matrix width: {path}")
    labels = (rows[0], [row[0] for row in rows[1:]])
    values = [None if value == "NA" else float(value)
              for row in rows[1:] for value in row[1:]]
    if any(value is not None and not math.isfinite(value) for value in values):
        raise ValueError(f"nonfinite probability: {path}")
    return labels, values


def compare_pairs(reference, current):
    if reference[0] != current[0] or len(reference[1]) != len(current[1]):
        raise ValueError("probability matrix labels or dimensions differ")
    maximum_absolute = maximum_relative = 0.0
    for before, after in zip(reference[1], current[1]):
        if before is None or after is None:
            if before is not after:
                raise ValueError("probability matrix NA cells differ")
            continue
        absolute = abs(before - after)
        relative = absolute / max(abs(before), abs(after)) if before or after else 0.0
        maximum_absolute = max(maximum_absolute, absolute)
        maximum_relative = max(maximum_relative, relative)
    return maximum_absolute, maximum_relative


def summarize(samples):
    summary = []
    for case in dict.fromkeys(row["case"] for row in samples):
        item = {"case": case}
        for variant in ("base", "candidate"):
            rows = [row for row in samples if row["case"] == case and row["variant"] == variant]
            if not rows:
                continue
            item[variant] = {"samples": len(rows)}
            for metric in ("elapsed_seconds", "user_seconds", "system_seconds", "max_rss_kib"):
                values = [row[metric] for row in rows]
                item[variant][metric] = {
                    "median": statistics.median(values), "min": min(values), "max": max(values)}
        if "base" in item and "candidate" in item:
            base_time = item["base"]["elapsed_seconds"]["median"]
            item["candidate_over_base_wall_time"] = (
                item["candidate"]["elapsed_seconds"]["median"] / base_time if base_time else None)
        summary.append(item)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="baseline IntaRNA executable")
    parser.add_argument("--candidate", type=Path, required=True, help="candidate IntaRNA executable")
    parser.add_argument("--output-dir", type=Path, required=True, help="new or empty results directory")
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--pair-output", action="store_true", help="request and numerically compare bpProb CSV")
    parser.add_argument("--time", type=Path, default=Path("/usr/bin/time"), help="GNU time executable")
    args = parser.parse_args()
    if args.repetitions < 1 or len(set(args.cases)) != len(args.cases):
        parser.error("repetitions must be positive and cases must not contain duplicates")
    output = args.output_dir.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("--output-dir must be new or empty; raw results are never overwritten")
    binaries = {"base": args.binary.resolve(), "candidate": args.candidate.resolve()}
    timer = args.time.resolve()
    for path in (*binaries.values(), timer):
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f"not an executable file: {path}")
    env = dict(os.environ, OMP_NUM_THREADS="1", LC_ALL="C", LANG="C")
    time_version = subprocess.check_output([str(timer), "--version"], text=True, env=env)
    if "GNU" not in time_version:
        parser.error("--time must identify itself as GNU time")
    output.mkdir(parents=True, exist_ok=True)
    cpu_info = Path("/proc/cpuinfo")
    cpu_model = next((line.split(":", 1)[1].strip() for line in cpu_info.read_text().splitlines()
                      if line.startswith("model name")), None) if cpu_info.exists() else None
    manifest = {
        "schema_version": 1, "started_at": stamp(), "invocation": [sys.executable, *sys.argv],
        "harness_sha256": digest(Path(__file__).resolve()), "gnu_time": time_version.splitlines()[0],
        "workload": {"target": TARGET, "query": QUERY, "target_id": "fhlA", "query_id": "OxyS",
                     "model": "P", "seed_bp": 7, "threads": 1, "cases": args.cases,
                     "pair_output": args.pair_output, "repetitions": args.repetitions, "warmups": 0},
        "environment": {"platform": platform.platform(), "cpu_model": cpu_model,
                        "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
                        "variables": {key: env[key] for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                      "LD_LIBRARY_PATH", "LC_ALL", "LANG") if key in env},
                        "timing_scope": "fresh complete process including accessibility and output",
                        "workstation_activity_and_cpu_frequency": "uncontrolled"},
        "verification": {"stdout": "byte-identical to first baseline sample of the same case",
                         "pair_csv": "exact labels/NA cells; every finite value compared",
                         "relative_tolerance": RELATIVE_TOLERANCE, "absolute_tolerance": 0},
        "raw_file_prefix": "{variant}-{case}-{repetition}",
        "timing_format": "%e\\t%U\\t%S\\t%M", "binaries": {}, "commands": {}}
    for variant, binary in binaries.items():
        manifest["binaries"][variant] = {
            "path": str(binary), "sha256": digest(binary),
            "version": subprocess.check_output([str(binary), "--version"], text=True, env=env).strip()}
        manifest["commands"][variant] = {}
        for case in args.cases:
            mode, seed_state = case.split("-")
            command = [str(binary), f"--query={QUERY}", f"--target={TARGET}", "--qId=OxyS", "--tId=fhlA",
                       "--model=P", f"--mode={mode}", "--threads=1", "--seedBP=7",
                       f"--noSeed={'false' if seed_state == 'seeded' else 'true'}"]
            if args.pair_output:
                command.append(f"--out=bpProb:{output / (variant + '-' + case + '.bpProb.csv')}")
            manifest["commands"][variant][case] = command
    save(output / "manifest.json", manifest)
    samples, reference_stdout, reference_pairs = [], {}, {}
    result = {"complete": False, "output_checks_passed": False}
    save(output / "results.json", result)
    fields = ("variant", "case", "repetition", "elapsed_seconds", "user_seconds", "system_seconds",
              "max_rss_kib", "stdout_sha256", "pair_output_sha256", "stdout_equal",
              "pair_max_absolute_error", "pair_max_relative_error")
    with (output / "samples.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for case in args.cases:
            for repetition in range(1, args.repetitions + 1):
                variants = ("base", "candidate") if repetition % 2 else ("candidate", "base")
                for variant in variants:
                    prefix = f"{variant}-{case}-{repetition}"
                    command = manifest["commands"][variant][case]
                    stdout_path, stderr_path = output / f"{prefix}.stdout", output / f"{prefix}.stderr"
                    timing_path = output / f"{prefix}.time.tsv"
                    print(f"Running {prefix}", flush=True)
                    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                        subprocess.run([str(timer), "-f", "%e\t%U\t%S\t%M", "-o", str(timing_path), *command],
                                       env=env, stdout=stdout, stderr=stderr, check=True)
                    elapsed, user, system, rss = timing_path.read_text().strip().split("\t")
                    stdout_hash = digest(stdout_path)
                    reference_stdout.setdefault(case, stdout_hash)
                    sample = {"variant": variant, "case": case, "repetition": repetition,
                              "elapsed_seconds": float(elapsed), "user_seconds": float(user),
                              "system_seconds": float(system), "max_rss_kib": int(rss),
                              "stdout_sha256": stdout_hash, "pair_output_sha256": "",
                              "stdout_equal": stdout_hash == reference_stdout[case],
                              "pair_max_absolute_error": 0.0, "pair_max_relative_error": 0.0}
                    if args.pair_output:
                        pair_path = output / f"{prefix}.bpProb.csv"
                        (output / f"{variant}-{case}.bpProb.csv").rename(pair_path)
                        pairs = read_pairs(pair_path)
                        reference_pairs.setdefault(case, pairs)
                        absolute, relative = compare_pairs(reference_pairs[case], pairs)
                        sample.update(pair_output_sha256=digest(pair_path),
                                      pair_max_absolute_error=absolute, pair_max_relative_error=relative)
                    writer.writerow(sample)
                    stream.flush()
                    samples.append(sample)
                    if not sample["stdout_equal"] or sample["pair_max_relative_error"] > RELATIVE_TOLERANCE:
                        save(output / "results.json", {**result, "failed_sample": prefix, "summary": summarize(samples)})
                        raise RuntimeError(f"output comparison failed: {prefix}; see raw outputs and samples.tsv")
    result.update(complete=True, output_checks_passed=True, ended_at=stamp(),
                  samples=len(samples), summary=summarize(samples),
                  pair_max_absolute_error=max(row["pair_max_absolute_error"] for row in samples),
                  pair_max_relative_error=max(row["pair_max_relative_error"] for row in samples))
    save(output / "results.json", result)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
