#!/usr/bin/env python3
"""Compare complete CLI runs for issue #262; Linux/GNU time, Python 3.8+."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense", required=True, type=Path)
    parser.add_argument("--hash", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--cpu", type=int, default=0)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE")
    binaries = {key: str(getattr(args, key).resolve()) for key in ("dense", "hash")}
    source = args.source.resolve()
    sparse_t = "A" * 200 + "G" * 8 + "A" * 192
    sparse_q = "A" * 200 + "C" * 8 + "A" * 192
    cases = {
        "OxyS-fhlA": (str(source / "doc/handson/fhlA.fasta"),
                      str(source / "doc/handson/OxyS.fasta"), "C", 4, 30),
        "GcvB-phoB": (str(source / "doc/handson/phoB.fasta"),
                      str(source / "doc/handson/GcvB.fasta"), "C", 7, 60),
        "sparse400": (sparse_t, sparse_q, "N", 4, 30),
        "dense40": ("G" * 40, "C" * 40, "N", 4, 30),
    }
    samples, verified = [], []
    for name, (target, query, acc, seed, span) in cases.items():
        for pairs in (False, True):
            expected = None
            # One warmup per variant; alternate which variant runs first thereafter.
            for repeat in range(-1, args.repeats):
                order = ("dense", "hash") if repeat % 2 else ("hash", "dense")
                for variant in order:
                    stem = output / f"{name}-{int(pairs)}-{variant}"
                    pair_file = stem.with_suffix(".pairs.csv")
                    rss_file = stem.with_suffix(".rss")
                    command = [binaries[variant], "--target=" + target, "--query=" + query,
                               "--model=P", "--mode=M", "--seedBP=" + str(seed),
                               "--intLenMax=" + str(span),
                               "--qIntLoopMax=3", "--tIntLoopMax=3", "--threads=1",
                               "--outNoLP=false", "--outNoGUend=false", "--outNumber=0",
                               "--outMode=E", "--acc=" + acc,
                               "--default-log-file=/dev/null"]
                    if pairs:
                        command.append("--out=bpProb:" + str(pair_file))
                    wrapped = ["taskset", "-c", str(args.cpu), "/usr/bin/time", "-f", "%M",
                               "-o", str(rss_file)] + command
                    start = time.perf_counter()
                    result = subprocess.run(wrapped, env=env, check=True, capture_output=True)
                    elapsed = time.perf_counter() - start
                    digest = hashlib.sha256(result.stdout).hexdigest()
                    pair_digest = hashlib.sha256(pair_file.read_bytes()).hexdigest() if pairs else ""
                    actual = (digest, pair_digest)
                    if expected is None:
                        expected = actual
                    if actual != expected:
                        raise RuntimeError(f"Output mismatch: {name}, pairs={pairs}, {variant}")
                    stem.with_suffix(".stdout").write_bytes(result.stdout)
                    stem.with_suffix(".stderr").write_bytes(result.stderr)
                    if repeat >= 0:
                        samples.append(dict(case=name, pairs=int(pairs), variant=variant,
                                            repeat=repeat, wall_s=elapsed,
                                            rss_kib=int(rss_file.read_text().strip())))
            verified.append(dict(case=name, pairs=int(pairs), stdout_sha256=expected[0],
                                 pairs_sha256=expected[1]))
            print(f"verified {name} pairs={int(pairs)}", flush=True)
    with (output / "native-samples.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)
    summary = []
    for item in verified:
        row = dict(item)
        for variant in binaries:
            group = [s for s in samples if s["case"] == item["case"]
                     and s["pairs"] == item["pairs"] and s["variant"] == variant]
            row[variant + "_median_s"] = statistics.median(s["wall_s"] for s in group)
            row[variant + "_min_s"] = min(s["wall_s"] for s in group)
            row[variant + "_max_s"] = max(s["wall_s"] for s in group)
            row[variant + "_median_rss_kib"] = statistics.median(s["rss_kib"] for s in group)
        row["hash_time_ratio"] = row["hash_median_s"] / row["dense_median_s"]
        summary.append(row)
    (output / "native-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "native-config.json").write_text(json.dumps(
        dict(cpu=args.cpu, repeats=args.repeats, binaries=binaries,
             source=str(source), cases=cases), indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
