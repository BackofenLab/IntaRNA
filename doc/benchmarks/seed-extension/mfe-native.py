#!/usr/bin/env python3
"""Compare MFE seed extensions with identical commands and byte-exact output."""

import argparse
import csv
import importlib.util
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

# Share the recorded biological sequences and summary format with the ensemble
# benchmark; this runner always selects model X and includes traceback and Eall.
spec = importlib.util.spec_from_file_location(
    "native", Path(__file__).with_name("anchored-native.py"))
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)

WORKLOADS = {
    "oxys": (native.QUERY, native.TARGET, []),
    "oxys-noacc": (native.QUERY, native.TARGET, ["--acc=N"]),
    "dense": ("G" * 48, "C" * 48,
              ["--acc=N", "--qIntLenMax=24", "--tIntLenMax=24"]),
    "dense-nolp": ("G" * 48, "C" * 48,
                   ["--acc=N", "--qIntLenMax=24", "--tIntLenMax=24", "--outNoLP=true"]),
    "bulged": (native.QUERY, native.TARGET,
               ["--acc=N", "--seedBP=5", "--seedMaxUP=2"]),
    "no-seeds": ("A" * 109, "C" * 113, ["--acc=N"]),
}
CASES = tuple(f"{name}-{mode}" for name in WORKLOADS for mode in ("M", "H"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    args = parser.parse_args()
    if args.repetitions < 1 or len(set(args.cases)) != len(args.cases):
        parser.error("repetitions must be positive and cases unique")
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", LC_ALL="C")
    binaries = {"base": args.binary.resolve(), "candidate": args.candidate.resolve()}
    manifest = {
        "started_at": native.stamp(), "invocation": sys.argv,
        "platform": platform.platform(), "cpu": next(
            (line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
             if line.startswith("model name")), None),
        "environment": {key: env[key] for key in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "LC_ALL")},
        "timing": "fresh complete process; perf_counter elapsed; GNU time user/system/RSS",
        "workstation_activity_and_cpu_frequency": "uncontrolled",
        "repetitions": args.repetitions, "warmups": 0,
        "binaries": {name: {"path": str(path), "sha256": native.digest(path),
                           "version": subprocess.check_output(
                               [str(path), "--version"], text=True, env=env).strip()}
                     for name, path in binaries.items()},
        "harness_sha256": native.digest(Path(__file__)), "commands": {},
    }
    for variant, binary in binaries.items():
        manifest["commands"][variant] = {}
        for case in args.cases:
            workload, mode = case.rsplit("-", 1)
            query, target, extra = WORKLOADS[workload]
            manifest["commands"][variant][case] = [
                str(binary), f"--query={query}", f"--target={target}",
                "--qId=query", "--tId=target", "--model=X", f"--mode={mode}",
                "--threads=1", "--outMode=C", "--outNumber=5", "--outOverlap=B",
                "--outCsvCols=start1,end1,start2,end2,E,Eall,hybridDP,seedStart1,seedStart2",
                *extra,
            ]
    native.save(output / "manifest.json", manifest)
    samples, reference = [], {}
    fields = ("variant", "case", "repetition", "elapsed_seconds", "user_seconds",
              "system_seconds", "max_rss_kib", "stdout_sha256", "stdout_equal")
    with (output / "samples.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for case in args.cases:
            for repetition in range(1, args.repetitions + 1):
                variants = ("base", "candidate") if repetition % 2 else ("candidate", "base")
                for variant in variants:
                    prefix = output / f"{variant}-{case}-{repetition}"
                    print(prefix.name, flush=True)
                    timing = prefix.with_suffix(".time.tsv")
                    with prefix.with_suffix(".stdout").open("wb") as stdout, \
                            prefix.with_suffix(".stderr").open("wb") as stderr:
                        start = time.perf_counter()
                        subprocess.run(["/usr/bin/time", "-f", "%U\t%S\t%M", "-o", str(timing),
                                        *manifest["commands"][variant][case]],
                                       env=env, stdout=stdout, stderr=stderr, check=True)
                        elapsed = time.perf_counter() - start
                    user, system, rss = timing.read_text().strip().split("\t")
                    digest = native.digest(prefix.with_suffix(".stdout"))
                    reference.setdefault(case, digest)
                    row = dict(variant=variant, case=case, repetition=repetition,
                               elapsed_seconds=elapsed, user_seconds=float(user),
                               system_seconds=float(system), max_rss_kib=int(rss),
                               stdout_sha256=digest, stdout_equal=digest == reference[case])
                    writer.writerow(row)
                    stream.flush()
                    samples.append(row)
                    if not row["stdout_equal"]:
                        raise RuntimeError(f"output differs: {prefix}; inspect retained raw files")
    native.save(output / "results.json", {
        "complete": True, "output_checks_passed": True, "samples": len(samples),
        "ended_at": native.stamp(), "summary": native.summarize(samples),
    })


if __name__ == "__main__":
    main()
