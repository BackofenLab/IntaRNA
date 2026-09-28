#!/usr/bin/env python3
"""Compare two identically built IntaRNA binaries; save every sample and output.

Requires Python 3 and Linux /usr/bin/time. No optional Python packages.
Run after builds/tests finish, on an otherwise idle machine.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import signal
import statistics
import subprocess
import time


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run_timed(command, env, timeout=180):
    """Bound the complete process group, including GNU time's child process."""
    with subprocess.Popen(command, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, start_new_session=True) as process:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except BaseException as error:
            # The child is in a separate session, so terminal interrupts must
            # also trigger explicit cleanup of the complete process group.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
            if isinstance(error, subprocess.TimeoutExpired):
                raise subprocess.TimeoutExpired(command, timeout, output=stdout,
                                                stderr=stderr) from error
            raise
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repetitions", type=int, default=7)
    parser.add_argument("--cpu", type=int)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("repetitions must be positive")
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binaries = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
    env = dict(os.environ, LC_ALL="C", OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    root = Path(__file__).resolve().parents[2]
    rng = random.Random(246)
    target = "".join(rng.choices("ACGU", k=1000))
    query = "".join(rng.choices("ACGU", k=120))
    inputs = {}
    for name, seq in [("t1000", target), ("t300", target[:300]), ("t80", target[:80]),
                      ("q120", query), ("q80", query[:80])]:
        path = out / (name + ".fa")
        path.write_text(">" + name + "\n" + seq + "\n")
        inputs[name] = path
    for name in ("fhlA", "OxyS"):
        path = out / (name + ".fa")
        path.write_bytes((root / "doc/handson" / (name + ".fasta")).read_bytes())
        inputs[name] = path
    common = ["--threads=1", "--outMode=C", "--outCsvCols=id1,id2,start1,end1,start2,end2,E,bpList",
              "--default-log-file=/dev/null"]
    cases = []

    def case(name, target_name, query_name, *options):
        cases.append({"name": name, "args": ["--target=" + str(inputs[target_name]),
                      "--query=" + str(inputs[query_name]), *common, *options]})

    case("biological-default", "fhlA", "OxyS")
    case("banded-default", "t1000", "q120")
    case("banded-narrow", "t1000", "q120", "--accL=30", "--accW=100")
    case("dense-no-accessibility", "t1000", "q120", "--acc=N")
    case("seed-bulges", "t300", "q120", "--seedMaxUP=2")
    case("helix-block", "t300", "q120", "--model=B")
    case("exact-mfe", "t80", "q80", "--mode=M", "--model=S", "--noSeed", "--intLenMax=30")
    case("exact-ensemble", "t80", "q80", "--mode=M", "--model=P", "--noSeed", "--intLenMax=30")
    case("seed-extension", "t300", "q120", "--mode=H", "--model=P", "--intLenMax=40")
    case("triangular-base-pair", "t300", "q120", "--energy=B")
    metadata = {"platform": platform.platform(), "affinity": sorted(os.sched_getaffinity(0)),
                "repetitions": args.repetitions, "warmups": 1, "seed": 246,
                "binaries": {name: {"path": str(path), "sha256": sha(path.read_bytes()),
                             "version": subprocess.check_output([str(path), "--version"], env=env).decode()}
                             for name, path in binaries.items()},
                "inputs": {name: sha(path.read_bytes()) for name, path in inputs.items()},
                "cases": cases, "environment": {key: env[key] for key in ("LC_ALL", "OPENBLAS_NUM_THREADS", "OMP_DYNAMIC")}}
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    expected = {}
    samples = []
    for repetition in range(args.repetitions + 1):
        schedule = [(case, variant) for case in cases for variant in binaries]
        rng.shuffle(schedule)
        for current, variant in schedule:
            stem = out / f"{current['name']}.{variant}.{repetition}"
            command = ["/usr/bin/time", "-f", "%U %S %M", "-o", str(stem) + ".time",
                       str(binaries[variant]), *current["args"]]
            started = time.perf_counter()
            try:
                result = run_timed(command, env=env)
            except subprocess.TimeoutExpired as error:
                Path(str(stem) + ".stdout").write_bytes(error.output or b"")
                Path(str(stem) + ".stderr").write_bytes(error.stderr or b"")
                raise RuntimeError(f"{stem.name} timed out after {error.timeout}s") from error
            wall = time.perf_counter() - started
            Path(str(stem) + ".stdout").write_bytes(result.stdout)
            Path(str(stem) + ".stderr").write_bytes(result.stderr)
            if result.returncode:
                raise RuntimeError(f"{stem.name} failed: {result.stderr.decode()}")
            digest = sha(result.stdout)
            if digest != expected.setdefault(current["name"], digest):
                raise RuntimeError(f"output mismatch: {stem.name}")
            user, system, rss = map(float, Path(str(stem) + ".time").read_text().split())
            row = {"case": current["name"], "variant": variant, "repetition": repetition,
                   "warmup": repetition == 0, "wall_s": wall, "user_s": user, "system_s": system,
                   "max_rss_kib": int(rss), "sha256": digest}
            samples.append(row)
            with (out / "samples.jsonl").open("a") as handle:
                handle.write(json.dumps(row) + "\n")
        print(f"Completed round {repetition}/{args.repetitions}", flush=True)
    summary = []
    for current in cases:
        row = {"case": current["name"], "output_identical": True}
        for variant in binaries:
            selected = [x for x in samples if not x["warmup"] and x["case"] == current["name"] and x["variant"] == variant]
            row[variant] = {key: statistics.median(x[key] for x in selected) for key in ("wall_s", "user_s", "max_rss_kib")}
            row[variant]["wall_min_s"] = min(x["wall_s"] for x in selected)
            row[variant]["wall_max_s"] = max(x["wall_s"] for x in selected)
        row["wall_change_percent"] = 100 * (row["candidate"]["wall_s"] / row["baseline"]["wall_s"] - 1)
        summary.append(row)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
