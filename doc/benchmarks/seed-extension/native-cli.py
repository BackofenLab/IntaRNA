#!/usr/bin/env python3
"""Measure the native OxyS/fhlA CLI workload from issue #270 with GNU time."""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import statistics
import subprocess
import sys


QUERY = "GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGC"
TARGET = "AGUUAGUCAAUGACCUUUUGCACCGCUUUGCGGUGCUUUCCUGGAACAACAAAAUGUCAUAUACACCGAUGAGUGAUCUCGGACAACAAGGGUUGUUCGACAUCACUCGGACA"


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def capture(command, env):
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    return {"command": command, "exit_code": result.returncode,
            "stdout": result.stdout, "stderr": result.stderr}


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def environment(env):
    cpu_info = Path("/proc/cpuinfo")
    model = next((line.split(":", 1)[1].strip()
                  for line in cpu_info.read_text().splitlines()
                  if line.startswith("model name")), None) if cpu_info.exists() else None
    keys = ("OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
            "LD_LIBRARY_PATH", "LC_ALL", "LANG", "CONDA_PREFIX")
    return {"platform": platform.platform(), "machine": platform.machine(),
            "cpu_model": model, "logical_cpu_count": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "load_average_at_start": os.getloadavg() if hasattr(os, "getloadavg") else None,
            "python": sys.version, "variables": {key: env[key] for key in keys if key in env},
            "working_directory": str(Path.cwd()),
            "timing_scope": "fresh complete CLI process, including accessibility and output",
            "cpu_frequency_and_other_workstation_activity": "uncontrolled"}


def summarize(samples):
    result = []
    for label, case in sorted({(row["variant"], row["case"]) for row in samples}):
        rows = [row for row in samples if (row["variant"], row["case"]) == (label, case)]
        item = {"variant": label, "case": case, "samples": len(rows)}
        for metric in ("elapsed_seconds", "user_seconds", "system_seconds", "max_rss_kib"):
            values = [row[metric] for row in rows]
            item[metric] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
        result.append(item)
    return result


def verify(records):
    """Compare each case only with itself, never seeded with unseeded or M with H."""
    outputs = {}
    checks = []
    for record in records:
        for sample in record["samples"]:
            hashes = (sample["stdout_sha256"], sample.get("pair_output_sha256"))
            case = sample["case"]
            if case not in outputs:
                outputs[case] = hashes
            checks.append({"variant": sample["variant"], "case": case,
                           "repetition": sample["repetition"],
                           "equal_to_first_sample_for_case": hashes == outputs[case]})
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="native IntaRNA executable")
    parser.add_argument("--label", choices=("base", "candidate"), default="base",
                        help="label for --binary; useful for separate candidate runs")
    parser.add_argument("--candidate", type=Path, help="optional second executable, alternated with --binary")
    parser.add_argument("--compare-results", type=Path,
                        help="compare output hashes with an earlier results.json; does not rerun it")
    parser.add_argument("--output-dir", type=Path, required=True, help="new or empty results directory")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--modes", nargs="+", choices=("M", "H"), default=["M", "H"])
    parser.add_argument("--pair-output", action="store_true", help="also request and compare bpProb CSV output")
    parser.add_argument("--time", type=Path, default=Path("/usr/bin/time"), help="GNU time executable")
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    if args.candidate and args.label != "base":
        parser.error("--candidate requires --label base")
    if len(set(args.modes)) != len(args.modes):
        parser.error("--modes must not contain duplicates")

    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("--output-dir must be empty (existing raw results are never overwritten)")
    binaries = {args.label: args.binary.resolve()}
    if args.candidate:
        binaries["candidate"] = args.candidate.resolve()
    for binary in binaries.values():
        if not binary.is_file() or not os.access(binary, os.X_OK):
            parser.error(f"not an executable file: {binary}")
    time_binary = args.time.resolve()
    if not time_binary.is_file() or not os.access(time_binary, os.X_OK):
        parser.error(f"GNU time is not executable: {time_binary}")
    env = dict(os.environ, OMP_NUM_THREADS="1", LC_ALL="C", LANG="C")
    time_version = capture([str(time_binary), "--version"], env)
    if time_version["exit_code"] or "GNU" not in time_version["stdout"]:
        parser.error("--time must identify itself as GNU time")

    workload = {"query_name": "OxyS", "query": QUERY, "target_name": "fhlA", "target": TARGET,
                "model": "P", "modes": sorted(args.modes), "seed_states": ["seeded", "unseeded"],
                "seed_bp": 7, "threads": 1, "pair_output": args.pair_output,
                "other_options": "CLI defaults"}
    previous = None
    if args.compare_results:
        previous = json.loads(args.compare_results.read_text(encoding="utf-8"))
        if previous["workload"] != workload:
            parser.error("saved results have a different workload or options")
        if not previous.get("complete"):
            parser.error("saved results are incomplete")

    output.mkdir(parents=True, exist_ok=True)
    record = {"schema_version": 1, "started_at": timestamp(), "complete": False,
              "workload": workload, "environment": environment(env), "gnu_time": time_version,
              "harness_sha256": digest(Path(__file__).resolve()),
              "invocation": [sys.executable, *sys.argv], "binaries": {}, "samples": []}
    if args.compare_results:
        record["comparison_source"] = {"path": str(args.compare_results.resolve()),
                                       "sha256": digest(args.compare_results)}
    for label, binary in binaries.items():
        record["binaries"][label] = {"path": str(binary), "sha256": digest(binary),
                                     "version": capture([str(binary), "--version"], env)}
    results_file = output / "results.json"
    save(results_file, record)

    for mode in args.modes:
        for seeded in (True, False):
            case = f"P-{mode}-{'seeded' if seeded else 'unseeded'}"
            for repetition in range(1, args.repetitions + 1):
                # Alternate executable order to reduce a consistent order bias.
                labels = list(binaries)
                if repetition % 2 == 0:
                    labels.reverse()
                for label in labels:
                    prefix = f"{label}-{case}-{repetition}"
                    stdout_path = output / f"{prefix}.stdout"
                    stderr_path = output / f"{prefix}.stderr"
                    timing_path = output / f"{prefix}.time.tsv"
                    command = [str(binaries[label]), f"--query={QUERY}", f"--target={TARGET}",
                               "--qId=OxyS", "--tId=fhlA",
                               "--model=P", f"--mode={mode}", "--threads=1", "--seedBP=7",
                               f"--noSeed={'false' if seeded else 'true'}"]
                    pair_path = output / f"{prefix}.bpProb.csv"
                    if args.pair_output:
                        command.append(f"--out=bpProb:{pair_path}")
                    timed_command = [str(time_binary), "-f", "%e\t%U\t%S\t%M", "-o", str(timing_path), *command]
                    sample = {"variant": label, "case": case, "repetition": repetition,
                              "command": command, "shell_command": shlex.join(command),
                              "timed_command": timed_command, "started_at": timestamp(),
                              "stdout": stdout_path.name, "stderr": stderr_path.name,
                              "gnu_time_output": timing_path.name}
                    print(f"Running {prefix}", flush=True)
                    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                        process = subprocess.run(timed_command, env=env, stdout=stdout, stderr=stderr)
                    sample.update(ended_at=timestamp(), exit_code=process.returncode,
                                  stdout_sha256=digest(stdout_path), stderr_sha256=digest(stderr_path))
                    if process.returncode:
                        record["failed_sample"] = sample
                        save(results_file, record)
                        raise RuntimeError(f"{prefix} exited {process.returncode}; see {stderr_path}")
                    elapsed, user, system, rss = timing_path.read_text().strip().split("\t")
                    sample.update(elapsed_seconds=float(elapsed), user_seconds=float(user),
                                  system_seconds=float(system), max_rss_kib=int(rss))
                    if args.pair_output:
                        sample.update(pair_output=pair_path.name, pair_output_sha256=digest(pair_path))
                    record["samples"].append(sample)
                    save(results_file, record)

    record["complete"] = True
    record["ended_at"] = timestamp()
    record["summary"] = summarize(record["samples"])
    record["output_checks"] = verify(([previous] if previous else []) + [record])
    record["output_equality_passed"] = all(check["equal_to_first_sample_for_case"]
                                           for check in record["output_checks"])
    if previous:
        record["comparison_summary"] = summarize(previous["samples"] + record["samples"])
    save(results_file, record)
    print(json.dumps(record["summary"], indent=2))
    print(f"Output equality: {'PASS' if record['output_equality_passed'] else 'FAIL'}; {results_file}")
    return 0 if record["output_equality_passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        sys.exit(str(error))
