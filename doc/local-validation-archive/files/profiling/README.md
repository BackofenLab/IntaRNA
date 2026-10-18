# IntaRNA profiling artifacts

Profiling target: commit `1b8377989ba79a0f1b8a450163b71bd8ed7c3e45` from the main `refactor-to-c++23` checkout. This directory contains an isolated source archive and fresh optimized executable. The working checkout and its pre-existing executables are unchanged.

Read [REPORT.md](REPORT.md) for conclusions, measurements, limitations and optimization priorities. [results/overview.png](results/overview.png) is the visual overview. All timings refer to this machine and the recorded build, not a cross-version comparison.

`intarna-profile-evidence.zip` contains the report, figures, scripts, inputs, text profile outputs, audited timing evidence and C/C++ sources for inspecting findings. Large build artifacts and raw gprofng experiment directories remain in this local folder. The archive includes a SHA-256 manifest of its evidence files.

## Evidence

- `environment.json`: CPU, compiler, linked libraries, revision and executable hash, relevant environment, affinity, governor and cgroup limits.
- `cases.json`: exact benchmark argument arrays, workload categories, thread counts and supplementary environment overrides.
- `inputs/manifest.json`: deterministic input hashes, actual lengths and GC fractions. Biological inputs come from the repository's handson examples. Synthetic generator seed is `20260919`.
- `results/runs.jsonl`: original per-run commands, monotonic elapsed time, GNU time resource metrics, timestamps, exit codes and output hashes. Individual `.stdout`, `.stderr` and `.time` files are preserved.
- `results/validated-runs.jsonl`: derived schema and hash validation after filtering IntaRNA's `# INFO` lines. Original run records are retained unchanged.
- `results/summary.csv` and `.json`: medians, min/max, CPU time, peak RSS and output stability for measured repetitions only. `failures.json` lists unsuccessful/invalid runs.
- `profiles/`: independent CPU sampling, heap tracing and verbose phase measurements, plus original capability probes. `probe*` uses the old executable solely to test tooling; these files do not supply final performance claims.
- `council/`: three agents' planning discussion, methodological and source reviews, and final audits.
- `logs/`: complete build/test logs, benchmark progress, cache-generation cost, and profile orchestration.

## Reproduce

Keep the current evidence intact. Copy `scripts/` into a new writable directory, then run the commands below there. The build requires the recorded conda dependency environment; another toolchain constitutes a different experiment. `INTARNA_DEPS` can point to an equivalent local dependency environment.

```bash
INTARNA_REPO=/home/alex/CodexProjects/IntaRNA/intarna_optimization_supervised bash scripts/build.sh
python3 scripts/benchmark.py prepare
python3 scripts/benchmark.py run
python3 scripts/supplement.py
python3 scripts/structural_checks.py
python3 scripts/profile.py
python3 scripts/profile_audit.py
python3 scripts/phase_summary.py
bash scripts/stream_cleanup_probe.sh
python3 scripts/benchmark.py summarize
python3 scripts/audit.py
python3 scripts/plots.py
```

Run the stages sequentially and avoid concurrent compute-intensive work. The harness refuses to overwrite a run's output. Main benchmarks use one excluded warmup and three sequential, randomized measured rounds. Supplementary checks use a separately seeded randomized order. Compilation, testing, CPU sampling, verbose logs and heap tracing do not overlap clean timing.

The controlled environment sets `LC_ALL=C`, `OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores` and `OPENBLAS_NUM_THREADS=1`. Supplementary BLAS cases explicitly vary the last setting. The default host CPU affinity covers logical CPUs 0–11; it is not an isolated machine. Timings measure complete process latency including startup and output. RSS is process peak resident memory, in MiB (1024² bytes); it is distinct from traced heap allocations.

CPU profiling uses gprofng's nominal 100 ms user-space sampling, with its timer warning preserved and discussed in the report. Hardware counters and syscall tracing were unavailable under the current host permissions. Heap tracing runs independently from the timing suite.

The structural-check stage first runs the two explicit 149/99-cap cache controls, then checks complete base-pair lists. Original cache-generation files must retain their 150/100 columns. The reader's one-column limit reduction and the confirmed input-stream cleanup defect are documented in the report. The isolated stream probe links the unmodified freshly built library; its process exits after reporting descriptor counts.
