# Final methodology review: clean timing and correctness evidence

The independent evidence audit passed after all clean timing and structural checks finished. CPU/heap profiling was still running when this review was written; those profiles are outside this review's completed scope.

## Audit outcome

- **56 cases and 224 run records:** each case has one excluded warmup and three measured repetitions. All commands completed successfully, with no timeout and no nonempty stderr.
- All preserved stdout hashes, input/cache manifest hashes, per-run executable hashes, GNU time resource fields, and command definitions match their recorded provenance.
- Independently recomputed medians, ranges, CPU statistics, RSS statistics, and output hashes match both JSON and CSV summaries. There are no missing or duplicate records.
- Every case has stable output across its warmup and measured repetitions. Output equality holds across all thread counts, the tested window comparisons, BLAS environment variants, and the matched accessibility-cache pair.
- Additional output containing complete `bpList` values is identical between 1 and 12 threads and between the matched computed/cached accessibility cases. Other case equality checks use the recorded seven output columns (IDs, coordinates, energy), not complete structures.
- Exactly four raw records needed a parsing correction: `algorithm_no_seed` printed one `# INFO` line in its warmup and each repetition. The independent audit verifies the corrected row count, CSV validity, and canonical hash against preserved stdout. Original raw records and raw-output hashes remain unchanged; `validated-runs.jsonl` is the corrected derived view.
- Repository test logs report **4,281 assertions across 36 API test cases**, plus the CLI regression script passing all 20 parameter-file cases. Automake reports both test executables/scripts passed with zero failures.

Machine-readable audit: `results/audit.json`; readable coverage and scatter tables: `results/audit.md`. Reproduce with `python3 profiling/scripts/audit.py` after data collection has stopped.

## Findings supported by the clean data

### Parallel throughput

The workload is 48 real 300 nt target records against ChiX (84 nt). Forty-six interactions are reported, but the throughput denominator is all 48 attempted target–query predictions.

| Threads | Median seconds | Speedup over 1 thread | Predictions/s | Median RSS MiB |
|---:|---:|---:|---:|---:|
| 1 | 4.944 | 1.00× | 9.71 | 18.71 |
| 2 | 2.556 | 1.93× | 18.78 | 22.98 |
| 4 | 1.457 | 3.39× | 32.94 | 31.36 |
| 6 | 1.104 | 4.48× | 43.50 | 39.95 |
| 12 | 0.974 | 5.08× | 49.29 | 65.25 |

Six threads reach 74.7% speedup efficiency; twelve reach 42.3%. On this six-core/twelve-thread CPU, moving from six to twelve workers improves throughput by about 13.3% while increasing RSS by about 63%. These are batch-throughput results, not evidence that a single unwindowed prediction gains the same speedup.

### Reusing accessibility calculations

The valid matched comparison explicitly caps target/query interaction lengths at **149/99**. Computed accessibility takes **0.805 s** median, while loading saved ED values takes **0.492 s**: **1.64× speedup**, or **39.0% lower latency**. Median RSS falls from **19.89 to 14.63 MiB**. Coordinates, energy, and base-pair lists match for this pair.

The original computed-versus-cache comparison must remain observational. The ED reader subtracts one from the available header length for dangling-end treatment (`AccessibilityFromStream.cpp:90`), so its original effective caps were 150/100 versus 149/99. Equal output on this input does not remove that search-space confound. The later matched pair corrects it.

One-time ED generation took 0.847 s and included a complete prediction plus serialization. This cost is excluded from the cached-call latency. The experiment reuses both input sequences' accessibility files across separate process invocations; it does not show the same saving for every existing batch workflow. IntaRNA already reuses per-sequence accessibility within its target/query loops.

### Windowing trades memory for repeated computation

For the 3,000 × 800 nt pair with model S, a 300 nt prediction window reduces median RSS from **127.50 to 24.14 MiB** (81.1% lower, 5.28× smaller), while runtime increases from **21.38 to 59.71 s** (2.79×).

For default model X on the same pair, memory is already about **24.2 MiB**. Windowing gives no meaningful RSS reduction and increases runtime from **41.34 to 71.17 s** (1.72×). Windowing is therefore workload- and model-dependent; it is not a blanket performance recommendation.

### Threaded BLAS startup is a separate environmental cost

The matched biological control with `OPENBLAS_NUM_THREADS=1` takes **0.0522 s** median. Removing that environment setting gives **0.1472 s** median; explicitly requesting twelve BLAS threads gives **0.0861 s**. Median total CPU time rises from approximately **0.04 s** to **1.48 s** and **0.89 s**, respectively, despite identical output.

These BLAS cases have much larger scatter than the controlled suite: their full measured ranges are **45.0%** and **56.6%** of the median. Report their individual ranges and the qualitative overhead rather than presenting either ratio as a precise universal speedup. The controlled suite explicitly limits BLAS to one thread while allowing IntaRNA's requested OpenMP parallelism.

### Algorithm choices materially change cost

For the selected 100 × 100 nt pair, exact model X takes **0.187 s** versus its heuristic default's **0.0436 s** (4.29× longer). Within model S, exact mode takes **6.42 s** versus heuristic mode's **0.0955 s** (67.2× longer). The recorded seven output columns match for this input, but one input does not establish equivalence between exact and heuristic algorithms.

## Statistical and scope limits

- Across the 54 cases excluding the two deliberate BLAS multithreading probes, the largest full measured range is **4.33% of the median**. Most main comparisons have a difference much larger than observed scatter. Three repetitions still support descriptive statistics only, not tight confidence intervals or universal claims.
- The CPU is one laptop-class AMD Ryzen 5 7530U. Power management, SMT, affinity, libraries, compiler, input sequences, and model choices all affect applicability elsewhere. The source revision and linked dependencies are fixed in the environment record.
- Timing measures complete CLI invocations, including startup, dynamic loading, input parsing, and output. Tiny cases are not pure algorithm-kernel benchmarks. GNU time CPU fields have centisecond granularity; Python monotonic wall timing avoids that granularity for latency.
- Peak RSS measures resident high water, not total allocations, retained heap, or leaks. Allocation and sampling claims require the separate profiler evidence and its own limitations.
- Disabling accessibility changes energies and predictor behavior. The accessibility-off case is actually slower on the measured 1,000 × 100 pair, so computed-minus-disabled runtime cannot estimate accessibility cost. Use phase instrumentation for attribution.
- The 75 nt accessibility comparison changes both folding-window width and base-pair span. Sequence-composition probes contain one pair each and do not establish population averages. Nested-prefix size sweeps describe this measured range rather than proving asymptotic complexity.
- Exact base-pair equality is established only for the explicit structural checks. Output equality for selected cases does not establish correctness for every possible input, algorithm, or parameter combination.
- Hardware performance counters remain unavailable under current perf permissions. No instruction, cache-miss, branch-miss, or hardware-stall conclusions should be presented as measured data.

The clean timing evidence is internally consistent and supports the bounded conclusions above. No reruns are required by this audit.
