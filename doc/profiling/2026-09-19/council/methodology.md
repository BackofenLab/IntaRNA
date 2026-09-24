# Profiling council: experimental rigor

## Initial evidence

- Reviewed source checkout `intarna_optimization_supervised`, HEAD `1b8377989ba79a0f1b8a450163b71bd8ed7c3e45`. Tracked files were clean; `configure~`, `src/bin/IntaRNA`, and `tests/runApiTests` were untracked. Existing Makefiles reference an older external Conda dependency directory; their flags alone do not prove the executable was built from current source.
- Host: AMD Ryzen 5 7530U, six physical cores and twelve hardware threads, one NUMA node, 30 GiB RAM, performance frequency governor, CPUs 0–11 permitted. Ancestor cgroups inspected reported unlimited CPU quota. Initial system load was low, but this is a shared interactive host, not a controlled dedicated benchmark machine.
- `perf`, `gprof`, `strace`, `gdb`, and `bpftrace` are installed. `perf_event_paranoid=4`; hotspot council confirmed even task-clock recording is denied. Valgrind and heaptrack are absent from PATH. Report unavailable hardware-counter/cache/allocation profiling explicitly, without presenting guesses as measurements.

## Recommended protocol

1. Profile a fresh optimized build from a recorded source revision, preferably with debug symbols and frame pointers. Save configure invocation, compiler version, effective flags, dependency versions, executable hash, and linked-library information. Keep source changes separate from profiling artifacts.
2. Record exact commands, workload hashes and sequence lengths, environment, thread count, output settings, start times, exit status, and resource limits. Synthetic inputs require a fixed generation seed; keep real biological examples as the primary anchors.
3. Measure clean runs independently of logging, tracing, and instrumentation. Each case receives one excluded warmup and at least three measured runs, preferably interleaved/randomized across comparisons to reduce temperature and order effects. Retain individual results; summarize median, minimum, and maximum. Three measurements do not support tight confidence intervals or significance claims.
4. Collect high-resolution monotonic wall time plus user/system CPU time and peak RSS. GNU time reports only centisecond wall time, so use Python `perf_counter` or a similar monotonic timer for short cases. Peak RSS is resident memory at process high water, not total allocations or a leak measurement. CPU/wall ratios can indicate aggregate parallel use, not per-core efficiency on their own.
5. Avoid competing CPU-intensive jobs during clean measurements, including compilation and other agents' benchmarks. Rerun a case only for a stated failure, uncertainty, or clearly demonstrated interference. Do not silently delete outliers.
6. Use process timeouts and bounded input sizes for expensive exact/ensemble modes. Report timeouts as censored results with the bound, not as measured completed runtime. Preserve failing stderr and exit codes.

## Workload and comparison controls

- Threads parallelize target sequences first, then query sequences if only one target, then window pairs if there is only one target/query pair. Query accessibility preprocessing also has a parallel loop. Use at least 24–48 independent targets for a 1, 2, 4, 6, and 12 thread sweep. A single pair without multiple windows cannot demonstrate throughput scaling. Six and twelve threads distinguish physical-core scaling from SMT on this host.
- Workload council identified `doc/handson/NC_000913.fa` as 4,319 separate 300 nt records, not a single contiguous genome. A recorded prefix of these targets with ChiX (84 nt) or GcvB (201 nt) is a reproducible throughput case.
- Hold prediction model, mode, seed rule, interaction lengths, accessibility windows, and output fields fixed across ablations. In particular, disabling accessibility removes an implicit interaction-length restriction: both accessibility on/off cases must explicitly set `--intLenMax=150` (or another shared value). Otherwise the result conflates accessibility work with an enlarged prediction search.
- Algorithm modes and window widths can alter answers. Their timing difference is not a drop-in optimization speedup unless outputs are shown equivalent for the tested inputs. Report scientific/accuracy choices alongside resource effects.
- Inspect memory growth using several lengths and a thread sweep. Any fitted exponent is empirical over the measured range; a small grid cannot establish asymptotic complexity or universal scaling.
- For windowed comparisons, choose overlap consistent with the documented interaction-length requirements and preserve output parity checks. Heuristic boundaries and ties may affect predictions even when overlap is valid.

## Correctness and attribution

- Build and run repository API tests and all 20 existing CLI regression parameter cases in the isolated profiling tree. The regression script writes `.testout` files, so do not run it in the user's checkout. Record any expected or unexpected failures rather than rewriting baselines.
- Every timed command must exit successfully and produce nonempty, parseable output appropriate to its mode. Preserve input and output hashes. For thread sweeps normalize/sort CSV data rows before comparison because order may vary; verify the schema/header separately.
- Equal output between an instrumented build and the clean build validates the profiled workload's behavior; it does not prove all possible inputs equivalent.
- Built-in `--verbose=9` timings cover query/target accessibility and prediction phases. Hotspot council found that `predict` includes seed filling, matrix setup, hybrid DP, and reporting. Seed time is nested; never add it to predictor time as a disjoint category. Single-thread logs permit clearer phase attribution.
- A gprof fallback measures instrumented code and perturbs small frequent calls. Shared-library work such as ViennaRNA may be invisible. Label its observed function rankings and percentages with that coverage; do not interpret them as percentages of clean end-to-end runtime. Match source symbols and recognize inlined/optimized functions before recommending changes.

## Required limitations in the final interpretation

This is a broad profile of one recorded build and a bounded workload suite on one laptop-class CPU. It is not exhaustive across every IntaRNA option, dataset, platform, or deployment. Keep clean timing, phase instrumentation, function instrumentation, syscall tracing, and unavailable hardware/allocation metrics distinct. Prioritize optimization opportunities only where measured evidence and source inspection agree.

## Council exchanges

- Challenged throughput design: independent pairs are required; workload council proposed real 300 nt genomic-record batches and biological small RNAs.
- Challenged profiler coverage: hotspot council verified permission restrictions, nested timers, and the need for a fresh executable; clean timing and logging runs remain separate.
- Highlighted accessibility-disable confounding and output-equivalence checks to both workload and coordinating agents.
- Coordinating agent confirmed a fresh isolated archive build with `-O3 -g -fno-omit-frame-pointer`, one excluded warmup, three randomized measured rounds, Python monotonic timing plus GNU time resource statistics, and canonical CSV hashes. Compilation/testing/profiling will finish before clean measurements.
- Hotspot council found working user-space gprofng CPU/heap profiling and detected an independent OpenBLAS worker pool dominating an old executable's tiny-input sample. Set `OPENBLAS_NUM_THREADS=1` explicitly for controlled runs; report that choice. Hold OpenMP dynamic scheduling/affinity environment constant (`OMP_DYNAMIC=FALSE`, `OMP_PLACES=cores`, `OMP_PROC_BIND=close` where supported). Separate CPU and allocation instrumentation where practical because allocation interception can perturb timings. Quantify any default-environment overhead with clean matched commands and output equivalence, not from instrumented CPU share alone.

## Actual-result audit plan

1. Inspect harness commands and input manifests for seed/length/window/algorithm confounds; verify explicit thread and BLAS settings and that the recorded executable is the tested build.
2. Check run counts, one excluded warmup per case, randomized order, return codes, timeouts, stderr, CSV headers, nonempty expected outputs, and consistent output hashes within each case.
3. Recompute medians and ranges from raw records; verify speedup is median one-thread wall time divided by the comparable multithread median and throughput uses the actual prediction count.
4. Check normalized output equality across thread counts, instrumentation versus clean builds, and any claimed semantics-preserving change. Keep algorithm/parameter changes that change outputs labeled separately.
5. Compare CPU/wall ratio, RSS, variation, and timing order for interference or suspicious results; retain anomalies and explain reruns transparently if needed.
6. Audit phase totals for nesting/double counting and sampler tables for own-code versus library coverage, sample count, interval, thread behavior, and profiling overhead. Do not promote weak rankings into precise percentages.
7. Check every final conclusion against its actual data scope. Explicitly list omitted instrumentation dimensions and identify measurements requiring a different host/tool setup.

Actual-result audit is pending until measurements are available.
