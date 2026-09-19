# Hotspot council: method and source hypotheses

## Agreed measurement approach

- Root owns clean timing/build execution; other agents avoid competing benchmarks.
- Profile a fresh isolated optimized build (`-O3 -g -fno-omit-frame-pointer`) of a recorded revision. The preexisting executable was older than some source files and lacked debug sections, so probes below validate tools only.
- Use clean GNU time/Python runs for performance results. Run phase logging, sampled CPU, and heap tracing separately. Track output equivalence where options preserve semantics.
- Set `OPENBLAS_NUM_THREADS=1`: the linked GSL/OpenBLAS dependency starts independent workers even with IntaRNA `--threads=1`. A tool probe of a tiny input without this setting recorded most CPU samples in `blas_thread_server`. Quantify real impact only with separate clean before/after timings.

## Available tools and validated commands

`perf` is installed but `perf_event_paranoid=4` blocks even `perf stat -e task-clock true`. `strace` fails `PTRACE_TRACEME: Operation not permitted`. Valgrind/heaptrack are absent from PATH. gprofng 2.42 successfully collects user-space CPU samples and heap allocations.

Example CPU collection (use a longer representative workload and a fresh unique experiment path):

```sh
OPENBLAS_NUM_THREADS=1 gprofng collect app -p 100 -o cpu.er ./IntaRNA --target=targets.fa --query=query.fa --threads=1 --outMode=C

gprofng display text -metrics e.%totalcpu:i.%totalcpu:name -sort e.totalcpu -limit 40 -functions cpu.er
gprofng display text -metrics e.%totalcpu:i.%totalcpu:name -sort i.totalcpu -limit 40 -functions cpu.er
gprofng display text -metrics e.%totalcpu:i.%totalcpu:name -limit 40 -lines cpu.er
gprofng display text -source 'FUNCTION_NAME' cpu.er
gprofng display text -header -overview -statistics cpu.er
```

Heap collection and extraction:

```sh
OPENBLAS_NUM_THREADS=1 gprofng collect app -p off -H on -o heap.er ./IntaRNA --target=targets.fa --query=query.fa --threads=1 --outMode=C
gprofng display text -heapstat -limit 30 -allocs -leaks heap.er
gprofng display text -metrics i.heapallocbytes:i.heapalloccnt:name -sort i.heapallocbytes -limit 40 -functions heap.er
```

Heap's `bytes leaked` means allocations outstanding when profiling ends; it is not proof of an application leak. Inspect callers, runtime caches and collector frames. Peak tracked heap is different from RSS.

Do not use `-i on`: it immediately fails here with `iotrace_init COL_ERROR_IOINIT llseek`.

### Sampling limitation discovered by calibration probes

All tested gprofng CPU experiments warn `Collection interval timer period was changed (... -> 0); profile data may be unreliable`. In old-binary phoB/GcvB probes lasting about 0.9 seconds CPU, both `-p hi` (997 us) and `-p on` (10007 us) collected only about 9 samples, reporting 0.009 and 0.090 CPU seconds respectively. `-p 100` collected 8 samples and reported 0.800 seconds for a 0.868-second experiment. Effective delivery appears roughly 100 ms in this execution environment. Therefore use 100 ms nominal period, a sufficiently long workload, and disclose the warning. Treat proportions as coarse sampled evidence and cross-check against clean CPU time and built-in phases. Do not interpret sub-percent differences or old probe percentages as final results.

Probe artifacts are under `profiling/profiles/probe*`. Some intentional unsuccessful capability probes have no useful profile data (invalid CLI forms or I/O failure).

## Built-in phase measurement

Use plain `--verbose` (NOT `--verbose=9` or `--v=9`, which Boost rejects). EasyLogging++ maps the flag to maximum verbosity 9 and enables `TIMED_FUNC_IF`.

- `AccessibilityVrna::fillByRNAplfold` covers fold-compound setup, constraints, `vrna_probs_window`, and RAII destruction. It excludes allocation of `edValues` in the AccessibilityVrna constructor initializer.
- Predictor `predict` time includes seed filling, matrix resizing, dynamic programming, traceback/reporting.
- Seed handler `fillSeed` is nested inside predictor time. Subtract it only as an approximate exclusive remainder; do not add it to predictor time again.
- Logs preceding each accessibility call identify query versus target. Use one thread for phase attribution; concurrent phase totals are not elapsed wall time.
- ELPP phase clocks are not an independent clean wall-clock benchmark. Short phases have millisecond resolution; zero milliseconds means below resolution, not no work. `DateTime::formatTime` (`src/easylogging++.cc:1195`) integer-divides durations at 1900 ms and above into whole seconds, and subsequently into whole minutes/hours. Long phase logs are therefore truncated coarse measurements, not millisecond-precise data. Do not derive precise phase percentages or narrow differences from these strings.

Old-binary illustrative probe (not final results): phoB299/GcvB201 gave query accessibility 62 ms, target accessibility 90 ms, predictor 848 ms inclusive, seed 0 ms (below resolution).

## Distinct code paths to exercise

1. Default X/H seed extension (`PredictorMfe2dHeuristicSeedExtension`). Main cost depends on number and location of valid seeds as well as sequence lengths and interaction/loop bounds.
2. X/M exact seed extension: computes left and right extension matrices, then combines four boundary loops. Use small inputs and bound interaction length.
3. S/H seed-aware complete matrix DP (`PredictorMfe2dHeuristicSeed`). Two matrices of `BestInteractionE` containing energy and right endpoints. Different memory behavior than X/H extension matrices of integer energies.
4. `--noSeed` switches X to S and removes seed gating; not a pure subtraction of seed time.
5. P/H or P/M ensemble path, B/H helix-block path, `--seedMaxUP=2` general seed DP, and seed-only mode each cover distinct algorithms.
6. `--acc=N` changes semantics and by default can remove an implicit interaction-length restriction. Explicitly fix `--intLenMax` for fair ablation.
7. Target batches exercise OpenMP. The main program parallelizes target loop if multiple targets, otherwise query loop if multiple queries, otherwise windows; ordinary single-pair DP is not internally parallelized. Query accessibility is precomputed/reused, target accessibility is reused across queries.

## Source hypotheses to test against profiles

- `src/IntaRNA/PredictorMfe2dHeuristicSeedExtension.cpp:71`: iterates valid seeds; resizes left/right matrices and repeats extension DP per seed. Seed count and extension dimensions can make same-length sequences dramatically different.
- `src/IntaRNA/PredictorMfe2dHeuristicSeedExtension.cpp:112` and `:209`: right/left extension loops visit complementary cells and inspect two internal-loop dimensions. Larger `--tLoopMax`/`--qLoopMax` increases candidate predecessor work.
- `src/IntaRNA/InteractionEnergyVrna.h:381`: `getE_interLeft` repeatedly checks validity, accesses sequence codes and calls inline ViennaRNA `E_IntLoop`; expected self/inlined hotspot in prediction-dominated runs. Optimize only after attribution establishes impact.
- `src/IntaRNA/PredictorMfe2dSeedExtension.cpp:105`: exact mode's four boundary-combination loops and repeated `updateOptima` may dominate even when extension DP is similar.
- `src/IntaRNA/AccessibilityVrna.cpp:241`: `vrna_probs_window` expected to dominate accessibility-specific cases; local folding window and span control work independently from pair DP.
- `src/bin/IntaRNA.cpp:166`: choice of target/query/window parallelism predicts poor speedup for a single pair without windows; batch throughput and memory per concurrent worker must be measured separately.
- `src/IntaRNA/SeedHandlerNoBulge.cpp`: no-bulge seed search uses diagonal traversal, skips impossible seeds and a ring buffer; default seed search may be cheap despite seed-extension DP being expensive. General bulged seeds follow a different DP.

These are hypotheses and mechanistic explanations from source, not measured final hotspot claims.

## Additional hypotheses prepared for final-profile audit

- X/H extension matrices resize at `PredictorMfe2dHeuristicSeedExtension.cpp:94` and `:100` for each seed. In the exact linked Boost headers, `matrix.hpp:217` defaults `preserve=true`, allocates a temporary and copies overlapping elements. The DP fill methods initialize every visited matrix cell. If heap data confirms large allocation churn, nonpreserving resize or retained-capacity storage is a targeted candidate; assess correctness and measured benefit before recommending implementation. This may reduce allocation traffic without materially changing CPU bottlenecks.
- No-seed P/M (`PredictorMfeEns2d`) loops over possible right boundaries, fills a left-boundary matrix for each, and computes Boltzmann weights inside internal-loop recurrences. Anticipate energy evaluation and exponential/Boltzmann conversion alongside matrix DP; inclusive rows overlap and cannot be added.

## Final-profile audit checklist

1. All reference/profile commands exit successfully; canonical prediction output matches.
2. Requested 100 ms period, actual sample count, collector warning, clean CPU time and profiled CPU estimate agree in scale. Report coarse CPU shares with enough samples; warn if low counts.
3. Separate exclusive self CPU from inclusive call-subtree CPU; do not sum overlapping parent and child rows. Use source to explain the algorithm path represented by each scenario.
4. Attribute libm/ViennaRNA functions as library work, accounting for inlining into application functions. Missing source lines or optimized/inlined boundaries limit fine attribution.
5. Heap allocated bytes are cumulative churn; peak tracked heap is concurrent tracked allocations; clean peak RSS is process residency. These are different metrics.
6. Outstanding heap allocations at process exit require caller inspection; collector/runtime buffers are not proven IntaRNA leaks.
7. Long ELPP phase strings truncate to whole seconds/minutes; present approximate phase ranges or short-case figures, never false precision.
