# Seed-extension kernel comparison

This Linux benchmark compares a selected seeded-kernel revision with the current
source checkout. The current shared-anchor comparison uses baseline
`3866c62d10846fab1adc38b937d0c9b96cc06e51`.
It runs the exact ensemble by default; `--heuristic` applies frozen endpoint
choices to the same anchored forward and outside implementation. True left
extensions retain the no-right endpoint and the selected right endpoint;
interactions beginning at the anchor retain every right endpoint.
It first checks every reported boundary hybrid weight and boundary coefficient,
the total partition, and every raw pair numerator. It then times partition-only
and pair-output calculations in separate processes, alternating baseline and
candidate order. Three repetitions are the default. Each timed invocation also
counts ordinary C++ allocation calls and cumulative requested bytes during
`compute()`, including callbacks but excluding seed construction and serialization.
The benchmark replaces ordinary `new`/`new[]`; all measured types use fundamental
alignment. These counters measure allocation traffic, not peak live memory.

Configure the source tree first so `src/IntaRNA/intarna_config.h` is available.
Requirements are a C++23 GCC toolchain, the project's headers and Boost headers,
Python 3, and GNU `/usr/bin/time`. Run the commands below from the repository
root.

For dependencies outside the compiler's search paths, add
`--include /path/to/dependencies/include` to the build command. An optional
`--runtime-library-dir /path/to/compiler/lib` embeds a runtime library search
path. Both options can be repeated. `--baseline` selects another Git revision
when a different comparison is intended; its historical default remains
`985ba44bdc16e92162e954b7baffa3f3923d4776`, so specify the intended baseline
explicitly. The current driver requires the first-acceptance kernel API
introduced in #266. To measure the heuristic ensemble, use `--heuristic` on
the runner. Use separate `--output-dir` directories for exact and
heuristic results, since the runner replaces its output files. The same built
executables can be used for both ensembles. `--verify-only` compares all dumps
without collecting timed samples, which permits checking correctness before a
quiet period reserved for timing.

The build script snapshots both versions of `SeededPartitionFunction` and
`StackSeedDomain`, including their headers; other project headers and configured
types are shared. It records the compiler, exact commands, baseline commit,
and SHA-256 source digests in `manifest.json`, including the harness and shared
numerical/configuration headers. The runner writes full verification
dumps, `samples.tsv`, and `results.json` in the build directory. Timed runs do not
write the full boundary dumps.

## Shared anchored heuristic comparison (#270)

Build both kernels against baseline `3866c62d10846fab1adc38b937d0c9b96cc06e51`.
Verify the exact and heuristic full dumps before collecting timings:

```sh
python3 doc/benchmarks/seed-extension/build.py \
  --source "$PWD" --output /tmp/intarna-anchored-benchmark --cxx g++-14 \
  --baseline 3866c62d10846fab1adc38b937d0c9b96cc06e51
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-anchored-benchmark \
  --output-dir /tmp/intarna-anchored-exact-check --verify-only
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-anchored-benchmark \
  --output-dir /tmp/intarna-anchored-heuristic-check --verify-only --heuristic
```

After both checks pass, reserve a quiet period without concurrent builds or
tests. Collect five alternating baseline/candidate samples for each fixture,
including partition-only and pair-output calculations:

```sh
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-anchored-benchmark \
  --output-dir /tmp/intarna-anchored-exact --repetitions 5
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-anchored-benchmark \
  --output-dir /tmp/intarna-anchored-heuristic --repetitions 5 --heuristic
```

The native CLI harness complements those synthetic fixtures with the OxyS/fhlA
reproduction from #270. Build complete baseline and candidate executables with
the same compiler, dependencies and release options. The default four cases
combine `--model=P` with modes M/H and seeded/unseeded prediction, using one
thread and five alternating samples per executable:

```sh
python3 doc/benchmarks/seed-extension/anchored-native.py \
  --binary /path/to/baseline/src/bin/IntaRNA \
  --candidate /path/to/candidate/src/bin/IntaRNA \
  --output-dir /tmp/intarna-anchored-native --repetitions 5
```

Use a separate invocation for seeded H/M pair output, without timing the
unrelated unseeded pair cases:

```sh
python3 doc/benchmarks/seed-extension/anchored-native.py \
  --binary /path/to/baseline/src/bin/IntaRNA \
  --candidate /path/to/candidate/src/bin/IntaRNA \
  --output-dir /tmp/intarna-anchored-native-pairs --repetitions 5 \
  --cases H-seeded M-seeded --pair-output
```

Native result directories must be new or empty. The harness requires identical
stdout for each case across both binaries and all repetitions. With
`--pair-output`, it also compares every finite CSV probability at relative
tolerance `2e-12`, with no absolute tolerance; matrix labels and `NA` cells must
match exactly. This permits harmless changes in floating-point summation order
without ignoring differences in small probabilities.

`manifest.json` records executable hashes, versions, environment, and one command
per case/variant. `samples.tsv` records individual timings, maximum process RSS,
output hashes and comparison errors; `results.json` contains the summary and
verification status. Raw stdout, stderr, GNU time output and requested pair CSVs
are retained for every sample. Native timings cover fresh complete CLI processes,
including thermodynamic/accessibility computation and output; they have no
warmup exclusions. They measure a different scope from the isolated kernel
runner. Record the resulting measurements and their interpretation in
[the anchored heuristic runtime analysis](../../analysis/heuristic-anchored-runtime.md).

## Workloads and measurement limits

The dense cases admit every vertex on a 40 by 40 grid. The banded cases use a
60 by 60 grid and admit only vertices with `abs(i-j) <= 3`. All cases use
inclusive interaction spans of 24 on each strand, loop limits of three unpaired
bases per strand, and stack seeds of length four. The many-seed cases admit all
contained seeds on valid diagonals: 1,369 dense or 387 banded occurrences. The
few-seed cases have three diagonal seeds starting at `n/4`, `n/2`, and `3n/4`.
`noLP` is disabled; separate oracle tests validate that restriction. Heuristic
runs select `(min(n-1,q.first+2+q.second%3), min(n-1,q.second+2+q.first%3))` for
anchor end `q`, creating both shared and distinct retained endpoints. Selection
is frozen and identical in both versions.

`driver.cpp` declares deterministic coordinate-dependent initiation, transition,
and full-boundary weights. It measures the complete kernel invocation, including
a small boundary callback. Seed construction, native thermodynamic/accessibility
computation, and serialization are outside this timing. These fixtures exercise
the recurrence, rather than representing complete biological workloads.

## Historical factorization measurements

The `recorded-*` files retain the original exact-kernel comparison between
`cb1c9d419f96c7d2a10bae864d922dcb0234f4d4` and `392357c`, before workspace reuse.
Those runs used the original driver saved by its source hash, without allocation
instrumentation. Their timing does not measure this allocation revision.

Those recorded runs used GCC 14.4.0, `-std=c++23 -O3 -DNDEBUG`, one thread, and an
AMD Ryzen 5 7530U on 2026-10-14. Samples were not CPU-pinned and other build/test
work could run on the same workstation. Treat small differences as noise.
`recorded-results.json` contains correctness differences and timing medians with
minimum/maximum values; `recorded-samples.tsv` contains the individual samples;
`recorded-source-hashes.json` identifies the measured source snapshots.

GNU time records whole-process maximum RSS, including runtime libraries and
allocator overhead. These small fixtures do not isolate the recurrence's peak
allocation. Timings and RSS are local observations, not asymptotic bounds.
See the [algorithm and measured tradeoff](../../analysis/seed-extension-partition.md).

## Workspace reuse measurements (PR #266 review)

On 2026-10-15, five alternating fresh-process samples per fixture/mode used
the same GCC 14.4.0 flags and workstation as above. No project build or test
was running during these samples; processes were not CPU-pinned and other
workstation activity was uncontrolled. Both variants include the allocation
counters. `reuse-exact-*` and `reuse-heuristic-*` contain results and samples;
`reuse-source-hashes.json` identifies the measured sources.

All eight full dumps are byte-identical: 599,242 boundary/partition/pair entries.
The allocation counts below exclude dump collection (timed runs pass `-`).
The remaining calls include the frozen selector map, eligible-start traversal
and initial workspace/result allocation.

| Ensemble | Fixture | Partition allocations, base → reuse | Pair allocations, base → reuse |
| --- | --- | ---: | ---: |
| exact | dense40-many | 23,025 → 21 | 28,502 → 27 |
| exact | dense40-few | 7,376 → 13 | 10,217 → 19 |
| exact | band60-many | 6,384 → 19 | 7,933 → 25 |
| exact | band60-few | 3,139 → 12 | 4,374 → 18 |
| heuristic | dense40-many | 418,631 → 1,399 | 554,966 → 1,405 |
| heuristic | dense40-few | 8,774 → 18 | 12,167 → 24 |
| heuristic | band60-many | 84,243 → 413 | 114,136 → 419 |
| heuristic | band60-few | 4,080 → 17 | 5,658 → 23 |

Median candidate/base runtime ratios range from 0.943 to 1.058;
median process RSS is 3,840–4,352 KiB. These measurements show the
allocation reduction, with broadly similar runtime and process RSS on these
small fixtures. They do not establish a general speedup or a peak-memory
reduction. For example, dense heuristic pair output reduces cumulative requested
bytes from 817,421,572 to 281,640, while retaining maximum-sized workspaces.

To reproduce the recorded sample count, pass `--repetitions 5` to the runner.
Run once without `--heuristic` and once with it, using separate build directories.

## Left-most exact-anchor comparison (#268)

On 2026-10-16, the exact anchor rewrite was compared with
`985ba44bdc16e92162e954b7baffa3f3923d4776`, using GCC 14.4.0,
`-std=c++23 -O3 -DNDEBUG`, one thread, and the same AMD Ryzen 5 7530U.
Five alternating fresh-process samples were collected per fixture/output mode.
Exact and heuristic measurements ran sequentially; no project build, test or
installation ran during the timed samples. CPU affinity, frequency and other
workstation activity were uncontrolled. `leftmost-*-results.json` and
`leftmost-*-samples.tsv` retain all timings, allocation counts, requested bytes
and process RSS; `leftmost-source-hashes.json` records the source, executable
and full-dump hashes, compiler, commands and measurement environment.

All 328,236 exact boundary/partition/pair entries agree within a maximum
relative difference of `5.96e-15` (the comparison tolerance is `1e-12`).
All 271,006 heuristic entries are identical. Each boundary dump includes both
the complete hybrid weight and the full-boundary coefficient; the comparison
also rejects duplicate boundary callbacks and nonfinite values.

The exact kernel uses at most 256 left boundaries per batch. Each reusable
boundary table targets at most 2 MiB, with one row as the minimum when a single
span rectangle exceeds that budget. All four fixtures use 256-row batches.
The table reports median kernel time in seconds; ratios are candidate/base.

| Exact fixture | Partition seconds, base → candidate | Ratio | Pair seconds, base → candidate | Ratio |
| --- | ---: | ---: | ---: | ---: |
| dense40-many | 0.093366 → 0.304718 | 3.264 | 0.238276 → 1.219707 | 5.119 |
| dense40-few | 0.055158 → 0.004673 | 0.085 | 0.146767 → 0.008607 | 0.059 |
| band60-many | 0.015683 → 0.047910 | 3.055 | 0.032126 → 0.176494 | 5.494 |
| band60-few | 0.012318 → 0.002428 | 0.197 | 0.026294 → 0.004654 | 0.177 |

The few-seed fixtures improve, while both many-seed fixtures regress
substantially. Recomputing bounded anchor tables and contracting many anchors
against complete boundaries remains expensive. These measurements do not
establish a general speedup.

| Exact fixture | Partition allocations, base → candidate | Pair allocations, base → candidate |
| --- | ---: | ---: |
| dense40-many | 21 → 1,393 | 27 → 1,398 |
| dense40-few | 13 → 19 | 19 → 24 |
| band60-many | 19 → 409 | 25 → 414 |
| band60-few | 12 → 18 | 18 → 23 |

Exact allocation traffic also increases: the anchor index allocates nodes per
retained occurrence, and the larger reusable boundary tables require more
storage. Partition-only requested bytes rise from 65,368–98,072 to
1,261,507–1,404,889; pair-output requested bytes rise from 124,312–156,952 to
2,518,467–2,661,849. Median whole-process RSS rises from 3,840–4,224 KiB to
4,480–6,144 KiB. These figures measure allocation traffic and process RSS,
respectively; neither is a bound on live kernel memory.

The unchanged heuristic ensemble retains identical allocation counts and
requested bytes, but its measured median runtime ratios are 1.046–1.225.
The largest observed regression is dense40-few partition-only
(0.048834 → 0.059830 seconds). The complete sample ranges are recorded rather
than assuming that an unchanged recurrence has unchanged performance.

To reproduce the final comparison after configuring the checkout:

```sh
python3 doc/benchmarks/seed-extension/build.py \
  --source "$PWD" --output /tmp/intarna-leftmost-benchmark --cxx g++-14 \
  --baseline 985ba44bdc16e92162e954b7baffa3f3923d4776
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-leftmost-benchmark \
  --output-dir /tmp/intarna-leftmost-exact --repetitions 5
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-leftmost-benchmark \
  --output-dir /tmp/intarna-leftmost-heuristic --repetitions 5 --heuristic
```

Add the dependency include/runtime options described above when necessary.
The recorded manifest contains the exact compiler and dependency paths used
for these measurements.

## MFE extensions (#274)

[`mfe-native.py`](mfe-native.py) compares native `--model=X --mode=M/H`
executables with byte-exact energy, traceback, seed and `Eall` output. It includes
OxyS/fhlA, dense, noLP, bulged-seed and no-seed cases. See the
[analysis and reproduction instructions](../../analysis/mfe-seed-extension-reuse.md)
and the `mfe-native-*`/`mfe-source-hashes.json` evidence files. Those measurements
describe the historical cached prototype at `4994edf`, before the PR #276
review removed local pairing and loop-energy caches; they do not measure the
current direct-evaluation implementation.
