# Seed-extension kernel comparison

This Linux benchmark compares the seeded kernel at PR #266 revision
`392357c00f0314659cd4c665446f3f6185f4814a` with the current source checkout.
It runs the exact ensemble by default; `--heuristic` uses frozen right-endpoint
groups to exercise their forward and outside workspace reuse.
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
Python 3, and GNU `/usr/bin/time`. Run from the repository root:

```sh
python3 doc/benchmarks/seed-extension/build.py \
  --source "$PWD" --output /tmp/intarna-seed-extension-benchmark --cxx g++-14
python3 doc/benchmarks/seed-extension/run.py \
  --build-dir /tmp/intarna-seed-extension-benchmark --repetitions 3
```

For dependencies outside the compiler's search paths, add
`--include /path/to/dependencies/include` to the build command. An optional
`--runtime-library-dir /path/to/compiler/lib` embeds a runtime library search
path. Both options can be repeated. `--baseline` selects another Git revision
when a different comparison is intended; the current driver requires the
first-acceptance kernel API introduced in #266. To measure heuristic groups, use
`--heuristic` on the runner. Keep separate build directories for exact and
heuristic results, since the runner replaces its output files.

The build script snapshots both versions of `SeededPartitionFunction` and
`StackSeedDomain`, including their headers; other project headers and configured
types are shared. It records the compiler, exact commands, baseline commit,
and SHA-256 source digests in `manifest.json`. The runner writes full verification
dumps, `samples.tsv`, and `results.json` in the build directory. Timed runs do not
write the full boundary dumps.

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
