# Seed-extension kernel comparison

This Linux benchmark compares the exact seeded kernel at
`cb1c9d419f96c7d2a10bae864d922dcb0234f4d4` with the current source checkout.
It first checks every reported boundary hybrid weight and boundary coefficient,
the total partition, and every raw pair numerator. It then times partition-only
and pair-output calculations in separate processes, alternating baseline and
candidate order. Three repetitions are the default.

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
when a different comparison is intended.

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
`noLP` is disabled; separate oracle tests validate that restriction.

`driver.cpp` declares deterministic coordinate-dependent initiation, transition,
and full-boundary weights. It measures the complete kernel invocation, including
a small boundary callback. Seed construction, native thermodynamic/accessibility
computation, and serialization are outside this timing. These fixtures exercise
the recurrence, rather than representing complete biological workloads.

Recorded runs used GCC 14.4.0, `-std=c++23 -O3 -DNDEBUG`, one thread, and an
AMD Ryzen 5 7530U on 2026-10-14. Samples were not CPU-pinned and other build/test
work could run on the same workstation. Treat small differences as noise.
`recorded-results.json` contains correctness differences and timing medians with
minimum/maximum values; `recorded-samples.tsv` contains the individual samples;
`recorded-source-hashes.json` identifies the measured source snapshots.

GNU time records whole-process maximum RSS, including runtime libraries and
allocator overhead. These small fixtures do not isolate the recurrence's peak
allocation. Timings and RSS are local observations, not asymptotic bounds.
See the [algorithm and measured tradeoff](../../analysis/seed-extension-partition.md).
