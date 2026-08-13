# Native performance comparison

`compare.cpp` is a standalone C++23 benchmark driver for fair black-box
comparisons between IntaRNA Legacy and IntaRNAnew. It does not import either
implementation and never inspects internal classes.

The driver first compares stdout byte-for-byte for every parameter file. A case
with different scientific/output behavior is rejected and never contributes to
the reported speedup. It then performs configurable warm-ups and alternating
timed subprocess runs at the same thread count. The CSV report contains median
wall time, median child peak RSS, per-case legacy/new speedup, and a geometric
mean over compatible cases.

Build and run from the repository root:

```sh
g++ -std=c++23 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wconversion -Wshadow \
  IntaRNAnew/benchmarks/compare.cpp -o IntaRNAnew/IntaRNAnewBenchmark

LD_LIBRARY_PATH="$PWD/.conda-env/lib" \
  IntaRNAnew/IntaRNAnewBenchmark \
  --legacy "$PWD/src/bin/IntaRNA" \
  --new "$PWD/IntaRNAnew/IntaRNAnew" \
  --cases "$PWD/IntaRNAnew/tests/fixtures/parameters" \
  --warmups 2 --repetitions 9 --threads 1
```

The legacy path can instead be supplied once through `INTARNA_LEGACY_BIN`. The
automated correctness gate is:

```sh
INTARNA_LEGACY_BIN="$PWD/src/bin/IntaRNA" \
  make -C IntaRNAnew legacy-check
```

It asserts that the fixture directory contains exactly 17 cases. When
`INTARNA_LEGACY_BIN` is not set, the gate reports `SKIP` rather than making the
normal native test suite depend on a separately installed legacy executable.

Each child run has a 300-second wall-time limit by default; use `--timeout N`
to select another positive limit in seconds.

The Makefile also provides `make -C IntaRNAnew benchmark`. CMake test builds
compile the driver for the parity gate. If tests are disabled, configure with
`-DBUILD_TESTING=OFF -DINTARNANEW_BUILD_BENCHMARKS=ON` to build it explicitly.

For a correctness-only executable gate, replace the warm-up/repetition options
with `--verify-only`. The driver routes legacy informational logs away from
stdout before comparing the scientific output.

The runner uses POSIX process and resource-accounting APIs and is therefore a
development benchmark, not an installed portable library component. Keep
correctness and performance separate: never suppress the parity gate merely to
produce a faster number.

`benchmarks/cases/` contains bounded biological workloads using the public
fhlA/OxyS, phoB/GcvB, and ilvE/GcvB-ST FASTA data shipped with this
repository. They select the table-independent base-pair model and compare
the first three exact-mode interactions byte-for-byte before timing. This
bounded cutoff avoids ending inside the workloads' large equal-energy tie groups.
The tiny unit-style
parameter corpus is useful as a compatibility gate, but its timings are mostly
process startup and should not be presented as algorithmic performance.
Run these biological cases from the repository root as in the command above.

## Reference result after the output-aware redesign

On 2026-08-13, using GCC 13.3.0 on an AMD Ryzen 5 7530U, the release build
produced these medians with two warm-ups, nine alternating repetitions, and one
thread:

| Case | Legacy | IntaRNAnew | Legacy/new speedup | Legacy RSS | New RSS |
|---|---:|---:|---:|---:|---:|
| fhlA/OxyS | 504.936 ms | 36.066 ms | 14.000x | 7,424 KiB | 5,120 KiB |
| ilvE/GcvB-ST | 4008.653 ms | 208.473 ms | 19.229x | 7,364 KiB | 5,120 KiB |
| phoB/GcvB | 3821.776 ms | 205.155 ms | 18.629x | 7,364 KiB | 5,120 KiB |
| geometric mean | | | **17.117x** | | |

All three cases passed the byte-for-byte gate before timing. These cases exercise
the output-aware exact/base-pair/model-S scalar kernel. Configurations requesting
partitions, complete sites, tracebacks, constraints, seeds, other models, or
multi-domain reduction deliberately use the general predictor and need separate
performance characterization.
