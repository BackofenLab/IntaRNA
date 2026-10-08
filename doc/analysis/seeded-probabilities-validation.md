# Exact seeded probabilities: implementation and validation

PR #258 implements the seven steps of the reviewed plan in separate commits.
The library ships the specialized suffix automaton; maximal stacks remain a
correctness/performance comparison in the tests. See the [kernel decision](seeded-kernel-comparison.md)
for measured tradeoffs and exact cell-storage counts.

## Delivered contract

- Stack-only capability is geometric and forwarded through index offsets. Active
  occurrences retain each handler's admission and per-start selection policies;
  both endpoints and inclusive span limits are checked.
- Seed-free suffix lengths and two seeded run-length states count ordered chains
  once. Whole-interaction noLP forbids departure/termination of singleton runs.
  End GU and full-span ED filters affect termination, not intermediate states.
- Initiation, successive-pair loop energies and the complete boundary coefficient
  own their contributions once. The complete coefficient is shared with outside
  initialization and virtual `updateZ` dispatch. Overrides can observe/delegate;
  reweighting belongs in `exactBoundaryWeight`, which serves both passes. Changing
  the hybrid mass through an old update override is rejected.
- The pair sink is non-owning. Each prediction creates private region state and
  converts internal pairs once through `getBasePair`. The accumulator rejects
  intersecting rectangles, validates every merge before mutation, and requires
  explicit finalization. CLI ownership spans all regions of one sequence pair.
  Cancellation/failure guards apply before output, including the OpenMP path.
- `bpProb` is independent of coverage trackers and ranked reporting. Pair-only
  output does not retain `Z_partition`. Simultaneous legacy trackers retain their
  existing callbacks/map and native energy rounding semantics.

## Numerical limits

Arithmetic uses `Z_type`, exact zero and checked nonnegative sums/products,
exponentials, division and region merges. Structural forbidden energies become
zero before exponentiation. NaN, infinity, positive underflow and positive
subnormals are errors. No scaling is implemented. Ordinary rounding loss in a
sum with widely different terms remains possible. Final row/column exclusivity
and pair mass bounds are checked; endpoint correction is limited to
`512 * epsilon(Z_type) * max(1, targetLength + queryLength)` above one.

Independent tests use a 2e-12 relative tolerance for raw double-precision totals,
boundaries and masses against long-double enumeration and 50-digit theory gold
values. The CLI base-pair-energy oracle also uses a 2e-14 absolute tolerance for
printed zeros/small values. Native integer ED is preserved and is not compared
against unrounded theoretical ED as though they were identical models.

## Memory

Beyond the active suffix grids documented in the kernel comparison, pair output
holds one full target-by-query raw matrix, the private regional outside matrix,
and briefly a regional matrix in original coordinates. Explicit finalization and
writing allocate a normalized full matrix; writing also buffers the complete CSV
block. These phases are sequential, not all simultaneous. Native accessibility
and retained seed preprocessing are additional costs. Legacy trackers add their
own coverage matrices and complete-boundary map. No global run tape or global
boundary-weight table is retained. Without pair output, no pair matrix or reverse
state is allocated.

## Reproducible checks

Use GCC 14 with C++23, Boost, ViennaRNA and the bundled Kokkos mdspan fallback:

```sh
bash autotools-init.sh
./configure --disable-debug # dependency prefixes as in AGENTS.md when needed
make -j2
make tests -j2
# Repeat in a clean source/build directory with --enable-debug.
```

The focused tags are `[SeededPartitionFunction]`, `[StackSeedDomain]`,
`[BasePairProbabilities]`, `[PredictorSeedOracle]` and `[PredictorMfeEns]`.
`tests/runBasePairProbabilities.sh` checks CLI behavior with a separate Python
chain enumerator. `tests/tools/export-seeded-theory.py` reproduces all ten
fixture extractions from the cited theory revision. No legacy CLI golden file
is regenerated for this change.

Final local verification used GCC 14.4.0, ViennaRNA 2.7.2, Boost 1.85,
Z_type=double and Kokkos mdspan on Linux x86-64. Release and debug each passed all
seven Autotools test suites; the API suite passed 91,361 assertions in 85 cases.
The initial incremental run exposed stale objects after the virtual interface
change; clean release/debug rebuilds passed. No failing golden files were replaced.

All 74 installed public headers compiled independently. The standalone probability
consumer compiled, linked through installed pkg-config metadata, and produced the
expected matrix. This check exposed a missing `-leasylogging` dependency in the
metadata; it is fixed, and CI now builds/runs this predictor consumer on Linux
and macOS. `make dist` includes every new source, fixture, test, example and report.
Apple Clang execution was not available locally; its existing CI job remains the
platform gate.

The final integration checks also cover too-short sequences as successful empty
matrices, dynamic OpenMP exception propagation and visible fatal diagnostics even
when a log file is selected. Seed-constraint lazy initialization is synchronized
for concurrent sequence pairs. Pure header/Matrix consumers did not previously
exercise the logging dependency; the new consumer calls a real predictor.

## Native measurements

The hidden `[SeededNativeBenchmark]` case uses the production predictor, native
ViennaRNA transition energies, inclusive spans of 30, loop limits three and a
four-pair computed stack seed. OxyS/fhlA accessibility is ViennaRNA-derived;
synthetic cases disable accessibility. Inputs and the timing adapter are committed
in `tests/SeededNativeBenchmark_test.cpp`. Set `INTARNA_BENCH_CASE` and
`INTARNA_BENCH_OUTSIDE=0|1` to measure an individual case in a fresh process.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/time -v \
  ./tests/runApiTests '[SeededNativeBenchmark]'
```

Milliseconds from one controlled single-thread run; "predictor" excludes measured
seed filling but includes active-domain construction, merging and finalization.
The added outside cost is the difference between the two predictor columns; it
includes pair-result allocation/normalization, not just differentiation.

| Input | Accessibility | Seed filling (Z / pairs) | Predictor Z | Predictor Z+pairs | Z |
| --- | ---: | ---: | ---: | ---: | ---: |

| none | 0.005 | 0.024 / 0.019 | 0.010 | 0.198 | 0 |
| sparse100 | 0.000 | 0.074 / 0.028 | 0.751 | 1.488 | 1.702e+15 |
| sparse400 | 0.000 | 0.308 / 0.288 | 0.712 | 5.621 | 1.702e+15 |
| dense | 0.000 | 0.296 / 0.304 | 413.301 | 972.852 | 2.00398e+68 |
| OxyS-fhlA | 27.627 | 0.180 / 0.178 | 53.376 | 134.278 | 37643.8 |

Raw Z was bitwise identical with pair output off/on for every measured case. No
representative measured case exceeded the checked numerical range. Increasing
unrelated sparse-sequence length from 100 to 400 leaves partition DP work nearly
constant; seed preprocessing and full-matrix pair output still grow. These costs
are reported separately rather than attributed to the recurrence.

Fresh-process peak RSS (KiB, includes the executable, native libraries,
accessibility and allocator overhead; small differences include measurement noise):

| Input | Z only | Z+pairs |
| --- | ---: | ---: |
| sparse100 | 15,436 | 15,376 |
| sparse400 | 15,436 | 18,688 |
| dense | 15,560 | 16,000 |
| OxyS-fhlA | 17,784 | 17,836 |

Full CLI measurement (including startup, accessibility, computation, normalization
and file output) used the repository's `doc/handson/fhlA.fasta` and `OxyS.fasta`:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/time -v src/bin/IntaRNA \
  --target=doc/handson/fhlA.fasta --query=doc/handson/OxyS.fasta \
  --model=P --mode=M --seedBP=4 --intLenMax=30 \
  --qIntLoopMax=3 --tIntLoopMax=3 --threads=1 \
  --outNoLP=false --outNoGUend=false --outNumber=0 --outMode=E \
  --default-log-file=/dev/null --out=ensemble.txt
# Repeat with --out=bpProb:pairs.csv.
```

- Without pair matrix: median 0.09 s over three runs, maximum observed RSS 16,916 KiB.
- With pair matrix: median 0.16 s over three runs, maximum observed RSS 16,788 KiB.

Both runs report Eall=-6.52 kcal/mol with the native presentation convention.
These measurements quantify added output work; they are not a speedup claim over
the legacy seeded ensemble, whose multi-anchor partition is a different result.
