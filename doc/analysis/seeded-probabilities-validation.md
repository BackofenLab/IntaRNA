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


## Steps 8–9: independent output and documentation audit (2026-10-09)

Separate agents reviewed code correctness, recurrence depictions, and performance
before the audit fixes were pushed to PR #259. The audit found and corrected:

- SVG accessibility frames used the interaction model's RT for imported Pu/ED.
  The writer now requires separate target/query accessibility RT values. CLI
  conversion selects RT=1 for computed base-pair accessibility and the configured
  physical temperature for ViennaRNA/imported accessibility. The same selection
  fixes the pre-existing `tPu`/`qPu` conversion of computed base-pair ED.
- SVG and CSV stream buffers were copied when returned; moving the completed
  strings removes that duplicate allocation without changing document bytes.
- The heuristic depiction omitted its initial right-end GU gate, and the inside
  depiction placed aggregation before incoming nonstack contributions. These
  were corrected, along with complete-site filtering and outside-index notation.

The added CLI regressions fail on the pre-audit binary and pass after the fix.
They cover Pu/ED imports, mixed computed/imported strands, 22/37°C, original
query orientation, independent input-Pu round trips (allowing ED quantization),
and computed-base-pair ED conversion at RT=1. API tests reject nonpositive or
nonfinite accessibility scales without writing output. GCC 14.4 release and
debug `make tests -j2` each pass all seven suites, with 91,394 API assertions in
86 cases. Installed public headers and an installed SVG consumer compile; the
consumer produces the expected matrix and seed annotations. Corrected SVG
recurrence documents parse, render, and were visually inspected. Apple Clang
was not available locally.

### Performance method and existing-mode comparison

Compared base `48efc541e23571ce1bab88e2754e1db7445cbc89` with steps 8–9 at
`fa2f84d6cdbb1e120a98eab9d2dfe3f8154ebad5`, then compared that saved binary with
the audit fixes. All builds use identical GCC 14.4 C++23 `-O3
-fno-strict-aliasing -fopenmp`, ViennaRNA 2.7.2, Boost 1.85 and Kokkos mdspan
settings on Linux x86-64 / Ryzen 5 7530U. Processes are pinned to CPU 0 with
`--threads=1`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and
`OMP_DYNAMIC=FALSE`. Builds and tests are stopped during measurements.

There is one warmup per configuration, then seven fresh-process runs in shuffled,
interleaved order for the base comparison. Time is wall time around the GNU time
subprocess; RSS is GNU time `%M`, including executable/libraries, prediction and
serialization. Output goes to the same local filesystem. Initial median seconds:

| Input | Base ensemble only | Steps 8–9 ensemble only | Base CSV | Steps 8–9 CSV | Steps 8–9 SVG | Steps 8–9 both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OxyS/fhlA | 0.0975 | 0.0978 | 0.1710 | 0.1703 | 0.1910 | 0.1947 |
| dense40 | 0.3876 | 0.3913 | 0.8958 | 0.9026 | 0.9090 | 0.9060 |
| sparse400 | 0.1191 | 0.1202 | 0.1616 | 0.1633 | 0.7326 | 0.8260 |
| sparse1000 | 1.5816 | 1.5801 | 1.8284 | 1.8372 | 4.0202 | 4.2578 |

Existing-mode medians differ by −0.4% to +1.1%, within observed timing variation.
Ensemble outputs and CSVs are byte-identical across base/candidate and modes;
all SVG pair values equal their CSV values. Sparse SVG timings have wider
variation than the compute-heavy cases; these are local workload measurements,
not a general throughput guarantee. Full SVG size scales with all target/query
pairs: the sparse1000 document is 220,685,460 bytes despite having few seed cells.

Reproduce each process with the following command from the corresponding build;
add either/both output options for the CSV/SVG configurations:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
  /usr/bin/time -f '%e %U %S %M' taskset -c 0 src/bin/IntaRNA \
  --target=doc/handson/fhlA.fasta --query=doc/handson/OxyS.fasta \
  --model=P --mode=M --seedBP=4 --intLenMax=30 \
  --qIntLoopMax=3 --tIntLoopMax=3 --threads=1 \
  --outNoLP=false --outNoGUend=false --outNumber=0 --outMode=E \
  --default-log-file=/dev/null --out=ensemble.txt
# Optional: --out=bpProb:pairs.csv --out=bpsvg:plot.svg
```

For dense40, replace the inputs with `G`×40 and `C`×40 and add `--tAcc=N
--qAcc=N`. For sparse length n (400 or 1000), use target `A`×(n/2) + `G`×8 +
`A`×(n/2−8), query with `C`×8 instead of `G`×8, and the same disabled-accessibility
options. All other flags remain identical.


### Audit-fix comparison

Five interleaved fresh-process runs per configuration compare the retained
`fa2f84d` binary with the rebuilt audit fix. All complete SVG, CSV and ensemble
outputs match byte-for-byte on these native-energy inputs. Time is median
[min, max] seconds; RSS is the maximum over the five runs, in KiB.

| Input/output | Before time | Fixed time | Before RSS | Fixed RSS |
| --- | ---: | ---: | ---: | ---: |
| OxyS SVG | 0.1938 [0.1932, 0.1968] | 0.1934 [0.1918, 0.1938] | 22,132 | 20,720 |
| OxyS both | 0.1965 [0.1945, 0.1976] | 0.1959 [0.1955, 0.1967] | 22,004 | 20,596 |
| sparse400 SVG | 0.5020 [0.5003, 0.7993] | 0.7491 [0.4806, 0.9672] | 87,876 | 84,548 |
| sparse400 both | 0.5553 [0.5441, 0.8527] | 0.7544 [0.5118, 0.8899] | 87,872 | 84,420 |
| sparse1000 SVG | 4.0260 [3.9535, 5.8240] | 3.8527 [3.8373, 5.5779] | 470,084 | 301,128 |
| sparse1000 both | 4.2407 [4.2148, 6.2242] | 4.0592 [4.0216, 4.1580] | 469,960 | 301,260 |

Transferring the document buffer reduces 1000×1000 SVG peak RSS by **35.9%**.
There is no general elapsed-time speedup claim: filesystem writes varied widely,
and the sparse400 median increased even though its median process CPU time fell
from 0.44 to 0.42 seconds (SVG) and 0.48 to 0.46 seconds (both). Dense40 controls
remain comparable: ensemble-only 0.3848→0.3902 seconds (+1.4%) and CSV
0.8982→0.8869 seconds (−1.3%), with unchanged output bytes. The output is still
buffered as a full document, so very large matrices remain memory intensive.

One comparison process received an unexplained SIGTERM after the shorter cases;
partial long-case warmups were discarded and the long case completed in a new
process. The table includes only complete five-run groups.

To separate serializer cost from file-write stalls, a final sparse400 control
sent SVG to `/dev/null` and ensemble output to STDOUT redirected to `/dev/null`.
Seven interleaved runs gave wall medians 0.412419→0.390979 seconds (ranges
0.411222–0.413401 and 0.388902–0.392267), CPU medians 0.40→0.38 seconds, and peak
RSS 87,744→84,544 KiB. The slowdown did not reproduce without filesystem writes;
the fixed serializer was 5.2% faster in this specific control. Use the same
sparse400 command with `--out=STDOUT --out=bpsvg:/dev/null > /dev/null` to repeat.
