# Dense versus hashed seeded-partition states (issue #262)

This investigation compares the dense implementation at
`e3479d693ae335780537c8bf3311aaa1b5ac6f76` on
`docs/257-base-pair-probability-plan` with a storage-only hash prototype.
The [issue](https://github.com/BackofenLab/IntaRNA/issues/262) requests time and
space measurements before choosing whether to replace the production matrix.


**Decision: retain dense production storage.** Hashing substantially helps large,
sparse rectangles, but the two measured natural RNA workloads are 19–23% slower
for partition-only prediction and 10–13% slower with pair output, with essentially
unchanged peak RSS. The experimental implementation is included as a patch.

## Representations

The original `Matrix<Z_type>(n*m, states)` already stores all states of a pair
contiguously. It is not physically a collection of separately allocated state
matrices. The proposed gain comes from omitting invalid pairs, rather than
changing the order of the state values.

The reproducible [prototype patch](../seed-partition-storage/hash-storage.patch)
replaces only `h` and `dh` with
`unordered_map<Interaction::BasePair, vector<Z_type>, Interaction::BasePair::Hash>`.
Each valid local pair owns one contiguous zero-initialized vector with `states`
forward values, followed by `states` reverse values when pair output is enabled.
`states = max(2, maximum seed length) + 2` is determined at runtime, so a single
fixed-size `std::array` cannot represent all supported seed lengths. The vector
implements the requested runtime-length array without imposing a seed limit.

The table reserves exactly the admitted-pair count. Current-cell and stacked
predecessor pointers are looked up once outside state loops. All entries are
inserted before recurrence traversal. The dense `total`, `dt`, `valid`, and
`ending` matrices remain common to both versions, as does the full output mass
matrix. Every valid pair is allocated, including unreachable/zero-weight cells;
there is no reachability pruning. Forward/reverse coordinate traversal, arithmetic
order, numerical checks, and complete-boundary callbacks are unchanged. No
recurrence iterates hash buckets.

## Memory accounting

Let `N = n*m` be one cropped start rectangle, `A` its admitted pairs,
`s = max(2, maximum seed length)+2`, and `k = 1` for partition-only or `2` for
partition plus outside. With eight-byte `Z_type` and `size_t`, the dense working
payload is

```
N * (8*k*(s+2) + 9) bytes.
```

This includes `h`/`dh`, two aggregate totals per direction, validity and seed
endings, but excludes container objects and allocator overhead. Pair mode also
allocates `8*d.n*d.m` bytes for the full raw output. Starts are processed
sequentially, so these working rectangles are not all live at once.

The [allocation probe](../seed-partition-storage/allocation.cpp) uses empty
counting allocators and checks that vector, map-value and map-object sizes equal
their native-allocator equivalents. GCC 14/libstdc++ requests 56 bytes per map
node plus eight bytes per bucket. Reserving 10,000 admitted pairs gives 10,273
buckets. The hashed **state-storage** allocation requests are therefore

```
A * (8*k*s + 56) + 8*bucket_count bytes,
```

compared with `N*8*k*s` for dense state storage. The common metadata/output terms
cancel when comparing these requests. Measured break-even occupancies for this
10,000-entry probe are:

| Seed length | States | Partition only | With outside |
| --- | ---: | ---: | ---: |
| 1 or 2 | 4 | 33.26% | 49.91% |
| 4 | 6 | 42.77% | 59.92% |
| 7 | 9 | 52.86% | 69.16% |
| 20 | 22 | 73.27% | 84.57% |

These are requested allocation bytes, not RSS break-even predictions. They
exclude malloc headers, rounding, fragmentation, and allocator retention. Bucket
growth also varies with the admitted count. The layout is library/platform
specific; [raw probe output](../seed-partition-storage/allocation.csv) is included.

## Measurement method

The host is Linux x86-64, AMD Ryzen 5 7530U, with GCC 14.4.0, ViennaRNA 2.7.2,
Boost 1.85, Kokkos mdspan, and `Z_type=double`. Both release variants use the same
Autotools configuration and compiler options, including C++23, `-O3`,
`-fno-strict-aliasing` and OpenMP. Processes are sequential and pinned to CPU 0,
with one OpenMP thread. Compilation and test runs finish before measurements.

The [kernel driver](../seed-partition-storage/benchmark.cpp) calls the actual
production `compute` translation unit or its patched counterpart. It covers
64/192-square domains, requested densities 1%, 10%, 37.5%, 100%, seed lengths
1/7/20, both pair-output modes and both noLP settings; it also includes
multiple-start domains and 1024-square sparse memory stress cases. All 116
workloads use deterministic masks and weights. Seed diagonals are forced valid;
the output records actual occupancy and admitted-start counts. The 37.5% mask
is a density control, not a model of RNA sequence correlations. Real RNA is
tested separately through the CLI.

Before timing, the runner byte-compares hexadecimal floating-point dumps of
every ordered complete-boundary callback, raw partition total and every output
mass cell, including zeros. Timed calls use only a callback counter, without
dump I/O. The core timer includes per-start allocation, initialization, both
requested recurrence passes, and destruction of their working state; it excludes
input-mask/seed-domain preparation and output serialization. Each workload gets
one warmup per variant, then repeated fresh processes in alternating variant
order. GNU time measures peak RSS including process startup, input setup and
the returned matrix. It also provides a small launcher process, avoiding an
inherited Python-parent RSS high-water mark in a directly launched child.
Very short timings and small RSS differences should not be treated as precise
speedups or memory savings.

The [native runner](../seed-partition-storage/native.py) measures complete CLI
wall time and GNU time peak RSS: startup, accessibility, seed preprocessing,
prediction, normalization and output. It compares OxyS/fhlA (112x108, seed 4,
span 30), GcvB/phoB (299x201, seed 7, span 60), dense G40/C40, and sparse400
with eight complementary bases per strand. Natural inputs compute ViennaRNA
accessibility; synthetic inputs disable it. Dimensions above are target by
query. Loop limits are three, noLP and
terminal-GU restrictions are disabled, and no ranked interactions are requested.
Both ensemble-only and pair-CSV modes run with one warmup and seven measured
repetitions per variant. Printed ensemble output (including rounded `Eall`) and
all pair-CSV bytes must match on every run; the kernel dumps provide the stronger
raw-value comparison.

## Results and decision

Retain the dense production kernel. The tested hash representation is useful for
large sparse rectangles and long seeds, but it slows the measured normal-length
native predictions without a material whole-process memory reduction. This PR
ships the experiment and evidence; it does not replace production storage or
introduce an automatic backend threshold.

All 116 kernel workloads have byte-identical raw results, every mass cell and
ordered callbacks. The final run contains 1,624 timed kernel samples (seven per
variant/workload); all eight native input/output configurations match too.

### Kernel: density at seed length seven

Single start, 192x192 rectangle, loop limits two, noLP disabled. Values are
medians of seven fresh processes; times are milliseconds. RSS columns are MiB
for pair-output mode (whole process, not isolated state allocations).

| Requested / actual valid density | Z dense | Z hash | Z+pairs dense | Z+pairs hash | Pair RSS dense / hash |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1% / 0.98% | 2.116 | 0.726 | 4.174 | 1.244 | 17.75 / 12.75 |
| 10% / 9.70% | 2.386 | 1.425 | 4.874 | 2.430 | 17.88 / 13.50 |
| 37.5% / 37.13% | 3.101 | 4.395 | 7.260 | 8.091 | 17.75 / 15.62 |
| 100% / 100.00% | 5.389 | 11.878 | 15.539 | 25.182 | 17.88 / 21.00 |

At the 37.5% density control the hash prototype is 42% slower for Z and 11%
slower with outside, despite somewhat lower memory. Fully valid rectangles make
it 2.20x and 1.62x slower respectively and increase RSS. Longer seeds change the
tradeoff: at the same 192-square/37.5% control, length 20 reduces Z time from
5.718 to 5.197 ms and pair time from 13.816 to 10.324 ms. Thus the allocation
crossover alone is insufficient to predict a runtime crossover. The full data
also includes noLP and multiple-start cases; the tables are representative
slices, not the entire test space.

### Kernel: large sparse memory stress

Single start, 1024x1024 rectangle, length 20, full spans, noLP disabled. These
intentionally large working rectangles expose savings hidden by normal bounded
spans. Times are milliseconds and peak RSS is MiB.

| Requested density | Output | Dense time | Hash time | Dense RSS | Hash RSS |
| --- | --- | ---: | ---: | ---: | ---: |
| 1% | Z | 121.246 | 23.332 | 212.88 | 39.62 |
| 1% | Z+pairs | 235.802 | 39.243 | 413.00 | 65.38 |
| 10% | Z | 134.938 | 73.862 | 212.88 | 63.25 |
| 10% | Z+pairs | 271.951 | 119.126 | 413.00 | 104.62 |

The 1% pair-output case is about 6x faster and reduces peak RSS from 413 MiB
to 65.4 MiB. This is a real benefit for this workload. It is not a prediction
for span-30 RNA runs: their per-start dense state payload with seed length four
is only 43,200 bytes for Z, or 86,400 with outside. The prototype still retains
all metadata and the global pair matrix, even when almost no pairs are valid.

### Native CLI

Seconds and KiB, medians of seven fresh processes. The ratio is hash time divided
by dense time. Full samples and min/max times accompany this report.

| Input | Output | Dense seconds | Hash seconds | Time ratio | Dense / hash RSS (KiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| OxyS-fhlA | Z | 0.0995 | 0.1183 | 1.189 | 17,180 / 17,184 |
| OxyS-fhlA | Z+pairs | 0.1715 | 0.1881 | 1.097 | 17,056 / 17,308 |
| GcvB-phoB | Z | 6.1391 | 7.5595 | 1.231 | 18,712 / 18,708 |
| GcvB-phoB | Z+pairs | 14.1424 | 16.0409 | 1.134 | 18,708 / 18,712 |
| sparse400 | Z | 0.1152 | 0.1146 | 0.994 | 17,920 / 17,920 |
| sparse400 | Z+pairs | 0.1546 | 0.1561 | 1.010 | 18,768 / 18,048 |
| dense40 | Z | 0.3644 | 0.4274 | 1.173 | 16,472 / 16,420 |
| dense40 | Z+pairs | 0.8414 | 0.9117 | 1.084 | 16,496 / 16,768 |

Hashing costs about 19–23% on the two natural ensemble-only inputs and 10–13%
with pair probabilities. Peak RSS is essentially unchanged. Sparse400 shows no
consistent end-to-end speed benefit: its active working rectangles are small,
while sequence preprocessing and full pair output remain. Small RSS differences
should not be extrapolated beyond these runs.

The recommendation applies to this `std::unordered_map` plus runtime-array
prototype. A pooled/compact representation or a workload-specific backend could
have different tradeoffs; neither is measured here. In particular, do not turn
the allocation-only break-even percentages into a production switch without
separate runtime/correctness validation.

Raw artifacts: [kernel summaries](../seed-partition-storage/results/kernel-summary.csv),
[all kernel samples](../seed-partition-storage/results/kernel-samples.csv),
[exact-output hashes](../seed-partition-storage/results/kernel-correctness.csv),
[kernel run metadata](../seed-partition-storage/results/kernel-run.json),
[native summaries and output hashes](../seed-partition-storage/results/native-summary.json),
[native samples](../seed-partition-storage/results/native-samples.csv), and
[environment/source hashes](../seed-partition-storage/results/environment.json).
All committed kernel RSS values are from the corrected GNU time runner. The
initial direct-launch run's inherited Python RSS floor was detected and that
entire kernel run was repeated; its measurements are not used in these tables.

## Reproduction

Use a clean checkout of this PR with C++23 dependencies installed; its production
kernel is unchanged from the recorded base.
Configure and build according to [AGENTS.md](../../AGENTS.md); use identical
options for both variants. From the configured source root, set a scratch output
directory outside the checkout, then:

```sh
storage_results=$(mktemp -d)
mkdir -p "$storage_results"
cp src/bin/IntaRNA "$storage_results/IntaRNA-dense"
cp src/IntaRNA/SeededPartitionFunction.cpp "$storage_results/dense.cpp"
python3 doc/seed-partition-storage/build.py --root . \
  --kernel "$storage_results/dense.cpp" --output "$storage_results/benchmark-dense"

git apply doc/seed-partition-storage/hash-storage.patch
make -j2
cp src/bin/IntaRNA "$storage_results/IntaRNA-hash"
cp src/IntaRNA/SeededPartitionFunction.cpp "$storage_results/hash.cpp"
python3 doc/seed-partition-storage/build.py --root . \
  --kernel "$storage_results/hash.cpp" --output "$storage_results/benchmark-hash"
# Run the validation described below before restoring the source.
git apply -R doc/seed-partition-storage/hash-storage.patch
make -j2

python3 doc/seed-partition-storage/run.py --cpu 0 --repeats 7 \
  --dense "$storage_results/benchmark-dense" \
  --hash "$storage_results/benchmark-hash" --output "$storage_results/kernel"
python3 doc/seed-partition-storage/native.py --cpu 0 --repeats 7 --source . \
  --dense "$storage_results/IntaRNA-dense" \
  --hash "$storage_results/IntaRNA-hash" --output "$storage_results/native"
```

`build.py` extracts the actual compiler, flags and libraries from the configured
Makefile. Its source translation unit resolves `compute` before the static
archive, so only that implementation is linked. Python 3.8+, Linux affinity
support, `taskset`, and GNU time are needed for these runners. Choose an available
CPU if CPU 0 is not permitted. Compile `allocation.cpp` with the same compiler
and include paths (no IntaRNA library is needed), then run `allocation 10000`.

## Correctness and baseline validation limitations

The dense and hash release kernels each pass the existing independent seeded
oracle tests:

```sh
make -C tests runApiTests
./tests/runApiTests '[SeededPartitionFunction]'
# 15,364 assertions in seven test cases, for each implementation.
```

This tag covers the heterogeneous-weight enumeration, native seeded oracle and
checked-arithmetic tests, beyond the deterministic benchmark fixtures. The hash
debug kernel passes the same 15,364 assertions. Full `make tests -j2` was also
attempted for both implementations in both configurations. The full suite is
**not green on the recorded base**, independently of this storage experiment:

| Configuration | Dense | Hash |
| --- | --- | --- |
| Release CLI suites | All six pass | All six pass |
| Release API suite | Aborts in `StackSeedDomain_test.cpp:26` | Same abort |
| Debug CLI suites, with syntax prerequisite below | All six pass | All six pass |
| Debug API suite, with syntax prerequisite below | 96/98 cases; 100,962/100,964 assertions | Identical counts and failures |

The release abort occurs when a test passes an unfilled bulged seed handler,
expecting rejection. The base moved the stack-capability rejection into a
debug-only guard, so release proceeds to access an empty seed matrix.
The conflicting-seed-length exception test likewise expects a release check
that the base now compiles out. An exploratory release run excluding the
stack-domain tag also found a missing overlap exception and an invalid-region
crash in `BasePairProbabilities_test.cpp`, whose checks were similarly moved
behind debug guards. These are not reported as successful full API runs.

The unmodified base does not compile in debug mode because
`PredictorMfeEns2d.cpp:60` has one excess closing parenthesis. Only the isolated
debug validation copy received this prerequisite correction:

```diff
-       || r2.from>=energy.getAccessibility2().getSequence().size()))
+       || r2.from>=energy.getAccessibility2().getSequence().size())
```

After that correction, dense and hash debug runs reproduce the same two failed
assertions: `StackSeedDomain_test.cpp:108` (explicit-handler stack capability)
and `Matrix_test.cpp:147` (nonzero lower-band rejection). The debug prerequisite
and these unrelated failures are not changes to production in this PR.
No golden output or existing test was altered to conceal them.

For selective Catch runs, explicitly exclude `[.]` if excluding tags: an
exclusion-only filter otherwise enables hidden benchmark cases. The exploratory
filtered debug run was cancelled and is not counted as validation.
