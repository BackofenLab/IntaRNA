# Reusing work in MFE seed extensions

This implements [issue #274](https://github.com/BackofenLab/IntaRNA/issues/274).
Following the [PR #276 review](https://github.com/BackofenLab/IntaRNA/pull/276#discussion_r4312235144),
the current implementation retains matrix-capacity reuse and bounded predecessor
scans, and evaluates pairing and physical loop energies directly through the
energy model. The local loop-energy cache and its supporting pair map have been
removed. Historical measurements below describe the earlier cached prototype.

The analysis and measurement baseline is `perf/270-anchored-heuristic` at
`21b9d8fc33aadd3c44d788c4da10a91cf6694ce1`. Both native MFE seed-extension
predictors are selected with `--model=X`: `--mode=M` is exact over the
handler-selected MFE seeds, and `--mode=H` chooses one optimal right extension
per seed for combination with left extensions. As before, neither considers
every possible seed structure at a given start.

## Comparison with ensemble prediction

The ensemble predictors (`--model=P`) now share left-most seed anchors and
bounded reusable left/right workspaces. See the
[ensemble recurrence](seed-extension-partition.md) and
[heuristic analysis](heuristic-anchored-runtime.md).
The following distinctions determine which techniques transfer safely.

| Technique in the ensemble implementation | MFE application |
| --- | --- |
| Reuse workspace; reshape without retaining obsolete cells | Applied to both MFE extension matrices, including traceback. Each recurrence initializes every active cell before reading it. |
| Evaluate a pair domain once and restrict recurrence transitions | Compute loop predecessor bounds once per cell. Preserve predecessor traversal order and evaluate complementary pairs directly. |
| Share reusable work across many seed anchors | Reuse scratch capacity across seeds. Evaluate physical loop energies directly in both left/right fills and traceback, as requested in review. |
| Left-most seed ownership and elimination of redundant containing anchors | Not applied. Ensemble ownership prevents double counting in sums; MFE predictions retain their existing seed enumeration, tie handling, tracker updates and approximate `Zall` accounting. |
| Batch complete boundaries and differentiate partition functions for pair probabilities | Not needed for MFE minimization. It would add state and alter the organization of candidate publication. |

An extension matrix cell is the optimal path **to or from a particular seed
boundary**. Overlapping seeds have different conditioned subproblems; copying a
matrix cell merely because the seeds overlap is incorrect. The ensemble kernel
also resets its per-anchor tables rather than copying numerical values between
arbitrary overlapping anchors. This implementation reuses their storage while
recomputing each anchor's MFE recurrence. Physical loop energies still come from
`energy.getE_interLeft()`.

## Implementation and limits

`PredictorMfe2dSeedExtension` owns the two scratch matrices used by both exact
and heuristic MFE extension prediction. The `resize(..., false)` calls
grow extension matrices on demand, reuse capacity, and avoid copying values that
are immediately overwritten. They do not eagerly reserve the full sequence
rectangle, which would waste space when only sparse or central seeds exist.

Pair checks use `energy.areComplementary()` and loop-energy queries use
`energy.getE_interLeft()` in recurrence fills and traceback. There is no local
pair map, loop table, known-value flag array, cache budget or cache lifecycle.
Removing these members changes the predictor object layout; C++ consumers must
rebuild with matching headers and library.

The predecessor ranges are the same inclusive ranges as before; only their
construction moves out of the innermost scans. Energy arithmetic, seed order,
candidate update order, right-end ties, output filtering, traceback and the
existing underestimated MFE partition remain unchanged. Worst-case recurrence
complexity is unchanged. The exact predictor's exhaustive left/right
combination can still dominate runtime, whereas the heuristic spends more of
its time building extension tables.

## Cache-removal review validation

GCC 14.4.0 release and debug `make tests -j2` each passed all eight test drivers
and 119 API cases. The independent recurrence, workspace, tie, maximal-loop,
range-change, traceback and seed-witness tests remain enabled. Both affected
installed public headers compile independently.

One alternating baseline/candidate pass across all 12 native harness cases
produced 24 byte-identical outputs, including energies, boundaries, traceback,
seed starts and `Eall`. This comparison used the cached PR #276 head
`cd193f7b58f3706882ce2d1f24bf0cb9e4f448cd` as baseline and the direct-evaluation
revision as candidate. It checks output equivalence; one sample per variant
does not establish a performance result. All 177 tracked `src` files in each
validation copy were checked against the intended source revision.

The shared base contains the separate [issue #282](https://github.com/BackofenLab/IntaRNA/issues/282)
Autotools bootstrap defect. Its two Makefile corrections were applied only to
disposable validation copies so these tests could run; they remain a separate
change from this cache-removal review. Linux was validated locally; macOS was
not run.

## Historical cached-prototype validation and measurements

The following results belong to the cached implementation at
`4994edf7bb754d305c00503f57bbb3d937bfad5b`; they do not measure the current
implementation after cache removal. Cache-only probes were removed along with
the cache. Independent recurrence, tie, maximal-loop, range-change and traceback
regressions remain applicable to the current implementation.

Release and debug `make tests -j2` each passed all seven test drivers. The API
suite passed 119 cases, with 462,526 assertions in release and 462,527 in debug.
Four new focused cases passed 4,691 assertions in each configuration. They
check independently enumerated MFE/valid heuristic traces and seed witnesses,
native and base-pair energy models, equal-energy ties, noLP/noGUend, asymmetric
and maximal loop bounds, shrinking/growing offset ranges, and seedless calls.
Cache probes verify finite/infinite energy reuse, invalidation, direct long-loop
fallback, and both memory-budget fallback paths. Installed public-header
compilation and source-snapshot checks also passed. A separate deterministic
120-case CLI comparison also matched byte-for-byte across M/H, native/base-pair
energies, noLP/noGUend, cropped ranges, bulged seeds and asymmetric loop limits.
Linux/GCC validation was
performed locally; macOS/Apple Clang was not run.

The measurement setup uses GCC 14.4.0, ViennaRNA 2.7.2, Boost 1.85, C++23 release
builds (`-O3`, `NDEBUG`), and an AMD Ryzen 5 7530U on Linux. No project builds or
tests ran during timing. CPU placement/frequency and other workstation activity
were uncontrolled. These are measurements of the listed inputs on this machine,
not a general speedup guarantee. Seven fresh-process samples per variant were
run in alternating order for each case; values below are median wall seconds.

| Case | Baseline s | Candidate s | Candidate / baseline | Median RSS KiB, baseline → candidate |
| --- | ---: | ---: | ---: | ---: |
| oxys-M | 0.2282 | 0.2106 | 0.923 | 16,992 → 20,244 |
| oxys-H | 0.0783 | 0.0634 | 0.810 | 16,944 → 20,572 |
| oxys-noacc-M | 0.6930 | 0.6169 | 0.890 | 14,416 → 18,176 |
| oxys-noacc-H | 0.1597 | 0.0835 | 0.523 | 14,652 → 18,176 |
| dense-M | 3.6327 | 2.7054 | 0.745 | 14,520 → 16,256 |
| dense-H | 1.1801 | 0.3434 | 0.291 | 14,644 → 16,384 |
| dense-nolp-M | 2.4051 | 1.7558 | 0.730 | 14,524 → 16,384 |
| dense-nolp-H | 0.9321 | 0.2791 | 0.299 | 14,656 → 16,256 |
| bulged-M | 20.4365 | 18.7458 | 0.917 | 15,872 → 20,096 |
| bulged-H | 3.1374 | 1.2525 | 0.399 | 16,000 → 20,224 |
| no-seeds-M | 0.0089 | 0.0089 | 0.999 | 14,276 → 14,464 |
| no-seeds-H | 0.0088 | 0.0090 | 1.020 | 14,432 → 14,464 |

All 168 measured predictions were byte-identical to the corresponding baseline,
including energies, boundaries, traceback, seed starts and `Eall`. The ordinary
OxyS/fhlA case improves by about 8% in exact mode and 19% in heuristic mode.
With accessibility disabled, the heuristic case improves by about 48%. Dense
heuristic cases take about 29–30% of baseline time (roughly 3.3–3.4 times faster).
Exact dense cases improve less because all left/right boundary combinations
still need to be published. The bulged exact case likewise remains dominated
by combination work. The table includes the no-seed controls; their small
startup-time differences do not measure the optimized recurrences.

Peak RSS increases on seeded cases, consistent with the bounded loop cache.
The recorded sample ranges and whole-process memory measurements should be
considered alongside these medians. The implementation does not claim universal
improvement, especially for very sparse seeds or cheap custom energy functions.

Raw timing samples, full command/binary manifests, summaries, and source/build
verification are retained in the benchmark directory as
[`mfe-native-samples.tsv`](../benchmarks/seed-extension/mfe-native-samples.tsv),
[`mfe-native-manifest.json`](../benchmarks/seed-extension/mfe-native-manifest.json),
[`mfe-native-results.json`](../benchmarks/seed-extension/mfe-native-results.json), and
[`mfe-source-hashes.json`](../benchmarks/seed-extension/mfe-source-hashes.json).

## Reproduction

To reproduce the historical cached-prototype measurements, build the baseline
revision above and candidate `4994edf7bb754d305c00503f57bbb3d937bfad5b` with the
same release compiler, dependencies and flags. The harness can also compare
any other two revisions, but the historical timings do not apply to them. Run:

```sh
python3 doc/benchmarks/seed-extension/mfe-native.py \
  --binary /path/to/baseline/src/bin/IntaRNA \
  --candidate /path/to/candidate/src/bin/IntaRNA \
  --output-dir /tmp/mfe-seed-extension --repetitions 7
```

The harness records full commands, sequences, binary hashes, versions and
environment, alternating baseline/candidate order. It compares byte-exact CSV
output including interaction boundaries, energies, `Eall`, hybrid traceback and
seed starts for up to five interactions. Every process uses `--threads=1` and
`OMP_NUM_THREADS=1`. Fresh-process wall time includes accessibility, prediction,
traceback and serialization; GNU time records CPU time and peak process RSS.
The OxyS/fhlA cases use the same full sequences as the ensemble benchmark.
Additional cases exercise disabled accessibility, dense overlapping seeds,
`noLP`, bulged seeds and a no-seed control. Process startup dominates the
no-seed case, which should not be interpreted as an algorithmic speedup.
