# Seed-extension ensemble memory audit (issue #260)

The additional memory in default `--model=P --mode=H` comes primarily from
retaining a partition for every considered interaction site. Sharing the same
recurrence dimensions with `--model=X` does not imply the same total storage.
The fix releases completed site partitions during seed traversal and stores a
smaller key for the partitions that remain active. Predictions and partition
values are preserved. The measurements below distinguish total process peak RSS
from the sizes of individual recurrence matrices.

## Source audit

The reference is `master` at `c5f1d51ce0f8d473b9a3e23ca2e0d410d92d77d7`.
The forward branch is `feat/258-bpsvg` at
`c0ef61f1603ffc41c2ae4f8bbcf98b56027880e9` (PR #259's October 12 revision).

The CLI routes seeded model X/mode H to
`PredictorMfe2dHeuristicSeedExtension`, and seeded model P/mode H to
`PredictorMfeEns2dHeuristicSeedExtension`. Their common costs include accessibility,
seed preprocessing, local left/right extension tables and ranked-output storage.
They differ in two ways:

1. Model X uses integer `E_type` energies. The default ensemble `Z_type` is a
   double, so equal-size partition matrices need twice the payload bytes on the
   measured platform. Optional multiprecision increases this difference.
2. `PredictorMfeEns::updateZ()` also aggregates hybrid partition sums by the
   four-coordinate boundary `(i1,j1,i2,j2)` in `Z_partition`. Each entry needs
   the boundary, a partition value, hash-node storage and bucket storage. The
   seed-extension predictor visits overlapping contributions to the same site;
   the site cannot be ranked correctly until those contributions have been
   summed. Model X can retain its current minimum and discard worse candidates.

For example, two seeds contributing to the same site require a sum of their
accepted terms, rather than selecting the better individual term. Reducing this
map to one candidate per left start would change the ensemble objective.
`--outNumber=1` limits reporting, not the number of site sums needed to select
that result. The legacy O(n²) space descriptions account for recurrence tables,
but omit the potentially much larger boundary map.

The exact **unseeded** 2D predictor already streams complete boundaries through
`updateCompleteZ()`. Seed-extension updates are partial, so applying that path to
each individual seed contribution would be incorrect.

PR #259 retains the same map on ordinary seeded heuristic runs. Its pair matrix,
outside calculations, signed seed-extension traces and SVG document are optional
and do not explain the issue's commands, which request no probability output.
The patched forward branch was built and tested separately, including its
probability regressions. The benchmark also compares normal CSV predictions
between the master and forward implementations.

The [tested forward-branch patch](ensemble-memory-probability-port.patch) can be
applied with `git apply --unidiff-zero` at the PR #259 revision above. It preserves that branch's
full heuristic seed-start ranges (`range_size1-1`, `range_size2-1`); master retains
its existing nominal seed-length limits. The exact seeded probability backend's
separate regional finalization path is unchanged. This PR targets master; the
port is supplied separately so the SVG review and memory change remain independent.

## Finalization without changing recurrences

The seeded predictors collect starts in their existing handler order and compute
a suffix minimum of the remaining target seed starts. Before processing a seed,
any buffered site with `j1` below this minimum is complete: no remaining seed
can fit inside that site's target boundary. Its final partition updates the
ranked candidates and is then released. Its contributions have already entered
`Zall`, so finalization must not increment `Zall` again.

Partitions are bucketed by target right boundary. This avoids repeatedly scanning
the whole active map, and makes `j1` redundant in each key. On the measured 64-bit
build, keys shrink from four 8-byte coordinates to three. The final bucket scan
consumes every remaining boundary at `reportOptima()`.

The recurrence, per-site addition order, global partition addition order, site
filters, seed-overlap corrections and virtual `updateZ()` call are preserved.
`initZ()` resets all storage for predictor reuse and empty subsequent regions.
Tracker-enabled predictions keep the previous deferred map, preserving callback
timing and ordering. Completed boundaries are sufficient for ranked and
nonoverlapping suboptimal selection; that selection uses the existing
`mfe4leftEnd` structure rather than the full partition map.

The suffix minimum is necessary because computed bulged-seed handlers can advance
query-first, whereas stack/explicit handlers advance target-first. It prevents
premature finalization for either traversal, without changing summation order.
Early release is most effective with target-ordered seeds and a target span much
shorter than the sequence. Query-first traversal, large allowed spans, very sparse
seeds, or requested trackers can reduce the benefit. The additional seed-start
list and right-boundary bucket descriptors are linear storage; the ensemble is
not expected to match model X's memory in every configuration.

## Measurements

132 timed runs plus 40 correctness/warmup runs, pinned to CPUs 0–3. Each row
reports the median of three repetitions (six for GcvB-32/model P). All final
runs produce identical CSV field values across the four binaries; parallel row
ordering is ignored. Normal model X controls are included in the JSON.

Peak RSS, MiB (KiB / 1024):

| Input | Master X | Master P before | Master P after | Forward P before | Forward P after |
| --- | ---: | ---: | ---: | ---: | ---: |
| OxyS-one | 18.67 | 18.68 | 18.63 | 18.68 | 18.45 |
| OxyS-eight | 23.16 | 25.86 | 25.12 | 26.35 | 25.15 |
| GcvB-32 | 24.39 | 87.34 | 78.49 | 87.79 | 78.42 |
| GcvB-32-noLP | 24.62 | 43.17 | 41.09 | 43.86 | 41.32 |
| dense256x24 | 14.00 | 17.49 | 18.12 | 17.61 | 18.38 |

Model P wall time, seconds:

| Input | Master before | Master after | Forward before | Forward after |
| --- | ---: | ---: | ---: | ---: |
| OxyS-one | 0.140 | 0.140 | 0.140 | 0.140 |
| OxyS-eight | 0.780 | 0.780 | 0.770 | 0.780 |
| GcvB-32 | 10.195 | 9.895 | 10.305 | 10.025 |
| GcvB-32-noLP | 3.610 | 3.620 | 3.640 | 3.660 |
| dense256x24 | 31.900 | 31.540 | 31.790 | 32.140 |

The GcvB-32 cases reproduce the large gap relative to model X on both source
branches. The reduction comes from the boundary storage changes, not from
changing the ensemble or skipping work. The short OxyS cases are dominated by
startup/accessibility and centisecond timing granularity. RSS differences of a
few hundred KiB at the baseline footprint also include allocator/process noise.
The dense stack-only case has little boundary-map storage to release; the seed
list and bucket descriptors can therefore slightly increase its peak RSS.
See the raw results for ranges; these measurements do not promise universal
speedups or equality with the model X footprint.

## Correctness and integration checks

- GCC 14.4, C++23, default double partitions and Kokkos mdspan: full release and
  debug `make tests -j2`, all six master suites pass; 76,009 API assertions in
  68 cases in both builds.
- Full tests after porting the fix to the current probability implementation:
  all seven suites pass, including 101,247 API assertions in 97 cases and the
  CSV/SVG probability regressions.
- A buffered tracker run is the independent storage-path reference for the
  streaming regression. Exact partition equality and complete ranked-output
  equality hold for exact/heuristic seed extension, LP/noLP, computed stack,
  computed bulged and explicit seed handlers, cropped regions and predictor
  reuse. Tracker callbacks remain active; completed maps are empty after reporting.
- All 69 installed master public headers compile independently. An installed
  pkg-config predictor consumer compiles, links and predicts with the documented
  Easylogging initialization and installed `libeasylogging`.
- `git diff --check`, UTF-8/LF checks and distribution inclusion are checked.

A preliminary four-thread CSV comparison failed once before differing output
was retained. The saved warmup CSVs matched. Eight targeted repeats and a further
28 interleaved checks then matched exactly; the cause of that initial observation
could not be established. The final measurements require output equality on
every run, and the supplied script saves any future mismatch for diagnosis.

Apple Clang, multiprecision builds and the full 4,319-target benchmark screen were
not run locally. The sampled benchmark results describe these inputs and this
machine; they do not establish a universal runtime or memory bound. Requested
pair-probability output and tracker-heavy workloads have additional storage costs.

## Reproduction

Build the before/after revisions with the same Autotools configuration and
dependencies. Preserve both binaries before modifying the checkout. Use
`--enable-debug` in a separate source/build tree for debug validation.

The fixtures come from [BackofenLab/IntaRNA-benchmark](https://github.com/BackofenLab/IntaRNA-benchmark),
revision `2b6a23547039f5a7d883706f12674109d046dec7`: the first 1, 8 and 32 E. coli
target records, OxyS and GcvB queries. An additional synthetic input has 256 G
target nucleotides and 24 C query nucleotides. The script writes the selected
FASTA fixtures, exact commands, raw CSVs, timings, summaries and binary hashes.
GNU time reports wall/user/system seconds and peak RSS in KiB.

```sh
taskset -c 0-3 python3 doc/benchmark-ensemble-memory.py \
  --before /path/to/master-before --after /path/to/master-after \
  --feat-before /path/to/probability-before --feat-after /path/to/probability-after \
  --benchmark-dir /path/to/IntaRNA-benchmark \
  --outdir /path/to/new-results --repeats 3
```

Select four available CPUs for the affinity mask on another machine. Each
configuration has a correctness/warmup run, followed by three deterministically
interleaved fresh-process repetitions (six for the GcvB-32 ensemble audit). The normal X/P output fields and `Zall`
must match across all four binaries; parallel CSV rows are sorted to ignore
thread completion order. Model X's reported `Zall` is only an output-parity check,
not an exact ensemble partition. No compilation runs during the final timing
pass. [The recorded measurements](ensemble-memory-benchmark.json) include all
model X controls and runtime ranges, rather than only the improved cases.
