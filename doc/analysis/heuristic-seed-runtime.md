# Seeded heuristic runtime on OxyS/fhlA

This investigation addresses the [PR #270 review request](https://github.com/BackofenLab/IntaRNA/pull/270#issuecomment-6215321389)
against `3866c62d10846fab1adc38b937d0c9b96cc06e51`, using ViennaRNA 2.7.2.
The review reports about 21 seconds of user CPU time for seeded model P/mode H
and 2.3 seconds for seeded mode M, despite H retaining fewer right endpoints.

## Reproduced timings and improvement

Five alternating fresh-process samples per case on 2026-10-16 used GCC 14.4.0
release (`-O3`, C++23, `NDEBUG`), ViennaRNA 2.7.2, Boost 1.85.0 and an AMD Ryzen
5 7530U running Linux. No project builds or tests ran during these samples.
The [full measurements](../benchmarks/seed-extension/native-cli-results.json)
retain individual times, commands, versions, output hashes and environment;
the [source manifest](../benchmarks/seed-extension/native-source-hashes.json)
identifies both binaries, compilation settings and correctness comparisons.

| Model P case | Baseline median user s | Revised median user s | Baseline → revised median wall s | Baseline → revised median RSS KiB |
| --- | ---: | ---: | ---: | ---: |
| M, seeded | 1.42 | 1.41 | 1.42 → 1.42 | 19,732 → 19,988 |
| M, unseeded | 12.06 | 12.12 | 12.07 → 12.13 | 16,788 → 16,664 |
| H, seeded | 14.67 | 6.38 | 14.68 → 6.39 | 17,300 → 17,280 |
| H, unseeded | 0.22 | 0.22 | 0.23 → 0.23 | 30,756 → 30,660 |

The baseline reproduces the review's ordering: seeded H is 10.3 times slower
than seeded M by median user time. Cropping reduces seeded H user time by
56.5% (2.30 times faster); its five samples range from 6.37 to 6.48 seconds,
compared with 14.50 to 14.99 seconds before. Prediction output is byte-identical
for every case across all repetitions and both executables.

The revised H is still 4.5 times slower than seeded M on this sparse-seed
workload. The per-source prefix formulation remains its dominant structural
cost; achieving the seed-anchor algorithm's reuse would require a further
change in evaluation strategy while preserving heuristic ownership/selection.
The change makes no universal H-versus-M speed claim. Small changes in the
unaffected modes and in process RSS are within the observed variation; these
RSS values do not isolate recurrence storage.

Final validation passes `make tests -j2` in release and `--enable-debug`:
all seven drivers, 114 API cases and 425,806 / 425,807 assertions. The eight
exact/heuristic synthetic full dumps are byte-identical across 599,242
boundary/partition/pair entries. Twelve additional native CLI comparisons,
covering dense/sparse/cropped inputs and `noLP` on/off, have identical prediction
output and actual-pair CSVs. New independent oracles cover cropped prefixes and
right endpoints beyond their bounds, in addition to the group/outside tests
described below. macOS/Apple Clang was not tested locally.

## Why the smaller ensemble can take longer

For stack-only seeds, both native predictors enter
`PredictorMfeEns2dSeedExtension::predictStackSeeds`, but they then use different
branches of `SeededPartitionFunction::compute`:

- Exact mode M uses the left-most anchor implementation introduced in #270.
  It builds left/right tables around each irredundant seed anchor and shares
  them across a bounded batch of interaction left boundaries.
- Heuristic mode H first freezes each admitted seed's selected right endpoint
  through `prepareStackSeedExtensions`. It then runs the first-completed-seed
  recurrence separately for each eligible interaction left boundary. Within
  that source, injections sharing a selected endpoint share a right table;
  different selected endpoints require separate right tables.

Thus H is not evaluated by running the new M algorithm and dropping endpoint
combinations. It preserves its own canonical owner and selected endpoint for
each chain. Restricting the ensemble does not by itself reduce the work of
these different implementations. The exact anchor rewrite explicitly left
the heuristic recurrence unchanged.

The dominant cost in this workload is the seed-free prefix recurrence, which
previously filled the entire span for each source, including positions past
every possible seed completion. A separate instrumented baseline sample spent
13.0132 of 14.3318 seconds (90.8%) inside prefix clearing/setup/filling, versus
1.28646 seconds (9.0%) in right-group clearing/filling. Freezing native endpoint
choices took 0.0113 seconds. This attributes the slowdown to repeated prefix
work, not to expensive endpoint selection or probability output (none was
requested).

The [attribution record](../benchmarks/seed-extension/native-profile-results.json)
contains raw stderr/time output, counters, source/binary hashes and timer
definitions. Its single instrumented wall-time sample is separate from the
uninstrumented repeated performance comparison above.

There are only three admitted seeds, with internal forward-coordinate starts
`(38,5)`, `(86,79)`, `(87,80)`, all length seven. The heuristic visits 2,568
valid eligible left boundaries and fills 12,463,452 prefix rectangle cells.
Bounding prefixes by the latest contained seed endings reduces that count to
6,195,379. Its 4,653 nonempty right groups visit 12,360,097 rectangle cells;
bounding them by their first-hit injection coordinates reduces this to 40,242.
These are loop-rectangle cell counts before validity/reachability checks, not
counts of biological structures or energy evaluations. Three groups retain all
right endpoints, and each nonempty group has one injection in this workload.

The heuristic right groups also contain avoidable scans. For interaction
source `p`, a group with local endpoint `e` previously cleared and visited all
`(e.first+1)*(e.second+1)` cells starting at `p`, although its injections start
at seed ends farther to the right. Each valid cell can additionally scan up to
`(loop1+1)*(loop2+1)` predecessor offsets. This work repeats for every source
and endpoint group. Pair output repeats forward groups and reverses them, so
it also pays for those empty prefixes.

## Change and correctness argument

For each source, bound the seed-free prefix recurrence by the largest local
end coordinates of seeds fully contained in its span. A prefix can contribute
only if it later reaches an admitted seed completion. Because both coordinates
increase along every edge, it cannot pass that bounding corner and return to
any such completion. Crop the forward prefix and its matching outside pass,
while preserving the full interaction span for retained right extensions.

Store the componentwise minimum `b` of each group's first-hit injection
coordinates. Every continuation increases both coordinates, so no path from
any injection can reach a cell with either coordinate below `b`. The active
rectangle is therefore `[b.first,e.first] x [b.second,e.second]`. Its visited
area becomes `(e.first-b.first+1)*(e.second-b.second+1)`.

Crop clearing, forward traversal, unrestricted endpoint accumulation and
outside traversal to that rectangle. Bound both stack and nonstack predecessor
lookups by the same lower corner: inactive storage can contain values from a
previous group and must never be read. The minima can come from different
anchors, so this is a safe bounding rectangle rather than a claim that all its
cells are reachable.

The outside pass additionally skips zero forward states and transitions into
zero predecessor states. With nonnegative weights, such states cannot carry
any mass from this source. This preserves derivatives of the partition with
respect to actual-pair fugacities, including pending/support states under
`noLP`. The independent enumeration tests exercise sparse zero edges, zero
boundary factors, crossed lower corners, shared and distinct selected endpoints,
mixed seed lengths, asymmetric spans/loop bounds, and reused rectangles.

The first-seed ownership rule, native right-endpoint selector, no-left and
no-right branches, exact backend, and legacy gapped-seed fallback are unchanged.
For every retained contribution, the forward arithmetic order is unchanged.
The optimization adds two coordinates per group, retains bounded reusable
matrices, and does not improve the worst-case asymptotic complexity. The
per-source prefix recurrence and repetition across different endpoints remain;
this change does not promise that H will always be faster than M.

## Reproduction

The [native CLI harness](../benchmarks/seed-extension/native-cli.py) embeds the
review's exact 109 nt OxyS and 113 nt fhlA sequences. It runs all four
`--model=P` combinations of `--mode=M|H` and enabled/disabled seeds. It makes
the default seven-pair seed length explicit, fixes one thread, and otherwise
uses CLI defaults, including computed accessibility. Default runs request
ordinary prediction output without a pair-probability file.

Build the baseline and candidate with the same release toolchain/dependencies,
then run while no other project build or test is active:

```sh
python3 doc/benchmarks/seed-extension/native-cli.py \
  --binary /path/to/baseline/IntaRNA \
  --candidate /path/to/candidate/IntaRNA \
  --repetitions 5 --output-dir /tmp/intarna-native-seed-runtime
```

The harness alternates executable order between repetitions and records the
complete commands, versions, executable/harness hashes, environment, all raw
GNU time samples, and prediction-output hashes. It rejects a changed result
within a case, across repetitions or between binaries. It never compares M's
prediction with H's prediction, since their ensembles differ. GNU time's
elapsed/user/system times and peak process RSS cover the whole process,
including accessibility, prediction, traceback, and output. Timing resolution
is 0.01 seconds. CPU affinity/frequency and unrelated workstation activity are
not controlled; small changes in the unaffected modes should be treated as
noise.

Optional `--pair-output` measures actual-pair CSV output and also checks its
bytes. Use a separate output directory for each experiment. The synthetic
[kernel harness](../benchmarks/seed-extension/README.md) separately checks full
hybrid coefficients, complete boundary factors, the partition and raw pair
numerators, including changes too small to appear in rounded CLI predictions.

For attribution, `native-profile.py` instruments a **disposable configured
source/build copy of the baseline** with steady-clock stage timers and counts
the old and cropped group rectangle areas. The ordinary performance comparison
must use uninstrumented executables. For example:

```sh
python3 doc/benchmarks/seed-extension/native-profile.py \
  --source /path/to/unmodified-3866c62-checkout \
  --build /path/to/disposable-configured-baseline-copy \
  --output /tmp/intarna-native-profile-source
make -C /path/to/disposable-configured-baseline-copy -j2
python3 doc/benchmarks/seed-extension/native-cli.py \
  --binary /path/to/disposable-configured-baseline-copy/src/bin/IntaRNA \
  --repetitions 1 --output-dir /tmp/intarna-native-profile
```

The script saves its exact patch and before/after source hashes. Stage timings
and counters appear in stderr (`#profile_predictor`, `#profile_kernel`).
These timers include instrumentation overhead and only attribute one workload;
they are separate from the repeated end-to-end performance samples. The prefix
timer includes source workspace clearing and validity/seed-end lookup as well
as first-hit recurrence filling. The right timer includes group clearing and
extension filling; neither timer includes all boundary publication overhead.
