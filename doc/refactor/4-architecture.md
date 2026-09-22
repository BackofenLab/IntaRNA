# Phase 4: exact ensemble boundary streaming

This document records the first deliberately architecture-changing performance
candidate from the roadmap in `2-c++23.md`. The change is narrow: it streams
complete four-boundary partitions produced by the exact, non-seed
`PredictorMfeEns2d` instead of retaining every one in `Z_partition` until
reporting. It does not change the recurrence, energy model, filters, command
line, or output schema.

The scientific parent executable is
`8bd1676d1e38ec6461c60e1b832ff3ace9ff77bc`. The final candidate source is
`91840835820d8e29b5c9dd9543104fec73d6f7a0` on
`refactor/phase4-ensemble-streaming`. Source and regression validation are
complete. The dedicated `phase4-stream` performance batch establishes a large
memory reduction and a repeatable runtime improvement with identical output,
so this candidate is accepted.

## Problem and bounded scope

For each accessible, complementary right boundary `(j1,j2)`, the exact 2D
ensemble predictor fills its two-dimensional `hybridZ` matrix for all eligible
left boundaries `(i1,i2)`. Before this change, every accepted result was also
stored under the four-coordinate key `(i1,j1,i2,j2)` in an unordered map. The
map was traversed only after the recurrence had finished, in order to convert
each boundary partition into an output candidate.

The number of valid four-boundary sites can grow as
`O(n1^2 * n2^2)`. Retaining one hash entry per site can therefore dominate the
`O(n1 * n2)` recurrence matrix even though the exact non-seed recurrence has
already finished a site's complete partition when the cell is reported. The
candidate introduces `updateCompleteZ()` for that one proven-complete call
site. With no prediction tracker, an accepted boundary is handed directly to
`PredictorMfe::updateOptima()` and never enters `Z_partition`.

The scope is intentionally limited:

- only `PredictorMfeEns2d::fillHybridZ()` marks its result as complete;
- exact seed extension and heuristic predictors continue to call the ordinary
  buffered `updateZ()` path because a boundary can receive multiple partial
  contributions there;
- tracker-enabled exact predictions retain the map and the established
  deferred tracker-callback behavior; and
- recurrence matrix allocation, traceback, top-k/overlap handling, and output
  ownership are unchanged.

Thus the tracker-free exact non-seed path removes the four-boundary map from
its live storage. This is a storage-ownership change, not a new recurrence or a
general streaming framework for all predictor families.

## Once-only boundary argument

The direct handoff is valid because a complete exact non-seed boundary is
finalized exactly once:

1. `predict()` visits each admissible right-boundary pair `(j1,j2)` once in its
   two nested loops and invokes `fillHybridZ(j1,j2,0,0,true)` once for that
   pair.
2. For fixed `(j1,j2)`, `fillHybridZ()` visits every eligible `(i1,i2)` cell
   once. It initializes `hybridZ(i1,i2)`, adds every valid internal-loop
   decomposition in the existing order, and calls the complete-boundary hook
   only after those additions are finished.
3. The tuple `(i1,j1,i2,j2)` identifies exactly one combination of these loop
   iterations. No second exact non-seed call contributes another disjoint
   structure set to that same tuple.
4. Zero/invalid partitions and filtered sites are rejected before either
   buffering or output update. Every remaining tuple was formerly one unique
   map entry and is now one direct `updateOptima()` call with the same
   ED-free boundary partition.

The proof does not apply to seed-extension or heuristic call sites, where the
same boundary may be assembled from separate contributions. They deliberately
remain buffered so `updateZ()` can aggregate those contributions before
ranking or overlap selection.

## Scientific and numerical equivalence

The former prefix of `updateZ()` was factored into
`addPartitionContribution()`. Its evaluation sequence is preserved:

1. reject zero partitions or an already infinite `Zall`;
2. apply the terminal-GU filter;
3. apply the two accessibility-energy filters;
4. add or remove ED and other non-hybrid terms according to `isHybridZ`;
5. perform the existing debug overflow check; and
6. add the with-ED value to `Zall`.

The recurrence loops and their multiplication/addition order are untouched,
and `Zall` receives accepted values at the same point and in the same traversal
order as in the parent. No floating-point reassociation, alternate partition
representation, or changed infinity test is introduced. Since the exact
non-seed path produces one value per four-boundary key, its old map path did
not need a same-key partition addition; direct conversion uses the same
`energy.getE(partZ_noED)` input.

What changes is the time at which a completed boundary reaches top-k/overlap
selection: recurrence-discovery order replaces a later unordered-map
traversal. This does not alter the scientific value of a site, and the focused
tests compare the resulting selected interactions under every overlap policy.
CLI golden tests provide the executable-output gate, and the dedicated
benchmark confirms identical raw and canonical output for the measured exact
ensemble workload.

## Modularity, trackers, and exception safety

`PredictorMfeEns::updateZ()` remains the virtual extension hook. The final
implementation has `updateCompleteZ()` establish a scoped complete-value mode
and then invoke virtual `updateZ()`. A subclass therefore still observes or
alters every update; an override that delegates to the base implementation
inherits the streaming behavior. This restores the modular hook that an
earlier direct-call version of the candidate would have bypassed.

The scoped mode is managed by a small RAII guard. It saves the previous mode,
sets complete-value mode only around the virtual call, and restores the saved
state both on normal return and during exception unwinding. Consequently, an
overriding `updateZ()` that throws cannot leave later ordinary updates in the
streaming path.

When `predTracker` is non-null, base `updateZ()` deliberately ignores the
streaming shortcut and writes the ED-free value to `Z_partition`. The existing
`reportOptima()` / `updateOptimaUsingZ()` traversal then invokes the tracker in
its established deferred phase. This fallback preserves tracker callback
timing and aggregation rather than silently redefining the tracker interface.
Ordinary calls from seeded and heuristic predictors likewise enter the map
because they are made outside the complete-value scope.

## Validation and compatibility

The final candidate at `91840835820d8e29b5c9dd9543104fec73d6f7a0` has the
following recorded C++23 validation:

| Gate | Result |
| --- | --- |
| Focused exact-ensemble regression | 54 assertions passed |
| Complete API suite | 4,319 assertions in 36 test cases passed |
| Command-line golden suite | all 20 cases passed |

The focused regression establishes that:

- the tracker-free exact predictor leaves `Z_partition` empty;
- the tracker fallback retains nonempty boundary storage and invokes its
  callback once per buffered boundary;
- `Zall`, the output partition, report count, and reported interactions agree
  between streaming and buffered paths for all four overlap policies;
- an overriding virtual `updateZ()` still receives complete-boundary calls;
  and
- a deliberately throwing override cannot leak complete-value mode into a
  later ordinary buffered update.

No CLI option, default, accepted input, error contract, output column, or
format is changed. The 20 golden cases therefore remain the immediate binary
program stability check. The accepted benchmark comparison independently
confirms byte-identical output for its exact-ensemble biological case.

Adding the complete-value mode flag changes the C++ object layout, and can
change the object size, of `PredictorMfeEns` and its derived classes. Although
the source-level predictor interfaces and virtual hook are retained, this
branch is not binary-layout compatible with client objects built against the
parent headers. Downstream C++ users and derived predictors must perform a
clean rebuild and relink against the new library. This ABI rebuild requirement
does not affect the IntaRNA command-line interface or its output format.

## Performance result and decision: accepted

Batch `phase4-exact-ensemble-streaming-gcc14-v1` in the private
[benchmark repository](https://github.com/BackofenLab/IntaRNA-refactor-benchmark)
compares the exact parent and candidate commits with suite `phase4-stream`.
The harness revision is
[`cdfea210c8b6a991f8e577f253447ee27d925e99`](https://github.com/BackofenLab/IntaRNA-refactor-benchmark/commit/cdfea210c8b6a991f8e577f253447ee27d925e99),
and the generated report, result JSON, and compressed raw outputs are stored in
evidence commit
[`30730f0d0932d31019d8461087f921c7e5497ea5`](https://github.com/BackofenLab/IntaRNA-refactor-benchmark/commit/30730f0d0932d31019d8461087f921c7e5497ea5).
The benchmark repository's validation Action for this evidence is green.

The `oxys-fhla-ensemble-exact-80` case is an exact, unseeded 80-by-80
interaction ensemble run with one thread. Both revisions used the same pinned
GCC 14.4, Boost 1.85 and ViennaRNA 2.7.2 environment, fingerprint
`52fb1f3e9e0f6e3806593ed89ca02aea5ea8fa2f2fa9544e9dfd53f07bdbd68e`.
Measurements comprise one warm-up and three fresh timed processes per
revision. Both builds and test gates passed, all timed processes exited zero,
and each revision was stable across repetitions.

| Metric | Parent `8bd1676` | Candidate `9184083` | Candidate change |
| --- | ---: | ---: | ---: |
| Wall time (s) | 3.110000 +/- 0.098489 | 2.346667 +/- 0.023094 | 1.325x faster; 24.5% lower |
| User time (s) | 3.050000 +/- 0.088882 | 2.340000 +/- 0.026458 | 23.3% lower |
| System time (s) | 0.050000 +/- 0.010000 | 0.000000 +/- 0.000000 | zero at timer resolution |
| Peak RSS mean +/- SD; max (KiB) | 89,184.0 +/- 156.6; 89,348 | 7,594.7 +/- 73.9; 7,680 | 91.5% lower mean; 11.74x lower |

Raw and canonical stdout are identical between parent and candidate; both
SHA-256 values are
`463dccb9f0c8d9c3490f35b2fb62d207bb018f826a3fe0577ad465e57165534e`.
The 0.763333-second wall-time reduction is well beyond the observed dispersion,
while the mean peak-RSS reduction is 81,589.3 KiB. Together with the focused,
full-API, and CLI gates above, these results accept the candidate: it removes
the targeted memory growth and improves runtime without a scientific-output or
interface tradeoff in scope.

## Deferred architecture roadmap

All other architecture-changing roadmap items remain explicitly deferred.
In particular, this candidate does not attempt to:

- stream seeded or heuristic boundary partitions without first proving their
  contribution ownership;
- change tracker callback ordering or remove tracker-required buffering;
- introduce scaled or log-domain partition arithmetic;
- define global partition ownership across multiple regions or overlapping
  windows;
- replace mutable run configuration with an immutable task graph;
- specialize virtual predictor/energy kernels; or
- introduce a shared tiled, wavefront, SIMD, GPU, or other accelerator engine.

Those changes require independent scientific proofs, regressions, benchmark
batches, and review. They must not be folded into this narrowly attributable
memory candidate.

## Conclusion

Phase 4 is complete for the exact non-seed ensemble boundary-streaming scope.
The implementation at `91840835820d8e29b5c9dd9543104fec73d6f7a0` is accepted
in the open [IntaRNA pull request #237](https://github.com/BackofenLab/IntaRNA/pull/237).
The remaining architecture roadmap stays deferred to separately proved and
benchmarked changes.
