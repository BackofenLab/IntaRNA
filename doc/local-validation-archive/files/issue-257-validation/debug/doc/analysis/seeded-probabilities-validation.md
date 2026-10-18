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

The full release/debug results, installed-header checks and native measurements
are recorded below after the final validation run.
