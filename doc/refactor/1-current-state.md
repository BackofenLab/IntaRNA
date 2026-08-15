# Phase 1: current state

This document records the phase-1 audit requested in discussion #232. It is
about scientific behavior and runtime cost, not a style rewrite. The audited
baseline is IntaRNA 3.4.1, commit `3b14bc0`.

## Reproducible baseline

The source was built as an optimized OpenMP build with GCC 13.3, Boost 1.85 and
ViennaRNA 2.7.2. Before any source change:

- the API binary passed 3,859 assertions in 30 test cases;
- all 17 command-line golden cases passed; and
- the golden cases covered only the table-independent base-pair energy model
  with disabled accessibility.

Thus the old suite is a useful output-stability gate, but it is not broad enough
to establish correctness of the ensemble predictors, ViennaRNA paths,
parallel/window reduction, or several public boundary cases.

## Architecture and data flow

`src/bin/IntaRNA.cpp` owns one invocation. `CommandLineParsing` parses the
selected personality and works as a factory for sequence, accessibility,
energy, seed/helix, predictor, tracker and output objects. Query accessibility
is cached across targets. For each target/query pair the program computes or
loads accessibility, constructs an interaction-energy facade, decomposes
allowed ranges into optional overlapping windows, runs one predictor per
window, and merges interactions in an `OutputHandlerInteractionList`.

The library keeps the biological pieces modular:

| Family | Responsibility | Scientific invariant |
| --- | --- | --- |
| `RnaSequence`, `IndexRange*`, `Interaction*` | sequence and inclusive coordinate value types | external/internal and reversed-query coordinates remain bijective and in range |
| `Accessibility*` | opening energy ED for a sequence interval | `ED=-RT log(Pu)`; unavailable intervals return the documented upper bound |
| `InteractionEnergy*` | initiation, loops, ends, dangles, ED and Boltzmann conversion | total energy uses integer centi-kcal/mol without arithmetic on infinity sentinels |
| `SeedHandler*`, `HelixHandler*` | admissible seeds/helices and traceback | preprocessing bounds and traceback describe the same structure |
| `PredictorMfe*` | exact/heuristic MFE dynamic programs | every reported boundary has a valid recurrence and traceback |
| `PredictorMfeEns*` | interaction partition functions and ensemble representatives | every valid interaction contributes exactly once with multiplicative independent weights |
| `OutputHandler*`, `PredictionTracker*` | filtering, reduction and rendering | filters apply consistently to both reported sites and their requested ensemble |

The second sequence is exposed in biological orientation and wrapped by
`ReverseAccessibility` for the antiparallel dynamic programs. This convention
is performance-sensitive: offsets can be cached, but changing the convention
would invalidate almost every recurrence and traceback.

### Numeric model

MFE values use `E_type=int` in centi-kcal/mol. `E_INF` and `E_MAX` are finite
sentinels chosen to leave arithmetic headroom. Partition values use `double`
unless the optional quadmath build is enabled. A hybrid path with base pairs
`(i1,i2),...,(j1,j2)` has the schematic energy

```
E = E_init + sum(E_interLeft) + ED1(i1,j1) + ED2(i2,j2)
    + E_dangleLeft + E_dangleRight + E_endLeft + E_endRight + E_add
```

and weight `exp(-E/RT)`. Exact energy ties intentionally use integer equality;
loop order and tie-breaking must remain deterministic in optimized code.

## Runtime structure and bottlenecks

The default personality uses ViennaRNA accessibility, a bulge-free seed and
`PredictorMfe2dHeuristicSeedExtension`. Important costs are:

1. ViennaRNA accessibility preprocessing and storage of banded ED tables.
2. Two-dimensional predictor matrices, repeatedly resized and initialized for
   windows, seeds, and right boundaries.
3. Nested internal-loop scans in which virtual offset, complementarity, ED and
   energy calls are repeated for nearby cells.
4. Exact ensemble prediction stores a hash entry keyed by all four site
   boundaries. That can grow as `O(n1^2*n2^2)` even though the recurrence matrix
   itself is two-dimensional.
5. Output/top-k and partition updates enter global OpenMP critical sections;
   parallel work is selected at one outer dimension rather than represented as
   independent pair/window tasks.
6. Automatic accessibility range decomposition is performed from mutable
   parser-owned vectors and can be repeated for sequence pairs.

The first optimization work should therefore remove unnecessary state and
calls from existing recurrences before changing their mathematics.

## Confirmed defects and regressions

The following findings were traced in this source tree. Tests are added before
fixes, as required by phase 1.

| Finding | Evidence / consequence | Phase-1 regression |
| --- | --- | --- |
| `Interaction::operator=` is not self-assignment safe | it clears `basePairs` before reading the same object | preserve pairs, energy and seed metadata across `x=x` |
| interaction equality dereferences asymmetric null seed pointers | `seed == i.seed || *seed == *i.seed` dereferences when only one side is null | seeded and unseeded interactions compare unequal without crashing |
| zero-capacity interaction storage dereferences an empty reverse iterator | `maxToStore=0` reaches `*storage.rbegin()` | count a report but retain no interaction |
| `NussinovHandler::getQb` returns one for an out-of-range paired interval | an invalid pair obtains the multiplicative identity instead of zero weight | `getQb(i,n)==0` |
| base-pair `getES*` includes the empty monomer structure | it stores `-RT log(Q)`, while the API specifies structures containing at least one pair, i.e. `-RT log(Q-1)` | a four-base sequence with one possible pair has `ES=-1`, not `-log(1+e)` |
| noLP heuristic ensemble counts the direct stack twice | the explicit direct continuation and the `w1=w2=1` loop iteration denote the same path | on `GGGG/CCCC`, heuristic `Zall` must not exceed exact `Zall` |
| ensemble `updateZ` bypasses `noGUend` and `maxED` | its direct `Zall`/boundary-map update does not use the filters in `PredictorMfe::updateOptima` | a single terminal GU contributes zero when terminal GU is forbidden |
| target base-pair accessibility uses query limits | target factory reads `qIntLenMax/qAccW` instead of target parameters | asymmetric CLI target-accessibility table retains the target maximum length |

Additional high-confidence findings are not treated as tiny local fixes because
they change aggregation or recurrence ownership and need the benchmark gates of
later phases:

- exact ensemble four-boundary storage violates the advertised practical
  two-dimensional memory bound;
- seed-extension ensemble code has early-return state reuse, additive where
  independent partition factors must be multiplicative, noLP duplication and
  unsigned boundary risks;
- windowed/global partition and tracker aggregation can count overlapping
  domains more than once;
- query/target range factories mutate parser-owned vectors from `const`
  methods, repeat decomposition and can race across pair tasks;
- partition accumulation can overflow in release builds; debug-only warnings
  neither prevent nor repair the result; and
- CLI object ownership uses raw factories and `const_cast` cleanup, making
  exceptional and parallel paths difficult to reason about.

## Missing coverage

The original suite has no direct test for exact MFE, exact seeded MFE, any
ensemble predictor, seed-extension predictors, zero requested output, or
partition/filter agreement. It also lacks:

- a brute-force small-instance partition oracle;
- reuse tests for stateful predictors;
- asymmetric target/query option tests;
- multiple-region and overlapping-window aggregation tests;
- deterministic threaded-output tests;
- optimized/release invalid-input tests;
- sanitizer coverage; and
- a GCC 14 plus macOS Clang portability build.

Phase 1 adds focused deterministic cases for the confirmed local defects.
Phase 2 adds compiler/build coverage. The benchmark and differential corpus
introduced before phase 3 supplies executable-level output parity for every
performance change.

## Change record

- Created the required `refactoring` branch from `master` at `3b14bc0`.
- Added this architecture, correctness and performance audit.
- Established the clean 30-case / 3,859-assertion API and 17-case CLI baseline.
- Added table-independent regressions for ownership, boundary, partition and
  target/query isolation defects before their fixes.

The discussion once calls the phase-1 document `refactor-changelog.md`; no such
file exists and the same phase otherwise consistently requires
`doc/refactor/1-current-state.md`. This file is therefore both the current-state
analysis and the reproducibility/change log.
