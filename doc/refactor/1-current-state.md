# Phase 1: current state

This document records the phase-1 audit requested in discussion #232. It is
about scientific behavior and runtime cost, not a style rewrite. The audited
baseline is IntaRNA 3.4.1, commit `3b14bc0`. The correction and regression
record below is current through `87f45f1`.

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

After the phase-1 regressions and local corrections through `a975c4d`:

- the API binary passes 3,895 assertions in 31 test cases; and
- all 20 command-line golden cases pass.

The late heuristic-cell regressions and correction through `96c14b3` raise the
API gate to 3,927 assertions in 33 test cases. The subsequent output-hub
regression and corrections through `87f45f1` raise it to 3,934 assertions in
34 test cases. The complete 20-case CLI suite still passes.

The added CLI cases cover asymmetric target/query accessibility options,
seed-free `outMinPu` range filtering, and rejection of partition output from
overlapping windows. The API additions cover the local ownership, boundary,
partition-filter, seed-extension, range-splitting and ViennaRNA cases listed
below.

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

The first eight findings below were exposed by `0b20908` before their
individual source fixes. Later local findings carry focused API or CLI evidence
in the same change sequence. The corrected accessibility oracle is called out
separately because it was a test-data error, not a production defect.

| Finding | Evidence / consequence | Phase-1 regression evidence | Correction |
| --- | --- | --- | --- |
| `Interaction::operator=` was not self-assignment safe | it cleared `basePairs` before reading the same object | `x=x` preserves pairs, energy and owned seed metadata | `11ac391` returns immediately for self-assignment |
| interaction equality dereferenced asymmetric optional seeds | `seed == i.seed || *seed == *i.seed` dereferenced when only one pointer was null | seeded and unseeded interactions compare unequal in both operand orders without crashing | `f498659` dereferences only when both seeds are non-null |
| zero-capacity interaction storage dereferenced an empty reverse iterator | `maxToStore=0` reached `*storage.rbegin()` | adding an interaction increments the report count but stores nothing | `8db8b0c` gates all storage access on positive capacity |
| `NussinovHandler::getQb` returned one for an out-of-range paired interval | an invalid pair obtained the multiplicative identity instead of zero weight | `getQb(0,n)==0` on a four-base sequence | `9e2a39e` returns zero for `j>=n` |
| base-pair `getES*` included the empty monomer structure | it stored `-RT log(Q)`, while the API specifies structures containing at least one pair, i.e. `-RT log(Q-1)` | `ACGU` with one admissible pair has `ES=-1`; intervals without a pair have infinite ES | `0f901cc` removes the empty structure's unit weight |
| the noLP heuristic ensemble counted the direct stack twice | the explicit direct continuation and the `w1=w2=1` loop iteration denoted the same paths | on `GGGG/CCCC`, heuristic `Zall` does not exceed exact `Zall` | `302bd38` skips that duplicate loop term when noLP is active |
| ensemble `updateZ` bypassed site filters | its direct global and boundary-partition updates did not apply `noGUend` or `maxED` | a single terminal GU contributes zero when terminal GU is forbidden | `f9612c3` applies both filters before either partition update |
| target base-pair accessibility used query limits | the target factory read `qIntLenMax/qAccW` instead of target parameters | an asymmetric CLI case with target limit 4 and query limit 3 retains the target length-4 accessibility column | `eae6555` uses `tIntLenMax/tAccW` |
| the initial asymmetric-accessibility oracle assigned a nonzero ED to the full `ACGU` target | at the default minimum loop length this four-base target has no admissible intramolecular pair, so its unpaired probability is one | the golden target-accessibility table expects zero ED for every interval, including length 4 | `a129205` corrects the test oracle only |
| seed-free `outMinPu` filtering removed every range | `RnaSequence::lastPos` was used as the minimum resulting length when no seed was required | a `GGGG/CCCC` CLI case with `noSeed=true,outMinPu=0.5` retains and reports the full interaction | `18bdac7` uses minimum length one for seed-free prediction |
| overlapping windows could report a double-counted global partition | the same interaction can occur in more than one window, so `Zall/Eall` cannot be summed safely | a CLI case requesting CSV `Eall` in window mode expects a deterministic error | `eefbda5` rejects ensemble output and CSV columns needing `Zall` in window mode |
| seed ordering collapsed distinct equal-energy seed ranges | the set comparator used only energy and the first sequence-1 boundary | four seeds differing in any remaining boundary all survive; only an exact duplicate is rejected | `77e01d3` uses all four boundaries as lexicographic tie breakers |
| accessibility range splitting included the forbidden split position | the closed prefix ended at `i` and its length counted `i`, although singleton ED at `i` exceeded the threshold | blocked positions 3 and 7 split `[0,7]` into `[0,2]` and `[4,6]`; minimum length four retains neither | `b7a4ba1` measures `i-lastStart` and closes at `i-1` |
| ViennaRNA ES ignored the accessibility base-pair-span model | `computeES` configured `curModel.max_bp_span` but passed the unmodified model to `vrna_fold_compound` | unrestricted `GGGGAAAACCCC` has finite ES while span 3 makes full-range ES infinite | `c132a61` passes the configured model |
| ViennaRNA ensemble temporaries lacked complete scoped ownership | `computeIntraEall` leaked its allocated sequence and fold compound on normal return; `computeES` cleanup was not exception-safe | the API suite exercises span-sensitive ES plus both `getEall1/getEall2` paths; direct leak detection still belongs to sanitizer coverage | `595029f` gives Vienna allocations and fold compounds scoped deleters |
| exact noLP seed extension added independent partition factors | the stack weight and remaining left subensemble were added, violating the sum-product recurrence | with one allowed two-pair seed in a `3x3` complementary grid, the exact partition is `exp(2)+exp(3)` | `1ed5819` multiplies the stack and subensemble weights |
| empty seed-extension ranges retained prior boundary partitions | early return reset scalar optima but did not clear `Z_partition`, so predictor reuse replayed stale sites | exact and heuristic predictors run on a seeded range and then a singleton range; the second partition and boundary store are zero | `a975c4d` calls `initZ()` on both no-seed paths |
| heuristic cell incumbents leaked, and ensemble noLP direct control was over-scoped | `curCellEtotal` was not reset in `PredictorMfe2dHeuristic`, `PredictorMfeEns2dHeuristic` or `PredictorMfe2dHeuristicSeed`, so a stale strict-tie incumbent could suppress a valid cell; ensemble direct extension was also nested under singleton `noGUend` acceptance, while its empty/too-long guards continued the whole cell and skipped later loop extensions | `08d375e` records five pre-fix failures: unseeded and ensemble energies `-2` instead of `-3`, seeded energy `-2` instead of `-4`, and missing direct/bulged boundaries of `exp(3)` and `exp(4)` | `96c14b3` resets each incumbent per cell, independently guards the ensemble direct extension, and leaves later loop enumeration reachable |
| `OutputHandlerHub` could neither instantiate nor report child counts correctly | its inline `add` definition retained an obsolete two-argument signature, while `reported()` initialized zero and repeatedly took the minimum, so a nonempty hub always reported zero | `0c9cbc1` first fails to compile against the public declaration; after the forwarding repair, a two-child hub reports zero instead of the required maximum two | `ce0599b` aligns and forwards the one-argument `add`; `87f45f1` aggregates child counts with `max` |

Additional high-confidence findings are not treated as tiny local fixes because
they change aggregation or recurrence ownership and need the benchmark gates of
later phases:

- exact ensemble four-boundary storage violates the advertised practical
  two-dimensional memory bound;
- the two proven seed-extension algebra/state defects are corrected, but the
  family still needs a broader recurrence and traceback oracle for overlapping
  bulged seeds, noLP corrections and unsigned boundary arithmetic;
- window-mode `Zall/Eall` is now rejected, but aggregation across general
  multiple or overlapping input regions and tracker domains still needs an
  ownership/counting audit;
- query/target range factories mutate parser-owned vectors from `const`
  methods, repeat decomposition and can race across pair tasks;
- partition accumulation can overflow in release builds; debug-only warnings
  neither prevent nor repair the result; and
- local ViennaRNA ensemble allocations now have scoped ownership, but CLI
  factories still return raw objects and use `const_cast` cleanup, making
  exceptional and parallel paths difficult to reason about.

## Missing coverage

At the audited baseline the suite had no direct test for exact MFE, exact
seeded MFE, any ensemble predictor, seed-extension predictors, zero requested
output, or partition/filter agreement. Phase 1 now has focused evidence for
zero-capacity output, one exact-versus-heuristic noLP ensemble comparison,
terminal-GU partition filtering, seed-extension sum-product/reuse behavior,
asymmetric target/query options, seed-free range filtering, window-partition
rejection, accessibility splitting, the span-sensitive ViennaRNA ES path, and
per-cell state plus direct/bulged noLP continuation in the three affected
heuristic families, and output-hub forwarding/count aggregation.

Coverage still lacks:

- a brute-force small-instance partition oracle and direct `maxED` partition
  filter regression;
- systematic small-grid differential coverage for heuristic pruning and ties
  beyond the five focused counterexamples;
- direct exact-MFE and exact-seeded-MFE regressions, plus wider bulged-seed and
  heuristic seed-extension oracles;
- reuse tests for other stateful predictor families;
- general multiple-region and overlapping-domain aggregation tests;
- deterministic threaded-output tests;
- optimized/release invalid-input tests;
- sanitizer/leak coverage for the scoped ViennaRNA ownership change; and
- a GCC 14 plus macOS Clang portability build.

The local seed-extension sum-product and lifecycle audit is complete; broader
seed/traceback enumeration remains a phase-2 correctness gate. Phase 2 also
adds compiler/build coverage. The benchmark and differential corpus introduced
before phase 3 supplies executable-level output parity for every performance
change.

## Change record

- Created the required `refactoring` branch from `master` at `3b14bc0`.
- Added this architecture, correctness and performance audit.
- Established the clean 30-case / 3,859-assertion API and 17-case CLI baseline.
- `0b20908` added table-independent regressions for ownership, boundary,
  partition and target/query isolation defects before their fixes.
- `11ac391` and `f498659` made interaction assignment and optional-seed
  comparison safe.
- `8db8b0c`, `9e2a39e` and `0f901cc` corrected zero-capacity storage,
  invalid paired-interval weight and base-pair ES semantics.
- `302bd38` and `f9612c3` removed the noLP duplicate contribution and made
  ensemble accumulation honor site filters.
- `eae6555` isolated target accessibility limits; `a129205` then corrected
  the scientific golden oracle used by that regression.
- `18bdac7` retained seed-free ranges under `outMinPu`, and `eefbda5`
  rejected undefined global partition output from overlapping windows.
- `77e01d3` retained distinct equal-energy seeds, and `b7a4ba1` excluded
  inaccessible split positions from returned ranges.
- `c132a61` honored the configured ViennaRNA base-pair span, and `595029f`
  introduced scoped ownership for the ViennaRNA ensemble resources.
- `de135fd` exposed seed-extension algebra and reuse defects before `1ed5819`
  restored multiplication and `a975c4d` cleared stale partitions.
- Established the intermediate 31-case / 3,895-assertion API and 20-case CLI
  gate at `a975c4d`.
- `08d375e` exposed stale cell incumbents in all three affected heuristic
  families and the ensemble direct/loop control-flow defects in five focused
  sections; `96c14b3` reset the incumbents and corrected the ensemble guard.
- `0c9cbc1` exposed the unusable output-hub forwarder and zero report count;
  `ce0599b` aligned its public forwarding call and `87f45f1` returns the
  documented maximum child count.
- Established the current 34-case / 3,934-assertion API and 20-case CLI gate at
  `87f45f1`.

The discussion once calls the phase-1 document `refactor-changelog.md`; no such
file exists and the same phase otherwise consistently requires
`doc/refactor/1-current-state.md`. This file is therefore both the current-state
analysis and the reproducibility/change log.
