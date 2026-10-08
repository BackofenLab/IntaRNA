# Implementation plan: seeded base-pair probabilities

Status: design only. This document addresses [issue #257](https://github.com/BackofenLab/IntaRNA/issues/257);
it does not implement the recurrences, proposed APIs, or proposed output options.

The code assessment uses IntaRNA revision
[`5454030a04bbcde79b2e75a0614dadf6f1115e6f`](https://github.com/BackofenLab/IntaRNA/tree/5454030a04bbcde79b2e75a0614dadf6f1115e6f).
The theory assessment uses IntaRNA-probabilities branch
`issue-10-fixed-seed-recursions` at
[`f62778504c8abbec8fc9894da16a357b60707680`](https://github.com/BackofenLab/IntaRNA-probabilities/tree/f62778504c8abbec8fc9894da16a357b60707680),
particularly [Step 5: fixed-length seeds](https://github.com/BackofenLab/IntaRNA-probabilities/blob/f62778504c8abbec8fc9894da16a357b60707680/theory/fixed-length-seeds.md)
and [Step 4: seed extension](https://github.com/BackofenLab/IntaRNA-probabilities/blob/f62778504c8abbec8fc9894da16a357b60707680/theory/seed-extension.md).
Recheck these integration points against the implementation branch before coding.

## 1. Decision and scope

Use the maximal-stack formulation as the leading implementation candidate for
stack-only seeds, with the optimized outside pass specified below. It has a
direct uniqueness argument, scalar seeded/seed-free state families, and simple
whole-interaction `noLP` handling. Do not port the reference implementation's
extra loop over every marked pair in every run.

Compare it with the specialized strict-seed suffix automaton before selecting
the production default. Fewer explicit state families do not prove faster
execution: a long stack can exceed the seed length substantially. Both are
valid stack-only algorithms. A general seed-pattern automaton is a separate,
larger extension for bulged seeds; it is not needed for the initial scope.

Initially support actual-pair probabilities for the exact seeded
`--model=P --mode=M` path with a stack-only handler. Preserve other prediction
modes; reject requests for the new output when their predictor cannot provide
the stated ensemble. The normal exact seeded path should use the selected new
backend whether or not pair probabilities are requested, so enabling output
does not change its partition function.

Mixed-length and singleton **stacked** explicit seeds are included. Their run
predicate follows the same uniqueness proof, but extends the fixed-length
Python API and therefore needs independent tests. Computed seeds retain their
existing allowed lengths. Do not extend the computed-seed CLI to length one.

## 2. Ensemble and output contract

Let an interaction be its complete ordered intermolecular pair chain. Include
it if it contains at least one admitted seed, satisfies interaction constraints,
and lies wholly in one searched target/query region pair. Count the interaction
once, irrespective of the number or overlap of its seed witnesses.

The searched ensemble is the union over the selected or automatically derived
nonoverlapping region pairs. `qRegionLenMax`, `tRegionLenMax`, and `outMinPu`
can restrict this union before prediction. Interactions spanning excluded
region boundaries are not included. Preserve model P's rejection of every
nonzero `windowWidth`; overlapping windows cannot be summed as disjoint sets.

For the full IntaRNA Boltzmann weight $w(I)$, compute

$$
Z=\sum_{I\in\mathcal I}w(I),\qquad
M_k=\sum_{I\in\mathcal I:\,k\in I}w(I),\qquad P_k=M_k/Z.
$$

These are conditional probabilities given an allowed seeded interaction;
there is no added unbound-state weight. Merge raw region numerators and
denominators before normalizing. Ranked-output controls (`outNumber`, energy
reporting cutoffs, overlap selection, and `outPerRegion`) do not select the
probability ensemble. Structural `noLP`, endpoint `noGUend`, accessibility
limits, seed admission, loop limits, and both interaction span limits do.

Actual pairing is different from site coverage. Preserve `spotProb`,
`qSpotProb`, and `tSpotProb`: their existing callbacks describe boundaries and
their writers accumulate coverage, not the interior pair chain. No selected-pair
"none of these pairs" probability is implied by $1-\sum_k P_k$, since pair
events can coexist. Per-nucleotide pair profiles, if added later, are row or
column sums of the actual-pair matrix.

## 3. Establish a numerical policy before the DP

Use `Z_type` for partition values, adjoints, raw pair masses, and normalization.
Use existing energy units and native energy functions. "Exact" refers to the
enumerated ensemble and disjoint recurrence; ordinary floating-point rounding
still occurs and must be measured against an independent higher-precision
reference. It does not mean exact real-number arithmetic.

The new nonnegative backend must not use `Z_equal(x,0)` to prune positive DP
or boundary contributions. At the reviewed revision, that macro uses the
absolute threshold `sqrt(float::min())`, approximately $1.08\times10^{-19}$.
`PredictorMfeEns::addPartitionContribution()` currently applies that test.
Calling it unchanged while reversing all positive contributions would make
the denominator inconsistent with the pair numerators.

Introduce a scoped exact-zero accumulation policy for the new path. Preserve
the virtual `updateZ()` dispatch from `updateCompleteZ()`, legacy tracker
callback behavior, and exception-safe restoration of scoped state. Retain the
existing policy for unrelated predictors. Observing/delegating overrides remain
supported; an override that changes boundary weights must also provide the same
objective coefficients to the reverse calculation. Arbitrary reweighting inside
an old override is not automatically a supported pair-probability backend.

Use one boundary-admission/weight contract for forward aggregation and outside
initialization. Carry its raw denominator through predictor, region collector,
and final output exactly once per disjoint contribution. There must not be an
independent energy-to-weight reconstruction for the new probabilities.

Check energy sentinels before exponentiation: a forbidden transition or boundary
is a structural zero. For finite admissible energies and nonnegative arithmetic,
check NaN and infinity, exponential and positive-product underflow to zero,
overflow in sums and region merges, and underflow/overflow in final $M/Z$.
Do not rely on `Z_isINF` alone; its comparison does not identify NaN. A computed
positive subnormal value may need scaling or a numerical-range failure if its
precision cannot meet the declared validation tolerance. Do not claim that
adding a small positive term to a large sum avoids ordinary rounding loss.

Distinguish successful nonempty, successful empty, and failed computations.
Only a successfully empty ensemble produces undefined probabilities (`NA`).
Numerical failure must be reported as an error, not as emptiness or a valid
matrix of zeros. Do not clamp substantial probability or mass inconsistencies.
Define any roundoff-only endpoint correction through the tested error tolerance.

Checked arithmetic with explicit range failure is the initial strategy. If
representative supported workloads fail, introduce a tested scaling scheme
before release, with common scales for forward/outside quantities and raw
region merging. A global replacement of all IntaRNA arithmetic is out of scope.

Follow the exact-zero policy through serialization too. Raw
`OutputHandler::incrementZ()` currently adds small values, but `Eall`,
`EallTotal`, `Zall`, and `P_E` output in `OutputHandlerCsv`, ensemble output in
`OutputHandlerEnsemble`, and ensemble-energy reporting in `OutputHandlerText`
use approximate-zero tests. Carry the new backend's partition status/policy
through collector/hub forwarding so these outputs recognize a small positive
total. Keep unrelated legacy behavior scoped. Existing `P_E` numerators retain
their documented native energy representation; do not use that field as the
independent oracle for actual-pair marginals.

## 4. SeedHandler capability and active occurrence mask

Add a proposed virtual `guaranteesStackOnlySeeds() const` method. Its contract
is geometric: every exposed seed is a consecutive diagonal pair chain in
internal coordinates. It does not assert that a seed exists or that all seeds
have the same length. Default to false for an unknown handler.

| Handler | Planned implementation |
| --- | --- |
| `SeedHandlerNoBulge` | Return true. |
| `SeedHandlerMfe` | Return true when both effective, normalized per-strand unpaired allowances are zero; otherwise return false conservatively. |
| `SeedHandlerExplicit` | Inspect retained accepted patterns after parsing and minimum-energy-per-start selection; both dot-bar strings must contain only paired positions. |
| `SeedHandlerIdxOffset` | Forward to the wrapped handler; offsets do not change geometry. |

After `fillSeed()`, compile the admitted occurrences in the prediction range.
Check **both complete endpoints**, since explicit `fillSeed()` counts contained
seeds without pruning its stored family and iteration bounds constrain starts.
Do not infer seed acceptance from complementarity or `seedMaxUP` alone.

Preserve existing policies: explicit patterns ignore computed-seed filters;
competing explicit seeds with the same start retain only the selected minimum-
energy pattern; NoBulge and Mfe differ on equality at their energy thresholds.
Do not replace these contracts with one generic admission test. Explicit
`getConstraint().getBasePairs()` is a minimum, not a common pattern length.

For a run from $a$ to $q$, define $\chi(a,q)$ as a Boolean: at least one admitted
occurrence is wholly contained in the run. In a backward scan, adding start $a$
introduces a hit exactly when its retained seed ends no later than $q$ on the
same run. OR this with the flag for the shorter run. This handles variable
lengths and singleton occurrences without counting witnesses. Do not multiply
the seed's standalone accessibility or full energy into the final interaction;
those values determine admission only.

Move the exact predictor's constructor assertion requiring more than one seed
pair into the legacy paths that require it. Check all other singleton assumptions
in seed annotation, range bounds, and reporting before enabling explicit
singletons in the new exact path. Leave the heuristic's restrictions explicit.

## 5. Sparse search domain and weight ownership

Use zero-based internal coordinates, with sequence 2 reversed so both pair
coordinates increase along a path. A stacked successor is $(i_1+1,i_2+1)$.
Pair vertices must pass accessibility/complementarity checks. Transitions must
pass the native internal-loop checks, including internal GU restrictions.

Preserve the immediate empty-seed return. For an admitted occurrence beginning
at $a$ and ending at $b$, every possible left boundary $p$ obeys

$$
\max(0,b_k-W_k+1)\le p_k\le a_k\quad(k=1,2),
$$

where $W_k$ is the effective inclusive span limit, clipped to the active region.
Discard occurrences longer than either limit. Enumerate the deduplicated union
of these rectangles intersected with valid pair vertices. Use signed or
saturating arithmetic for lower bounds; avoid unsigned underflow. This is a
necessary-condition filter, not a heuristic seed/path selection.

For each surviving $p$, compute only states within its two inclusive span
limits. Use on-demand transitions and run scans rather than materializing all
boundary/run incidences. Preserve intermediate states that cannot currently
terminate due to endpoint GU or ED filters; later extension may make their
terminal boundary admissible. Whole-sequence seed preprocessing remains a
separate cost.

Set

$$
d(p)=\operatorname{BW}(E_{\rm init}),\quad
t(r,a)=\operatorname{BW}(E_{\rm interLeft}(r,a)),\quad
T(a,q)=\prod_{b\to b^+\text{ in the stack }a\ldots q}t(b,b^+).
$$

The empty stack product is $T(q,q)=1$. For admissible complete boundaries use

$$
B(p,q)=\operatorname{BW}(\texttt{energy.getE}(p_1,q_1,p_2,q_2,0));
$$

otherwise $B(p,q)=0$. This includes complete-interval ED, native accessibility-
weighted dangling energies, terminal penalties, and `energyAdd` once. Keep
IntaRNA's existing dangle convention; do not replace it with a different
microstate model. Initiation occurs in $d$ once, and each pair transition owns
one local interaction energy.

## 6. Maximal-stack forward recurrence

Fix left boundary $p$. Let $N_p(q)$ and $S_p(q)$ be hybrid weights of seed-free
and seeded paths through $q$. Let $J^0_p(a)$ and $J^1_p(a)$ be contributions
entering a new maximal run at $a$. For **nonstack** transitions $r\rightsquigarrow a$:

$$
J^0_p(a)=\mathbf1_{a=p}d(p)+\sum_{r\rightsquigarrow a}N_p(r)t(r,a),
\qquad J^1_p(a)=\sum_{r\rightsquigarrow a}S_p(r)t(r,a).
$$

For complete allowed runs with $p\preceq a\preceq q$:

$$
N_p(q)=\sum_{a:\chi(a,q)=0}J^0_p(a)T(a,q),\qquad
S_p(q)=\sum_a[J^1_p(a)+\chi(a,q)J^0_p(a)]T(a,q).
$$

Visit $q$ in increasing pair order, first computing its entries and then its
completed-run sums. Unreachable states are zero. Under `noLP`, omit completed
runs of length one; keep pending $J$ entries so such runs may later extend.
Without `noLP`, source singletons follow the same equations, including an
admitted singleton seed. Enumerate only runs connected by allowed stack edges.

Every path has a unique decomposition at its nonstack edges. Each run either
contains an admitted seed or does not; the seeded-before-run and newly-seeded
cases are disjoint. This proves uniqueness even for overlapping, disjoint,
mixed-length, and singleton seed witnesses. Allowing stack edges in the entry
sum would introduce artificial cuts and duplicate paths.

Accumulate $Z=\sum_{p,q}B(p,q)S_p(q)$. Submit each complete seeded boundary
partition once through `updateCompleteZ(..., S_p(q), true)` with the numerical
policy from section 3. Use the same terminal weight/acceptance for that update
and the outside objective. Cache or recompute coefficients deterministically
without storing every boundary globally. Keep the heuristic's legacy protected
extension methods; its own `predict()` explicitly calls them.

## 7. Outside pass for every actual-pair numerator

This differentiates the forward computation after attaching a marker to the
initial pair and every transition's destination pair. Seed masks and boundary
weights are fixed while differentiating. Write $H^0=N$, $H^1=S$ and
$f(c,a,q)=c\lor\chi(a,q)$.

For each $p$, initialize $\bar H^0_p(q)=0$, $\bar H^1_p(q)=B(p,q)$ for every
terminal $q$, and all $\bar J$ to zero. Reverse $q$ in pair order. For every
forward run term, first accumulate

$$
\bar J^c_p(a)\mathrel{+}=\bar H^{f(c,a,q)}_p(q)T(a,q),\qquad
g(a)=\sum_{c=0}^1\bar H^{f(c,a,q)}_p(q)J^c_p(a).
$$

Avoid the reference's explicit marking of every pair in every run. For this
$p,q$, reconstruct $T(a,q)=t(a,a^+)T(a^+,q)$ backward from $T(q,q)=1$, with
$O(L)$ scratch. Reverse that product chain from its earliest start toward $q$:

$$
M_{a^+}\mathrel{+}=g(a)T(a,q),\qquad
g(a^+)\mathrel{+}=g(a)t(a,a^+).
$$

Use the updated, propagated $g$. Stop before $a=q$, since the empty product
has no pair marker. Set the direct $g$ contribution to zero for excluded forward
terms, but retain necessary product-chain nodes. In particular, `noLP` excludes
singleton completed paths, not the $T(q,q)=1$ node needed for longer products.

After all runs ending at $q$, mark its entry pair and reverse nonstack entries:

$$
M_q\mathrel{+}=\sum_{c=0}^1\bar J^c_p(q)J^c_p(q),\qquad
\bar H^c_p(r)\mathrel{+}=\bar J^c_p(q)t(r,q)
\quad(r\rightsquigarrow q).
$$

Sum $M$ across left boundaries. Each entry owns a run's first pair; product
edges own the rest. One reverse pass per $p$ covers all right boundaries and
all marked pairs, without division, subtraction, or a global product tape.
Skip reverse storage and computation when pair output is not requested.

## 8. Algorithm comparison before production integration

Build correctness-checked numerical kernels for the comparison as an early
implementation milestone, before replacing the production default. Do not
infer production performance from exact-polynomial Python timings.

Let $V$ be the number of allowed pair vertices, $E$ the number of nonstack
edges, $L$ the maximum run length, and $s$ a common fixed seed length. Bounds
below are conservative before span and eligible-start restrictions:

| Candidate | Forward plus all-pair outside work | Logical DP values per left boundary |
| --- | --- | --- |
| Literal Step 5 reference | $O(VE+V^2L^2)$ | Depends on retained reference outputs |
| Maximal stacks with section 7 | $O(VE+V^2L)$ | $O(V)$ plus $O(L)$ run scratch |
| Specialized strict-seed suffix automaton | $O(VE+sV^2)$ | $O(sV)$ |

These are graph-based bounds assuming iteration/storage indexed by allowed
vertices. A dense `Matrix` for an active span rectangle also stores invalid
pair positions: if that rectangle contains $A_p$ cells, physical DP storage is
$O(A_p+L)$ for maximal stacks or $O(sA_p)$ for the suffix candidate. Charge dense
initialization and scanning their actual cell visits: $O(\sum_p A_p)$ for
scalar tables and $O(s\sum_p A_p)$ when every suffix-state table is visited.
Do not substitute the number
of complementary vertices for rectangular cells in reported memory or timing.

For the suffix candidate, store seed-free suffix run length and one absorbing
seeded family. A stack edge increments the capped run length; a nonstack edge
resets it to one. Accept only if an admitted occurrence ends at the new vertex
and the suffix contains its full length. For fixed length, aggregate seed-free
states before processing nonstack edges, avoiding an $s$ factor on those edges.
Initialize source-pair acceptance explicitly for singleton seeds.

The same construction can compare mixed lengths by capping at the maximum
admitted length and using the shortest admitted length ending at each vertex;
longer occurrences ending there are redundant for existential acceptance.
With `noLP`, retain whether the current run has length one or at least two,
including after seed acceptance. Allow nonstack departure and final termination
only from runs of length at least two. This adds constant state. Reverse the
same DP transitions, initialize all accepted terminal states with $B$, and mark
each transition's destination and the source once. For mixed lengths, substitute
the maximum admitted length for $s$ in the cost estimate.

A general prefix-pattern DFA can recognize retained bulged seed paths, but
needs state construction, failure/acceptance transitions, and a decision about
the seed family. Current handlers expose a selected MFE seed per start, not all
bulged structures satisfying constraints. A general DFA does not remove that
semantic limitation automatically. Its state count can be larger than $s$;
do not apply the strict-seed bound to it.

Benchmark no seeds, one sparse seed amid increasing unrelated sequence,
overlapping/disjoint seeds, dense admission, long stacks with $L\gg s$, small
loop limits, mixed explicit lengths, asymmetric spans, and realistic inputs.
Measure seed preprocessing, partition-only work, added outside work, complete
runtime, and peak memory separately. Use identical inputs, models, numerical
policy, compiler settings, and thread count. Validate both candidates first.

Keep maximal stacks if its lower state/storage cost and measured runtime meet
the representative workload requirements. Choose the suffix candidate if its
runtime advantage justifies its state/storage cost. Record the cases and
tradeoffs before selecting one production backend; shipping both or an adaptive
selector is not required. Neither the theory nor this plan claims a universal
speedup. Existing legacy results are correctness references only on independently
verified ensembles; elsewhere use them only as performance baselines.

Report total memory separately: active span grids and run scratch; a full
target-by-query output matrix if requested; optional transition caches; and
legacy boundary-map storage when old trackers are also enabled. Do not retain
all run lists or boundary weights and still claim only $O(V)$ working memory.

## 9. Raw result API, lifecycle, and output

Introduce a proposed `BasePairProbabilities` result with raw `Z_type` numerator
storage, raw denominator, region identity, and explicit completion status. Pass
an optional non-owning result sink to the supported predictor; document lifetime
and keep existing constructor calls source-compatible. The pair-level owner
lives outside predictor/range scopes in `IntaRNA.cpp`.

Do not implement this by attaching an ordinary `PredictionTracker`. Its energy
callbacks lose partition precision and cannot identify interior pairs; moreover,
a non-null tracker currently disables complete-boundary streaming. Retain those
callbacks for old outputs. Pair output alone must not populate `Z_partition`.

Compute each region result privately and commit its raw masses only after
successful forward/outside completion. Allocate/reset region state on every
`predict()` call. Create one accumulator per target/query sequence pair, merge
disjoint completed regions once, and finalize explicitly after all requested
regions succeed. If one region fails or required work is cancelled, mark the
pair result failed and do not write a valid-looking normalized matrix. A writer
destructor must only release resources. Guard publication even on the existing
OpenMP exception path, where final output code can run before the saved exception
is rethrown. Check finite merged totals before forwarding ordinary exact-path
partition diagnostics as well.

Store published matrix coordinates as original target/query positions. Convert
local DP pairs once through `InteractionEnergyIdxOffset::getBasePair()`, then
use each original sequence's `getInOutIndex()` only when printing labels.
Rows are target positions in original 5'-to-3' order; columns are query positions
in original 5'-to-3' order. Do not reverse these columns again with the old spot
writer, which expects internally reversed query storage.

Add proposed `--out=bpProb:FILE`: a CSV-style matrix with the existing separator
and nucleotide/index-label conventions, but its own `bpProb` label and writer.
Use enough significant digits to preserve `Z_type`-appropriate validated
probability precision, including scientific notation for small values. Successful
empty ensembles use `NA`; successful nonempty ensembles use zero for pairs absent
from the searched ensemble. Validate the entire normalized result before writing
its first row. An output I/O failure is an error, not a successful completion.

Register the output enum, parser lookup, help, validation, filename handling,
and writer together. Require a supported exact stack-seed predictor and
`needZall`; do not turn on `needBPs` merely to compute ensemble marginals.
Reject unsupported heuristic, kinetic, evaluation, unseeded, seed-only, or
bulged-seed requests rather than silently changing model/mode. Reuse multi-FASTA
filenames and buffer complete output blocks for synchronized shared-stream
emission. Keep computation local to each sequence pair; no new within-pair
parallelism is necessary. Check API callers' ranges for disjoint merging too.

## 10. Validation gates

### A. Mathematical kernel, independent of native energy rounding

Use independent exhaustive enumeration of ordered pair chains, not the new
recurrence or the handler's own mask as its expected answer. Compare each
boundary partition, total $Z$, and every $M_k$ before normalization. Include
arbitrary positive heterogeneous weights, zero/forbidden boundaries, overlaps,
disjoint seeds, filtered windows, singleton/mixed lengths, missing edges,
asymmetric spans, `noLP`, and the eligible-start optimization.

Use all ten fixed-length theory fixtures through a test-weight adapter where
their model matches. In particular, `random03` has two seeded interactions;
its four seed pairs have probability one and extension pair `(10,4)` has
$e/(1+2e)\simeq0.4223187983$ with **unrounded** accessibility. Map its one-based
antiparallel coordinates explicitly. Compare optimized outside marking with
literal run marking or independent forward derivatives on tiny cases.

### B. Native IntaRNA model

Extend `PredictorSeedOracle_test.cpp` and reuse appropriate fixtures from
`PredictorTinyOracle_test.cpp`. Existing seeded ensemble assertions cover only
a restricted unique-anchor domain; do not assume legacy multi-anchor output is
the oracle. Independently enumerate the admitted family with each handler's
specified threshold semantics, including exact energy/ED equalities, GU at
seed versus interaction ends, explicit selection, and range-crossing seeds.

Native `AccessibilityBasePair` converts ED values to integer hundredths via
`Z_2_E`; the theory's decimal is therefore not an exact native end-to-end golden
value. Native Nussinov and nearest-neighbor oracles must use independently
enumerated structures with the production energy/rounding contract, including
initiation, loops, terminal/dangling terms, and full-span accessibility. Do not
loosen tolerances or change energy precision merely to match the theory decimal.

Check $H_{\rm all}=N+S$ on the same allowed domain, and seeded plus seed-free
partition equality with the independently enumerated unconstrained ensemble.
Production skips starts that cannot contain a seed, so its seed-free tables do
not describe the whole unconstrained ensemble. For this identity either disable
eligible-start pruning in a test-only calculation or restrict both sides of
the comparison to the identical retained-start domain.
For nonempty ensembles verify $0\le P_{ij}\le1$, each row/column sum at most one,
and $\sum_{ij}P_{ij}$ equal to the independently computed expected pair count.
For marker tests, check $\sum_k M_k$ against the derivative that marks only
intermolecular pairs, holding accessibility and seed admission fixed.

### C. Numerical and integration regressions

- Tiny positive totals and positive boundary masses below the old epsilon;
  both raw APIs and exact-path serialized diagnostics must remain nonempty.
- Exponential/product/normalization underflow, overflowing intermediate or
  merged totals, NaN, and successful structural emptiness as distinct cases.
- Invariance of normalized probabilities under a representable constant
  `energyAdd` shift with a **fixed admitted seed family**; a shift can otherwise
  alter computed-seed admission. Test out-of-range shifts as explicit failures.
- Equal raw $Z$ with pair output off/on, legacy trackers off/on, and no duplicate
  denominator increment through predictor, collector, hub, or final writer.
- Virtual update interception and scoped-state restoration after exceptions;
  legacy heuristic paths and existing spot semantics remain covered.
- Asymmetric nonpalindromic sequences, nonzero region offsets, positive/negative
  display shifts, and known interior pairs to expose double query reversal.
- Several disjoint regions, including a structure crossing a region boundary
  that the region-union oracle excludes; `outPerRegion` changes reporting only.
- Injected failure after an earlier successful region: no partial probability
  matrix is finalized. Repeated predictions reset all per-region state.
- Multi-FASTA file naming, shared stdout/stderr blocks, supported thread counts,
  compressed streams where supported, and simultaneous spot/pair output.
- A bulge whose rectangle covers a coordinate pair that never actually pairs:
  spot coverage may be positive while the actual-pair probability is zero.
- Explicit unsupported-mode/window errors and no hidden `Z_partition` retention
  from pair output alone. Account separately for legacy trackers' map cost.

Use justified absolute/relative tolerances with a higher-precision reference;
do not demand bitwise equality after changing sum order or hide model differences
inside a large tolerance. Benchmark after correctness passes, with the same
seed-locality and memory cases used in the early algorithm comparison.

## 11. Implementation sequence and affected components

Each milestone has an exit condition; none is implemented by this document.

| Order | Work | Main files/components | Exit condition |
| --- | --- | --- | --- |
| 1 | Freeze ensemble, weights, numerical outcomes, and independent fixtures | This specification; `tests/PredictorSeedOracle_test.cpp`, `tests/PredictorTinyOracle_test.cpp` | Mathematical and native-model expected results distinguished |
| 2 | Add stack capability and active seed/domain compilation | `SeedHandler*`; proposed seed-domain helper | Handler-specific admission, mixed/singleton, offsets, and sparse boxes tested |
| 3 | Compare forward/outside numerical kernels | Proposed `SeededPartitionFunctionStack` or suffix helper; focused tests | Both compared kernels match oracle; measured default decision recorded |
| 4 | Integrate selected partition backend and scoped numerical policy | `PredictorMfeEns2dSeedExtension`, `PredictorMfeEns`, relevant output status/forwarding | Correct complete boundary sums with output off/on; virtual hooks preserved |
| 5 | Add raw pair result and success-only region aggregation | Proposed `BasePairProbabilities`; `src/bin/IntaRNA.cpp` | Pair numerators, disjoint merges, numerical failures and lifecycle validated |
| 6 | Add writer and CLI support, including tiny-total presentation | `CommandLineParsing.*`, proposed pair writer, `OutputHandlerCsv`, `OutputHandlerEnsemble`, `OutputHandlerText` | Correct indices, formats, statuses and compatibility checks |
| 7 | Finish performance and regression validation, documentation and registration | Library/test `Makefile.am`, README/help, Doxygen, ChangeLog | Required checks pass and intentional output changes are independently explained |

Use the repository's `Matrix<Z_type>` and existing conventions rather than
introducing a new general matrix library. New helper and method names above
are proposals, not existing public API. Preserve existing predictor constructor
ownership, tracker lifetimes, and legacy subclass dependencies.

For the implementation, register new headers/sources and tests, bootstrap after
Autotools source changes, and run `make tests -j2` on supported release and debug
builds. Run installed public-header/consumer checks for new public interfaces.
Record toolchain, options, inputs, timings, peak memory, and numerical differences.
An intentionally corrected partition can change ensemble rankings; justify those
changes against the oracle instead of regenerating expected files blindly.

For this documentation-only plan, validate references, equations, pseudocode
contracts, Markdown structure, distribution registration, and `git diff --check`.
No C++ implementation build or performance result is claimed by this plan.
