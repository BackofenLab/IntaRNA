# Seed-extension partition functions and actual-pair probabilities

This document describes the nonnegative stack-seed extension formulation for
[issue #264](https://github.com/BackofenLab/IntaRNA/issues/264). Its reference is
the IntaRNA-probabilities branch
[`theory-seed-extension`](https://github.com/BackofenLab/IntaRNA-probabilities/tree/theory-seed-extension),
reviewed at revision
[`e904b65101294f782b0938590335f3e23631a1aa`](https://github.com/BackofenLab/IntaRNA-probabilities/tree/e904b65101294f782b0938590335f3e23631a1aa),
especially the
[seed-extension chapter](https://github.com/BackofenLab/IntaRNA-probabilities/blob/e904b65101294f782b0938590335f3e23631a1aa/theory/seed-extension.md)
and its `SE-LEFT-REC`, `SE-FIRST-HIT`, and `SE-PAIR-FIRST-HIT` equations.
The historical subtraction formulas discussed in that chapter are not the
recurrences used here. The exact backend now evaluates anchor-local left/right
extensions as requested in the [PR #268 review](https://github.com/BackofenLab/IntaRNA/pull/268#issuecomment-6212869312).
The heuristic backend retains first-completed-seed ownership and its frozen
endpoint selector. Both compute their respective ensembles with nonnegative
weights and differentiate the same computation for actual-pair probabilities.

## Ensemble and weight ownership

An interaction is an ordered nonempty chain of complementary intermolecular
pairs. Coordinates below are local, zero-based, and increase on both strands;
the second strand has already been reversed. A legal edge joins successive
pairs and respects the interaction loop limits. An admitted seed is a concrete
stack occurrence retained by the active seed handler, after its energy,
accessibility, position, and seed-specific pairing filters.

For an interaction with outer pairs `p` and `r`, the hybrid weight is its
initiation weight multiplied by one transition weight for each successive pair.
The complete boundary coefficient contains accessibility and the native
boundary energy corrections. It is evaluated for the full outer intervals and
multiplied once. A seed's stored energy is used by seed admission; it is not an
additional factor in the exact chain weight. This avoids counting initiation,
anchor pairs, or their internal loops twice.

The exact ensemble contains every legal chain with an admitted seed, once,
regardless of how many overlapping or disjoint seed occurrences it contains.
Interaction span limits are inclusive on both strands. Whole-interaction
restrictions such as `noLP` are part of the recurrence; endpoint restrictions
such as `noGUend` and complete-span accessibility are applied to the full
boundary. Neither ranked-result limits nor display cutoffs select the ensemble.

## Exact left-most seed-anchored evaluation

The exact backend visits seed anchors and computes left and right extension
tables around each anchor. Its seed family is first reduced to an equivalent
irredundant family: remove an occurrence if it contains another admitted stack
occurrence. Every chain containing the removed seed also contains the retained
inner seed, so the union of interactions, every boundary partition, and every
pair numerator remain unchanged. This reduction is private to exact evaluation;
the admitted seed family used for annotations and heuristic choices is unchanged.

This matters for variable-length explicit seeds. An outer seed may start before
an inner singleton but finish later. In the irredundant family, such nesting is
absent. Seeds on a common chain are either disjoint or overlap on a strict stack;
their starts and ends therefore have the same order. Every seeded chain has a
unique left-most retained seed. That seed owns the chain, even if later seeds
overlap it or follow it after a loop.

For an anchor from `s` to `t`, the left table enumerates paths preceding the
anchor while rejecting earlier retained seeds. Its stack-run state detects
seeds completed by prepending a pair, including earlier seeds that overlap the
anchor. The right table extends from `t` without any seed restriction. Both
respect the active loop limits and pair validity. Their contraction includes
one initiation factor at the complete left boundary, one factor for each
transition, and one complete boundary factor at the outer endpoints. Seed
admission energies are not multiplied into the interaction a second time.

Writing `S[a]` for the product of internal anchor transitions, `L[a,p,h]` for
the left extension, and `R[a,r,k]` for the right extension gives

```text
H[p,r] = initiation(p) * sum_a S[a] * sum_(h,k compatible) L[a,p,h] * R[a,r,k]
Z      = sum_(p,r) H[p,r] * boundary(p,r).
```

The support flags `h,k` matter only for singleton anchors under `noLP`:
compatibility requires `h || k`. Longer seeds already support all anchor pairs.
The right table also records whether its current outer pair is supported;
the left table obtains that information from its initial stack-run length.
Only complete-boundary contractions enforce initiation and total span. Table
rectangles are additionally cropped to the sources present in the current batch.

The left and right offsets include zero, representing no extension. An anchor
of length `k` leaves at most `span-k` extension positions on each strand, shared
between its left and right sides. Each side's table must therefore contain
`span-k+1` cells before cropping at the sequence ends. Full-boundary contraction
checks the sum of the left and right offsets, rather than allowing each side
to consume the whole span independently. Shared buffers are allocated once and
reset only over each anchor's active rectangle.

The `noLP` states retain whether a terminal pair already has a stacked neighbor.
Only a supported outer end may terminate, and only supported pairs may precede
a nonstack edge. A singleton anchor is supported when its left or right
extension supplies a stack; support from both sides is one interaction. An
isolated singleton and a singleton connected only through nonstack edges are
excluded. The outside pass uses the same compatible states.

### Bounded boundary aggregation and outside differentiation

Different anchors can contribute to the same complete boundary. Reporting each
anchor contribution separately would change the complete hybrid coefficient
seen by native ensemble-energy reporting and extension hooks. Instead, the
exact backend accumulates batches of at most 256 left boundaries, then publishes
each positive complete hybrid sum once in row-major boundary order. The same
hybrid coefficient and full-boundary factor determine the total partition and
outside derivatives, even when the boundary factor is zero.

The anchor-local forward tables and their adjoints are shared across anchors.
For pair output, the backend revisits the anchors of a batch, contracts the
saved boundary coefficients into left/right adjoints, and reverses the same
extension recurrences. Extension pairs are marked by forward weight times
outside derivative; every pair in the anchor receives its complete anchor-class
mass exactly once. No pair is counted merely because it lies inside an outer
boundary rectangle. Partition-only mode omits the outside storage and pass.

Batching avoids storing all source/endpoint combinations across the sequences.
The batch is additionally limited to a 2 MiB target for its hybrid table, with
at least one left boundary even when a single span rectangle exceeds that
budget. Pair output saves one equally sized table of boundary coefficients.
The bound is independent of whether pair output is requested, preserving the
forward computation in both modes. Larger batches reuse anchor factors across
more sources but retain more boundary coefficients.

With `W1=min(span1,n)`, `W2=min(span2,m)`, maximum irredundant seed length
`L`, and bounded batch size `B<=256`, workspace storage is
`O(W1*W2*(max(2,L)+B) + n*m + S)` for `S` anchors. The `n*m` term includes
a byte-valued validity cache and, for pair output, the returned mass matrix.
Partition-only mode does not allocate adjoints or saved boundary factors.
Its cost is repeated anchor-table work between batches. This is an explicit
memory/runtime tradeoff, not a claim of a general speedup. The
[benchmark records](../benchmarks/seed-extension/README.md) compare all complete
boundary coefficients, partitions and pair numerators against `985ba44`, and
report allocation calls, requested bytes, runtime and process RSS separately.

Five-sample kernel measurements show the tradeoff: with three seeds, revised
runtime is `0.059–0.197` times the previous runtime; with many seeds it is
`3.055–5.494` times the previous runtime. Exact median process RSS increases
from `3,840–4,224 KiB` to `4,480–6,144 KiB` on these small fixtures. The anchor
index and boundary batches also increase allocation traffic. These results
characterize the requested factorization, not a universal optimization; consult
the benchmark records for per-case timings, allocation counts and source hashes.

### Exact rewrite validation

GCC 14.4.0 release and `--enable-debug` builds both pass `make tests -j2`:
all seven API/CLI drivers pass, with 112 API cases and 399,390 / 399,391
assertions. Independent complete-chain enumeration checks nested and overlapping
mixed-length seeds, singleton support from either side, cropped/native reversed
coordinates, asymmetric spans and loop limits, full hybrid coefficients with
zero boundary factors, and more than 256 eligible source boundaries. Native
checks also verify the reported optimal boundary and rounded ensemble energy.

The two changed public headers compile independently after installation, and
the installed pkg-config probability example compiles, links and runs.
The final kernel snapshots agree across 599,242 boundary/partition/pair entries:
exact maximum relative difference is `5.96e-15`; heuristic dumps are identical.
All 12 native CLI predictions match the base byte-for-byte. The six heuristic
pair-probability files are identical; the six exact files agree to maximum
relative difference `1.53e-15`. Floating-point order changes are expected from
the different nonnegative factorization; the scientific ensemble is unchanged.

## First acceptance in the heuristic backend

The heuristic and exact implementations use different traversal strategies.
See the [OxyS/fhlA runtime investigation](heuristic-seed-runtime.md) for the
cost of repeated heuristic endpoint groups and the lower-corner bounds that
exclude empty right-extension prefixes without changing this ensemble.

Scan a chain from its left boundary. Its canonical anchor is the first admitted
seed to finish. Before that completion, the prefix contains no admitted seed;
after it, the remaining chain has no further seed requirement. This event is
unique even when later seeds overlap the anchor.

Strict stack occurrences ending at the same pair can have different lengths
when explicit seeds are used. Choose the shortest admitted occurrence ending
there as canonical owner. Every chain containing a longer occurrence ending at
that pair also contains this shorter suffix, so these simultaneous matches form
one acceptance event. Ownership is chosen before applying the heuristic
selector described below.

To detect an event, the nonaccepting state retains the terminal strict-stack
run length, capped at the longest relevant seed length. A new pair accepts when
its terminal run is at least the shortest admitted seed length ending there.
The run can exceed that length because earlier potential seed windows may have
been filtered out. A filtered-out window neither accepts nor resets the run.
A nonstack transition resets the run to one pair.

The empty declared seed family produces an empty seeded ensemble. A singleton
explicit seed can first accept at the initial pair or at a later appended pair;
these alternatives are disjoint.

## Seed-free left prefixes, first hits, and right extensions

Fix the exact left boundary `p`. Let `N[p,h,a]` be the hybrid weight of prefixes
ending at pair `h` in nonaccepting state `a`. The initial prefix has the sole
initiation factor. Its weight is assigned to `N` only when the initial pair is
not itself a declared singleton seed. Each legal transition appends its weight
to `N` only if the resulting state remains nonaccepting.

Let `F[p,q,c]` hold first-hit weights ending at pair `q`, with structural state
`c` needed by the right extension. It receives:

- the initiation weight for a singleton source hit at `p = q`; or
- `N[p,h,a] * transition(h,q)` for an accepting edge from a seed-free prefix.

An accepting contribution leaves `N`; it is not retained in both families.
The left prefix owns initiation and all pairs before `q`. The accepting edge
owns the new pair `q`. Thus `F` already includes the complete anchor and must
not be multiplied by another standalone seed weight.

The theory writes a normalized right continuation `B` with value one for an
empty continuation. In the absence of additional structural state, its
first-hit contraction is

```text
Hseed[p,r] = sum_q F[p,q] * B[q,r].
```

The heuristic backend shares the work of compatible anchors through an inhomogeneous
right-extension recurrence:

```text
T[p,q,c] = F[p,q,c]
           + sum_(h,b -> q,c) T[p,h,b] * transition(h,q).
```

Here the edge relation includes any structural-state transition. The right
recurrence has no seed state, seed test, or second initiation factor. Expanding
it gives exactly the sum of first-hit weights multiplied by their unrestricted
right continuations. Combining injections before extending is valid by
distributivity: each complete chain still has its unique first-hit event.
This shares continuation work without enumerating complete interactions or
subtracting overlapping seed classes.

Sum terminally admissible `T` states to obtain the seeded hybrid partition for
each `(p,r)`. Apply the complete boundary coefficient and publish that boundary
once, after all anchor contributions have been combined. Keeping `p` throughout
the calculation is necessary because accessibility and span restrictions can
couple the two outer ends.

## The `noLP` structural state

Seed acceptance does not remove the requirement that every pair have a stacked
neighbor. Before acceptance, a terminal run of length one needs support from a
future stack. A run of length at least two is already supported. Nonstack edges
may leave only supported states.

After acceptance, the right extension retains these two structural states:

- **Pending:** the terminal pair has no stacked neighbor yet.
- **Supported:** the terminal pair already has a stacked neighbor.

A stack from either state enters the supported state. A nonstack edge can leave
only the supported state and enters the pending state. Only the supported
state may terminate. A singleton source hit therefore starts pending and can
become a valid interaction through a right stack. The same rule handles a
singleton first hit following a nonstack edge.

Consequently the right extension is unrestricted with respect to seeds, but
requires structural state when `noLP` is enabled. Using a scalar unconstrained
continuation in that case would admit lonely pairs.

## A defined heuristic subset

The stack-seed heuristic retains the original canonical first-hit ownership.
Unlike exact evaluation, it does not remove redundant containing occurrences:
that could change which anchor selects a right endpoint. For canonical anchor `S`, let `start(S)` and `end(S)` be its
outer pairs, and let `selectedRight(S)` be the deterministic endpoint selected
by the heuristic's right-extension policy. Select the endpoint before the
partition and outside calculations; marking a pair must not change it.

The default predictor preserves its existing right-extension selection hooks.
It enumerates strict right extensions in increasing first-strand endpoint,
then increasing reversed-second-strand endpoint. Candidates pass the existing
right-partition numerical-zero test. Their scores are the full seed-plus-right
ensemble energies, including initiation and the complete boundary energy for
`start(S)` and the candidate endpoint. `noGUend` excludes GU right endpoints.
An endpoint replaces the current choice only for a strictly lower score, so
equal scores retain the first visited endpoint. If no candidate qualifies,
the choice remains `end(S)`. Custom overrides of the right-extension hooks
remain responsible for deterministic choices. The default heuristic predictor
retains its existing rejection of singleton seeds, although the shared kernel
and exact predictor support them.

Keep a complete chain with outer boundaries `(p,r)` precisely when it satisfies
the ordinary interaction constraints and at least one of these equalities:

```text
p == start(S) OR r == end(S) OR r == selectedRight(S).
```

The first alternative retains every right extension with no left extension.
The second retains every left extension with no right extension. The third
retains both-sided extensions ending at the selected right endpoint. The
selector is a union: a chain satisfying several alternatives is counted once.
All legal right paths to an admitted endpoint are retained, rather than only
one representative path. The canonical owner is fixed independently of whether
a later seed would offer a different selected endpoint.

For fixed `p`, anchors with `p == start(S)` can share an unrestricted right pass.
For the remaining anchors, add the zero-right-extension contribution directly
at `end(S)`, and group the other first-hit injections by their selected right
endpoint. Each group shares one right recurrence and contributes only at that
group's terminal endpoint. If the selected endpoint equals the anchor end,
the direct contribution already covers it. Anchors admitted by the unrestricted
case must not also contribute through the restricted cases.

Grouping is valid because the right recurrence is linear and independent of
seed history. Keeping separate selected-endpoint groups preserves the selector
for each canonical anchor; merging all groups into one unconstrained pass would
incorrectly admit other endpoints. Sum the retained contributions for each
complete boundary before reporting it or initializing its outside coefficient.

With unchanged nonnegative weights, this defines a subset of the exact seeded
ensemble, hence `0 <= Zheuristic <= Zexact`. The heuristic probabilities are
conditional on that retained subset. They need not equal the exact seeded
probabilities, and this definition does not promise numerical identity with
the historical signed-overlap heuristic.

## Actual-pair probabilities by the outside pass

The raw numerator for pair `k` is the sum of complete interaction weights whose
chains actually contain `k`. It is not the sum of weights of boundary rectangles
covering `k`; the latter is a spot probability.

For the heuristic first-acceptance recurrence, initialize outside derivatives
from the same complete boundary coefficients used by forward accumulation.
Reverse the right-extension passes, transfer
their derivatives into the corresponding first-hit injections, and then
reverse the accepting edges and seed-free prefix recurrence. A transition
owns its destination pair; the initialization owns the first pair. Accumulating
each state's forward weight times its outside derivative, `N * dN + T * dT`,
counts that pair exactly once in every retained complete chain, including
anchor pairs. First-hit states are not marked separately because their final
pair already belongs to the right-extension state. A directly retained
zero-right-extension branch accounts for its accepting pair in the same way.

For heuristic calculations, reverse the same grouped passes and frozen
endpoint selector. The outside calculation must not independently recompute a
different ensemble or normalize each anchor separately. Merge raw region pair
numerators and denominators before final normalization. An empty ensemble has
undefined conditional probabilities; a numerical failure is an error, not an
empty ensemble.

## Scope, checks, and resource considerations

The run-length automaton supports admitted stack occurrences, including
different-length and singleton explicit seeds. A general gapped seed requires
a different prefix-matching automaton. The existing legacy gapped-seed fallback
is retained; the stack-seed proof and nonnegative formulation do not establish
its correctness or expand its supported probability modes.

The production weights use native IntaRNA energies and `Z_type` arithmetic.
"Exact" describes the ensemble and recurrence, not exact real-number
arithmetic. Checked partition operations must report nonfinite values and
unsupported numerical ranges rather than silently discard positive terms.

The independent tests should compare complete boundary partitions and raw
pair numerators, not only rounded final probabilities. Small path-enumeration
oracles should cover overlapping and disjoint seeds, filtered windows,
co-ending different-length seeds, singleton source and edge hits, both span
limits, and `noLP` states on both sides of acceptance. Heuristic checks must use
the same frozen selector, verify each retained chain once, and compare its
partition with the exact ensemble.

The heuristic seed-free table still retains run-length state for each active
left boundary. Sharing right extensions avoids repeating identical continuation
work for every anchor, but does not imply a lower asymptotic bound for the
whole calculation. Heuristic cost also depends on the number of distinct
selected endpoints. Pair output adds outside storage and work. Performance or
memory improvements require measurements against the previous implementation
with identical inputs, options, build configuration, and thread count; the
factorization alone is not a performance claim.

## First-acceptance workspace reuse (PR #266 review)

The allocation revision requested in [PR #266](https://github.com/BackofenLab/IntaRNA/pull/266#issuecomment-6176325864)
kept the then-shared exact/heuristic recurrence and its arithmetic order unchanged.
The heuristic still uses these buffers; exact evaluation now uses the anchor-local
batched algorithm above. One first-acceptance invocation allocates
its prefix, first-hit, right-extension and (when requested) outside matrices at
`min(span1,n)` by `min(span2,m)`, with the prefix run-length dimension retained.
The inclusive maximum extension offset is `span-1`; storing offset zero requires
`span` cells. This bound covers singleton and mixed-length explicit seeds without
assuming the configured nominal seed length. These shared contractions include
the source and accepting pairs, so reducing the cell count to `span-1` would omit
valid boundaries.

Each left boundary or endpoint group uses its own upper index bounds within the
same fixed-size buffers. Accumulation buffers are cleared over their active
rectangle; validity and seed-ending tables are fully overwritten there. This also
resets cells when a rectangle shrinks at a sequence end and grows at the next
row. Right contractions and their adjoints share one reusable buffer each across
all endpoint groups. Partition-only calls omit the boundary/outside buffers.

Groups retain their first-encounter order. A reusable flat anchor list and dense
local endpoint index replace per-group vectors and hash nodes. The native
heuristic selector also preallocates its right-extension table once and overwrites
each active rectangle, preserving the extension hooks and their logical bounds.
Workspaces are local to the invocation (or predictor for native selection), with
no global shared mutable storage. Gapped-seed fallback behavior is unchanged.

Preallocation retains the worst-case span capacity for the duration of the call;
it does not reduce the asymptotic live-memory bound and can retain more memory
than a particular small active rectangle requires. The goal is to remove repeated
allocation/destruction. The [benchmark](../benchmarks/seed-extension/README.md)
records allocation calls and cumulative requested bytes separately from process
RSS and runtime.

The reuse revision passes `make tests -j2` in GCC 14.4 release and debug builds:
all seven test drivers pass, with 108 API cases and 188,838 / 188,839 assertions.
The added regressions compare changing active rectangles with independent path
enumeration and mixed-length native selector choices with isolated seeds.
Exact/heuristic benchmark dumps match `392357c` byte-for-byte across 599,242
entries; 12 native CLI predictions and pair-probability matrices also match.

## Original factorization measurements

These measurements describe revision `392357c`, before the workspace reuse
revision above; the benchmark README records the newer allocation comparison.

The [reproducible kernel benchmark](../benchmarks/seed-extension/README.md)
compares the exact backend with revision
`cb1c9d419f96c7d2a10bae864d922dcb0234f4d4`, using GCC 14.4.0 release flags
on an AMD Ryzen 5 7530U. Every reported boundary hybrid value, boundary
coefficient, total partition, and pair numerator was compared before timing:
328,236 entries agree with maximum relative difference `1.30e-14`.

Median kernel seconds over three alternating fresh-process samples are below.
Both versions use the same deterministic weights, spans, loop limits, and seed
families. The benchmark directory contains individual samples and source hashes.

| Fixture | Seeds | Partition only, previous / revised | With pair numerators, previous / revised |
| --- | ---: | ---: | ---: |
| Dense 40 by 40 | 1,369 | 0.0618 / 0.0951 | 0.1911 / 0.2509 |
| Dense 40 by 40 | 3 | 0.0446 / 0.0606 | 0.1349 / 0.1580 |
| Banded 60 by 60 | 387 | 0.00960 / 0.0169 | 0.0240 / 0.0335 |
| Banded 60 by 60 | 3 | 0.00811 / 0.0130 | 0.0205 / 0.0277 |

The revised exact kernel is slower on these fixtures: `1.36–1.76` times the
previous partition-only cost and `1.17–1.40` times the cost with pair numerators.
Whole-process median maximum RSS is `3,840–4,096 KiB` for both versions; these
small cases do not resolve the recurrence's allocation cost separately from
runtime and allocator overhead. Samples were not CPU-pinned and other work
could run on the workstation. These measurements characterize the cost of this
factorization on the stated fixtures and do not establish a general runtime
or memory bound.

## Original factorization validation

GCC 14.4 release and `--enable-debug` builds both pass `make tests -j2`: all
seven API/CLI suites pass in each build. The API suite reports 106 cases with
175,007 assertions in release and 175,008 in debug. Coverage includes independent
complete-chain enumeration for exact and frozen-selector ensembles, randomized
tiny selectors, native predictor integration, singleton and mixed-length
anchors, endpoint/region constraints, `noLP`, and outside pair numerators.
The installed pkg-config example and standalone compilation of the changed
public headers also pass. The benchmark's portable build and runner were
executed from the committed scripts, and `git diff --check` passes.

Validation also exposed existing boundary-admission problems that are repaired
with this change. Seed-handler iteration now starts from the remaining
sequence ends, which are invalid cursor positions even when an explicit seed
lies exactly at a cropped region's excluded end. This prevents an out-of-range
seed from causing the first valid seed in the region to be skipped. Conflicting
retained seed lengths at one start are rejected in release as well as debug.
The public probability collector validates region dimensions, annotation-mask
dimensions, and overlapping searched regions in both build modes. A nonfinite
partition produced by a custom predictor override is explicitly rejected before
comparison with the kernel result.
