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
recurrences used here.

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

## Unique first acceptance

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

The implementation can share the work of all anchors through an inhomogeneous
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

The stack-seed heuristic uses the same canonical first-hit ownership as the
exact computation. For canonical anchor `S`, let `start(S)` and `end(S)` be its
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

Initialize outside derivatives from the same complete boundary coefficients
used by forward accumulation. Reverse the right-extension passes, transfer
their derivatives into the corresponding first-hit injections, and then
reverse the accepting edges and seed-free prefix recurrence. A transition
owns its destination pair; the initialization owns the first pair. Accumulating
each owned pair's forward weight times its outside coefficient counts that
pair exactly once in every retained complete chain, including anchor pairs.

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

The seed-free table still retains run-length state for each active left
boundary. Sharing right extensions avoids repeating identical continuation
work for every anchor, but does not imply a lower asymptotic bound for the
whole calculation. Heuristic cost also depends on the number of distinct
selected endpoints. Pair output adds outside storage and work. Performance or
memory improvements require measurements against the previous implementation
with identical inputs, options, build configuration, and thread count; the
factorization alone is not a performance claim.
