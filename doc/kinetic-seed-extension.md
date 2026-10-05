# Deterministic kinetic seed extension

## Scope and modes

`--model=X --mode=K` follows a deterministic downhill path from every seed
provided by the selected seed handler. It compares both ends and commits the
best strictly favorable complete move. This is a zippering-inspired heuristic:
it has no calibrated transition rates or time axis, does not cross barriers
between committed states, and does not guarantee a global minimum. A favorable
two-pair move does not establish a barrier-free physical reaction pathway.

`--mode=L` selects the experimental subclass
`PredictorSeedExtensionKineticPruned`. It applies local-energy and accessibility
pruning **before** checking complementarity and evaluating full energies.
Its additional assumptions can change the selected path; K remains the
reference for measuring those changes. Both modes accept scores A/B/C and
support ordinary energy trackers. They reject seedless operation, other models,
and requests for equilibrium partition functions or probabilities.

These semantics incorporate the October 5 review of
[PR #254](https://github.com/BackofenLab/IntaRNA/pull/254#issuecomment-5991039020).

## Seeds, states and allowed extensions

The predictor trusts the seed handler's structure and hybridization energy.
It does not repeat seed complementarity, noLP or GU-loop checks. Explicit seeds
may contain lonely pairs, including at their ends. Finite energy, prediction
ranges and per-strand span limits still apply. Seed annotations come from the
seed handler without an additional predictor-specific whitelist.

A state stores its complete ordered base-pair chain, inclusive boundaries
`(i1,j1,i2,j2)`, hybridization energy `H`, and full interaction energy `E`.
Sequence 2 uses reversed energy indices; existing wrappers handle prediction
range offsets and conversion to original coordinates. Initially,
`H = seedHandler.getSeedE(i1,i2) + energy.getE_init()`.

Extensions **always** use the no-lonely-pair strategy, independent of the API
output constraint. The CLI sets `--outNoLP=true` when absent or false and emits
an INFO message using the normal logging destination. This applies to new
extensions, not to revalidation of handler-provided seeds. Allowed moves are:

- One pair stacked directly onto the current boundary (`|SEED`).
- Two successive stacked pairs (`||SEED`), evaluated atomically.
- A loop-closing pair followed immediately by its outward stack (`//.SEED`).

For each strand, `sk` denotes the number of skipped bases between the current
boundary and the closing pair. A single-pair move has `s1=s2=0`. Two-pair moves
include all gap pairs from zero through the separate strand loop limits,
subject to available sequence range and maximum interaction span. They advance
each boundary by `sk+2`. Neither an isolated closing pair nor an intermediate
state of a two-pair move is separately committed or reported.

Under `--outNoGUend`, a newly formed nonstacking loop must have non-GU closing
pairs. Stacking can temporarily expose a GU outer endpoint; a state is only
reported if its outer endpoints satisfy the flag. The energy model's own
internal-loop GU policy also applies. Earlier reportable prefixes remain
available if a trajectory stops at an unreportable endpoint.

## Complete energies and deterministic scores

Every candidate surviving geometric and optional pruning checks is evaluated
with the active energy model:

```
H_next = H_current + E_loop + E_stack_if_two_pairs
E_next = energy.getE(i1_next, j1_next, i2_next, j2_next, H_next)
delta  = E_next - E_current
```

`getE()` includes both accessibility penalties, terminal terms, both weighted
dangling ends, and the configured additive term. Updating one end can change
the opposite end's dangling weight, so complete energies must be refreshed
even when local loop energies are reused. Infinite values are excluded before
subtraction. Only `delta < 0` is accepted. Output energy and accessibility
thresholds filter reports, without imposing additional trajectory barriers.

| Score | Quantity minimized |
| --- | --- |
| A (default) | `delta` |
| B | `delta / (1+s1+s2)` |
| C | `delta / (1+2*max(s1,s2))` |

B/C are optional uncalibrated distance preferences. Their denominators count
one move, including two-pair moves. Equal scores prefer left, smaller total
gap, smaller first-strand gap, then a single-pair move. Wide integer cross
products avoid rounding during comparison. Output remains ranked by full `E`.
Two-pair stacking can skip a prefix that the previous implementation visited;
only states actually committed by the revised walk are retained.

## Candidate storage and reuse

Each end has a contiguous rectangular table of two-pair candidates plus its
single-stack candidate. A shared position-pair table records unknown, possible
or impossible complementarity. The closing pair of one candidate may be the
outer pair of another; each such check is resolved once per unchanged end.
The first pass resolves pairing and GU feasibility before the energy pass.
The second pass caches local loop-plus-stack energies and updates the best
candidate as it evaluates full energies, without a separate selection pass.

After committing a move, only that end's tables are rebuilt. The opposite
end's pair checks, feasibility and local energies remain valid. Its total
energy, span eligibility and optional pruning decision are refreshed. An
initially uphill candidate is retained because it may become downhill when
the opposite end changes. Experimentally pruned entries remain unresolved and
can be reconsidered in a later state. Tables use O((m1+2)(m2+2)) space per end.

## Experimental pruning in L

For ViennaRNA, precompute local loop-plus-stack minima for each of the six
oriented root base-pair types, both extension sides and each gap pair. Use the
**active temperature-scaled parameter set**, including custom parameter files.
Minimize over the closing and outer pair types and all four adjacent nucleotide
identities, including unknown nucleotide code 0 within loops. Relaxing consistency between these identities can only make this
local estimate more optimistic. For the base-pair model, two added pairs have
twice its configured base-pair energy. Unknown energy subclasses and API loop
limits above the CLI maximum of 30 fall back to unpruned K enumeration.

Apply componentwise suffix minima to the gap tables. Each cell then bounds
the local loop-plus-stack term for that gap and every larger gap pair. Before
checking candidate pairs, compare:

```
local_suffix_bound + ED1_next + ED2_next - ED1_current - ED2_current >= 0
```

If true, skip that candidate and all componentwise larger gaps. Within the
ordered rectangular traversal this rejects suffixes without more ED, pairing
or loop-energy lookups. A single stacked pair is always evaluated exactly.
Surviving moves still require strictly downhill **complete** energy changes.

This is deliberately an experiment, not a certified bound on the complete
energy change. It assumes ED is monotone as intervals grow and omits changes in
terminal and dangling contributions. Imported/nonmonotone accessibility and
favorable endpoint changes can invalidate the filter; tests include a concrete
case where L stops while K continues. The local table includes the mandatory
stack, so it avoids the original bare-loop error: a Turner2004 loop of +0.50
kcal/mol can be rescued by a -3.30 kcal/mol stack.

Tables require O(12(m1+1)(m2+1)) space and parameter enumeration at construction.
Setup and extra ED lookups may outweigh pruning benefits on small or short-path
inputs. Thus L remains a separate subclass/mode for later real-world evaluation.
See [the reproducible benchmark](benchmark-kinetic.py) and measurements below.

## Reporting and validation

Retain each reportable visited prefix and its actual pair chain. For duplicate
boundaries, keep the lowest full energy, then lexicographically smallest chain.
Reduce duplicates before the normal optimum collector. Overlap-constrained
output can select shorter retained prefixes; traceback restores the selected
path directly. Memory for retained paths is proportional to their total length,
not constant per seed. Repeated predictions reset trajectory and ED caches.

The tests compare K with an independent absolute-endpoint oracle that rebuilds
and reevaluates whole chains. They cover scores/ties, single and double stacks,
positive-loop rescue, strict stopping, separate spans and regions, nonmonotone
ED, GU restrictions, retained prefixes, explicit seeds and annotations, cache
reuse and repeated calls. L has differential and limitation tests. CLI tests
exercise both modes, automatic noLP INFO logging, incompatible requests, and
independent reevaluation of predicted structures through `--rri`.

## Benchmark record (2026-10-05)

Linux x86-64, AMD Ryzen 5 7530U, GCC 14.4 release (`-O3`), ViennaRNA 2.7.2,
Boost 1.85, Kokkos mdspan. These are synthetic measurements, not a real-world
screening benchmark. Each cell is the median of five single-threaded process
runs after one warm-up. Run order rotates. Timing includes startup, folding,
seed enumeration, pruning-table setup and prediction. The script records each
sample, peak RSS and output hashes in
[the raw results](kinetic-benchmark-20261005.json).

The comparison binary uses the same revised move rules, scoring and seed
semantics. Its only change is rebuilding **both** end tables after each move.
It still shares complementarity checks within an update. To reproduce it in
a separate build, replace this line in `extendSeed()`:

```cpp
buildCandidates(sides[best.left ? 0 : 1], bounds, best.left, last1, last2);
```

with:

```cpp
buildCandidates(sides[0], bounds, true, last1, last2);
buildCandidates(sides[1], bounds, false, last1, last2);
```

Build both versions with identical release flags, then run:

```sh
python3 doc/benchmark-kinetic.py --cached /path/to/revised/IntaRNA \
  --uncached /path/to/rebuild-both/IntaRNA --repeat 5 > measurements.json
```

The random cases use deterministic Python seed 254, a 600-nt target and 80-nt
query; the folded case uses `accW=150`, `accL=100`. The stack-rich case uses
100 Gs against 30 Cs without accessibility costs. All use the ViennaRNA energy
model, seven-pair seeds, `intLenMax=60`, `intLoopMax=10`, score A and ten reports.

| Input | Rebuild both ends K (s) | Cached K (s) | Experimental L (s) |
| --- | ---: | ---: | ---: |
| random-no-ED | 0.0173 | 0.0165 | 0.0690 |
| random-folded | 0.2624 | 0.2539 | 0.2807 |
| stack-rich | 0.6640 | 0.5823 | 0.6505 |

Every variant produced the same ten reported structures and energies for these
inputs. Candidate reuse reduced the stack-rich median by about 12%; the short
random/folded runs showed only small gains. Peak RSS was about 15.5 MiB for
random/no-ED, 20.5–20.7 MiB for folded input and 22.3–22.6 MiB for stack-rich
input, without a meaningful memory improvement. Runtime gains vary with the
input and host load; this does not establish a general speedup over other
IntaRNA predictors.

L was slower than cached K on all three samples. Precomputation and extra ED
lookups outweighed any saved candidate work. Together with the endpoint and
monotonicity limitations, this supports keeping L as an explicit experimental
subclass for later real-world benchmarking, rather than enabling it by default.
