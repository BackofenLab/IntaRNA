# Deterministic kinetic seed extension

## Council decision and scientific scope

`PredictorSeedExtensionKinetic`, selected exclusively with `--model=X --mode=K`,
implements a deterministic, seed-conditioned, downhill extension heuristic.
It starts from each seed provided by the selected seed handler, compares moves
at both ends of the current duplex, and commits one move at a time. It does not
sample transition rates, simulate elapsed time, cross uphill barriers, remove
base pairs, or guarantee the global minimum free energy interaction. A
macro-step is a coarse move of this heuristic: favorable combined energy does
not establish a barrier-free physical reaction pathway.

This scope follows the distinction between gradient walks to local minima and
kinetic dynamics in the ViennaRNA ecosystem. [RNAlocmin](https://www.tbi.univie.ac.at/RNA/BHG/RNAlocmin.html)
uses gradient walks; [Kinfold and treekin](https://www.tbi.univie.ac.at/software/)
use stochastic moves or transition-rate matrices. The name "kinetic" identifies
the proposed extension strategy, not a calibrated kinetic prediction. The
choice of a single seed structure per seed start inherits the selected seed
handler's behavior; this algorithm does not enumerate every seed conformation.

The second council resolved the first review's open questions as follows.
These decisions supersede contradictory optimization and biological claims in
the supplied `PredSeedExtKinetic.md` proposal.

## State, energy and accepted moves

A state stores the complete, ordered chain of intermolecular base pairs, its
inclusive boundaries `(i1,j1,i2,j2)`, its hybridization energy `H`, and its full
interaction energy `E`. Indices use IntaRNA's internal coordinate system:
sequence 2 is reversed and prediction-range offsets are handled by the existing
wrappers. The initial `H` includes the seed's loop energies and `getE_init()`
exactly once. Every accepted step retains its actual pairs for traceback.

Evaluate every candidate with the active `InteractionEnergy` instance:

```
E_current = energy.getE(i1, j1, i2, j2, H_current)
E_next    = energy.getE(i1_next, j1_next, i2_next, j2_next, H_next)
delta     = E_next - E_current
```

This includes both accessibility penalties, both accessibility-weighted dangling
ends, terminal penalties, the configured temperature and parameters, and the
configured additive energy term. Changing one boundary can alter the opposite
end's dangling contribution through its accessibility weight. Consequently,
loop-plus-stack-plus-accessibility differences alone are insufficient. Infinite
states are rejected before subtraction. A constant additive term cancels in the
difference while remaining part of reported energy.

Only moves with **strictly negative full energy difference** are eligible.
Zero and uphill moves are rejected. Stop when no eligible move remains; the
finite, growing span also guarantees termination. Output energy and accessibility
thresholds filter reported states, rather than introducing additional barriers
in the trajectory. Model-inaccessible states and maximum span limits still
make a move infeasible.

For each side, enumerate all `0 <= s1 <= m1` and `0 <= s2 <= m2`, where `mk`
is the active energy model's maximum unpaired loop size for strand `k`.
`(s1,s2)=(0,0)` is a stack. Other combinations include bulges and internal loops.
The maximum total skipped length is `m1+m2`, equal to `2*m` only when the
per-strand limits are equal. The proposal's fixed 10 is not a separate limit.

A normal move adds the complementary loop-closing pair and advances each
boundary by `sk+1`. With `--outNoLP`, a move across any nonzero loop additionally
requires the immediately following outward stack. Evaluate and commit these
two pairs **atomically**, advancing each boundary by `sk+2`. Neither the
isolated closing pair nor its energy is a separately accepted or reported state.
Check both pairs, the intervening loop, and both complete strand spans before
acceptance. Each span must fit that strand's accessibility maximum length and
the requested prediction range; no single ambiguous shared `W` is introduced.

The complete initial seed must have ordered complementary pairs and finite
per-loop energies. Under `--outNoLP`, every seed pair must already have a direct
stack neighbor. Incompatible explicit seeds are skipped, including seeds with
lonely terminal pairs; this mode does not perform a preliminary seed repair.
Together with atomic macro-steps, this preserves the no-lonely-pair invariant.

Under `--outNoGUend`, both boundaries of every nonstacking extension loop must
be non-GU. Direct stacking may temporarily expose a GU outer endpoint, as in
existing IntaRNA extension recurrences. A state is reportable only if both outer
endpoints satisfy the flag. The active energy model's internal-loop GU policy
also remains authoritative. A valid earlier state remains available if descent
later stops at an unreportable GU endpoint.

## Scores and deterministic choice

`--kineticScore=A|B|C` selects the score; `C` denotes the proposal's C1:

| Option | Score minimized | Interpretation |
| --- | --- | --- |
| `A` (default) | `delta` | Steepest decrease of the actual modeled interaction energy |
| `B` | `delta/(1+s1+s2)` | Heuristic preference per total skipped length |
| `C` | `delta/(1+2*max(s1,s2))` | Heuristic preference penalizing the longer skipped strand |

A is the default because it requires no uncalibrated length-to-time assumption.
B and C are retained as explicit alternatives for exploring the supplied
proposal; their denominators are not experimentally calibrated rates or times.
The written `+1` denominator is retained even for a two-pair macro-step: it
counts one candidate move, not the number of pairs added. Scores choose the next
move only. Reported interactions remain ranked by their modeled total energy `E`.

Compare candidates from **both** sides together. Equal scores are resolved by
left before right, then smaller `s1+s2`, then smaller `s1`. Use sufficiently wide
arithmetic for score comparison, avoiding integer division truncation. Seed
iteration and output comparison follow deterministic existing IntaRNA ordering.

## Why the proposed pruning is removed

A loop-only lower bound is not a lower bound for a loop-plus-stack macro-step:
a positive loop may be rescued by the negative following stack. The initial
review reproduced a Turner2004 example at 37 degrees Celsius with a `+0.50`
kcal/mol loop and `-3.30` kcal/mol following stack, totaling `-2.80` kcal/mol.
The full delta also contains terminal and dangling changes absent from the
proposal's filters. Accessibility at the farthest candidate endpoint is not a
certified lower bound for nearer candidates; imported accessibility values need
not obey monotonicity assumptions.

Therefore the initial implementation exhaustively evaluates all feasible moves
for the current state. It uses neither the six-base-pair-type precomputed Turner
tables nor their `mdspan`/suffix-minimum representation. The suffix-minimum
operation itself is valid, but cannot repair an invalid underlying bound.
Populating an exact full-delta table first would add storage and selection work
without avoiding those energy evaluations. Omitting that cosmetic optimization
is deliberate, rather than replacing one unsafe bound with another.

There are at most `2*(m1+1)*(m2+1)` candidate shapes per state, with constant-size
local energy updates plus full boundary-energy evaluation for each. No runtime
improvement over other predictors is claimed. Any future pruning must supply a
certified lower bound for the **complete** move under the active energy model
and pass differential tests against this exhaustive implementation.

## Reporting, traceback and compatibility

The seed and every atomically committed, reportable prefix are candidates for
normal IntaRNA output. All earlier valid prefixes remain available for suboptimal
and overlap-constrained reporting. Keeping only the last state would lose valid
GU-end prefixes; keeping only the best state for each left boundary would lose
shorter candidates needed after excluding overlap with another report.

Use existing MFE ordering, output filters, seed annotations and index conversion.
Cache the actual selected pair chain for each retained candidate, including all
macro-step pairs. Traceback must recover that chain rather than run an unrelated
minimum-energy recurrence between its boundaries. Resolve duplicate boundaries
by retaining the lowest total energy, then the lexicographically smallest full
base-pair chain on an energy tie. Reduce duplicates before feeding the ordinary
optimum collector, so no stale energy can select a different cached path. Seed
annotations include only starting seeds that passed validation. Output energies
must agree with an
independent sum over the reported chain and complete boundary terms.

The cache stores full paths: its memory cost is proportional to the sum of the
retained paths' lengths, in addition to the seed handler's storage. Retaining
prefixes is therefore not a constant-memory walk. It supports reproducible
traceback and shorter alternatives for overlap-constrained output.

Only `model=X` accepts `mode=K`; seedless operation and `kineticScore` outside K
are rejected. Existing defaults and other models remain unchanged. Ordinary
energy and minimum-energy tracker output are supported. Equilibrium partition
functions and normalized equilibrium probabilities are not defined by this
selected collection of greedy trajectories. Requests needing `Zall`, including
ensemble output and probability trackers, are rejected in CLI validation, with
an API constructor guard for `needZall`. The algorithm does not manufacture an
ensemble by summing repeated prefixes from multiple seeds.

## Implementation plan and acceptance checklist

1. Add the public predictor header and implementation, derived from
   `PredictorMfe`, with owned seed handler, the existing index-offset wrappers,
   validated A/B/C selection and `needZall` rejection. Register both files for
   building and installation.
2. Initialize and trace each handler-provided seed; validate its full chain,
   seed range, energy, strand spans and active structural constraints. Record
   reportable seed states.
3. Enumerate both sides and all feasible loop shapes, form complete normal or
   atomic noLP moves, recompute full candidate energy, and select the strictly
   downhill winner using the chosen score and the specified tie order. Repeat
   until stalling or range exhaustion.
4. Preserve complete committed paths and valid prefixes; integrate ordinary MFE
   output and custom candidate lookup for overlap-constrained suboptimals.
   Reconstruct seed annotations without changing the retained greedy chain.
5. Wire `--model=X --mode=K` and `--kineticScore=A|B|C` into parsing, help and
   factory construction. Reject incompatible model, seedless and ensemble or
   probability requests with clear diagnostics. Update README and ChangeLog.
6. Add an independent tiny-sequence reference that enumerates absolute candidate
   endpoints and recomputes the entire chain energy. Compare reported energy,
   coordinates and traceback against it across scores, constraints and offsets.
   Include targeted regressions for a positive-loop/negative-stack rescue,
   strict stopping, tie order, nonmonotone accessibility, full boundary-energy
   changes, invalid explicit seeds and GU-prefix retention.
7. Run focused API tests, CLI mode/flag compatibility checks, full `make tests`,
   debug validation of bounds and ownership, installed standalone-header checks
   where supported, and `git diff --check`. Inspect failures; never regenerate
   expected outputs merely to hide a change.
8. Review the complete diff and open a pull request documenting this scientific
   scope, the deliberate pruning correction, implemented behavior, validation
   results and any remaining environmental limitations.

Test and build results belong in the pull request and development record; this
checklist specifies the required work and does not imply an unrun check passed.
