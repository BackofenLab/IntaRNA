# Output overlap analysis for issue 212

This report examines `--outOverlap=B,N,T,Q` at master revision
`afede020fe7d1c44cf686684dc7eb40927f91583` (2026-10-02). It answers the
[request for analysis and consultation](https://github.com/BackofenLab/IntaRNA/issues/212#issuecomment-5927769865).
It proposes changes but implements none of them.

There are reproducible correctness bugs beyond the documented enumeration
heuristic. Region merging can violate the requested overlap rule, and some
heuristic predictors bypass output validation or use invalid energies for later
results. Separately, even exact prediction retains only one right extension per
left boundary for non-overlapping output. Consequently, a missing suboptimal
result does not establish that no compatible interaction exists.

## Findings and proposed priority

| Finding | Scope | Assessment | Proposed priority |
| --- | --- | --- | --- |
| Independent region results can overlap on a forbidden sequence | `N`, `T`, `Q`, with multiple region combinations | Correctness bug under the default output per sequence pair | High |
| Exhausted seed-free heuristic ensemble enumeration emits an invalid energy and repeats a site | `P/H --noSeed`, `N/T/Q`, requesting more results than available | Undefined numeric conversion | High |
| Later seed-free heuristic ensemble results omit accessibility contributions | `P/H --noSeed`, `N/T/Q` | Incorrect reported energy, ranking and energy filtering | High |
| Later heuristic MFE results can ignore `outNoGUend` | Confirmed for `S/H`, seeded and seed-free, with `T` | Output constraint violation | High |
| Region merging does not apply a global `outDeltaE` threshold | All four overlap modes, multiple regions | Same local versus global scope problem | Medium |
| Valid non-overlapping alternatives can be discarded before selection | `N/T/Q`, including `mode=M` | Documented heuristic with a substantial completeness limitation | Consult on semantics and cost |

These priorities are recommendations, not implementation decisions. The
ensemble-specific findings also explain a single-region overlap violation;
regional merging is not the only way to obtain one.

## What the four settings mean

Overlap concerns inclusive interaction **intervals**, including unpaired bases
between the outermost intermolecular pairs. It does not mean shared base pairs
or only shared paired nucleotides. Coordinates below use the normal one-based
forward target and query indices.

| Setting | Target overlap allowed | Query overlap allowed | Reason to reject a later interaction |
| --- | --- | --- | --- |
| `B` | Yes | Yes | Neither interval is excluded by this option |
| `T` | Yes | No | Its query interval intersects an already selected query interval |
| `Q` | No | Yes | Its target interval intersects an already selected target interval |
| `N` | No | No | Either interval intersects a selected interval on that sequence |

The CLI mapping and `PredictorMfe::reportOptima()` implement this interpretation
correctly within their local selection state. The internally reversed query
coordinates are converted through the energy interface. There is no evidence
here of a simple swapped `T`/`Q` flag or reversed-query indexing error.

The [README](../../README.md#subopts) promises **up to** `outNumber` results,
including the optimum. `B` does not enumerate every structure: exact prediction
normally retains the best structure for each boundary pair, and heuristic
predictors can consider fewer candidates. Allowing an overlap does not require
one. `N` does not search for the best *set* of simultaneously occupied sites;
the current algorithm chooses an optimum first and then compatible alternatives.

## Reproduction setup

Build IntaRNA using the [repository instructions](../../AGENTS.md), then run
from the repository root:

```bash
python3 doc/analysis/out-overlap/reproduce.py ./src/bin/IntaRNA \
  --output /tmp/intarna-overlap-observations.json
```

The script captures the arguments, exit status, output, and forbidden interval
overlaps for 118 invocations. It checks 64 non-exhausting controls: all four
overlap settings for 16 CLI predictor configurations (including RIblast mode).
It independently enumerates every antiparallel structure for three tiny
seed-free examples and verifies their exact `B` site energies. Known defects
are recorded as observations, not installed as regression-test expectations.

All commands below use a Bash helper to keep the examples short:

```bash
rna() {
  ./src/bin/IntaRNA --threads=1 --outMode=C \
    --outCsvCols=start1,end1,start2,end2,E "$@"
}
```

Except where stated, examples use `energy=B` (minus one kcal/mol per pair) and
`acc=N` (zero accessibility penalty) to make energies directly inspectable. These
are demonstrations of program behavior, not biological predictions. Explicit
small seeds and seed-free examples are intentional. `--outDeltaE=100` makes
the energy window non-limiting in these examples.

## 1 Regional merging permits forbidden overlaps

```bash
rna -t CCAACC -q GG --energy=B --acc=N --seedBP=2 \
  --tRegion=1-2,5-6 --outOverlap=N -n 10 --outDeltaE=100
```

Actual output:

```text
start1;end1;start2;end2;E
1;2;1;2;-2
5;6;1;2;-2
```

The query interval `1..2` is reused. With the default `outPerRegion=false`,
`N` should permit at most one of these two interactions. `T` also produces this
violation; `Q` and `B` permit these rows. Conversely:

```bash
rna -t CC -q GGAAGG --energy=B --acc=N --seedBP=2 \
  --qRegion=1-2,5-6 --outOverlap=Q -n 10 --outDeltaE=100
```

reports target `1..2` with query `1..2` and `5..6`, violating `Q` and, if selected,
`N`. The permitted cases are `T` and `B`.

Automatic decomposition reaches the same path:

```bash
rna -t UAUCGGCC -q GG --energy=B --acc=C --seedBP=2 \
  --tRegionLenMax=4 --outOverlap=N -n 10 --outDeltaE=100
```

This reports `(3..4,1..2,-2)` and `(7..8,1..2,-2)`. With `-v`, the prediction
messages show the two target regions. The problem also occurs with ViennaRNA
energies: replace the first example's target/query with `CCCACCC`/`GGG`, use
`--tRegion=1-3,5-7`, and omit `--energy=B`. It reports `-4.2` and `-3` kcal/mol
while reusing query `1..3`.

**Cause.** [IntaRNA.cpp](../../src/bin/IntaRNA.cpp) creates a fresh predictor
for every region/window combination. Each predictor clears its own
`reportedInteractions` lists.
[OutputHandlerInteractionList::add](../../src/IntaRNA/OutputHandlerInteractionList.cpp)
merges, sorts, deduplicates and truncates their results but never enforces
`reportOverlap` across combinations. Disjoint target regions do not prevent
their interactions from sharing a query region, and vice versa.

**Proposed change.** For default output per sequence pair, apply overlap
selection in a common sequence-pair context using original sequence
coordinates. Preserve independent selection when users explicitly request
`outPerRegion=true`, and document that scope. The reproduction with
`outPerRegion=true -n 1` intentionally produces one result per region, so its
combined overlaps are not evidence of a violation of that independent scope.

A final overlap filter alone can ensure validity, but cannot ensure the best
available compatible output or fill the requested count. A region's local
winner may suppress a local alternative and later be rejected against a better
winner from another region. That suppressed alternative can now be admissible.
Selection therefore needs unpruned candidate streams, or a way to request or
recompute alternatives after global exclusions. Applying a fixed top-k cap
before compatibility filtering has the same shortfall.

Window mode is explicitly different: the CLI rejects `outNumber>1` with
`N/T/Q` when `windowWidth` is set. The reproduction script checks that guard.
Manual and automatic regions currently lack equivalent protection or global
selection. `--rri` evaluation intentionally ignores overlap and other prediction
filters and is outside this report's prediction contract.

## 2 Exhausted ensemble enumeration converts infinity to an integer

```bash
rna -t CC -q GG --energy=B --acc=N --model=P --mode=H --noSeed \
  --outOverlap=N -n 2 --outDeltaE=100
```

Actual output on the checked GCC build:

```text
start1;end1;start2;end2;E
1;2;1;2;-2.14748e+07
1;2;1;2;-2
```

Only the `-2` row is valid. The corrupted row is sorted ahead of the real
optimum by the collector. `T` and `Q` reproduce the problem; exact `mode=M`
returns only the valid row. `B` follows a different enumeration path and does
not encounter this exhaustion conversion.

**Cause.** In
[PredictorMfeEns2dHeuristic::getNextBest](../../src/IntaRNA/PredictorMfeEns2dHeuristic.cpp),
`curBestCellE` has type `Z_type` and starts at `Z_INF`. With no eligible cell,
the function assigns floating-point infinity to the integer `curBest.energy`.
It leaves the previous coordinates intact. The invalid converted energy can
pass the next report-loop checks. This is undefined behavior; the exact printed
number and even the failure mode are not portable.

A focused rebuild of this unmodified translation unit with GCC 14.4 and
`-fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all` confirmed:

```text
PredictorMfeEns2dHeuristic.cpp:291:19: runtime error:
inf is outside the range of representable values of type 'int'
```

**Proposed change.** Represent candidate energies with `E_type` and its
`E_INF` sentinel, and explicitly return exhaustion before updating coordinates
or emitting another result. Test zero, one and several remaining candidates
for all restricted overlap modes, with both boundary-only and traced output.
Reusing the final validated candidate path described next would remove the
need for this independent raw-matrix enumeration.

## 3 Later ensemble results use the wrong energy

This is distinct from exhaustion; `-n 2` below stops while a real second
candidate exists.

```bash
rna -t AGAGC -q GAUUC --energy=B --acc=C --model=P --mode=H --noSeed \
  --outOverlap=N -n 2 --outDeltaE=100 --outMaxE=0
```

Actual output:

```text
start1;end1;start2;end2;E
1;2;3;4;-2
5;5;1;1;-1
```

The second site has `ED1=0`, `ED2=1.31`, and total energy `0.31`. The same
predictor reports `0.31` for that site with `--outOverlap=B -n 100
--outMaxE=100`; include `ED1,ED2` in `outCsvCols` to inspect the penalties.
There is only one pair in the site, so no alternative interior structure
explains this difference. It should be excluded by `outMaxE=0`.

**Cause.** The heuristic ensemble `getNextBest()` uses
`energy.getE(curCell->val)`, the hybrid score from its recursion matrix.
The first-result path instead calls `updateOptimaUsingZ()` on finalized site
partition values and `updateOptima(..., isHybridE=true)`, which adds the site
energy contributions. Thus changing the overlap option changes the energy
definition for later rows. This also invalidates their ranking and their
comparison with `outDeltaE`/`outMaxE`.

**Proposed change.** Enumerate later results from the same finalized, filtered
site energies used for the first result. Adding ED to a raw matrix score alone
would not establish equivalence with finalized per-site partition aggregation.
Test the energy of identical boundaries across overlap modes, including
nonzero ED and energy cutoffs. Seeded `P/H` uses a different predictor class;
this reproducer and the raw-matrix diagnosis concern `P/H --noSeed`.

## 4 Later heuristic results bypass terminal pair filtering

```bash
rna -t UUGA -q CAUU --energy=B --acc=N --model=S --mode=H --noSeed \
  --outNoGUend --outNoLP=false --outOverlap=T -n 10 --outDeltaE=100
```

Actual output:

```text
start1;end1;start2;end2;E
1;3;1;2;-2
3;4;3;4;-2
```

The second row pairs target G at position 3 with query U at position 4. This
terminal GU pair violates `outNoGUend`. The row is absent with `B` and with
exact `S/M`. Replacing `--noSeed` by `--seedBP=2 --seedMaxUP=2` reproduces the
same violation in seeded `S/H`. The corresponding seeded `X/H` control does
not reproduce it.

**Cause.** [PredictorMfe::updateOptima](../../src/IntaRNA/PredictorMfe.cpp)
checks terminal GU and maximum ED before recording output candidates. The
specialized `getNextBest()` functions in
[PredictorMfe2dHeuristic](../../src/IntaRNA/PredictorMfe2dHeuristic.cpp) and
[PredictorMfe2dHeuristicSeed](../../src/IntaRNA/PredictorMfe2dHeuristicSeed.cpp)
read recursion cells directly and bypass that validation. The raw state can
legitimately contain intermediate extensions that are unsuitable as complete
reported interactions.

**Proposed change.** Share complete-candidate validation and the final energy
calculation between initial and subsequent selection. Audit the analogous
helix and ensemble overrides for all output constraints, including maximum ED;
the GU reproducer above establishes the two `S/H` paths, not every possible
constraint failure in every override. Merely rejecting the bad final row may
still hide a valid alternative at its starting position, so candidate retention
must also be considered.

## 5 Energy windows are local to regions

```bash
rna -t CCAAC -q GG --energy=B --acc=N --model=S --noSeed \
  --tRegion=1-2,5-5 --outOverlap=B -n 10 --outDeltaE=0
```

Actual output:

```text
start1;end1;start2;end2;E
1;2;1;2;-2
5;5;1;1;-1
5;5;2;2;-1
```

The pair's best energy is `-2`. The `-1` rows are outside a zero-width energy
window around it. They survive because each predictor uses its own local MFE
as the reference and the collector does not apply `deltaE` again. The script
reproduces this in all four overlap settings; restricted modes change how many
of the local rows survive, not the scope of the threshold.

**Proposed change.** Under default output per sequence pair, evaluate the
energy window relative to that pair's best eligible candidate when combining
regions. Under explicit output per region, retain and document local windows.
This is coupled to the global selector in finding 1; `B` needs the global
energy-window check even though it needs no overlap rejection.

## 6 A valid right extension can be discarded before selection

```bash
rna -t CCCAAU -q GGGUU --energy=B --acc=N --seedBP=2 \
  --model=X --mode=M --outOverlap=N -n 100 --outDeltaE=100
```

Only `(1..3,1..3,-3)` is returned. Yet `(4..5,4..5,-2)`, a two-pair AU seed,
is disjoint on both sequences. It appears with `--outOverlap=B`. It is also
returned by rerunning this particular example with
`--tAccConstr=b:1-3 --qAccConstr=b:1-3`. `mode=H` shows the same omission;
changing `N` to `T` does too.

**Cause.** [PredictorMfe::updateMfe4leftEnd](../../src/IntaRNA/PredictorMfe.cpp)
stores just one best right end for each `(i1,i2)`. At target start 4 and forward
query end 5, it retains the better `(4..6,1..5,-3)` extension. That extension
overlaps the first result on the query. The shorter `(4..5,4..5,-2)` extension
has already been discarded, so `getNextBest()` cannot recover it. Several
heuristic predictors instead scan matrices with the same one-extension limit.

The following still smaller seed-free examples establish the limitation for
each restricted mode. The script verifies their complete exact `B` site
energies against an independent structure enumerator. Use `model=S`, `mode=M`,
`noSeed`, `outNoLP=false`, `energy=B`, `acc=N`, `outDeltaE=100`, and `-n 1000`.

| Mode | Target | Query | Only reported site and energy | Omitted compatible site and energy |
| --- | --- | --- | --- | --- |
| `N` | `CCAU` | `CGGU` | `(1..2,2..3,-2)` | `(3..3,4..4,-1)` |
| `T` | `GUA` | `GUA` | `(1..2,1..2,-2)` | `(2..2,3..3,-1)` |
| `Q` | `GCGC` | `AGUG` | `(2..4,2..4,-3)` | `(1..1,3..3,-1)` |

This behavior is already acknowledged in the README's warning about
non-overlapping output. It is therefore a documented completeness limitation,
not evidence that every missing suboptimal row is a newly introduced bug.
Calling `mode=M` exact can nevertheless mislead users unless the distinction
between exact candidate calculation and heuristic restricted enumeration is
explicit. These examples are small witnesses; no global minimality claim is
made.

**Proposed choices for consultation.** Keep and clearly document the heuristic,
or support complete greedy selection by retaining more endpoints or recomputing
the best candidate after excluding previously selected intervals. For `N`,
recomputation must cover all combinations of remaining target and query
intervals; for `T` or `Q`, only the forbidden-overlap sequence is split. Avoid
changing the underlying accessibility model when excluding interaction spans.
Retaining all boundaries can require quartic storage, whereas repeated
prediction trades memory for time. Do not promise exact completeness for a
heuristic predictor merely because output selection has improved.

Blocking paired nucleotides is not a general replacement for interval exclusion.
For example, `-t CCGG -q CCGG --model=S --mode=M --noSeed --energy=B --acc=N
--tAccConstr=b:2-3 --outOverlap=N -n 2 --outDeltaE=100` still produces target
interval `1..4`, with the blocked positions inside a loop. A span-based selector
must reject crossings over an excluded interval, not just forbid pairs at its
positions.

## Decisions and verification needed before implementation

The recommended contract is: choose the best eligible interaction first, then
greedily choose the best eligible interaction compatible with every previously
selected interval, within a documented sequence-pair or region scope. This
preserves the current interpretation of ranked alternatives. Optimizing the
number of interactions or the total energy of a compatible set is a different
objective that can discard the single best interaction; it should not be
introduced implicitly as an overlap bug fix.

Agree on that scope, the desired completeness/runtime tradeoff, and tie-breaking
before changing enumeration. `Interaction::operator<`, the base map scan, and
specialized reverse matrix scans currently use different tie rules. Equal-energy
choices can change which later intervals remain available. This report does
not classify that unspecified ordering as a separate bug, but a common selector
needs an explicit rule in original sequence coordinates.

Suggested implementation order:

1. Repair the invalid ensemble sentinel and use common validated final energies
   for all emitted candidates. Add exhaustion, ED and terminal-pair regressions.
2. Enforce sequence-pair overlap and energy-window scope across regions, with
   explicit independent behavior for `outPerRegion`. Include alternatives that
   become available when a local winner is globally rejected.
3. Decide whether to retain the documented one-extension heuristic or introduce
   more complete restricted enumeration, with measured memory/runtime costs.

Future regressions should check interval validity, energy identity, output
constraints, maximum count and completeness separately. Cover `B/N/T/Q`,
multiple regions on either sequence, automatic decomposition, nonzero region
offsets, shifted output indices, tied energies, seeds, terminal GU constraints,
nonzero accessibility and exhaustion. Existing overlap fixtures cap interaction
lengths at two bases, so they cannot detect a discarded longer/shorter extension
like finding 6. The independent toy enumerator is a useful oracle for those
regressions; the known-bug observations themselves are not correct expected
output.

## Validation and limits

The analysis used a fresh release build of the stated master revision with
GCC 14.4, bundled Kokkos mdspan, ViennaRNA 2.7.2 and Boost 1.85.0 on Linux.
The 118-run script passed all 64 control checks and three independent site
oracles. All 20 existing CLI fixtures passed via `tests/runIntaRNA.sh`, run
from `tests/` with `INTARNABINPATH` set to the repository root. The full API
suite was not rerun for this analysis-only change. Documentation distribution
was checked after regenerating the Autotools files. The focused sanitizer
probe instrumented the ensemble heuristic
translation unit, not the entire application. No production source, CLI
behavior, or existing test expectations were modified.

These controls establish the reported examples and shared code paths, not
exhaustive correctness of every predictor, energy model or constraint
combination. In particular, the successful `B` controls do not prove exhaustive
structure enumeration, and the platform-specific numeric manifestation of
undefined behavior must not be used as a portable regression expectation.
