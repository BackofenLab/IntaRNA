# Phase 3: architecture-preserving performance decisions

This document is the durable decision record for the architecture-preserving
experiments requested in
[discussion #232](https://github.com/BackofenLab/IntaRNA/discussions/232).
Phase 3 is deliberately a set of small scientific experiments, not a general
rewrite: each change has its own branch, focused regression, full API/CLI gate
and same-host benchmark against its exact scientific parent. A candidate is
merged only when it preserves the result and produces a useful end-to-end
improvement on a workload that exercises the changed path.

The exact executable parent is
`8bd1676d1e38ec6461c60e1b832ff3ace9ff77bc`. The later commits through
`321364c56e914e9eff56ddba4846852297fe5c55` change only phase-2 documentation
and CI branch selection, so they do not replace that scientific or benchmark
parent. No result in this record is compared across environment fingerprints.

## Candidate ledger

Implementation and focused regression are intentionally committed together on
each isolated candidate branch. Thus the code and test SHA are the same unless
the row explicitly says otherwise.

| Experiment | Branch | Code SHA | Focused-test SHA | Benchmark suite | Decision |
| --- | --- | --- | --- | --- | --- |
| Remove duplicate feasibility/accessibility predicates | `perf/ap1-redundant-feasibility` | `56ac5eb2771db9c59235f21cc7b329a1660a87f8` | `56ac5eb2771db9c59235f21cc7b329a1660a87f8` | `ap1` | **accept** |
| Cache the two ED values inside `InteractionEnergy::getE` | `perf/ap2-cache-ed-lookups` | `3af257542aa6c7d0e783438884d6b0b2766c3755` | `3af257542aa6c7d0e783438884d6b0b2766c3755` | `core` | **reject** |
| Discard, rather than preserve, scratch matrices on resize | `perf/ap3-nonpreserving-scratch-resize` | `7eae5f475e43ad432e266823c6f341ebdaa39fa3` | `7eae5f475e43ad432e266823c6f341ebdaa39fa3` | `seed-resize` | **reject** |
| Skip partition merges when no output requests `Zall` | `perf/ap4-skip-unused-z-merge` | `0f0a5edb4374a0eb80fe43373a345a1eecf2a6a4` | `0f0a5edb4374a0eb80fe43373a345a1eecf2a6a4` | `ap5` | **reject** |

Rejected branches remain as negative experimental evidence; they are not
stacked under accepted work. All candidates were rooted at the same executable
parent, so their effects remain attributable.

## Scientific and performance gates

Every candidate must pass all of the following before acceptance:

1. The optimized strict-C++23 GCC 14 build and the complete API and CLI suites
   pass at the exact candidate commit.
2. A focused regression proves that the intended operation was removed or
   changed and covers its short-circuit or reuse boundary.
3. Energies, interaction boundaries, base pairs, seeds, rank/filter decisions,
   error behavior and deterministic output fields match the scientific parent.
   Partition values must satisfy the predeclared phase-2 numeric contract.
4. Raw output hashes match where output order is deterministic. For the
   multithreaded throughput case, raw row order may vary, but the canonical
   hash must match exactly.
5. Same-host release measurements use one warm-up and five fresh-process
   repetitions. An accepted gain must be larger than timing dispersion on a
   workload that reaches the changed path; a micro-call-count reduction alone
   is insufficient.
6. Peak RSS must not regress materially unless a larger, scientifically neutral
   speed gain makes the tradeoff worthwhile.

The focused and full-suite evidence recorded before benchmarking is:

| Candidate | Focused evidence | Full GCC 14 API gate | CLI gate |
| --- | --- | ---: | ---: |
| Parent `8bd1676` | phase-1 regressions plus exact unseeded and seeded oracles | 4,275 assertions / 36 cases | 20 / 20 |
| Feasibility `56ac5eb` | 9 assertions: one accessibility delegation per sequence for a feasible pair, and none after complementarity/range failure | 4,284 assertions / 37 cases | 20 / 20 |
| ED cache `3af2575` | 12 assertions: one ED1/ED2 read on the finite path and unchanged left-to-right infinity short-circuiting | 4,287 assertions / 36 cases | 20 / 20 |
| Scratch resize `7eae5f4` | reused `SeedHandlerMfe` and RIblast prediction equal fresh objects over same, smaller and equal-product/different-shape ranges | 4,539 assertions / 36 cases | 20 / 20 |
| Unused-Z merge `0f0a5ed` | 9 assertions / 1 focused case: an unrequested exact ensemble retains positive internal `Zall` with zero output increments; the requested path is positive and increments exactly once | 4,284 assertions / 36 cases | 20 / 20 |

These tests establish correctness preconditions. They do not decide whether an
optimization is worthwhile; that decision comes from the corresponding
benchmark batch.

## Reproducible benchmark provenance

The private evidence repository is
[BackofenLab/IntaRNA-refactor-benchmark](https://github.com/BackofenLab/IntaRNA-refactor-benchmark).
All completed phase-3 batches below use GCC 14.4.0, Boost 1.85.0,
ViennaRNA 2.7.2, libgomp, `-O3 -fno-strict-aliasing`, an AMD Ryzen 5 7530U,
and environment fingerprint
`52fb1f3e9e0f6e3806593ed89ca02aea5ea8fa2f2fa9544e9dfd53f07bdbd68e`.
The result JSON contains the complete compiler, flags, dependency lock,
commands, repetitions, hashes and RSS values.

| Batch | Suite | Parent -> candidate | Harness revision used to run | Evidence commit |
| --- | --- | --- | --- | --- |
| `phase3-ap1-redundant-feasibility-gcc14-v1` | `ap1` | `8bd1676` -> `56ac5eb` | `c2077a54a503c1b1e9a86e0e87b92e6140a50f3e` | `5dfd52dc12ff35f5a6bf31f447537eb7460987af` |
| `phase3-ap2-ed-cache-gcc14-v2` | `core` | `8bd1676` -> `3af2575` | `bbbcc1452afef98faae8f371b9b8553e73549c3e` | `c2077a54a503c1b1e9a86e0e87b92e6140a50f3e` |
| `phase3-ap3-scratch-resize-gcc14-v1` | `seed-resize` | `8bd1676` -> `7eae5f4` | `5dfd52dc12ff35f5a6bf31f447537eb7460987af` | `8c4f649939334a48b31a8200b1db54a87439b1cf` |
| `phase3-ap4-unused-z-merge-gcc14-v1` | `ap5` | `8bd1676` -> `0f0a5ed` | `8c4f649939334a48b31a8200b1db54a87439b1cf` | `cdfea210c8b6a991f8e577f253447ee27d925e99` |

The evidence commit is necessarily newer than the harness revision embedded
in its JSON: it is the commit that adds those generated result files and the
regenerated report.

## Accepted: redundant feasibility checks

`56ac5eb` removes 34 predicates that repeated feasibility work already owned
by the delegated energy/seed check. It removes 68 direct accessibility calls
from 14 predictor/seed source files without changing a recurrence, loop bound,
arithmetic expression or traversal order. The focused `GA`/`AC` regression is
failing-first with respect to call count: a complementary in-range base pair
delegates accessibility exactly once per sequence, while an incompatible or
out-of-range pair stops before accessibility.

Batch `phase3-ap1-redundant-feasibility-gcc14-v1` produced:

| Case | Parent wall mean +/- SD (s), RSS (KiB) | Candidate wall mean +/- SD (s), RSS (KiB) | Parent/candidate | Scientific output |
| --- | ---: | ---: | ---: | --- |
| `nc000913-stratified64-xh-nacc-lp-t1` | 3.602000 +/- 0.068702, 8,678.4 | 3.502000 +/- 0.053572, 8,577.6 | 1.029x | raw and canonical hashes identical |
| `nc000913-stratified64-xh-nacc-nolp-t1` | 0.738000 +/- 0.008367, 8,850.4 | 0.672000 +/- 0.024900, 8,883.2 | 1.098x | raw and canonical hashes identical |

Decision: **accept**. The noLP workload improves by 9.8%, the LP workload by
2.9%, both changes have the same direction, and the strongest gain is larger
than the measured dispersion. The 32.8 KiB noLP mean-RSS increase is about
0.4% and lies inside process-level variation; the LP case decreases by
100.8 KiB. Exact output identity and the full tests show that the removed
checks were redundant rather than an altered feasibility rule.

## Rejected: cache ED inside total-energy evaluation

`3af2575` stores `getED1(i1,j1)` and `getED2(i2,j2)` once each in
`InteractionEnergy::getE`. It retains the literal guard and addition order,
including the rule that an infinite hybrid energy reads neither ED value, an
infinite ED1 reads no ED2, and finite ED1 is checked before ED2. The focused
counter regression proves that the finite path changes from two virtual reads
per ED value to one.

Batch `phase3-ap2-ed-cache-gcc14-v2` produced the following complete core-suite
comparison. `Parent/candidate` above one is faster; below one is slower.

| Case | Parent wall mean +/- SD (s), RSS (KiB) | Candidate wall mean +/- SD (s), RSS (KiB) | Parent/candidate | Scientific output |
| --- | ---: | ---: | ---: | --- |
| `gcvb-phob-exact` | 1.048000 +/- 0.025884, 15,790.4 | 1.046000 +/- 0.011402, 15,686.4 | 1.002x | raw and canonical hashes identical |
| `gcvb-st-ilve-seed-exact` | 0.654000 +/- 0.008944, 15,653.6 | 0.650000 +/- 0.012247, 15,892.0 | 1.006x | raw and canonical hashes identical |
| `nc000913-stratified64-t1` | 5.996000 +/- 0.075697, 13,007.2 | 6.064000 +/- 0.157099, 12,860.8 | 0.989x | raw and canonical hashes identical |
| `nc000913-stratified64-t4` | 1.888000 +/- 0.034205, 25,078.4 | 1.880000 +/- 0.058737, 25,044.0 | 1.004x | raw order unstable in both; canonical hash identical |
| `oxys-fhla-default` | 0.032000 +/- 0.004472, 10,863.2 | 0.034000 +/- 0.005477, 10,682.4 | 0.941x | raw and canonical hashes identical |
| `tiny-bp-exact` | 0.000000 +/- 0.000000, 7,449.6 | 0.000000 +/- 0.000000, 7,449.6 | n/a | raw and canonical hashes identical |
| `tiny-ensemble-nolp` | 0.000000 +/- 0.000000, 7,526.4 | 0.000000 +/- 0.000000, 7,500.8 | n/a | raw and canonical hashes identical |
| `tiny-seedfree-outminpu` | 0.000000 +/- 0.000000, 7,449.6 | 0.000000 +/- 0.000000, 7,500.8 | n/a | raw and canonical hashes identical |

Decision: **reject**. Scientific neutrality is established, but the two
longest workloads range from a 1.1% slowdown to a 0.4% speedup, the exact
cases improve by only 0.2--0.6%, and all are within their observed dispersion.
The very short default case is below useful timer resolution. RSS changes are
small and inconsistent. Reducing a provable virtual-call count therefore does
not produce an end-to-end gain worth carrying in the scientific code.

The branch and benchmark are retained so this tempting micro-optimization is
not repeated without new profiling evidence; it is not an accepted dependency
for any later candidate.

## Rejected: non-preserving scratch resize

`7eae5f4` changes the 27 remaining default-preserving uBLAS scratch-matrix
resizes to `preserve=false`; the one pre-existing non-preserving resize remains
unchanged. The matrices are workspaces whose active cells are recomputed, but
the known predictor-reuse defects from phase 1 make shape transitions the
critical scientific boundary. The focused tests compare reused objects with
fresh objects across `6x5`, shifted same-shape, smaller `4x3`, and equal-cell-
count/different-shape `3x4` ranges for `SeedHandlerMfe` and RIblast prediction.

The specialized `seed-resize` suite contains
`ilve-gcvbst-xh-dense-seed-resize`. Both source builds and their test gates
passed in batch `phase3-ap3-scratch-resize-gcc14-v1`.

| Case | Parent wall mean +/- SD (s), RSS mean +/- SD; max (KiB) | Candidate wall mean +/- SD (s), RSS mean +/- SD; max (KiB) | Parent/candidate | Scientific output |
| --- | ---: | ---: | ---: | --- |
| `ilve-gcvbst-xh-dense-seed-resize` | 0.846000 +/- 0.027019, 8,192.0 +/- 90.5; max 8,320 | 0.850000 +/- 0.029155, 8,038.4 +/- 140.2; max 8,192 | 0.995x | raw and canonical hashes identical |

Decision: **reject**. The candidate is slightly slower and the timing delta is
well inside dispersion. Its lower mean RSS does not establish a memory win:
the samples overlap at the measurement tool's 128 KiB granularity, and the
candidate has greater RSS dispersion. Scientific neutrality is established,
but avoiding a library-preservation copy produces no useful program-level
benefit on the range-shape workload.

## Rejected: skip unused partition merges

Candidate `0f0a5edb4374a0eb80fe43373a345a1eecf2a6a4` avoids calling
`incrementZ` in `PredictorMfe::reportOptima` and in the final CLI pair merge
when `OutputConstraint::needZall` is false. This targets needless entry into a
global OpenMP critical section for ordinary output that never requests
`Zall/Eall`; it does not change partition computation in an ensemble predictor
or the value merged when partition output is requested.

The focused gate passes 9 assertions in one case. In particular, exact
ensemble prediction without requested partition output still computes a
positive internal `Zall` but makes zero output `incrementZ` calls. With
partition output requested, the value is positive, exactly one increment is
made, and the output value equals the predictor's `Zall`. The optimized GCC
14.4 strict-C++23 release build passes the full 4,284-assertion / 36-case API
suite and all 20 CLI cases at the candidate SHA.

The specialized `ap5` suite contains the four-thread window/output-heavy case
`ilve-gcvbst-window-output-heavy-t4`, which exercises 2,016 window pairs. Both
source builds and test gates passed in batch
`phase3-ap4-unused-z-merge-gcc14-v1`.

| Case | Parent wall/user mean +/- SD (s), RSS mean +/- SD; max (KiB) | Candidate wall/user mean +/- SD (s), RSS mean +/- SD; max (KiB) | Parent/candidate | Scientific output |
| --- | ---: | ---: | ---: | --- |
| `ilve-gcvbst-window-output-heavy-t4` | wall 0.130000 +/- 0.000000; user 0.540000 +/- 0.010000; RSS 7,526.4 +/- 140.2; max 7,680 | wall 0.134000 +/- 0.005477; user 0.542000 +/- 0.021679; RSS 7,449.6 +/- 140.2; max 7,680 | 0.970x | raw and canonical hashes identical |

Decision: **reject**. The candidate is slightly slower in both wall and user
time; the four-millisecond wall difference is below useful timer/dispersion
resolution. Mean RSS shifts by less than one 128 KiB measurement step and the
maximum is unchanged. Scientific neutrality is established, but avoiding the
unused merges does not produce a program-level contention or memory benefit.

## Roadmap items not selected in phase 3

The phase-2 roadmap is a profile-driven candidate list, not a requirement to
rewrite every subsystem. The following items are deliberately deferred:

- General workspace retention, active-band clearing and cross-window reuse are
  deferred beyond the isolated resize experiment. Phase 1 found stale
  seed-extension partition state, so a broader lifetime change needs per-family
  fresh-versus-reused oracles and a measured allocation hotspot first.
- Flat/banded replacement of private predictor matrices is deferred. It
  duplicates index and traceback risk across many predictor families, while no
  completed profile yet attributes enough phase-3 runtime to justify that
  surface area. Exact-ensemble boundary storage is handled separately in phase
  4 because it changes recurrence-result ownership.
- Immutable pair-feasibility, reversed-index and accessibility-range
  preprocessing is deferred until profiling shows repeated preprocessing is
  material. It also needs asymmetric target/query, offset, constraint and
  multi-pair concurrency coverage before sharing state.
- Thread/window-local top-k buffering and deterministic bulk merge are deferred
  because the narrower unused-Z experiment found no measurable output-locking
  benefit. A bulk merge can change arrival order, tie handling and overlap
  reduction even when the final list appears similar.
- Replacing the global output critical section with per-instance mutexes is not
  selected: the output-heavy four-thread case showed no contention win. Such a
  change also alters handler layout and copy/move properties and is broader
  than skipping provably unused merges.
- ViennaRNA window APIs and narrower ES work are deferred because they can
  select different library algorithms or numeric paths. They require a full
  ED/ES-table differential across temperature, constraints and base-pair span,
  and profiling must first show preprocessing dominance.
- Boost/string conversion replacements and C++23 syntax-only rewrites are
  deferred because no measured hotspot points to them. Portability and byte-
  stable parsing/formatting costs are not justified by modernity alone.
- Fast-math, reassociation and changed partition accumulation order remain out
  of scope: they weaken the scientific numeric contract rather than optimize
  its current implementation.

These deferrals keep phase 3 proportional to measured benefit. A future branch
may reopen one item with new profile data, but it must use the same exact-parent
scientific and benchmark gates.

## Phase-3 conclusion

Phase 3 is complete: all four isolated candidates have exact code/test SHAs,
full GCC 14 and CLI gates, committed same-parent benchmark evidence, output
classification and an explicit decision. Only redundant feasibility removal
met the merge gate. It is proposed independently against
`refactor-to-c++23` in
[pull request #235](https://github.com/BackofenLab/IntaRNA/pull/235).

The ED cache, non-preserving scratch resize and unused-Z merge are scientifically
neutral negative results: each passed its tests and output comparison, but none
showed a useful end-to-end time or memory improvement. They remain documented
on isolated branches and are not merged merely because a local operation count
or allocation policy looks cheaper. Phase-4 architectural work therefore
remains independently attributable and must carry its own scientific parent
and benchmark decision.
