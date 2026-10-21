# Per-seed output selection (issue #285)

Status: implementation checkpoint awaiting maintainer guidance. The option is
not implemented by this change. This records the source investigation and a
proposed implementation plan for [issue #285](https://github.com/BackofenLab/IntaRNA/issues/285),
which requests stopping and asking in the PR when implementation raises questions.

## Requested behavior

- Add `--outNumberPerSeed`, default `0`, accepting nonnegative values using the
  same numeric validation convention as `--outNumber`.
- `0` preserves current behavior.
- A positive value `K` makes each seed's best `K` interactions eligible for
  output. In particular, `1` includes only the best interaction per seed.
- Merge eligible interactions, report duplicates once, and apply `--outNumber`
  to the final output.
- Cover standard prediction modes, excluding evaluation such as `IntaRNAeval`.

The issue's reference to values greater than one is followed by an explicit
example for one; the proposed interpretation is that every positive value
enables the limit.

## Findings in the current implementation

| Area | Current behavior and consequence |
| --- | --- |
| [CLI dispatch](../../src/bin/CommandLineParsing.cpp) | `getPredictor()` selects seeded and unseeded variants of models S, X, B and P, plus seed-only and kinetic modes. There is no shared seed-aware selection interface. |
| [Shared ranking](../../src/IntaRNA/PredictorMfe.cpp) | `initOptima()` reserves the global report count for overlapping output, or one optimum for restricted overlaps. `updateOptima()` prunes candidates before `reportOptima()` performs traceback. An output-handler filter cannot recover candidates already discarded here. |
| [Seed annotation](../../src/IntaRNA/SeedHandler.h) | The handler exposes an MFE seed at each start pair. `addSeeds()` annotates occurrences covered by the actual traced base pairs, rather than every occurrence within the site's intervals. |
| [Seed extension](../../src/IntaRNA/PredictorMfe2dSeedExtension.cpp) | Exact X enumerates candidates while a particular seed is active. That seed identity is not passed to the shared ranker. S instead minimizes seed-containing conformations in its dynamic program. |
| [Ensemble ranking](../../src/IntaRNA/PredictorMfeEns.cpp) | P ranks boundary sites by their partition-derived free energy. That ranking need not agree with either one conformation's MFE or a partition conditioned on a particular seed. |
| [Stack-seed ensembles](../../src/IntaRNA/PredictorMfeEns2dSeedExtension.cpp) | `predictStackSeeds()` uses disjoint first-completed-seed ownership to count chains once. An ownership contribution excludes chains owned by an earlier seed, even when they also contain the current seed. It is not the partition of all chains containing a specified seed. |
| [Regional merging](../../src/bin/IntaRNA.cpp) | Region/window results can be collected before final sequence-pair reporting. A per-seed quota needs a defined scope across these calls. |
| [Existing output contract](../../README.md#subopts) | Ranked predictions have distinct site boundaries, not distinct internal structures. Heuristic predictors retain fewer candidates; restricted-overlap enumeration also discards some alternatives. |

## Decisions needed before implementation

### 1. Which interactions belong to a seed?

Should the limit rank extensions constrained to contain each seed, or filter
the seeds found in each site's existing MFE traceback?

These produce different results. Suppose site A has a best conformation of
energy -10 containing seed t, and another conformation of energy -8 containing
seed s. Site B's best conformation has energy -7 and contains s. Ranking
extensions of s selects A; filtering only existing site-MFE tracebacks selects B.
These are illustrative energies in kcal/mol, not measured predictions.

Proposed direction: rank candidates constrained to the generating seed. Please
also confirm the deduplication rule when different seeds select different
conformations of the same site: retain the lowest-energy eligible conformation
and report the site once, consistent with the current boundary-based contract?

### 2. What is the per-seed objective for model P?

Should P select a seed's sites using:

- the existing full site's ensemble free energy;
- the ensemble free energy conditioned on containing that seed; or
- the energy of the best individual conformation containing that seed?

For example, let a seed's conditional partition weights at sites C and D be 2
and 3, while the full seeded-site weights are 10 and 4. Full-site ensemble
ranking selects C; seed-conditioned ensemble ranking selects D. A canonical
first-seed ownership contribution is a third quantity and cannot substitute
for a containing-seed partition.

Proposed direction: preserve ensemble ranking for P, with maintainer guidance
on whether to condition it on the seed. Keep `Zall`, probabilities and prediction
trackers independent of ranked-output selection, matching the current contract
for reporting limits. If conditioning is selected, also specify whether the
published site energy remains the full-site value or the conditional value.

### 3. Which special modes and scopes are intended?

Please confirm these proposed rules:

- Reject a positive value with `--noSeed`, since that prediction has no seed
  identity. Keep evaluation outside this feature.
- Support every currently dispatched seeded mode: S/H, S/M, X/H, X/M, X/R,
  X/K, B/H and seed-only S/X/P. Include P/H and P/M using decision 2.
- Apply quotas per sequence pair across merged regions/windows; reset per
  region combination with `--outPerRegion=true`.
- Select each seed's best `K` candidates, take their union, deduplicate sites,
  then apply overlap, energy-window and global-count restrictions. Do not
  replace a seed's rejected candidate with its next-ranked candidate after
  final overlap selection. Preserve each mode's documented heuristic scope.

## Implementation plan after review

1. Add CLI parsing/help and a trailing, defaulted per-seed member to
   [OutputConstraint](../../src/IntaRNA/OutputConstraint.h), preserving existing
   constructor calls. Include parameter-file and evaluation-option handling.
2. Add a shared per-seed candidate collector before global pruning. Use
   `Interaction::BasePair` keys in original sequence coordinates for the
   handler's seed starts, retaining the seed geometry needed for traceback.
   Keep the zero setting on the existing path without new candidate storage.
3. Supply the seed membership/objective established by decisions 1 and 2 from
   each predictor family. Preserve the qualifying candidate's energy and
   traceback together. Extend ensemble bookkeeping if conditional partitions
   are required; do not use first-seed ownership as containing-seed membership.
4. Merge retained candidates at the scope agreed in decision 3. Deduplicate
   with the agreed site/structure rule, use deterministic existing tie ordering,
   and apply the final output restrictions without losing eligible candidates
   to an earlier `outNumber` cutoff.
5. Document the option in the README, CLI help and ChangeLog, including model P,
   heuristic, overlap and region semantics.

## Planned validation

- Default and explicit zero: unchanged existing API and CLI results, including
  ensemble statistics, trackers and probability outputs.
- Limits one and two: independent short-sequence enumeration grouped by seed;
  a seed's best site hidden behind another seed's many globally better sites;
  multiple seeds selecting the same site; and the competing-conformation case
  from decision 1.
- Model P: a fixture where conditional, full-site ensemble and individual-MFE
  rankings differ, plus a chain containing several seeds, to enforce decision 2.
- Mode matrix: all agreed seeded modes, explicit/bulged seeds where supported,
  boundary-only output, no seeds, no valid candidates and `outNumber=0`.
- Restrictions: global count, energy window, all overlap policies, deterministic
  ties, disjoint regions, overlapping windows where supported, `outPerRegion`,
  reversed query coordinates and nonzero sequence/range offsets.
- CLI validation: negative/out-of-range values, parameter files, and the agreed
  unseeded/evaluation behavior.
- Run `make tests -j2` in release and debug builds after implementation; extend
  the seed and tiny-sequence oracles for the selected scientific semantics.

This checkpoint changes documentation only. Source references and Markdown
links can be checked without a C++ rebuild. No prediction or performance
results are claimed here.
