# Per-seed output selection (#285, review #288)

The [maintainer's answers](https://github.com/BackofenLab/IntaRNA/pull/288#issuecomment-6319284289)
resolve the planning checkpoint in PR #288:

- Rank each seed's own candidates before forwarding them to overall selection.
- If different seeds select different conformations at identical site boundaries,
  retain the lowest-energy eligible conformation and report the site once.
- For model P, use the ensemble assigned to the current seed as its left-most
  anchor. The first-seed ownership contribution is the requested objective.
- Apply quotas across merged regions/windows, or independently per region
  combination with `--outPerRegion`. Deduplicate, then apply final restrictions
  without refilling quotas. Reject positive quotas with `--noSeed`; ignore them
  for evaluation.

## Implementation

`OutputConstraint::reportMaxPerSeed` is a trailing constructor argument defaulting
to zero. CLI `--outNumberPerSeed` accepts integers in [0,1000]. Zero follows the
existing ranking, traceback, ensemble and output path with no candidate storage.

`OutputHandlerSeedCandidates` keeps at most K distinct sites for each generating
seed, keyed by its start pair in original sequence coordinates. Rediscovering a
site in an overlapping window spends no additional quota. The selected union is
deduplicated by its outer pairs; energy and conformation are retained together.
Deterministic energy/site ordering precedes final energy-window, overlap and
count filtering. No rejected candidate is replaced by the seed's next choice.
The CLI owns the collector across calls; direct predictor API users receive the
same selection within each `predict()` call. Independent regions reset quotas
only after all their windows have completed.

MFE prediction reuses each selected mode's own implementation with a non-owning
single-seed adapter and an unrestricted intermediate output constraint. This
preserves seed-constrained conformations, traceback and the mode's heuristic
recurrences, including S/H, S/M, X/H, X/M, X/R, X/K, B/H and seed-only modes.
The adapter reuses the admitted seed table and forwards full seed annotations.
Enabling the option costs one prediction per admitted seed. When statistics or
trackers are requested, a separate ordinary prediction supplies them exactly
once; seed-specific passes never contribute. This favors explicit, verifiable
semantics over adding seed-indexed DP state to every predictor family. Retained
candidate storage is O(K times the number of seed starts), excluding the reused
prediction workspace and retained traceback lengths. Custom MFE predictor
subclasses can override `predictPerSeed` to supply a factory preserving their
recurrence/state; incompatible inherited factories raise an explicit error when
a positive quota is requested. Complete-site spans are validated before each
seed-specific rank update, including legacy DP cells reserved for extensions.
Runtime can increase
substantially on long sequences with many seeds.

Model P consumes the existing anchor decomposition directly. The stack-seed
kernel's optional `anchorComplete` callback supplies each canonical left-most
anchor's whole hybrid partition at each complete site before full-site
aggregation. It sums noLP support states before publishing the candidate.
Legacy bulged-seed extension uses its existing anchor-specific left/right
contributions, retaining that mode's approximation and numerical limitations.
Model P reports the selected anchor ensemble's energy, including the site's
ordinary boundary/accessibility terms; it does not substitute the full site's
energy. Ensemble reports continue to represent site boundaries and annotate
the generating seed. The stack backend's canonicalization of nested explicit
seeds is unchanged: a containing seed is not an additional owner.

`Zall`, tracker callbacks and base-pair probability matrices remain the ordinary
prediction's quantities. Neither per-seed quotas nor final overlap/count/energy
selection changes the forward/outside objective. Exact/heuristic partition flags
are forwarded with those statistics. Probability SVG structure rendering keeps
its existing site-MFE representation.

## Validation

Regression coverage includes independent exhaustive chain enumeration grouped
by each contained seed for exact MFE models and by the left-most seed for exact
ensembles; quotas 1, 2 and effectively unlimited; cropped/reversed coordinates;
noLP; qualifying traceback energy; seed-only prediction; duplicate sites and
windows; final overlap/energy/count selection without quota refills; and all
currently dispatched seeded CLI modes. CLI checks also cover explicit zero,
parameter files, invalid integers, unseeded/evaluation behavior, merged versus
independent regions, deterministic multithreaded windows and unchanged actual
base-pair probabilities.

GCC 14 / ViennaRNA 2.7.2 release and debug `make tests -j2` both passed all nine
suites (123 API cases; 529,850 release / 529,851 debug assertions). The installed
pkg-config consumer and standalone compilation of the five affected public
headers passed. Out-of-tree builds used the source `src/` include path, the
generated source-side `IntaRNA.pc` expected by current build rules, and copies of
the unchanged CLI fixtures in each build's `tests/data/`. No expected prediction
outputs were regenerated.
