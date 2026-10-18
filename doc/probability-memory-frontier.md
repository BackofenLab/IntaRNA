# Seed-extension memory frontier on the probability branch

This is the forward port of [PR #261](https://github.com/BackofenLab/IntaRNA/pull/261)
to the base-pair probability branch in PR #258. It groups partial ensemble
partitions by their target right boundary. Once no remaining seed starts
at or before that boundary, the predictor can finalize the site's energy
and release the bucket. The key stores the other three coordinates.

The seed list retains the seed handler's traversal order and records a suffix
minimum of target starts. This keeps explicit seeds safe even when their
iteration order is not monotonic. Partition contributions and Zall keep
their existing addition order. Tracker-enabled runs retain the deferred map
and callback order. Each prediction resets the frontier, including empty
regions and predictor reuse.

The change applies to the legacy exact seeded extension path and the
heuristic seeded extension path, including its optional probability trace.
The canonical exact stack-seed kernel already emits complete boundaries;
its recurrence and reverse probability calculation are unchanged. Region
bounds and explicit/singleton seed admission remain those of PR #258.

## Why the code remained uncommitted

The recovered implementation lived in `issue-260-validation/feat-fixed` as
an integration-validation port. PR #261 intentionally targeted master and
included this port as `doc/ensemble-memory-probability-port.patch` rather
than committing it to the probability branch. The old validation worktree
also contained SVG changes that were subsequently committed in PR #259;
those changes are already present in this PR's base and are not duplicated.
This PR promotes the six-file memory port to a reviewable commit and adds
the direct public-header includes used by its containers and hash helper.

The saved October 12 API run passed 101,247 assertions in 97 cases. That is
historical evidence; current validation results are recorded in this PR's
description. Historical measurements and limitations remain in PR #261.
This port makes no new universal memory or runtime improvement claim.

## Regression coverage

The recovered tests compare streamed predictions with a buffered subclass
for partition totals and ranked interaction sites. They cover exact and
heuristic paths, stack/bulged/explicit seeds, LP/noLP, cropped regions,
reused predictors and tracker callbacks. The probability CLI suites also
check that CSV/SVG output still agrees with the supported ensembles.
