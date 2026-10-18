# Recovered local validation artifacts

This archive preserves local research and validation source recovered before
removing the IntaRNA workspace. It is a historical corpus for review, rather
than an installed tool or a new prediction implementation. The files are not
wired into the application build, its tests, or CI.

## Contents and original purpose

The `files/` tree preserves the original workspace-relative layout. It includes
standalone public-header/installed-library consumers, benchmark drivers,
configure/dependency probes, archived configure and CI fixtures, native-output comparisons, plotting/report
helpers, instrumentation patches and alternative source implementations used
in local measurements. Source snapshots include the abandoned C++17 work,
a matrix prototype, hash/dense partition experiments, an uncached kinetic
control, and profiling instrumentation. They are archived rather than
substituted for current production files.

The source distribution carries these files in `sources.zip`, whose members
use the original workspace-relative paths. This keeps long historical names
intact despite the existing distribution tar format's 99-character limit.
The individual files remain directly reviewable in Git.

`manifest.json` lists every stored file with its original path, byte count,
SHA-256 and Git blob hash. Copies already reachable from published repository
refs are omitted. The C++23 audit is proposed separately, as is the six-file
probability-memory integration port. Delivery payloads, generated build files,
downloaded dependencies, private correspondence, binaries and large raw data
are outside this source archive. A separate local recovery archive preserves
the remaining unique workspace material.

## Why these files were not committed before

These files lived outside the tracked application checkouts, in issue-specific
validation and profiling directories. They supported earlier PR verification,
benchmarking and exploratory comparisons, while selected production changes
were committed separately. Several are deliberately alternate controls or
instrumented implementations. No file-specific decision record exists for
every helper; their location and content establish their validation role,
not that each was an overlooked production fix. They are committed now so
Martin can inspect that material and the original folder can be removed.

The copied `issue-260-validation/debug-master` header is a pre-finalization
validation snapshot that omits two includes present in PR #261. Its apparent
missing documentation files came from sharing a copied `.git` pointer with
the real worktree; those deletions are not proposed here. The preserved source
snapshot is inactive and is not a header fix for current master.

## Using the archive

Scripts retain their original paths, assumptions and revision-specific APIs.
Reconstruct the relevant source revision and dependency environment before
running them, and use a separate workspace for scripts that edit source files.
The manifest preserves original bytes, including historical formatting. Its
scoped Git attribute avoids treating unified-diff context and copied formatting
as new production whitespace changes.

Python and shell files were checked for syntax when archived. This does not
certify old experiments against current IntaRNA or reproduce their recorded
performance. The active memory port has its own release/debug validation PR.
Historical reports retain their original measurement dates and limitations.

## Original directories

| Directory | Preserved files |
| --- | ---: |
| `issue-184-validation` | 4 |
| `issue-212-validation` | 11 |
| `issue-244-validation` | 10 |
| `issue-245-validation` | 6 |
| `issue-246-validation` | 17 |
| `issue-257-validation` | 4 |
| `issue-258-review-validation` | 3 |
| `issue-258-validation` | 13 |
| `issue-259-revision-validation` | 8 |
| `issue-260-validation` | 3 |
| `issue-262-validation` | 4 |
| `issue-264-validation` | 5 |
| `issue-265-validation` | 8 |
| `issue-267-validation` | 13 |
| `issue-268-validation` | 8 |
| `issue-270-anchored-validation` | 4 |
| `issue-270-validation` | 10 |
| `issue-274-validation` | 2 |
| `kinetic-validation` | 9 |
| `profiling` | 17 |
| `profiling-default-20260921` | 14 |
