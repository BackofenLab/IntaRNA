# Step 3: browser controls (#286, PR #287 review)

## Design and layout

The interface follows [vaRRI main-v2](https://github.com/BackofenLab/vaRRI/tree/main-v2):
header indigo `#344a9a`, action blue `#145a96`, dark blue `#0f4676`, pale blue
`#eaf2f8`, white `#ffffff` and focus gold `#f4bb44`. Segoe UI/Tahoma/Geneva serve
labels and body text; Consolas/Liberation Mono serve sequence and structure data.
The official group and university logos share the header with the IntaRNA mark.
Their local copies and source attribution are in `gui/assets/`.

The layout puts the scientific task first, then mutually exclusive parameter
groups, followed by the compute button, results, additional files and versions.
All working content is left aligned; two equally sized strand editors collapse
into one column on mobile. The shared institutional header is the strongest
visual element. The wrapping tab bar identifies parameter groups. Numeric inputs
are at most 9 rem wide within bounded, left-aligned row grids; select/text inputs
use compact type-specific tracks and text/data fields retain useful width. Rows
keep their scientific grouping without stretching across the complete panel,
use 16 px horizontal gaps and bottom-align labels above their controls.
The Seed tab places the constraint selector (12 rem) beside Explicit seeds on
desktop; the CSV sort selector uses a 19 rem track. Both controls stack at mobile
widths.
The institutional link reads “Bioinformatics - University of Freiburg” and opens
the group website. Footer and table notes use the available width, wrapping only
when the viewport requires it.

```text
Group logo       IntaRNA package / RNA–RNA interaction analysis        University logo
RNA–RNA interactions                             Computations stay on your device
Task: Prediction / Evaluation
[ Sequences | Seed | Output | CSV | Model & Energy | Query | Target | CLI-only ]
[ Active panel: e.g. query + target, start indices, evaluation structures      ]
Compute     Input Parameters: Examples / Store / Load / Reset
Results / Show-Hide / name / ZIP download / delete
Interactions / search / show-hide columns / CSV download
Selected interaction / embedded vaRRI
Additional output files / diagnostics / compiled versions
```

Every CLI group has a tab with one visible panel. Sequences is selected
initially. Tabs expose `aria-selected`, `aria-controls` and roving focus; Left/Right
move between available tabs, and Home/End select the first/last. Tab names wrap
on narrow screens. Prediction-only tabs and controls are disabled/hidden during
evaluation; if the selected tab becomes unavailable, Sequences is selected.
Changing task resets parameters to their task defaults while preserving sequence text; Load example also restores the
complementary pair. The browser defaults to 10 interactions per pair (the native default remains
one), and accepts 0–100. Zero enables accessibility-only output.

## Inputs, defaults and validation

`options.js` defines labels, types, defaults and ranges from
`src/bin/CommandLineParsing.cpp`. `protocol.js` validates them on both sides of
the worker boundary. Query/target retain their distinct interior-loop defaults
(16 and 10). Meta-options such as `acc`, `accW`, `accL`, `intLenMax` and
`intLoopMax` are represented by per-strand controls to avoid conflicting values.
FASTA headers supply sequence names instead of `qId`/`tId`. The dedicated CSV
tab contains `outCsvSort` and the interaction-table column choices; Output keeps
prediction filters and additional file formats. `outCsvSort` selects a native
download sort column; table sorting is independent.
CSV uses a fixed semicolon delimiter.
Process settings (threads, file/log paths, parameter files, help/version commands)
are replaced by the worker lifecycle, local data inputs, downloads and displayed
compiled versions. The base-pair maximization energy model is intentionally
excluded. The helix-block interaction model B is supported, with helix controls at the end
of Model & Energy, enabled only for model B. Seed constraint offers Disable,
Explicit seed and Identify seeds, enabling only the relevant controls. Query and
Target include their own accessibility/SHAPE columns and synchronized copies of
the two shared accessibility flags. Ignore disables accessibility controls;
import modes disable window/span and enable the data input. Disabled settings
remain editable when re-enabled but do not enter the request. Explicit row grids
align inputs even when their labels wrap.

CSV columns are grouped as Energy, Output, Sequence, Interaction and Partition
Function, with per-group selection actions and name-only labels with description
tooltips. Prediction adds ED1, ED2, seedE, seedStart1 and seedStart2 to the native
CSV columns; evaluation selects all energy outputs. JSON records these browser
defaults explicitly so native replay retains the same output.
Every additional output checkbox has a descriptive tooltip covering its content,
format and intended use, including the model restrictions shown in its label.

Prediction accepts plain sequence or up to 50 FASTA records per strand.
Evaluation requires one query and target and 1–100 predefined structures.
Wrapped/lowercase sequence text and DNA T normalize to uppercase RNA U. A/C/G/U/N
are accepted. Each record is 2–5,000 nt; seed predictions additionally require
at least `seedBP` nucleotides. FASTA IDs end at the first whitespace, matching the CLI. Duplicate IDs (including collisions after
native CSV semicolon replacement) are rejected. Ignored header descriptions are reported below the input. The longest selected
query length × longest selected target length is capped at 1,000,000, including
pairwise mode; the cap is not summed over records.
The per-strand start index applies to all records and defaults to one. Negative
indices skip zero, matching `RnaSequence`; output coordinates remain inclusive.

Evaluation encodes target then query in 5′ to 3′ dot-bar notation, such as
`1|||&1|||`, separated by colons. Bars pair antiparallel. Validation checks bounds,
equal nonzero pair counts and AU/GC/GU complementarity. Evaluation retains energy,
accessibility and SHAPE controls while excluding prediction/seed/helix/region
and output-filter constraints, matching the native personality.

Validation also covers numeric ranges, accessibility windows/spans, selected
FASTA subsets, regions, structural constraints, seed/helix compatibility,
probability-output requirements and window restrictions. SHAPE and RNAplfold
accessibility text can be pasted or loaded from a local file (up to 1 MB). Data
are written to fixed MEMFS paths, never used as user-provided CLI arguments.
Native diagnostics remain available for model-specific failures.

Edits update the form while completed sessions retain their inputs and results;
the form disables during execution. Abort and Escape terminate the worker during startup or C++ execution. Status, errors,
empty output and diagnostic text have dedicated regions. Keyboard focus,
reduced motion and mobile overflow are covered by real browser tests.
