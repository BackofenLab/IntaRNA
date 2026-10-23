# JSON settings and reproducible runs

Load and Store sit to the right of Reset parameters. Store downloads
`parameter.json` containing non-default settings and embedded input text in the
[native JSON format](../../doc/json-parameters.md). `_comment` records compiled
tool/library versions and executable reproduction commands. It is JSON metadata,
not a nonstandard comment syntax. Settings can be saved before sequences are
complete; Compute still validates the complete job before loading WASM.

Load validates the complete document before replacing the form, including the
task, signed start indices, input text, CSV selection/order and additional
outputs. Malformed documents and settings the browser cannot represent leave
the form unchanged with an error. Native-only settings and external file paths
must be used in the CLI; the browser requires embedded text. JSON uploads are
limited to 5 MB. The application never uploads these files to a server.

ZIP uses the same serializer on the completed run, including normalized FASTA.
It contains **only** `parameter.json` and output files. Native replay uses:

```sh
IntaRNA --parameterFile=parameter.json --out=intarna-interactions.csv
```

Additional report commands are listed in `_comment`. They override `outMode`
and CSV-only settings while retaining the same inputs and analysis settings.
Reports run independently, matching the worker. The WASM integration test
replays the JSON in fresh engines and compares every output byte.

## Computed sessions

Each successful computation, including a run with no interaction rows, creates a
new session with a unique ID within the page and a default name such as
`Prediction 1` or `Evaluation 2`. The session retains the computed inputs in JSON,
CSV output, additional files and diagnostics. Editing parameters, Load, examples,
Reset, invalid inputs, failed computations and cancellation preserve prior
sessions. Failed or aborted computations do not create a session.

The newest session opens automatically. Use **Show results** (or **Restore
inputs** after editing), or select its output area, to restore that session’s
computed inputs and fold the other sessions. The task, sequence text, signed
indices, SHAPE/accessibility text, CSV selection/order and additional output
selection are restored together. **Hide results** folds the output. Each session
keeps its own filter, sorting and visible columns; folded vaRRI frames are unloaded
and restored when their session is reopened.

Edit **Session name** directly in the header. Names are plain text; an empty name
returns to the default. **Download ZIP** remains available while folded and always
archives that session’s original JSON and complete outputs, regardless of the
current form or table filters. The session name supplies a safe download filename;
ZIP members remain `parameter.json` and output files. The garbage-bin button
removes only its session and releases its output; it does not change form inputs.

Sessions exist only in the current page’s browser memory. They are not written to
browser storage or uploaded. Reloading or closing the page discards them; use ZIP
exports to retain results. Memory use grows with the runs kept open, so delete
sessions that are no longer needed. Store continues to save the editable form.
If compiled metadata cannot be fetched, a completed run is still retained and its
JSON explicitly records that version information is unavailable.

## Reset and examples

Reset clears the editable sequences, interactions and uploaded files and restores
prediction defaults. Previously computed sessions remain available. Empty FASTA
fields show the short demonstration sequences as placeholders only. The Example disclosure menu sits to the left of Reset.

Maintainers edit `gui/examples/index.json` (an array of local JSON filenames).
Each listed file contains `name`, `description`, and `configuration`, where the
configuration is the same native JSON document used by Load. Two replaceable
examples cover prediction and evaluation. Loading validates their settings and
sequences through the same importer and computation validator. Missing or invalid
examples show an error in the menu; they never overwrite the form.

## Local FASTA files

Each strand has an Upload button and accepts a text-file drop into its textarea.
Files are read locally (at most 60 KB), validated with the same FASTA checks as
manual entry, and replace only that strand after successful validation. Prediction
accepts multi-FASTA; evaluation requires one record. Invalid files leave the
previous input intact. A later edit, reset, import or computation discards any
pending file read, so delayed reads cannot replace newer input.

## Validation of the PR #287 revision (2026-10-23)

- GCC 14 release and debug: `make tests -j2` passed all nine API/CLI groups,
  including `runParameterFile.sh`. The existing out-of-tree CLI fixture harness
  required source data files linked into each build's `tests/data` directory.
- `npm test --prefix gui`: all 18 tests passed.
- Emscripten 4.0.23 / ViennaRNA 2.7.2 / Boost 1.87.0: actual WASM build,
  `node scripts/wasm/smoke.mjs` and `node scripts/wasm/smoke-protocol.mjs` passed.
  Exported JSON reproduced selected CSV, signed/multiple sequences, additional
  outputs and reports byte for byte in fresh engines.
- `scripts/wasm/test-browser.py` passed in Chromium and Firefox at a nested URL,
  covering blank startup/reset, both JSON examples, Load/Store round trips,
  FASTA file upload/drop, stale-read protection, task keyboard controls, JSON ZIP
  contents, result exploration, cancellation/recovery and mobile overflow.
  The desktop screenshot was inspected.
- Embedded ViennaRNA energy-parameter text matched native file input. Existing
  dependency-download and ViennaRNA-release tests, shell syntax, workflow YAML,
  source-distribution contents, vendor checksums and `git diff --check` passed.

## Session follow-up validation (2026-10-23)

- `npm test --prefix gui`: all 18 Node tests passed.
- `scripts/wasm/test-browser.py` passed in Chromium and Firefox with real WASM
  at nested URLs. Session checks cover input/JSON restoration, task changes,
  signed indices, SHAPE text, CSV order, reports, independent filters/columns,
  edits/reset/imports, failed/aborted runs, keyboard selection, safe renaming,
  folded-session archives, deletion, stale file reads and page reload.
- Desktop screenshots were inspected; session controls remain accessible at a
  390 px viewport without document overflow.
- `node scripts/wasm/smoke-protocol.mjs` passed native-reference checks and
  byte-identical archive replay. Browser checks reused the parent branch’s
  compiled Emscripten 4.0.23 / ViennaRNA 2.7.2 assets with their original metadata;
  this follow-up changes no C++ code or worker protocol.
- A fresh GCC 14 release build passed `make tests -j2`: all nine native API/CLI
  groups passed. Tracked fixtures were linked into the out-of-tree build’s
  `tests/data` directory for the existing CLI harness.
- Static staging (including the new session module and engine digest checks),
  shell syntax and `git diff --check` passed.
