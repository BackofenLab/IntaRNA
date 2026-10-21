# Step 2: worker and I/O architecture (#286, PR #287 review)

`PredictionRunner` creates one module Web Worker per task and sends sequence
text, task, start indices, validated parameter values, CSV column IDs, optional
output IDs and evaluation structures. It receives status, result or error
messages. `prepareRequest` validates both the UI and worker boundary; argument
construction only accepts schema keys, enumerated output formats and fixed paths.

The worker imports `wasm/intarna.js`, maps `locateFile` to the lowercase WASM
asset and creates an isolated `createIntaRNA` instance. It writes normalized
multi-FASTA to `/query.fa` and `/target.fa`; local SHAPE/accessibility text uses
fixed `q/tShape.txt` or `q/tAccFile.txt` paths. CSV is written to `/results.csv`;
additional native output prefixes write files under `/output/`. Multi-FASTA
filenames use the CLI's native record-number suffixes.

Main is called once per WASM instance, keeping C++ globals, RNAlib state and files
isolated. Optional ASCII/detailed/ensemble reports each use another fresh
instance with the same inputs/settings. This costs an additional computation per
report and is explained next to those controls. Files are read only after main
returns and streams close. The worker sends CSV, selected columns, additional
file text and bounded diagnostics (25 lines, 2,000 characters each). Nonzero
exits and JS/runtime failures settle as errors, without partial successful files.

Both success and failure terminate the worker. Abort uses `Worker.terminate()`
on the main thread and immediately rejects with `AbortError`, including during
synchronous C++ execution. Stale messages are ignored and the next task starts
fresh. Inputs are never uploaded or persisted. No backend, shared memory or
cross-origin isolation is needed. A single WASM instance can grow up to 1 GiB.

Run `npm test --prefix webinterface`, then the native-reference WASM smoke tests
and Chromium/Firefox suites described in step 4.
