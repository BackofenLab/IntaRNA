# Step 2: worker and I/O architecture (#286)

`PredictionRunner` creates an ES module Web Worker on demand. The main thread
sends `{query, target, count}` as text and receives `status`, `result` or `error`
messages. `protocol.js` validates both the UI and worker boundary and constructs
fixed CLI options; users cannot inject CLI arguments or file paths.

The worker lazily imports `wasm/intarna.js`, awaits `createIntaRNA`, writes
normalized single-record FASTA to `/query.fa` and `/target.fa` in MEMFS, and calls
`callMain` once. It selects semicolon-delimited CSV output in `/results.csv` and
reads that file after main returns and the C++ output streams close. Logs are
separate from CSV, bounded to the most recent 25 lines of 2,000 characters.
Nonzero exit codes, missing assets, allocation failures and JS errors reject
the prediction. Warnings travel with successful results as diagnostics.

Every computation gets fresh C++ globals, RNAlib state and files. Both success
and failure terminate the worker. Abort uses `Worker.terminate()` from the main
thread, so it works during download, initialization and synchronous C++ execution.
Posting an abort message to a busy worker would not interrupt C++. Cancellation
rejects with `AbortError`, releases references and ignores queued stale messages;
the next computation creates a fresh worker. No shared memory, cross-origin
isolation, server, sequence upload or persistent sequence storage is needed.

Run protocol and lifecycle tests with `cd webinterface && npm test`. Browser
integration tests exercise the real worker and WASM assets after staging.
