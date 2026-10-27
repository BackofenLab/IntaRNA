# CopomuS Mutation Design in the browser

The **Mutation Design** task integrates the existing CopomuS scientific workflow
without a server or native subprocess. It accepts exactly one query and one
target FASTA record, normalizes DNA `T` to RNA `U`, and rejects ambiguous `N`
bases. CopomuS-specific controls live below the two sequence inputs:

- ranking measures, with arrow controls for priority order;
- MFE or MFE-region/suboptimal candidate selection;
- composable base-pair, lonely-pair and helix-end filters;
- flipped-pair or all-pairing mutation generation;
- the existing alpha/beta separation thresholds; and
- optional query-first explicit encodings such as `G1C&C10G`.

An explicit encoding replaces candidate selection, filters and generation, as in
the CLI. The advanced Seed, Model & Energy, Query and Target tabs configure the
IntaRNA calls. IntaRNA output-selection tabs are hidden because CopomuS chooses
the columns and position profiles required by its measures.

## Runtime boundary

`copomus-worker.js` owns the complete run so aborting terminates both runtimes.
It loads the checksum-pinned, same-origin Pyodide distribution and the staged
Python modules under `/python`. `copomus.browser.run()` reuses the CLI mutation,
filter, measure, comparison, ranking and TSV helpers. Only the process-launch
adapter changes: Python sends a bounded JSON request to the worker and awaits a
fresh IntaRNA WebAssembly instance. The worker validates this private protocol,
embeds the mutated sequences in an `IntaRNA-parameters-v1` document, and returns
CSV/profile text as JSON. No shell command, upload, CDN or external service is
used for computation.

Wildtype output is reused across candidates, matching the CLI optimization.
Every candidate evaluates wildtype/wildtype, wildtype/target-mutant,
query-mutant/wildtype and compensatory-mutant pairs. Browser runs are capped at
100 generated mutations because each candidate can require several complete
folding runs; larger screens remain available through `python/CopomuS.py`.

## Results and reproducibility

The sortable/filterable table contains the complete extended CopomuS TSV. Select
**View** to show the wild-type interaction once in vaRRI, with both designed
substitutions annotated at their original signed sequence coordinates in dark
green. A compact Plotly bar chart beside it compares the `ww`, `wm`, `mw` and
`mm` interaction energies and marks `ww` with a dotted baseline; unavailable
interactions remain empty. Session ZIPs contain `parameter.json` and
`copomus-mutations.tsv`. The JSON keeps native IntaRNA settings and records the
task and CopomuS settings under `gui`, so Load restores the complete editable
run. As with other GUI sessions, all data stays in the tab and disappears on
reload unless downloaded.

The WASM workflow vendors no generated runtime files in Git. Build staging
downloads Pyodide 314.0.7 from its npm release, verifies the pinned SHA-256,
stages only the core runtime/standard library and includes the MPL-2.0 notice.
`build-info.json` records its version and every staged Pyodide asset digest.
