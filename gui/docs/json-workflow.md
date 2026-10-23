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

## Reset and examples

Reset clears sequences, interactions, uploaded files and results and restores
prediction defaults. Empty FASTA fields show the short demonstration sequences
as placeholders only. The Example disclosure menu sits to the left of Reset.

Maintainers edit `gui/examples/index.json` (an array of local JSON filenames).
Each listed file contains `name`, `description`, and `configuration`, where the
configuration is the same native JSON document used by Load. Two replaceable
examples cover prediction and evaluation. Loading validates their settings and
sequences through the same importer and computation validator. Missing or invalid
examples show an error in the menu; they never overwrite the form.
