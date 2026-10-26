# JSON parameter files

`IntaRNA --parameterFile=parameter.json` accepts the same CLI option names and
values as a legacy `key = value` file. Detection uses the first non-whitespace
character, after an optional UTF-8 BOM, independent of the filename suffix.
Legacy files remain supported. Explicit scalar CLI options override file values;
composing options such as `out` accumulate as before. A file's `personality`
selects defaults before parsing; an explicit CLI personality takes precedence.

```json
{
  "schema": "IntaRNA-parameters-v1",
  "_comment": ["Reproduce: IntaRNA --parameterFile=parameter.json"],
  "parameters": {
    "query": {"text": ">query\nGGGAGGAGGG\n"},
    "target": {"text": ">target\nCCCUCCUCCC\n"},
    "temperature": 25,
    "outMode": "C",
    "out": ["interactions.csv", "qPu:query-pu.txt"]
  }
}
```

`parameters` contains strings, JSON numbers and booleans. Repeated/composing
options accept nonempty arrays. Flags without values accept `true` or `false`.
Unknown options, duplicate keys, null values, unsupported nested values and
recursive `parameterFile` are errors. Existing native type, range and cross-option
validation still applies. JSON uses normal JSON escapes; `#`, whitespace and
newlines in string values are not interpreted as legacy configuration syntax.

The optional `_comment` field holds provenance and reproduction instructions as
ordinary JSON data, so files work with standard JSON parsers. Optional `gui`
metadata describes browser-only output selections and is ignored by the CLI.
Other document fields and unsupported schema versions are rejected.

## Embedded inputs

For `query`, `target`, `qAccFile`, `tAccFile`, `qShape`, `tShape` and `energyVRNA`,
use `{"text": "file contents"}` to embed UTF-8 text. FASTA may contain multiple
records. Accessibility uses RNAplfold text; SHAPE and energy parameters retain
their native formats. Ordinary strings retain their original meaning (sequence,
preset or filesystem path). Relative paths remain relative to the working
directory, as with legacy configuration files. Binary `.agz` caches can still be
referenced by path; they are not a text encoding.

For accessibility of multiple selected records, `{"texts": ["record 1 data",
"record 2 data"]}` materializes files with the existing `-sN` suffixes, indexed by
original FASTA record number. Use `text` when selecting a single record. No
user-specified paths are used for embedded files. A private temporary directory
keeps generated files alive for the entire prediction and removes them on normal
exit or exception unwinding. WASM uses its private in-memory filesystem. Like
other temporary resources, abrupt process termination can leave native files.

## Dependency choice

We evaluated [Boost.JSON](https://www.boost.org/doc/libs/1_87_0/libs/json/doc/html/json/overview.html)
and [nlohmann/json](https://json.nlohmann.me/integration/). Both provide typed JSON
and can avoid a separately linked JSON library. Boost.JSON would raise our
current Boost 1.50 minimum to a release containing JSON, or require another
vendored Boost component and its dependency integration. Boost.PropertyTree is
already available but loses JSON scalar/array type distinctions needed here.

We chose the unmodified, MIT-licensed nlohmann/json 3.12.0 single header (953,436
bytes), compiled in one CLI translation unit. It adds no system package, linker
flag, public library API or WASM build dependency. The tradeoff is a vendored
header to maintain and extra compile/binary size; the parser is outside prediction
loops. Checksums and update instructions live in `src/bin/nlohmann/README.md`.
The source distribution and staged GUI include the license.

### Observed WASM footprint

With Emscripten 4.0.23, ViennaRNA 2.7.2 and Boost 1.87.0 using
`scripts/wasm/build.sh`, the cached pre-JSON engine at `d44d7a6` is 3,217,224 bytes
(781,737 bytes with Python `gzip.compress(..., mtime=0)`). The JSON-enabled engine
is 3,362,458 bytes (823,241 gzip): an increase of 145,234 bytes raw and 41,504
compressed. These are local artifact-size observations with matching dependency
versions, not runtime benchmarks or guarantees for other toolchains. Reproduce
by building both revisions, comparing `site/wasm/intarna.wasm`, and checking
`site/build-info.json` for compiler/library provenance.
