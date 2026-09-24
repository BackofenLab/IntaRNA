# ED reader: source and reproduced evidence for Martin

This addresses [the request for evidence and source links](https://github.com/BackofenLab/IntaRNA/pull/240#issuecomment-5763250539). The observation is reproduced at current master `0a4568ad6e3da52f935221ff0da98a4a171e9833`. It is separate from the default-mode performance matrix; cached ED is not offered as the optimization target.

## Observation and source path

1. The [writer](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/Accessibility.cpp#L68-L104) labels lengths through `min(sequence size, getMaxLength())` and writes valid values through that limit; it does not append an extra length column.
2. The [reader](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/AccessibilityFromStream.cpp#L87-L105) obtains the last header label, subtracts one, and assigns that reduced value to `availMaxLength` if it is below the requested maximum. The adjacent comment explicitly says the reserve is for dangling-end treatment.
3. [`getMaxLength()`](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/AccessibilityFromStream.h#L158-L164) returns `availMaxLength`, and [`getED()`](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/AccessibilityFromStream.h#L135-L153) returns its upper-bound sentinel for regions longer than that effective maximum.
4. The [computed-accessibility constructor](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L1955-L1965) caps default maximum length using sequence length and the accessibility window. On the 1,000/100-nt inputs here, those computed limits are 150/100. Reading the self-generated files exposes 149/99.

**This is a reproduced writer/reader limit asymmetry, not a claim that the documented dangling-end policy is necessarily wrong.** The open design question is whether self-generated ED reload is intended to reduce the effective interaction search limit, or whether the export/read contract should preserve it.

## Saved files contain their last advertised columns

| File | Last header length | Data rows | Fields per row, including row index | Numeric entries in last column | First valid last-column value (kcal/mol) |
| --- | ---: | ---: | --- | ---: | --- |
| `target.ed` | 150 | 1000 | [151] | 851 | row 150: 5.054000e+01 |
| `query.ed` | 100 | 100 | [101] | 1 | row 100: 2.235000e+01 |

Thus the observation is not caused by a missing or truncated final ED column. The complete files and their hashes are in the evidence archive under `default-mode/diagnostics/ed/`; extracted header/value excerpts are [available directly](diagnostics/ed-header-excerpts.txt).

## Reader messages from the fresh reproduction

Exact messages in [ed_reload.stdout](diagnostics/ed_reload.stdout):

```text
# INFO : initializing ED data for sequence 'synthetic_query_100 : available maximal window length 99 is smaller than maximal interaction length 100 : reducing maximal interaction length to 99
# INFO : initializing ED data for sequence 'synthetic_target_1000 : available maximal window length 149 is smaller than maximal interaction length 1000 : reducing maximal interaction length to 149
```

The uncapped CLI diagnostic requests default interaction limits (`--intLenMax=0`); the cached reader initializes from sequence lengths before applying the header-derived limits. Its diagnostic can therefore mention the 1,000-nt target even though computed local accessibility is capped at 150. The source paths above distinguish those two stages.

On this input, **coordinates and base-pair lists match, but energy does not**, including after both variants are explicitly capped at target/query **149/99**: computed **-14.64**, loaded **-14.66 kcal/mol**. Matching the limits removes the search-space confound but does not establish a numerically transparent round trip. No cache-latency claim is made in this follow-up.

## Separate numeric round-trip observation

[Energy-component output](diagnostics/energy-components.json) localizes the difference: `ED1=4.76` and `ED2=9.8` agree, while left/right dangling contributions change from `-0.21/-0.49` to `-0.22/-0.50` kcal/mol. The target ED entries for the one-base extensions at positions 649–677 and 650–678 are written as `4.890000e+00` and `4.770000e+00`.

The reader [parses into `double curVal` and calls `Ekcal_2_E`](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/AccessibilityFromStream.cpp#L135-L157); the [conversion macro](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/general.h#L103-L113) multiplies by 100 and casts to the integer energy type. A standalone C++ probe using the actual header and the reader's `double` type produces **488** and **476** for those two text values, rather than 489 and 477. The [probe output](diagnostics/ed-conversion.tsv) and source/build script are in the archive. The [conditional dangling probabilities](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/InteractionEnergy.h#L883-L949) use precisely these one-base-extended ED regions. This is an additional numeric conversion problem to review separately from the header-limit policy; no application fix is included.

## Reproduce

Use the pinned executable and the committed [target](inputs/t1000.fa) and [query](inputs/q100.fa) inputs. From the repository root, set `DATA="$PWD/doc/profiling/2026-09-21"`; alternatively, use the extracted archive's `default-mode` directory. From a writable directory, with `BIN` pointing to that executable and `DATA` retaining that absolute path:

```bash
"$BIN" --target="$DATA/inputs/t1000.fa" --query="$DATA/inputs/q100.fa" \
  --model=X --mode=H --acc=C --accW=150 --accL=100 --windowWidth=0 \
  --intLenMax=0 --threads=1 --out=tAcc:target.ed --out=qAcc:query.ed
"$BIN" --target="$DATA/inputs/t1000.fa" --query="$DATA/inputs/q100.fa" \
  --model=X --mode=H --accW=150 --accL=100 --windowWidth=0 \
  --intLenMax=0 --threads=1 --tAcc=E --tAccFile=target.ed \
  --qAcc=E --qAccFile=query.ed --verbose
```

The exact collection commands, including the machine-readable output settings, are [recorded here](diagnostics/commands.jsonl). `scripts/diagnostics.py` in the archive automates file inspection, checks the two reduction messages, and compares full output before/after explicit limit matching. Its JSON result is [diagnostics/checks.json](diagnostics/checks.json). No files were manually extended or truncated, and no application source was modified.
