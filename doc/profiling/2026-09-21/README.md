# Default-mode profiling requested by Martin

This is the 21 September 2026 follow-up to [Martin's requested optimization mode](https://github.com/BackofenLab/IntaRNA/pull/240#issuecomment-5762928155) and [his request for ED evidence and source links](https://github.com/BackofenLab/IntaRNA/pull/240#issuecomment-5763250539).

Read [REPORT.md](REPORT.md) for the new measurements and priorities, and [ED-READER-EVIDENCE.md](ED-READER-EVIDENCE.md) for the separate writer/reader reproduction. The [19 September study](../2026-09-19/README.md) remains a historical broad survey. Its different revision, inputs, and output configuration are not a before/after comparison with this study.

**Output checks found issues worth reviewing before optimization:** [intermittent energy changes at 12 threads](OUTPUT-STABILITY.md), and a separate ED reload energy difference even after matching the effective length limits. Evidence integrity passes; scientific output equality does not pass for those comparisons. Both findings and their raw outputs are retained.

## Requested mode and recorded source

Profiled current upstream `master` at commit [`0a4568ad6e3da52f935221ff0da98a4a171e9833`](https://github.com/BackofenLab/IntaRNA/commit/0a4568ad6e3da52f935221ff0da98a4a171e9833), using a clean archive and a fresh optimized build. No application patch is part of this study.

| Requirement | Measured configuration | Source at the profiled revision |
| --- | --- | --- |
| Seed-extension heuristic | `--model=X --mode=H`, seed required, 7 base pairs, no seed unpaired bases | [Default values](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L164-L185); [predictor selection](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L2458-L2466) |
| No prediction window mode | `--windowWidth=0` for every timed case | [Unlimited prediction range when width is zero](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L1880-L1887) |
| Local accessibility | Computed ED, `--acc=C --accW=150 --accL=100`; default interaction-limit derivation (`--intLenMax=0`) | [Defaults](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L148-L153); [local accessibility construction](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/CommandLineParsing.cpp#L1955-L1965) |

The 150-nt accessibility folding window is distinct from prediction tiling. Prediction tiling is disabled throughout; local accessibility still uses its normal folding window and 100-nt base-pair span.

All scientific flags above are default values. Only input sizes, sequences, and thread counts vary. Output uses CSV with full `bpList`, retaining traceback; formatting differs from the CLI's default text output. Three biological examples additionally check that explicit scientific flags and omitted/default flags produce identical coordinates, energies, and base-pair lists. Exact, ensemble, alternative predictor models, disabled accessibility, cached ED, and prediction windows are excluded from the performance matrix. ED reload is an untimed diagnostic responding to the second review comment.

## Study design

- 18 configurations, one excluded warmup and five measured repetitions each: **90 measured runs plus 18 warmups**.
- Three biological examples; nested synthetic target/query length sweeps; one 3,000 x 800 nt pair; two 48-target biological batches at 1, 6, and 12 threads.
- The 48 targets are selected at evenly spaced record indices across the 4,319-record repository FASTA. This fixed bounded sample differs from the previous first-48 sample. [Selection indices](inputs/selection.json) and [input hashes](inputs/manifest.json) are recorded. Synthetic generator seed: `20260921`.
- Fresh processes, sequential execution, deterministic randomized measured rounds. The recorded medians and min/max exclude warmups. Full output hashes are checked, including base-pair lists.
- Four CPU profiles (including an independent repeat of the GcvB profile), two heap traces, and three verbose phase runs execute after clean timing. Instrumentation does not supply the timing medians.
- Environment: `OPENBLAS_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, and `LC_ALL=C`. Exact host, compiler, libraries, binary hash, affinity, limits, and governor are in [environment.json](environment.json).

## Evidence and verification

The report links its timing, CPU, heap, phase, and review-response evidence. The [immutable evidence archive](default-mode-evidence.zip) includes raw stdout/stderr and resource records, sequence inputs, scripts, profile text exports, and build/test logs. Large compiled artifacts and raw gprofng experiment directories remain local. Every archived file has a SHA-256 manifest entry; [SHA256SUMS](SHA256SUMS) protects the archive as a whole.

From the repository root:

```bash
report_dir="$PWD/doc/profiling/2026-09-21"
(cd "$report_dir" && sha256sum -c SHA256SUMS)
review_dir=$(mktemp -d)
unzip -q "$report_dir/default-mode-evidence.zip" -d "$review_dir"
python3 "$review_dir/default-mode/scripts/audit.py"
```

The audit checks scientific settings, 108 run records, input/output and binary provenance, full CSV structure, independently recomputed statistics, sequential timing, thread-count comparisons, implicit-default equivalence, all nine instrumented output comparisons, and 28 additional untimed repeats. Expected result: `integrity_status: pass`, `scientific_output_status: differences_observed`. This separates intact evidence from failed scientific equality checks. A missing-executable warning is expected when checking an extracted archive: compiled binaries are excluded, while their recorded per-run hashes remain checked. Verification recomputes derived audit output in the extraction directory, preserving the committed archive.

## Run a fresh experiment

Requires the recorded or an equivalent dependency prefix (GCC 16.1.0, Boost 1.85.0, ViennaRNA 2.7.2), autotools, make, GNU time, GNU gprofng, Python 3, and Matplotlib for figures. A different machine or toolchain is a new experiment. Start at the repository root after extracting the evidence above:

```bash
export INTARNA_REPO="$PWD"
export INTARNA_REVISION=0a4568ad6e3da52f935221ff0da98a4a171e9833
export INTARNA_DEPS=/path/to/the/recorded-or-equivalent/conda-prefix
experiment_dir=$(mktemp -d)
cp -R "$review_dir/default-mode/scripts" "$experiment_dir/scripts"
cd "$experiment_dir"
bash scripts/build.sh
python3 scripts/benchmark.py prepare
python3 scripts/benchmark.py run
python3 scripts/diagnostics.py
python3 scripts/profile.py
python3 scripts/phase_summary.py
python3 scripts/summarize_profiles.py
python3 scripts/ed_numeric.py
python3 scripts/parallel_check.py
bash scripts/ed_conversion_probe.sh
python3 scripts/audit.py
python3 scripts/plots.py
```

Run stages sequentially and avoid other compute-heavy work. Scripts refuse to overwrite timing or profile outputs. Five repetitions describe observed medians/ranges on this host; they are not a population confidence interval, a release comparison, or a biological-accuracy evaluation.
