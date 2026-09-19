# IntaRNA profiling baseline: 19 September 2026

This study profiles `refactor-to-c++23` at commit [`1b8377989ba79a0f1b8a450163b71bd8ed7c3e45`](https://github.com/BackofenLab/IntaRNA/commit/1b8377989ba79a0f1b8a450163b71bd8ed7c3e45). It provides a measured baseline and proposed follow-up work. It does not compare releases or measure the effects of the other open optimization PRs.

Start with the [two-page briefing](BRIEF.pdf) or the [full report](REPORT.md). The report contains all 56 timing results, CPU and heap findings, experimental controls, and limitations.

## Findings and review decisions

| Finding | Evidence | Proposed action |
| --- | --- | --- |
| Confirmed input-stream cleanup defect | 100 open/read/cleanup cycles increased descriptors from 4 to 104. `general.cpp:136` tests for `filtering_ostream`, although the object is a `filtering_istream`. | Make an isolated cleanup fix and add repeated plain/gzip stream regression checks. Introduction history has not been investigated. |
| Default X/H CPU hotspot | Internal-loop energy evaluation accounts for 57.6% of 368 samples in the biological GcvB batch. Exact and ensemble modes have different hotspots. | Measure focused scoring/energy changes separately by mode; sample share is not a promised speedup. |
| Batch parallelism trades memory for throughput | The 48-target ChiX batch takes 4.944 s at one thread and 0.974 s at twelve (5.08x); peak RSS rises from 18.7 to 65.3 MiB. | Choose a thread count using both latency and memory requirements. Six threads reach 88% of twelve-thread throughput with 61% of its RSS on this host. |
| Saved accessibility can reduce repeated-call cost | With matched target/query limits of 149/99, loading ED takes 0.492 s versus 0.805 s computed, with equal base-pair lists. One-time generation takes 0.847 s and is excluded. | Confirm whether the reader's one-nucleotide limit reduction is intended before treating ED reload as transparent. |
| Windowing is model dependent | On a 3,000 x 800 nt pair, S/H RSS falls from 127.5 to 24.1 MiB, while runtime grows from 21.38 to 59.71 s. Default X/H gains no useful RSS reduction and takes 1.72x longer. | Use windowing when the measured memory benefit justifies its runtime cost; validate outputs for the chosen workload. |

Martin's review can resolve two decisions: whether the ED-file limit reduction is intended, and which production workload should serve as the acceptance benchmark for the first CPU optimization. The cleanup defect has a separate, reproducible resource-retention probe. No proposed application fix is included here.

## Evidence layout

- [REPORT.md](REPORT.md): complete findings, all timing results, and limitations.
- [results/overview.png](results/overview.png): visual overview; the original SVG is included in the evidence archive.
- [results/summary.csv](results/summary.csv) and [JSON](results/summary.json): medians, observed ranges, CPU time, and process peak RSS for the three measured repetitions per case.
- [environment.json](environment.json), [cases.json](cases.json), and [input manifest](inputs/manifest.json): recorded build/host provenance, exact original argument arrays, and sequence hashes.
- [results/audit.md](results/audit.md) and [JSON](results/audit.json): original audit of coverage, hashes, statistics, and output comparisons.
- [profiles/validated-outputs.json](profiles/validated-outputs.json): four CPU-profile and four heap-profile output comparisons.
- [results/stream-cleanup-probe.json](results/stream-cleanup-probe.json): the descriptor-retention result.
- [council/](council/): workload, methodology, and source reviews. Early planning notes describe proposed work; the report and final reviews state the completed scope.
- [intarna-profile-evidence.zip](intarna-profile-evidence.zip): 1,029 evidence files plus their SHA-256 manifest. This is the original immutable evidence archive, including inputs, raw stdout/stderr/resource records, scripts, text profile exports, build/test logs, and C/C++ sources. [bundle-info.json](bundle-info.json) records its size and digest.

Large gprofng experiment directories and compiled binaries remain local and are not in the archive. Original commands retain the capture machine's absolute paths for provenance. The original report inside the archive is unchanged; this PR's Markdown copy clarifies measured-versus-warmup counts and adjusts links for GitHub. The PDF references paths inside the extracted archive.

## Verify the saved evidence without rerunning IntaRNA

From the repository root, using Python 3, `unzip`, and `sha256sum`:

```bash
report_dir="$PWD/doc/profiling/2026-09-19"
(cd "$report_dir" && sha256sum -c SHA256SUMS)
review_dir=$(mktemp -d)
unzip -q "$report_dir/intarna-profile-evidence.zip" -d "$review_dir"
python3 "$review_dir/profiling/scripts/audit.py"
python3 "$review_dir/profiling/scripts/profile_audit.py"
```

The first audit must report `pass`, 56 cases, 224 records, zero integrity errors, and zero incomplete-coverage findings. One expected warning says the executable is unavailable: binaries are deliberately excluded, while the recorded per-run executable hashes remain auditable. The second command prints eight comparisons, all with `output_equal: true`. The audit writes its derived results only inside the temporary extraction directory; the committed archive remains unchanged. These checks validate saved evidence, not a new execution or a cross-version performance result.

Before publication, the archive's full manifest was also checked: all 1,029 file sizes and SHA-256 values matched. The original fresh build passed 4,281 assertions in 36 API cases and all 20 CLI golden cases; the archive retains those logs. Application tests were not rerun for this documentation-only PR.

## Reproduce the experiment in a new directory

Preserve the existing measurements. Extract the archive as above, then copy only its scripts into a fresh directory. The scripts require the recorded conda compiler/dependency environment, GNU time, GNU gprofng, autotools, make, and Python 3; plot generation also requires Matplotlib. `INTARNA_DEPS` must point to an equivalent dependency prefix containing the compiler executables and libraries. Different toolchains or machines constitute a new experiment.

```bash
# Start at the repository root, after setting review_dir above.
export INTARNA_REPO="$PWD"
export INTARNA_REVISION=1b8377989ba79a0f1b8a450163b71bd8ed7c3e45
export INTARNA_DEPS=/path/to/the/recorded-or-equivalent/conda-prefix
experiment_dir=$(mktemp -d)
cp -R "$review_dir/profiling/scripts" "$experiment_dir/scripts"
cd "$experiment_dir"
bash scripts/build.sh
python3 scripts/benchmark.py prepare
python3 scripts/benchmark.py run
python3 scripts/supplement.py
python3 scripts/structural_checks.py
python3 scripts/profile.py
python3 scripts/profile_audit.py
python3 scripts/phase_summary.py
bash scripts/stream_cleanup_probe.sh
python3 scripts/benchmark.py summarize
python3 scripts/audit.py
python3 scripts/plots.py
```

Run stages sequentially on an otherwise idle host. Each timed case has one excluded warmup and three measured repetitions in deterministic randomized rounds. Profiling, build/tests, and diagnostic logging are separate from clean timing. The baseline fixes `LC_ALL=C`, `OPENBLAS_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, and `OMP_PLACES=cores`; two supplementary BLAS probes intentionally vary the BLAS setting.

The scope is one laptop, one build, and a bounded workload suite. Three repetitions support descriptive medians and ranges, not tight confidence intervals. CPU sampling had a collector timer warning and incomplete optimized ancestry; hardware counters were unavailable. Window comparisons check reported coordinates and energy, while full base-pair equality was separately checked for thread scaling and the matched ED comparison. See the report for all limitations.
