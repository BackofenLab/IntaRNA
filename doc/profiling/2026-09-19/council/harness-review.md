# Harness review before timing

Reviewed `scripts/benchmark.py`, `cases.json`, `inputs/manifest.json`, and `environment.json` before the first benchmark invocation. This was a read-only audit apart from this note.

## Verdict

The timing method and workload matrix are suitable for the planned bounded profile. The blocking subprocess wait plus watchdog avoids timeout-polling quantization. No fundamental timing blocker was identified. The coordinating agent received the robustness improvements below before measurements began.

## Verified

- All 47 case names are unique; no benchmark results existed when reviewed.
- Every generated FASTA matches its manifest SHA-256.
- Fresh binary SHA-256 matches `environment.json`.
- Throughput input contains exactly 48 target records, each 300 nt; its query is ChiX, 84 nt.
- Environment controls BLAS threads, OpenMP dynamic behavior, placement, binding, and locale consistently. All 12 logical CPUs are allowed; source/compiler/library provenance is recorded.
- One warmup precedes three randomized measured rounds. Warmup failures/timeouts are retained in raw records and skipped for repeated measurements.
- Accessibility on/off cases explicitly share the interaction-length cap; windowed cases satisfy source validation with overlap equal to 150 and cap 150.
- Source validates the used model/mode letters, implicit boolean switches, and window-overlap equality. The default is model X, mode H.

## Recommended robustness changes

1. Prevent accidental reruns from appending new records under old identities and overwriting the corresponding stdout/time files. The original harness reused `case.phaseN` paths and appended to one `runs.jsonl`. Guard against nonempty results or give each execution a unique run identity.
2. Record a UTC timestamp per run and the executable hash for the invocation. A single environment snapshot plus appendable results can otherwise lose provenance.
3. Validate the expected seven-column CSV header and row shape, preserve the header outside row sorting, and allow a header-only result for a legitimate no-interaction case. The original canonicalization sorted every line, including the header, and did not validate output schema.
4. Emit failure/coverage records alongside successful summaries. The original summary omitted cases with no successful measured runs and did not flag fewer than three successful repetitions; present those explicitly to prevent a survivor-only report.

The coordinating agent owns any harness changes; this review did not edit the harness.

## Interpretation constraints

- Disabling accessibility changes interaction energies and can change prediction choices. Its runtime difference cannot be described as pure accessibility cost; use the separate phase timer for attribution.
- `accessibility_window75` changes both folding-window width and maximum base-pair span to 75, compared with the default 150/100. The interaction-cap control is useful but does not isolate window width from base-pair span.
- Defaults impose an effective maximum interaction length of 150 when accessibility is enabled. For the long default/windowed comparisons, this agrees with the explicit cap of 150 in the windowed cases.
- The synthetic size families use nested sequence prefixes, reducing input confounding. The three GC cases each use one synthetic sequence pair and characterize sensitivity, not population averages.
- Tiny cases measure end-to-end CLI latency including process launch and dynamic linking. Kernel-only efficiency must come from instrumented longer cases.
- Returned row count is not the number of attempted predictions: the throughput denominator is the 48 target–query combinations, even if some have no reportable interaction.
- A timeout during warmup means the case is censored at its limit; it cannot be reported as a completed runtime or silently omitted.

Actual measurement and final-report audit remain pending.
