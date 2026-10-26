# Heuristic ensemble prediction with shared seed-anchor tables

This change implements the [PR #270 request](https://github.com/BackofenLab/IntaRNA/pull/270#issuecomment-6216756696)
to reuse exact prediction's left-most seed-anchor extension tables for heuristic
prediction, combining one chosen right endpoint with the possible left
extensions. The implementation starts directly from parent revision
`3866c62d10846fab1adc38b937d0c9b96cc06e51`, which is also the measurement baseline.

## Recurrence and retained ensemble

Exact and heuristic stack-seed prediction now use the same left and right
tables, contraction, complete-boundary publication, and outside pass. For an
anchor starting at `s` and ending at `t`, heuristic mode retains a complete
path with boundaries `(p,r)` when

```text
p == s OR r == t OR r == selectedRight(anchor).
```

Thus a true left extension joins the seed ending and at most one selected right
endpoint. A source at the anchor start still retains every right endpoint.
The seed-ending branch is counted once when the selected endpoint equals it;
invalid or out-of-span choices add no two-sided extension. Each selected
endpoint retains all legal paths to that endpoint, not a single traceback.
These are the existing heuristic ensemble rules.

Removing containing seed occurrences preserves ownership as well as the exact
union: an inner seed completes earlier, or wins the shorter co-ending tie.
Among the remaining occurrences on a common path, start and end order agree.
The left-most irredundant anchor is therefore the original first-completed
owner, including for nested, overlapping, and mixed-length explicit seeds.

Choices are frozen once per canonical anchor. Native `fillHybridZ_right` and
`updateOptRightZ` customization hooks still choose endpoints using the original
admitted seed family. Their ranking, tie handling, terminal-GU filtering, and
native singleton-seed rejection remain in place. The kernel's singleton
`noLP` states preserve support supplied by either extension.

The shared workspace uses bounded batches of complete left boundaries. A batch
without the anchor's own start needs a right table only up to its selected
endpoint, or just the seed-ending cell if that endpoint is unavailable. Both
passes join the same retained roots and use the same complete boundary factors.
This reduces heuristic contraction work; one-sided paths can still require a
full right table. See the [recurrence description](seed-extension-partition.md).

For ordinary prediction without pair-probability output, each batch may use up
to 16 MiB for its only complete-boundary table, still capped at 256 sources.
This reduces how often anchors crossing batches are recomputed. Candidate
discovery scans admitted anchors directly when they are sparser than a span
rectangle and retains the original grid lookup for dense seed domains. Exact
and pair-output computations keep their 2 MiB-per-table target.

## Correctness checks

Release and debug validation each passed all seven test drivers. Both ran 115
API cases, with 457,827 release and 457,828 debug assertions. Added regressions
compare with enumeration of the full declared seed family, covering ownership,
one-sided paths, invalid selections, singleton support, native extension hooks,
and complete-boundary coefficients across more than 256 sources.

The isolated comparison checked every boundary hybrid, boundary coefficient,
partition, and raw pair numerator. All 328,236 exact entries were byte-identical.
The 271,006 heuristic entries had maximum relative difference
`6.622292098580306e-15`, below the `1e-12` comparison tolerance. The changed
heuristic factorization changes floating-point summation order.

The 12 native verification predictions were byte-identical between versions.
The largest relative difference in native pair-probability output was
`2.3073317770873946e-15`; matrix labels and `NA` cells matched exactly.
Installation, standalone compilation of the three changed public headers, and
the installed pkg-config probability example also passed.
Validation was performed on Linux; macOS/Apple Clang was not run locally.

## Measurements

Measurements on 2026-10-16 used GCC 14.4.0, ViennaRNA 2.7.2, Boost 1.85,
and an AMD Ryzen 5 7530U running Linux. Both binaries use the C++23 release
configuration with `-O3` and `NDEBUG`; the manifest records full build flags.
Native runs use `--threads=1` and `OMP_NUM_THREADS=1`.
No project builds or tests ran during timing. CPU placement, frequency, and
other workstation activity were uncontrolled. Results describe these inputs
and this machine; they do not establish a general speedup.

The native workload uses the complete 109-nt OxyS and 113-nt fhlA sequences
recorded in the harness. Five fresh-process samples per variant alternate
baseline and candidate order. Ordinary CLI output is measured separately from
the focused pair-output runs. GNU time includes accessibility, prediction, and
serialization, recording elapsed/user/system time and whole-process maximum RSS.
Below, times are median elapsed seconds; ratios are candidate/baseline.

| Native case | Baseline seconds | Candidate seconds | Ratio | Median RSS KiB, baseline → candidate |
| --- | ---: | ---: | ---: | ---: |
| P/M seeded | 1.46 | 1.45 | 0.993 | 19,888 → 19,804 |
| P/M unseeded | 12.44 | 12.47 | 1.002 | 16,860 → 16,988 |
| P/H seeded | 14.59 | 1.01 | 0.069 | 17,280 → 19,928 |
| P/H unseeded | 0.23 | 0.23 | 1.000 | 30,512 → 30,604 |
| P/M seeded, pair output | 3.83 | 3.80 | 0.992 | 23,580 → 23,568 |
| P/H seeded, pair output | 42.51 | 2.71 | 0.064 | 18,196 → 23,700 |

Seeded H drops from 14.59 to 1.01 seconds median wall time (14.4 times faster),
compared with 1.45 seconds for candidate seeded M. Median user time drops from
14.56 to 1.00 seconds. H samples range from 0.99 to 1.03 seconds, versus
14.29 to 14.76 seconds before. Its median process RSS increases from 17,280
to 19,928 KiB, similar to the exact predictor's 19,804 KiB. All 40 ordinary
CLI samples have identical prediction output.

The 20 pair-output samples also retain identical prediction output. Their
maximum probability difference is `8.59e-15` relative (`1.55e-15` absolute),
within the `2e-12` relative tolerance. Exact pair CSVs remain byte-identical.

The isolated driver uses the four documented dense/banded, few/many-seed
fixtures, each with partition-only and pair-output execution. Its eight
heuristic cases measure the kernel without native energy/accessibility setup.
Allocation counts and cumulative requested bytes describe allocation traffic;
GNU time RSS includes the whole process and is not a live-kernel memory bound.

| Heuristic fixture | Partition seconds, baseline → candidate | Ratio | Pair seconds, baseline → candidate | Ratio |
| --- | ---: | ---: | ---: | ---: |
| dense40-many | 0.424553 → 0.127418 | 0.300 | 2.665865 → 0.417940 | 0.157 |
| dense40-few | 0.059757 → 0.003792 | 0.063 | 0.131882 → 0.005355 | 0.041 |
| band60-many | 0.056363 → 0.017060 | 0.303 | 0.292022 → 0.053131 | 0.182 |
| band60-few | 0.011505 → 0.001943 | 0.169 | 0.024049 → 0.002943 | 0.122 |

The heuristic kernel improves in all eight measured cases, by 3.3–24.6 times.
Using exact mode's batched boundary tables increases cumulative requested bytes
from 88,656–222,760 to 1,261,555–1,426,793 for partition-only runs and from
147,600–281,640 to 2,518,515–2,683,753 for pair output. Median process RSS rises
from 3,840–4,352 to 4,480–6,144 KiB. This is the shared exact workspace footprint
plus frozen heuristic endpoints; it is a time/space tradeoff against the old
heuristic kernel, not a memory reduction.

Exact kernel allocation counts and requested bytes are unchanged. Its observed
candidate/baseline median runtime ratios range from 0.745 to 0.991, with
byte-identical numerical output. Small timing changes in the native exact and
unseeded controls should not be interpreted as improvements to those algorithms.

Raw samples and summaries are retained as `anchored-native-*`, `anchored-pairs-*`,
`anchored-exact-*`, and `anchored-heuristic-*` files in the
[benchmark directory](../benchmarks/seed-extension/).
[`anchored-source-hashes.json`](../benchmarks/seed-extension/anchored-source-hashes.json)
records source/binary hashes, commands, toolchain, environment, and native run
manifests. Sample ranges and unchanged exact/unseeded controls should be read
alongside the medians when assessing noise or regressions.

## Reproduction

Build native parent and candidate executables with the same release toolchain
and dependencies. From a configured candidate checkout, use new output folders:

```sh
python3 doc/benchmarks/seed-extension/anchored-native.py \
  --binary /path/to/parent/IntaRNA --candidate /path/to/candidate/IntaRNA \
  --output-dir /tmp/anchored-native --repetitions 5
python3 doc/benchmarks/seed-extension/anchored-native.py \
  --binary /path/to/parent/IntaRNA --candidate /path/to/candidate/IntaRNA \
  --output-dir /tmp/anchored-pairs --repetitions 5 \
  --cases M-seeded H-seeded --pair-output
python3 doc/benchmarks/seed-extension/build.py \
  --source "$PWD" --output /tmp/anchored-kernel --cxx g++-14 \
  --baseline 3866c62d10846fab1adc38b937d0c9b96cc06e51
python3 doc/benchmarks/seed-extension/run.py --build-dir /tmp/anchored-kernel \
  --output-dir /tmp/anchored-exact --repetitions 5
python3 doc/benchmarks/seed-extension/run.py --build-dir /tmp/anchored-kernel \
  --output-dir /tmp/anchored-heuristic --repetitions 5 --heuristic
```

Dependency include/runtime options and fixture definitions are documented in
the [benchmark guide](../benchmarks/seed-extension/README.md).
