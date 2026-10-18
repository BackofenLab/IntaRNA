#!/usr/bin/env python3
"""Render the report from audited measurements and reviewed profile interpretation."""
import json,re,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=json.loads((ROOT/'results/summary.json').read_text());D={r['case']:r for r in S}
ENV=json.loads((ROOT/'environment.json').read_text())
RUNS=[json.loads(x) for x in (ROOT/'results/validated-runs.jsonl').read_text().splitlines()]
def t(name):return D[name]['median_s']
def rss(name):return D[name]['median_rss_mib']
def table(names):
    lines=['| Case | Median wall s | Min–max s | Median user / system s | Median peak RSS MiB | Output rows |',
           '|---|---:|---:|---:|---:|---:|']
    for name in names:
        r=D[name]
        lines.append(f"| `{name}` | {r['median_s']:.3f} | {r['min_s']:.3f}–{r['max_s']:.3f} | {r['median_user_s']:.2f} / {r['median_system_s']:.2f} | {r['median_rss_mib']:.2f} | {r['rows']} |")
    return '\n'.join(lines)
thread=['| Threads | Median s | Speedup | Evaluated pairs/s | Peak RSS MiB |','|---:|---:|---:|---:|---:|']
for n in [1,2,4,6,12]:
    name=f'throughput_{n}t';thread.append(f'| {n} | {t(name):.3f} | {t("throughput_1t")/t(name):.2f}× | {48/t(name):.2f} | {rss(name):.2f} |')
cpu=['| Profile | Nominal 100 ms samples | Largest exclusive CPU samples |','|---|---:|---|']
for name in ['biological_batch','hybrid_noaccess','exact_x','ensemble_m']:
    text=(ROOT/'profiles'/f'{name}_self.stdout').read_text()
    parsed=[]
    for line in text.splitlines():
        m=re.match(r'\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)',line)
        if m:parsed.append(m.groups())
    total=float(parsed[0][0]);top=[]
    for own,pct,incl,ipct,fun in parsed[1:4]:
        fun=fun.split('(')[0].replace('IntaRNA::','')
        top.append(f'`{fun}` {float(pct):.1f}%')
    cpu.append(f'| {name} | {round(total/.1)} | '+ '; '.join(top)+' |')
review=(ROOT/'council/profile-results-review.md').read_text()
review=review.split('## Optimization priorities supported by the profiles')[0].strip()
review=review.replace('## ', '### ')
doc=f'''# IntaRNA profiling report

Profiled 2026-09-19. **56 workloads, 224 clean runs, four CPU profiles, four heap traces, and separate phase and structural checks.** The fresh build passed its regression tests; all timed commands and the independent evidence audit passed. The original checkout was left unchanged.

The strongest findings are a **{t('throughput_1t')/t('throughput_12t'):.2f}× batch speedup at 12 threads**, a default-mode CPU hotspot in internal-loop energy evaluation, and a substantial time/memory tradeoff from windowing. **An input-stream cleanup defect is confirmed: 100 open/cleanup calls leave 100 file descriptors open.** Reusing accessibility data saved **{100*(1-t('cache_matched_reused')/t('cache_matched_computed')):.1f}% latency** in a matched test, but the current ED reader lowers the available interaction limit by one nucleotide. Exact and ensemble modes have different hotspot distributions and need separate optimization work.

![Measured performance overview](results/overview.png)

## Scope and experimental controls

The council comprised workload, hotspot/tooling, and experimental-methodology reviewers. They exchanged proposals and challenged stale build provenance, thread scaling, accessibility-limit confounding, sampling validity, output equivalence and memory interpretation. Their recommendations are retained in [workloads](council/workloads.md), [hotspots](council/hotspots.md), and [methodology](council/methodology.md). An independent [final audit](council/final-review.md) checked the actual evidence.

Target: `refactor-to-c++23`, commit `1b8377989ba79a0f1b8a450163b71bd8ed7c3e45`, reporting IntaRNA 3.4.1. A clean Git archive was built in `profiling/source`; the pre-existing executable was used only for tooling probes. Build: GCC 16.1.0, C++23, release `-O3`, `-g -fno-omit-frame-pointer`, OpenMP, ViennaRNA 2.7.2, Boost 1.85.0. This is an optimized build with profiling-friendly stack frames, not a debug `-O0` build. Binary SHA-256: `{ENV['binary_sha256']}`.

Host: AMD Ryzen 5 7530U, six physical cores/twelve logical CPUs, approximately 30.7 GiB RAM, Linux 6.8.0-139, performance governor. All 12 logical CPUs were available; recorded cgroup ancestors imposed no CPU quota or memory maximum. Dynamic frequency and other desktop activity were not disabled. See [environment.json](environment.json) for exact provenance and linked dependencies.

Each case had one excluded warmup and three measured runs. Runs were sequential, with deterministic randomized order in each round; compilation, tests, CPU sampling, heap tracing and verbose logs were separate. The main matrix and later controls have separate randomization blocks. Wall time covers process startup, input, computation and output. GNU time records user/system CPU, process peak RSS, faults, context switches and filesystem I/O. Medians and min–max ranges are reported; three repetitions do not support strong tail-latency or population confidence claims.

The baseline fixes `OPENBLAS_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, and `LC_ALL=C`. The BLAS experiment explicitly varies that setting. The linked dependency stack otherwise starts extra BLAS workers independently of IntaRNA's thread option.

Inputs include repository biological examples, the first 48 real 300-nt target records in `NC_000913.fa`, and deterministic synthetic sequences (seed `20260919`). That FASTA is a collection of targets, not one contiguous genome. Length sweeps use nested sequence prefixes. GC-distribution and repeat cases are individual controlled examples, not a representative population of all RNAs. Actual lengths, composition and hashes are in [input manifest](inputs/manifest.json); every command is in [cases.json](cases.json).

## Runtime and scaling

### Biological examples

{table(['bio_fhla_oxys','bio_phob_gcvb','bio_ilve_gcvb'])}

These inputs are respectively 112×108, 299×201 and 299×200 nt. The last two have almost identical lengths yet different timing, illustrating sequence dependence.

With a fixed 100-nt query and default local accessibility settings, target lengths 100, 300, 1,000, 3,000 and 10,000 nt take {', '.join(f"{t('length_t'+str(n)+'_q100'):.3f}" for n in [100,300,1000,3000,10000])} seconds. The 1,000→10,000 range increases time {t('length_t10000_q100')/t('length_t1000_q100'):.2f}× for 10× target length. This describes the measured regime with bounded folding/interaction windows; it is not an asymptotic complexity proof.

At 300×100 nt, the composition cases span {t('no_pair'):.3f} seconds for nonpairing all-A sequences to {t('gc_repeat'):.3f} seconds for GC repeats. GC repeats take {t('gc_repeat')/t('composition_gc50'):.1f}× the random 50%-GC-distribution case. Seed density and placement matter in addition to sequence length. Query-length tests intentionally include queries longer than the target; no orientation-invariance claim is made.

### Algorithm choices

On the same 100×100 synthetic input, default X/H takes {t('algorithm_default'):.3f} s, exact X/M {t('algorithm_exact_x'):.3f} s, heuristic S/H {t('algorithm_heuristic_s'):.3f} s, exact S/M {t('algorithm_exact_s'):.3f} s, and helix B/H {t('algorithm_helix'):.3f} s. Exact S is {t('algorithm_exact_s')/t('algorithm_default'):.0f}× the default on this case. These modes implement different prediction guarantees and models; their timing differences are not evidence of interchangeable answers.

The suite also covers ensemble P/H and P/M on bounded 40×40 input, seed-only mode, gapped seeds, no lonely pairs, loop limits and suboptimal output. `--noSeed` changes X to S and is labelled accordingly in the council review. Requests for 10 and 100 outputs both return the same ten interactions here; their similar cost does not establish that arbitrary output enumeration is cheap. Seven fixed CSV columns avoid accidentally requesting extra ensemble computations through `--outCsvCols=*`.

### Parallel throughput

{chr(10).join(thread)}

The denominator is **48 evaluated target–query pairs**, with ChiX length84; 46 have reportable interactions. All five thread counts produce the same reported coordinates and energies after sorting. A separate comparison including complete base-pair lists confirms equality at one and twelve threads.

Six threads deliver {100*t('throughput_12t')/t('throughput_6t'):.1f}% of twelve-thread throughput using {100*rss('throughput_6t')/rss('throughput_12t'):.1f}% of its RSS. Twelve was fastest among the tested settings. This is throughput from independent targets; a single unwindowed pair does not have the same parallelism. Source scheduling prioritizes targets, then queries, then windows.

## Accessibility and memory

For the 1,000×100 input with an explicit 150-nt interaction cap, local accessibility takes {t('accessibility_computed'):.3f} s and {rss('accessibility_computed'):.1f} MiB. Global accessibility takes {t('accessibility_global'):.3f} s and {rss('accessibility_global'):.1f} MiB. Disabling accessibility takes {t('accessibility_none'):.3f} s, despite removing folding. Accessibility changes energies, admissible boundaries and prediction work; subtracting these runtimes would not isolate folding cost. Separate phase logs provide that attribution.

Reducing accessibility window/span to75/75 also changes the scientific model. A separate75-nt interaction-cap control retains the default150/100 folding settings. Both are recorded so their differences cannot be mistaken for a pure storage optimization.

### Reusing saved accessibility data

The matched target/query limits149/99 give {t('cache_matched_computed'):.3f} s when computing accessibility and {t('cache_matched_reused'):.3f} s when reloading ED: **{t('cache_matched_computed')/t('cache_matched_reused'):.2f}× speedup**, with RSS {rss('cache_matched_computed'):.1f}→{rss('cache_matched_reused'):.1f} MiB. Coordinates, energies and base-pair lists match. The separately recorded one-time generation took0.847 s including a complete prediction and serialization; amortize that cost when planning reuse.

**Cache semantic caveat:** the writer emits columns through its configured maximum, while `AccessibilityFromStream.cpp:90–100` advertises one less than the final header length, with a source comment about dangling ends. The saved files contain their last columns; they are not truncated. Original100/150-limit round-trip runs therefore become99/149. Their output happens to agree on this example, but only the explicitly matched controls support the cache timing claim. Do not assume this round trip preserves every requested interaction limit. See the [source and file review](council/accessibility_cache_review.md). No application fix was made.

### Prediction windowing

| Mode / 3,000×800 nt input | Whole input seconds | Window300 seconds | Whole peak MiB | Window300 peak MiB |
|---|---:|---:|---:|---:|
| X/H | {t('memory_t3000_q800'):.3f} | {t('memory_window300'):.3f} | {rss('memory_t3000_q800'):.2f} | {rss('memory_window300'):.2f} |
| S/H | {t('memory_s'):.3f} | {t('memory_s_window300'):.3f} | {rss('memory_s'):.2f} | {rss('memory_s_window300'):.2f} |

Overlap is150 and the effective interaction cap is150. X/H's extension matrices are already bounded by interaction length (at most144×144 per side with the seven-pair seed here). Whole-sequence accessibility still contributes to peak RSS. Windows add95 target/query window combinations and repeat work, explaining the {t('memory_window300')/t('memory_t3000_q800'):.2f}× time without useful RSS reduction.

S/H keeps two matrices covering the full target/query ranges, so windowing reduces RSS by **{100*(1-rss('memory_s_window300')/rss('memory_s')):.1f}%**, while taking {t('memory_s_window300')/t('memory_s'):.2f}× longer. All tested window comparisons preserve the seven reported CSV fields; full structural equivalence was not separately tested for these window cases. The [window review](council/workload_results_review.md) links these observations to source.

## CPU hotspots and allocations

CPU workloads are a48-target biological batch with GcvB; a3,000×400 default-X pair with accessibility disabled and cap150; 128 repeated100-nt targets for exact X; and32 repeated60-nt targets for exact ensemble P without seeds. Repeated copies make sampling long enough to resolve coarse rankings. These profiles characterize distinct paths, not mutually comparable biological workloads.

The user-space collector emits a timer-reset warning on this host. A nominal100ms interval gave sampled totals close to recorded process CPU; shorter nominal intervals failed that calibration in tooling probes. Percentages remain coarse sampled evidence. Optimized stack unwinding misses some ancestors, so inclusive rows must not be added or interpreted as exhaustive phase fractions. The complete tables and warning text remain in `profiles/`.

PROFILE_REVIEW_INSERT

## Recommended fixes and optimization priorities

1. **Fix input-stream cleanup.** Correct the stream-type handling in [general.cpp](source/src/IntaRNA/general.cpp) at line 136 and add a repeated-open/cleanup regression check. The supplied probe confirms retained descriptors against the freshly built library; this is a concrete resource defect, separate from speculative performance improvements.
2. **Default X/H: reduce repeated internal-loop energy work.** `InteractionEnergyVrna::getE_interLeft` is the strongest measured CPU hotspot. Inspect repeated sequence-code and base-pair-type lookups, invariant loop bounds, and validation in the nested extension loops before changing the algorithm. Keep exact energy and structural regression checks; no speedup from such a change has been measured here.
3. **Exact X: reduce repeated boundary-energy and dangling-probability work.** Repeated ED lookups and exponentials have a different cost profile from default X. Evaluate caching or hoisting only with a numerical-correctness check and a measured memory tradeoff; these calculations contribute to the thermodynamic score.
4. **Reuse work across predictions.** Batched targets demonstrably parallelize; matched precomputed ED reuse helps repeated runs. Address or explicitly control the reader/writer limit asymmetry before treating caching as transparent.
5. **Ensemble P: measure recurrence and partition-map changes separately.** `fillHybridZ`, loop energy evaluation, exponentials and boundary-map maintenance contribute distinct costs. Preserve the ensemble definition and numerical behavior when evaluating changes.
6. **Treat matrix storage and windowing by algorithm.** The measured RSS benefit is strong for S/H and negligible for this bounded X/H workload. Heap traces and source review identify allocation reuse candidates, but allocation traffic is not itself proof of a dominant runtime bottleneck.
7. **Keep the benchmark harness as a regression baseline.** Compare proposed patches against the same binary provenance, cases, environment, warmups, repetitions and output checks. Report algorithm changes separately from implementation speedups.

## Validation and limits

- Fresh regression tests:4,281 assertions in36 API cases, plus all20 CLI golden cases. No tracked source was changed. Original untracked executables and backup remain intact.
- All56 benchmark cases completed their warmup and three measured runs. Every case's output was stable. The independent [audit](results/audit.md) verifies224 records, original stdout/input/binary hashes, resource data, commands, recalculated summaries and structural comparisons.
- Four original no-seed records contain a normal informational log line on stdout. Derived validation excludes that line; the raw evidence is retained. Profiler output comparison similarly excludes only the collector's recognized experiment-creation banner and IntaRNA log lines.
- Excluding the two intentionally varied BLAS cases, the largest min–max range is4.33% of the median. BLAS-unset and BLAS-12 are noisy: their ranges are {D['blas_unset']['min_s']:.3f}–{D['blas_unset']['max_s']:.3f} s and {D['blas_twelve']['min_s']:.3f}–{D['blas_twelve']['max_s']:.3f} s, versus {t('blas_one_control'):.3f} s with BLAS1. They demonstrate startup sensitivity, not a precise universal penalty.
- Hardware cycles, IPC, cache misses and branch misses were unavailable: `perf_event_paranoid=4` blocked perf. Syscall tracing was denied; gprofng I/O tracing failed its initialization probe. GNU time resource counters and process CPU evidence are available, but no detailed disk/network/syscall bottleneck claim is made.
- Heap tracing reports allocated bytes, tracked high-water usage and allocations outstanding at process exit. It is not a memory-safety or leak proof. RSS includes additional process mappings and allocator state. Instrumented runs supply attribution, never the clean timing medians.
- This is one build on one laptop-class CPU, with a bounded workload suite and warm page caches. It excludes cross-PR/version comparisons, a full4,319-target scan, all scientific parameter combinations, NUMA/multi-host scaling, energy-consumption measurement and biological prediction-accuracy benchmarking. Large exact/ensemble runs were deliberately bounded. No universal worst-case or tail-latency claim is made.

## Complete clean timing results

Every row below has three measured repetitions; warmups are excluded. "Peak RSS" is the median of the process high-water marks; maximum observed values are also in [summary.csv](results/summary.csv). Exact arguments and experimental controls are in [cases.json](cases.json).

{table([r['case'] for r in S])}

## Reproduction and evidence

See [README.md](README.md) for the sequential build/benchmark/profile commands and artifact map. The machine-readable [summary](results/summary.json), [validated runs](results/validated-runs.jsonl), [raw runs](results/runs.jsonl), [CPU/heap output checks](profiles/validated-outputs.json), [phase duration bounds](profiles/phase-summary.json), and [council reviews](council/) preserve the basis for the findings. The application source and fresh binaries are retained locally under `source/`.
'''
doc=doc.replace('PROFILE_REVIEW_INSERT',review)
for old,new in {
    'ChiX length84':'ChiX length 84', 'to75/75':'to 75/75', 'separate75-nt':'separate 75-nt',
    'default150/100':'default 150/100', 'limits149/99':'limits 149/99', 'took0.847':'took 0.847',
    'Original100/150':'Original 100/150', 'become99/149':'become 99/149', 'Overlap is150':'Overlap is 150',
    'cap is150':'cap is 150', 'most144×144':'most 144×144', 'add95':'add 95',
    'a48-target':'a 48-target', 'a3,000':'a 3,000', 'cap150':'cap 150', 'and32':'and 32',
    'repeated100-nt':'repeated 100-nt', 'repeated60-nt':'repeated 60-nt', 'nominal100ms':'nominal 100 ms',
    'tests:4,281':'tests: 4,281', 'in36 API':'in 36 API', 'all20 CLI':'all 20 CLI',
    'All56':'All 56', 'verifies224':'verifies 224', 'is4.33%':'is 4.33%', 'BLAS1':'BLAS 1',
    'full4,319':'full 4,319',
    'Windows add 95 target/query window combinations':'Windowing processes 95 target/query window combinations',
}.items():doc=doc.replace(old,new)
(ROOT/'REPORT.md').write_text(doc)
print(ROOT/'REPORT.md')
