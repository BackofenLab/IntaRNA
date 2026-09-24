## CPU hotspots differ by algorithm

Four separate 100 ms sampling runs captured 1,182 samples. All four prediction outputs matched their clean references after removing the collector banner. Sampled CPU totals closely matched the profiler's independently recorded process CPU time.

| Profile | Samples | Sampled CPU | Dominant exclusive CPU functions |
|---|---:|---:|---|
| 48 biological promoters × GcvB, default X/H | 368 | 36.8 s | `getE_interLeft` 57.6%; left extension 12.2%; right extension 7.9% |
| Synthetic 3000 × 400, X/H, accessibility off, interaction cap 150 | 296 | 29.6 s | `getE_interLeft` 72.0%; left extension 8.4%; right extension 8.4% |
| 128 repeated 100 × 100 pairs, exact X/M | 211 | 21.1 s | `getBoltzmannWeight` 19.4%; accessibility lookup `getED` 18.5%; full-energy `getE` 9.5%; unnamed libm routine 9.5% |
| 32 repeated 60 × 60 pairs, ensemble P/M without seed | 307 | 30.7 s | `getE_interLeft` 31.6%; ensemble DP `fillHybridZ` 22.8%; unnamed libm routine 12.1%; `updateZ` 8.1% |

The default predictor spends most of its CPU evaluating internal-loop/stacking energies inside repeated seed-extension recurrences. Exact X/M instead spends most inclusive CPU in candidate scoring (`updateOptima` 83.4%, including its children). Its four boundary-combination loops repeatedly retrieve accessibility energies and compute probability-weighted dangling contributions. The exponential calls are part of the MFE energy model, not evidence of an unnecessary partition-function calculation. Ensemble mode adds Boltzmann-weighted DP and a boundary-keyed partition map; its insertion, reporting and destruction also appear in the profile.

These are exclusive function shares, not additive call-tree phase percentages. Optimized call ancestry is incomplete: for example, `main` contains only 87% of biological samples. The unnamed libm routines remain unnamed; attribution should not be more specific than the symbols support.

## Phase measurements

Built-in timers were collected in separate diagnostic runs. Predictor time includes nested seed filling and reporting.

| Diagnostic workload | Accessibility calls, summed | Predictor calls, summed | Nested seed filling |
|---|---:|---:|---:|
| Biological batch | 3.815 s, 49 calls | 32.848 s, 48 calls | 0.024 s |
| Exact X/M repeated batch | 1.253 s, 129 calls | 19.694 s, 128 calls | 0.015 s |
| Ensemble P/M repeated batch | 0.079 s, 33 calls | 29.244 s, 32 calls | No seed |
| Hybridization without accessibility | Disabled | Displayed “29 seconds”: approximately 29–30 s | 0.009 s |

The first three rows sum millisecond-formatted durations; each individual call loses less than a millisecond through formatting. Longer timers truncate to whole seconds/minutes, so the last row is a range. Among the biological batch's logged major phases, about 10% is accessibility and 90% prediction. Its equal-length promoters still take 435–1,145 ms per prediction with 26–63 valid seeds. Seed count and predictor duration correlate strongly in this sample (Pearson r = 0.899), consistent with repeated extension DP per seed.

## Heap allocation and memory representation

Each heap trace matched its corresponding clean prediction output. Byte counts below are profiler-tracked allocations, not resident memory.

| Heap workload | Allocation calls | Cumulative allocated bytes | Peak tracked heap bytes | Outstanding at shutdown |
|---|---:|---:|---:|---:|
| phoB/GcvB, default X/H | 10,997 | 13,796,710 | 4,554,589 | 14,174 bytes / 20 allocations |
| 100 × 100, exact X/M | 6,392 | 4,250,556 | 1,693,265 | 13,672 bytes / 18 allocations |
| 3000 × 800, default X/H | 64,100 | 238,886,646 | 10,750,986 | 26,554 bytes / 29 allocations |
| 3000 × 800, S/H | 61,660 | 181,619,663 | 118,003,656 | 27,802 bytes / 29 allocations |

The default 3000 × 800 case allocates about 172.5 MB cumulatively at two extension-matrix callsites, across 1,231 allocations at each site: 72% of total allocated bytes. Source inspection identifies per-seed `resize` calls with Boost's default `preserve=true`, which allocates a temporary matrix and copies overlapping contents. This confirms substantial allocation churn, although allocation itself is not a leading sampled CPU hotspot.

S/H instead makes two 57,600,000-byte matrix allocations. Each represents 3000 × 800 cells with a 24-byte energy/right-boundary record; together they explain most of its 118 MB peak tracked heap. This is a structural memory difference between algorithms, not a leak.

## Confirmed input-stream cleanup defect

All four traces retain 9,384 bytes in input-stream allocation stacks. Source `src/IntaRNA/general.cpp:136` casts an input stream to `boost::iostreams::filtering_ostream`, although `newInputStream` creates `filtering_istream`. The cast fails, cleanup is skipped, and the pointer is discarded.

An isolated probe linked against the freshly built IntaRNA library confirmed resource retention: 100 open/read/delete-input-stream iterations increased open file descriptors from **4 to 104**. The probe and result are saved in `scripts/stream_cleanup_probe.cpp`, `scripts/stream_cleanup_probe.sh`, and `results/stream-cleanup-probe.json`. The OS reclaims these resources when the CLI exits; repeated use within a long-lived process accumulates them. Application source was not changed.

Other outstanding allocations include a 4,096-byte stdout buffer, 192 bytes in OpenMP, and ViennaRNA allocation chains. Their presence at shutdown alone does not establish additional application leaks.

## Optimization priorities supported by the profiles

1. **Fix input-stream ownership/type cleanup.** This is a verified correctness/resource issue independent of throughput. Validate with repeated file and compressed-file reads and descriptor counts.
2. **Default X/H: optimize repeated internal-loop energy evaluation.** Start with `InteractionEnergyVrna.h:381` and extension recurrences in `PredictorMfe2dHeuristicSeedExtension.cpp`. Investigate repeated sequence-code access, validity checks and reusable energy context; retain the current thermodynamic model and compare outputs.
3. **Exact X/M: reuse interval scoring terms.** `InteractionEnergy.h:883`, `:920`, and `:957` repeatedly derive dangling probabilities from single-sequence accessibility differences. Caching/reusing those terms and ED lookups is a more relevant target here than seed search. Bound cache memory and preserve temperature, constraints and rounding behavior.
4. **Reduce matrix allocation churn and investigate S/H cell storage.** Reuse extension storage where correctness permits; inspect nonpreserving resize because the DP initializes cells. For S/H, the two full matrices dominate peak heap. Algorithm/window changes also alter work or predictions and must be presented with their measured tradeoffs.
5. **Ensemble mode: profile recurrence and partition-map changes separately.** `PredictorMfeEns2d.cpp` and `PredictorMfeEns.cpp` identify the DP, Boltzmann conversion and boundary storage costs. Do not assume default-mode optimizations address the full ensemble workload.
