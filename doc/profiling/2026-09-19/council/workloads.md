# Workload council: profiling coverage and interpretation

Scope: end-to-end runtime, process memory, algorithm/parameter sensitivity, throughput scaling, and phase attribution of the current source snapshot. This is a coverage plan, not a claim that every combination is measured; the final run manifest is authoritative. Source was inspected read-only; no production code changes are required.

## Existing biological inputs

| Input | Actual size | Use |
|---|---:|---|
| fhlA / OxyS | 112 / 108 nt | Documented interaction example and small exact-mode comparison |
| phoB / GcvB | 299 / 201 nt | Medium biological baseline |
| ilvE / GcvB.ST | 299 / 200 nt | Second medium biological pair, different organism sequence |
| ChiX | 84 nt | Query for a target-library screen |
| NC_000913.fa | 4,319 records, each 300 nt; 1,295,700 nt total | Real target library; use first 48 records for throughput |

These counts are parsed from the repository FASTA files. Despite its filename, NC_000913.fa is a library of target regions, **not** one contiguous chromosome. Do not report a first-48 subset as a whole-genome scan. Sequence provenance beyond the repository annotations is not independently verified; this study tests computational behavior rather than predictive biological accuracy.

## Agreed practical matrix

Use the default physical energy model and, for controlled comparisons, explicitly declare model X, mode H, seed BP 7, interaction length cap 150, local accessibility window 150/span 100, output count 1, CSV columns and thread count 1. Also keep true-default biological baselines if desired. State exceptions in every manifest row.

| Dimension | Suggested settings | Purpose |
|---|---|---|
| Biological anchors | Three single pairs listed above | Real-sequence end-to-end latency and phase split |
| Target-length scaling | Nested synthetic target prefixes 100,300,1000,3000 nt × fixed query100 | Isolate growth with target size |
| Query-length scaling | Fixed target300 × nested query50,100,200,400 | Test query size and input orientation; 400-nt query exceeds target length, so explicitly note it |
| Composition | Fixed lengths with GC20%,50%,80%, independently generated RNAs | Complementarity and seed-density effects |
| Sequence extremes | AU repeat, GC repeat, no-complement pair | Dense-candidate and sparse-candidate behavior; stress cases are not population averages |
| Predictors | X/H, X/M, S/H, S/M, B/H on small100 inputs | Distinct MFE DP implementations; exact methods restricted to small input |
| Ensemble prediction | P/H and P/M on 60–100 nt inputs | Distinct partition-function/ensemble DP paths; recommended addition to coverage |
| Seed behavior | noSeed, seedMaxUP2, seed-only mode S; optionally seedBP5/7/9 | Remove seed constraint, invoke gapped-seed DP, isolate seed enumeration |
| DP neighborhood | intLoopMax0,5,10,20 at fixed lengths; optionally outNoLP | Loop fanout/branch sensitivity, without changing input sizes |
| Accessibility | Computed vs disabled vs precached ED; local windows vs global | Folding contribution and reproducible precomputation tradeoff |
| Interaction extent | intLenMax30,60,150 | Bound extension work and compare user-tunable cost |
| Suboptimal/output | outNumber1,10,100 at fixed overlap B and columns | Additional predictions, traceback and serialization cost |
| Windows | Unwindowed vs widths300/600, overlap150, fixed cap150 | Runtime–memory tradeoff and overlapping-work overhead |
| Throughput | First48 real targets × ChiX, threads1,2,4,6,12 | Parallel independent predictions; aggregate throughput and RSS |

Do not run the full Cartesian product: vary one dimension against a defined baseline, plus a few explicit interactions (e.g. exact mode × length, seed constraint × composition). Run a small pilot before escalating expensive exact/global/dense-repeat cases. A 3000-nt global-accessibility run may be substantially more expensive than local folding; start at 300 or1000 with a timeout. Retain timed-out rows as censored observations and do not interpolate them as completed timings.

For synthetic length scaling, generate each longest sequence once using a fixed recorded RNG seed, then take nested prefixes. Record exact observed GC fractions, hashes and generation procedure. Composition comparisons confound secondary structure, complementarity and seed density by design and should be described as workload sensitivity rather than causal GC coefficients. A single synthetic sequence per fraction is not enough to generalize to sequence populations.

## Source-confirmed compatibility and confounders

- `--noSeed` changes model X to S internally (`CommandLineParsing.cpp`, seed option finalization). Label it S/no-seed, or explicitly supply S. This is a distinct algorithm, not merely skipping the seed prepass.
- Model B supports heuristic H, not exact M. Models X and S support H/M/S when seeded; P has its own H/M/S ensemble predictors. No-seed plus seed-only mode is invalid.
- Turning accessibility off otherwise removes the cap implied by accessibility window size. Keep an explicit `--intLenMax` identical in computed/disabled/cached comparisons. Accessibility-off changes energies and possibly winning interactions; the runtime difference is not a pure measurement of folding time.
- Larger local accessibility windows also change the folding ensemble and possibly interaction cap. For folding ablation, fix interaction length separately. Global `--accW=0 --accL=0` still needs an explicit interaction cap for a controlled comparison.
- Precached ED reuse must use the identical input, energy model, temperature, constraints and effective accessible length. Include precompute time separately; the steady-state improvement amortizes only over repeated reuse. Compare output energies and interaction coordinates within file serialization precision.
- Window CLI validation requires overlap >= max of the supplied per-sequence accW/intLenMax values (when accessibility is used), and width > overlap. The current CLI help permits equality even though the README says “larger.” Use overlap150 when all relevant limits <=150. Model P rejects window decomposition, and Zall/Eall output rejects windows because overlapping windows double-count interactions.
- Heuristic windowing can change results; record and explain differences rather than asserting scientific equivalence solely from overlap size.
- Parallelization is hierarchical over target records, then queries, then windows. A lone pair without windows cannot exercise thread scaling. The 48-record library is divisible by 1,2,4,6,12, but load varies by sequence; report CPU occupancy and speedup alongside wall time.
- Default target/query loop bounds differ (10 vs16); changing `--intLoopMax` sets both. Input swapping therefore is not a pure implementation symmetry test unless sequence-specific constraints are also swapped.
- Extra CSV fields such as Eall/Zall may request partition-function work. Keep ordinary output columns fixed and reserve expensive summary fields for an explicit experiment. `--outCsvCols=*` is not a neutral output formatting choice.
- Non-overlapping suboptimal output can select different paths/constraints; hold overlap B fixed in the basic outNumber study.

## Correctness and measurement guards

1. Fingerprint source tree, binary, compiler flags, linked ViennaRNA and Boost dependencies, input files and all command lines; rebuild current edits in an isolated source snapshot.
2. Run existing API/oracle tests plus CLI parameter/reference regressions. The CLI suite includes heuristic/exact, seed/no-seed, X/S, noLP, overlap and window behavior, but most fixtures use simplified energy B and disabled accessibility. Biological examples complement rather than replace these tests.
3. For every run record exit code, timeout, stderr and output hash. Preserve normalized CSV results; compare repeated/parallel runs after sorting rows so nondeterministic output order is not misclassified as changed predictions.
4. For equivalent execution options (thread count), require identical prediction rows. For scientific parameter changes (seed/accessibility/model/window), report differences rather than requiring equality. Exact-vs-heuristic differences are expected and informative.
5. Separate quiet wall-time/RSS runs from verbose9 phase timings and instrumented profiles. Instrumentation overhead and overlapping nested timers must not be interpreted as production runtime or summed naively.
6. Warm up once, then run at least three interleaved repetitions; report median, range and variation. Very short cases are startup dominated. Use sufficient independent predictions or repeat the entire process, retaining process-startup cost as explicitly part of CLI latency.
7. Serialize timed workloads across agents; do not overlap compilation, tests or profilers with measured runs. Record residual host noise and avoid overinterpreting small differences.

Council discussion: workload coverage exchanged with hotspot and methodology agents; consensus emphasizes biological plus controlled synthetic coverage, independent-pair thread scaling, fixed interaction caps in accessibility experiments, explicit output correctness checks, and source/phase evidence alongside wall-time measurements.
