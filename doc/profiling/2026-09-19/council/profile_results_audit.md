# Profile results audit

The final experiments use the isolated fresh optimized build, one IntaRNA thread, and `OPENBLAS_NUM_THREADS=1`. Timings below distinguish sampled CPU self attribution from clean reference elapsed time. The collector still warns that its interval timer changed to zero. The requested 100 ms interval produces CPU totals consistent with actual profiled process CPU time, unlike the initial shorter-period capability probes.

## Biological target batch

48 real 300-nt promoters against GcvB (201 nt), default model X/heuristic H.

- 368 samples imply 36.8 sampled CPU seconds. Collector process times are 36.780 user + 0.049 system seconds, agreeing closely in scale. Clean reference: 37.99 wall, 37.92 user seconds. Different runs are subject to normal host/thermal variance; the small lower profiled time is not a speedup caused by profiling.
- Exclusive CPU: `InteractionEnergyVrna::getE_interLeft` 57.61%; seed-extension left DP 12.23%; right DP 7.88%. These disjoint self rows account for about 77.7% of samples.
- Source interpretation: the extension recurrences repeatedly evaluate internal loops/stacking energies. `getE_interLeft` includes inlined ViennaRNA energy logic; the profile does not separate all inlined work into standalone library functions.
- Inclusive ancestry is incomplete: main/GOMP ancestors contain only 32.0/36.8 seconds. Therefore predictor inclusive 77.72% is not an exhaustive phase fraction. Do not call the remainder startup or use parent+child sums as independent costs.
- A separate `--verbose` run reports 49 accessibility calls totaling 3.815 seconds, and 48 predictor calls totaling 32.848 seconds. Nested seed filling totals 0.024 seconds and must not be added to predictor time. All these individual durations use milliseconds, avoiding the long-duration truncation caveat. Among these logged major phases, accessibility is approximately 10.4%, predictor 89.6%.
- Across the 48 same-length promoters, predictor times range 435–1145 ms while valid seed counts range 26–63. Descriptive Pearson correlation of seed count and predictor duration is 0.899. This supports the source mechanism of repeating extension DP per seed; it does not establish seed count as the sole predictor of work.

## Hybridization without accessibility

Synthetic target 3000 nt/query 400 nt, accessibility disabled, interaction length explicitly capped at 150.

- 296 samples imply 29.6 CPU seconds.
- Exclusive CPU: `getE_interLeft` 71.96%, left/right seed-extension DP 8.45% each. This reinforces energy evaluation as the main optimization candidate for this computation-heavy path.
- Inclusive predictor/main totals are 26.0/29.6 seconds, again demonstrating incomplete optimized ancestry. Use self attribution for ranking.

## Shared interpretation caveats

- Raw CPU stdout contains the gprofng `Creating experiment directory ...` banner. The initial raw comparison records `false`; a derived comparison must strip exactly that known banner (and intended `#` diagnostic lines) and compare prediction content. Preserve the raw artifacts.
- CPU profiles have hundreds, not millions, of samples. Report dominant percentages at sensible precision and do not overinterpret 0.1-second/sub-percent rows.
- Heap totals are cumulative allocated bytes, peak simultaneously tracked heap, and outstanding allocations at exit; none is equivalent to clean peak RSS. Outstanding allocations are not automatically application leaks.
- ELPP `DateTime::formatTime` integer-truncates durations at 1900 ms into whole seconds, then whole minutes/hours. Long phase strings are approximate; zero ms means below displayed resolution.
