# Workload council review of initial measurements

This note reviews warmup observations, not final repeated-run estimates. No profiling workloads were launched by this reviewer while the baseline run was active.

## Why prediction windows barely reduce default-X RSS

Observed initial warmups for target3000/query800:

| Setting | Wall seconds | Maximum RSS, KiB |
|---|---:|---:|
| Default X/H, no prediction windows | 41.21 | 24,872 |
| Window width300, overlap150 | 71.04 | 24,612 |

The windowed case completed before the 90-second timeout. It took approximately 1.724 times as long, with only approximately 1.0% lower peak RSS. Treat that small RSS difference cautiously; repeated results are authoritative.

Source explanation:

- `PredictorMfe2dHeuristicSeedExtension.cpp:81-100` bounds each per-seed left/right extension matrix by effective maximum interaction length minus seed length plus one. With the default150-nt cap and7-bp ungapped seed, each matrix is at most144×144 integer energies. The two cell arrays total approximately162KiB; this path does not allocate a full3000×800 DP matrix merely because the sequences have those lengths.
- `SeedHandlerNoBulge.h:29` and `SeedHandlerNoBulge.cpp:22` show a sparse hash of valid seeds, which is reset for each prediction. It does not store a dense seed matrix for every pair of sequence positions in this ungapped case.
- `IntaRNA.cpp:119` and `:190` compute whole-sequence query and target accessibilities before entering the prediction-window loop (`:239-281`). ED values are banded whole-sequence arrays (`AccessibilityVrna.h:89`; allocation at `AccessibilityVrna.cpp:86`). Prediction windows do not shrink this preprocessing or its peak. Maximum RSS records the whole-process high-water mark, including folding, even when later predictor buffers are smaller.
- `IndexRange.cpp:50-80` yields19 target windows and5 query windows at width300/overlap150. All95 combinations construct a fresh predictor and recompute their seed sets; overlapping candidates are revisited. Summed window lengths are5700 and1400, so the summed position-pair area is3.325 times the original3000×800 rectangle. This is a geometric work indicator, not a runtime prediction: boundary effects constrain seed extensions and accessibility is shared across windows. At width600 the corresponding count is7×2=14 combinations and pair-area multiplier1.54375.

Therefore these data support: "For this bounded seed-extension workload, prediction windows added repeated work without materially reducing peak process memory." They do not support a blanket claim that IntaRNA windowing never saves memory. Other predictors, gapped seeds, larger interaction limits, dense seed sets, and multiple threads can have different memory profiles. The case named `memory_stress` is a larger-input probe; observed24MiB peak means it was not a successful severe memory-pressure test on this machine.

## Result-count and equivalence interpretation

- Throughput is48 evaluated target/query pairs per invocation, not46. The48-record library produced46 reportable interactions under the chosen filtering and model. Report pairs/second as48/wall_seconds.
- All thread counts1,2,4,6,12 had identical canonical hashes in the reviewed warmups. The requested columns were `id1,id2,start1,end1,start2,end2,E`. This establishes identical reported identifiers, coordinates and energies after row sorting. It does **not** directly establish identical internal structures/base-pair lists, because these were omitted. If desired, a later untimed `bpList` comparison can strengthen the check.
- `outNumber=10` and `outNumber=100` produced the same10 output rows and same canonical hash in the initial warmups. This input saturates the reported heuristic results under the selected constraints. Similar time here does not show that producing100 interactions costs the same as producing10, and should not be generalized to other inputs.

## Remaining interpretation cautions

The matrix includes true-default baselines, where accessible interaction lengths can be limited by shorter input length, and explicit caps for controlled ablations. Describe the actual manifest settings rather than claiming that every matrix row used one universal150-nt cap. Composition and seed constraints can alter which candidates pass filters, so changes in runtime can reflect different explored work as well as arithmetic cost. Exact-vs-heuristic comparisons address computational behavior; no external biological truth set was used to assess predictive accuracy.
