# Intermittent energy differences in the default-mode ChiX batch

**Finding:** identical inputs, scientific parameters, and executable produced different energies at twelve threads, while coordinates and base-pair lists remained the same. This occurred in one original measured run and recurred once in 20 subsequent twelve-thread diagnostic runs. The root cause is not established. The affected timing row is an observed execution time, not a validated correctness-preserving speedup.

This uses the requested seed-extension heuristic, computed local accessibility, and no prediction windowing. It does **not** use cached ED; the separate ED reload observations cannot explain this as a cache-input mismatch.

## Preserved observations

| Collection | Target / query | One-thread reference energy | Observed twelve-thread energy | Other reported fields |
| --- | --- | ---: | ---: | --- |
| Original measured repetition 3 | `b3831` / `ChiX_NC_000913` | -8.9 | -6.6 | IDs, coordinates, and complete base-pair list agree |
| Additional repeat 08 | `b2410` / `ChiX_NC_000913` | -5.02 | -5.8 | IDs, coordinates, and complete base-pair list agree |

Energies are kcal/mol. Inspect the [one-thread reference](results/batch_chix_1t.measure1.stdout), [original differing run](results/batch_chix_12t.measure3.stdout), and [additional differing run](diagnostics/parallel/t12_repeat08.stdout). Output hashes and the exact binary hash are saved with each record; row order is canonicalized before comparison.

The original ChiX twelve-thread case has five measured repetitions and one excluded warmup; only measured repetition 3 differs. Original GcvB full output matches across all tested thread counts. All original one-thread measurements retain stable energies and structures.

## Focused follow-up

After all timing and profiler collection finished, 28 separate fresh-process checks used the same binary, inputs, and controlled environment. Their thread counts were interleaved with fixed shuffle seed `202609212`:

| Threads | Additional checks | Checks different from the one-thread reference |
| ---: | ---: | ---: |
| 1 | 3 | 0 |
| 6 | 5 | 0 |
| 12 | 20 | 1 |

These are diagnostic counts, not a reliable estimate of a population failure probability. No claim that six threads is universally safe follows from these few checks. The [repeat summary](diagnostics/parallel/summary.json) identifies the differing field; every additional stdout/stderr and command is in the evidence archive.

## Reproduce

Use the recorded executable at source revision `0a4568ad6e3da52f935221ff0da98a4a171e9833` and the archived inputs. Set `BIN` to that executable and `DATA` to the extracted `default-mode` directory. Repeat the following in fresh processes and compare sorted CSV rows against the same command with `--threads=1`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores LC_ALL=C \
"$BIN" --target="$DATA/inputs/promoters48_systematic.fa" \
  --query="$DATA/inputs/ChiX_NC_000913.fa" \
  --model=X --mode=H --acc=C --accW=150 --accL=100 --windowWidth=0 \
  --seedBP=7 --seedMaxUP=0 --intLenMax=0 --threads=12 \
  --outMode=C --outCsvCols=id1,id2,start1,end1,start2,end2,E,bpList
```

The issue is intermittent, so a single matching execution does not refute the saved observations. The archive's `scripts/parallel_check.py` automates the interleaved repeat check, records outputs without overwriting earlier evidence, and reports field-level differences.

## Interpretation and next step

The evidence establishes output instability associated with the tested parallel runs. It does not yet identify a particular race, library call, or source line as the cause. The CLI's [target loop](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/bin/IntaRNA.cpp#L171-L177) is parallelized; the CSV [energy field](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/src/IntaRNA/OutputHandlerCsv.cpp#L203-L205) prints the stored interaction energy. Those are investigation entry points, not a causal diagnosis.

Isolate and fix this in a separate application change with a repeated-output regression. Preserve this profiling revision and evidence so performance changes can be compared against an explicit parent. The current PR adds evidence and makes no application fix.
