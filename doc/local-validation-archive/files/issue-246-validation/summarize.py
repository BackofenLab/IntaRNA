from pathlib import Path
import csv, json, statistics
validation=Path(__file__).resolve().parent
repo=validation.parent/'pr-worktrees/issue-246'
results=validation/'benchmark'
summary=json.loads((results/'summary.json').read_text())
metadata=json.loads((results/'metadata.json').read_text())
rows=[json.loads(line) for line in (results/'samples.jsonl').read_text().splitlines()]
fields=['case','variant','repetition','warmup','wall_s','user_s','system_s','max_rss_kib','sha256']
with (repo/'doc/mdspan-storage-samples.tsv').open('w') as handle:
 writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
 writer.writeheader();writer.writerows(rows)
lines=['\n## Same-host performance results\n',
       'Baseline: upstream `c5823e0`. Candidate: the storage migration in this change.\n',
       'Both are optimized builds using conda-forge GCC/libstdc++ 16.2.0, Boost 1.85.0\n'
       'and ViennaRNA 2.7.2 on Linux x86-64, AMD Ryzen 5 7530U. Configure adds\n'
       '`-O3 -fno-strict-aliasing`; OpenMP is enabled with `--threads=1`.\n'
       'Runs are pinned to CPU 2 with `LC_ALL=C`, `OPENBLAS_NUM_THREADS=1`,\n'
       '`OMP_DYNAMIC=FALSE`. No build or sanitizer runs overlap these timings.\n',
       'One warmup plus seven measured executions per binary/case: 160 executions\n'
       'total, all stdout results byte-identical within and between variants.\n'
       'The table shows median wall seconds and median peak RSS in MiB. A negative\n'
       'time change means faster. RSS is a process-level measurement and includes\n'
       'the executable, shared libraries and ViennaRNA allocations.\n',
       '| Workload | uBLAS time (s) | mdspan time (s) | Change | uBLAS / mdspan RSS (MiB) |',
       '| --- | ---: | ---: | ---: | ---: |']
for row in summary:
 b,c=row['baseline'],row['candidate']
 lines.append(f"| {row['case']} | {b['wall_s']:.4f} | {c['wall_s']:.4f} | {row['wall_change_percent']:+.1f}% | {b['max_rss_kib']/1024:.2f} / {c['max_rss_kib']/1024:.2f} |")
lines += ['\nMeasured wall-time ranges (min–max across the seven measured samples):\n',
          '| Workload | uBLAS (s) | mdspan (s) |','| --- | ---: | ---: |']
for row in summary:
 b,c=row['baseline'],row['candidate']
 lines.append(f"| {row['case']} | {b['wall_min_s']:.4f}–{b['wall_max_s']:.4f} | {c['wall_min_s']:.4f}–{c['wall_max_s']:.4f} |")
lines += ['\nAll samples (including explicitly marked warmups) and output hashes are in\n'
          '[mdspan-storage-samples.tsv](mdspan-storage-samples.tsv). Full invocation\n'
          'arguments and deterministic input generation are in the benchmark script.\n',
          'Binary SHA-256 values:\n']
for name,item in metadata['binaries'].items(): lines.append(f"- {name}: `{item['sha256']}`")
lines += ['\nThese measurements describe the selected workloads on one host and one\n'
          'toolchain. They do not establish a universal speedup. Short biological\n'
          'runs include substantial startup overhead; differences near the observed\n'
          'spread should be treated as inconclusive. Clang/macOS and other CPU\n'
          'architectures still need their own measurements.\n']
with (repo/'doc/mdspan-storage.md').open('a') as handle: handle.write('\n'.join(lines))
print('\n'.join(lines))
