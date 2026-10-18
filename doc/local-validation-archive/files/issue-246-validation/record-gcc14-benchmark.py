from pathlib import Path
import csv
import json

repo=Path('/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-246')
root=Path(__file__).resolve().parent/'benchmark-gcc14-kokkos'
summary=json.loads((root/'summary.json').read_text())
metadata=json.loads((root/'metadata.json').read_text())
samples=[json.loads(line) for line in (root/'samples.jsonl').read_text().splitlines()]
assert len(summary)==10 and len(samples)==160 and all(row['output_identical'] for row in summary)
columns=['case','variant','repetition','warmup','wall_s','user_s','system_s','max_rss_kib','sha256']
with (repo/'doc/mdspan-storage-gcc14-samples.tsv').open('w') as out:
    writer=csv.DictWriter(out,fieldnames=columns,delimiter='\t',lineterminator='\n')
    writer.writeheader()
    writer.writerows(samples)
changes=[x['wall_change_percent'] for x in summary]
text='''\n## Bundled Kokkos performance results (GCC 14, 2026-09-29)\n\nThe parent `c5823e0` and current candidate were rebuilt cleanly with the same\nGCC 14.4.0 compiler and libstdc++ headers, Boost 1.85.0, ViennaRNA 2.7.2,\n`-O3 -fno-strict-aliasing`, and OpenMP. Both link to the same dependency-prefix\nlibstdc++ 16.2 runtime; mdspan availability is determined by the GCC 14 headers.\nThe candidate automatically selects the pinned Kokkos backend. The host, CPU 2\naffinity, workload arguments, random seed, environment controls, one warmup and\nseven measured runs are as described above. No builds, header compilation or\nsanitizer checks overlap these timings.\n\nAll 160 executions produce matching stdout for their respective workloads.\n'''
text+=f'Median wall-time changes range from {min(changes):+.1f}% to {max(changes):+.1f}%.\n'
text+='These observations apply to this host and workload set; they do not establish\na universal performance improvement.\n\n'
text+='| Workload | uBLAS time (s) | Kokkos time (s) | Change | uBLAS / Kokkos RSS (MiB) |\n| --- | ---: | ---: | ---: | ---: |\n'
for row in summary:
    b,c=row['baseline'],row['candidate']
    text+=f"| {row['case']} | {b['wall_s']:.4f} | {c['wall_s']:.4f} | {row['wall_change_percent']:+.1f}% | {b['max_rss_kib']/1024:.2f} / {c['max_rss_kib']/1024:.2f} |\n"
text+='\nMeasured wall-time ranges across the seven samples:\n\n'
text+='| Workload | uBLAS (s) | Kokkos (s) |\n| --- | ---: | ---: |\n'
for row in summary:
    b,c=row['baseline'],row['candidate']
    text+=f"| {row['case']} | {b['wall_min_s']:.4f}–{b['wall_max_s']:.4f} | {c['wall_min_s']:.4f}–{c['wall_max_s']:.4f} |\n"
text+='\nAll samples, warmups, and output hashes are in\n[mdspan-storage-gcc14-samples.tsv](mdspan-storage-gcc14-samples.tsv).\n\nBinary SHA-256 values:\n\n'
for name in ('baseline','candidate'):
    text+=f"- {name}: `{metadata['binaries'][name]['sha256']}`\n"
p=repo/'doc/mdspan-storage.md'
s=p.read_text()
assert '## Bundled Kokkos performance results' not in s
p.write_text(s+text)
p=repo/'doc/Makefile.am'
s=p.read_text().replace('mdspan-storage-samples.tsv \\\n','mdspan-storage-samples.tsv \\\n\t\t\tmdspan-storage-gcc14-samples.tsv \\\n')
p.write_text(s)
print(f'Added 160 samples; median wall changes {min(changes):+.1f}% to {max(changes):+.1f}%')
