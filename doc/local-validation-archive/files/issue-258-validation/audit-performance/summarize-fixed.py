import json,pathlib,statistics
P=pathlib.Path('/home/alex/CodexProjects/IntaRNA/issue-258-validation/audit-performance')
r=json.loads((P/'fixed/raw.json').read_text())+json.loads((P/'fixed-long/raw.json').read_text())
summary=[]
for case in ['OxyS-fhlA','sparse400','sparse1000']:
 for version in ['before','fixed']:
  for mode in ['svg','both']:
   a=[x for x in r if (x['case'],x['version'],x['mode'])==(case,version,mode)]
   assert len(a)==5,(case,version,mode,len(a))
   summary.append(dict(case=case,version=version,mode=mode,wall_median=statistics.median(x['wall_seconds'] for x in a),wall_min=min(x['wall_seconds'] for x in a),wall_max=max(x['wall_seconds'] for x in a),cpu_median=statistics.median(x['user_seconds']+x['system_seconds'] for x in a),rss_max_kib=max(x['rss_kib'] for x in a),svg_bytes=a[0]['svg_bytes']))
(P/'summary-fixed-all.json').write_text(json.dumps(summary,indent=2)+'\n')
controls=json.loads((P/'fixed/control-summary.json').read_text())
lines=['Performance audit after buffer-transfer and accessibility fixes (2026-10-09)', '', 'Direct comparison: retained fa2f84d binary vs final rebuilt candidate. Same environment/flags/input options as initial audit; five fresh-process timed repeats per configuration, deterministic shuffled/interleaved order, CPU0, one thread.', 'One warmup/check per configuration. Whole SVG/CSV/ensemble output matched byte-for-byte on native-energy cases.', 'The initial follow-up process received an unexplained SIGTERM after completing both shorter cases; sparse1000 was resumed in a separate process. Only complete five-run groups are reported; partial long-case warmups were discarded.', '', 'case | output | before/fixed wall median sec | before/fixed CPU median sec | before/fixed max RSS KiB | RSS reduction']
for case in ['OxyS-fhlA','sparse400','sparse1000']:
 for mode in ['svg','both']:
  before=next(x for x in summary if (x['case'],x['version'],x['mode'])==(case,'before',mode));fixed=next(x for x in summary if (x['case'],x['version'],x['mode'])==(case,'fixed',mode))
  lines.append(f"{case} | {mode} | {before['wall_median']:.4f}/{fixed['wall_median']:.4f} | {before['cpu_median']:.2f}/{fixed['cpu_median']:.2f} | {before['rss_max_kib']}/{fixed['rss_max_kib']} | {100*(1-fixed['rss_max_kib']/before['rss_max_kib']):.1f}%")
lines+=['', 'Dense40 control after fixes (same before/fixed binary pair):']
for mode in ['none','csv']:
 b=next(x for x in controls if (x['version'],x['mode'])==('before',mode));f=next(x for x in controls if (x['version'],x['mode'])==('fixed',mode));lines.append(f"{mode}: {b['wall_median']:.4f}s -> {f['wall_median']:.4f}s ({100*(f['wall_median']/b['wall_median']-1):+.2f}%), peak RSS {b['rss_max_kib']} -> {f['rss_max_kib']}KiB; ensemble/CSV byte-identical.")
lines+=['', 'Conclusion: transferring output buffers eliminates an avoidable full XML/CSV copy. On 1000x1000 output, SVG peak RSS drops by approximately 36% while preserving all bytes. No material existing-mode regression found. File-write elapsed times fluctuate on this nearly full filesystem, so avoid a speedup claim; process CPU and RSS support the improvement.', 'Full matrix SVG cost remains proportional to all target-by-query cells: 1000x1000 emits 220,685,460 bytes. This is expected from per-cell SVG elements and hover metadata.', '', 'Reproduce: compare-fixed.py (all SVG comparisons); compare-fixed-long.py (resumed sparse1000); control-fixed.py (dense control). Raw timed results: fixed/raw.json, fixed-long/raw.json, fixed/control-raw.json. Exact long-case commands and binary hashes: fixed-long/commands.json and fixed-long/binaries.json. Initial base-vs-feature audit: report-before.txt.', 'No tracked source files edited by the performance agent; root applied the verified fix.']
(P/'report-fixed.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
