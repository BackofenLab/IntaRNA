#!/usr/bin/env python3
"""Parse gprofng text exports without inventing missing symbols or call ancestry."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'profiles'
rows=[]
for name in ['gcvb','gcvb_repeat','chix384','long_pair']:
    text=(DEST/(name+'_self.stdout')).read_text();stats=(DEST/(name+'_statistics.stdout')).read_text();functions=[]
    for line in text.splitlines():
        match=re.match(r'\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)',line)
        if match:
            a,b,c,d,fun=match.groups();functions.append({'self_s':float(a),'self_pct':float(b),'inclusive_s':float(c),'inclusive_pct':float(d),'function':fun})
    total=next(f for f in functions if f['function']=='<Total>')['self_s']
    user=float(re.search(r'User CPU:\s+([\d.]+)',stats).group(1));system=float(re.search(r'System CPU:\s+([\d.]+)',stats).group(1))
    rows.append({'profile':name,'nominal_sampling_period_ms':100,'sampled_cpu_s':total,'samples':round(total/.1),'process_user_s':user,'process_system_s':system,'sampled_minus_process_cpu_s':total-user-system,'collector_warnings':[line for line in stats.splitlines() if 'Warning:' in line],'self_functions':functions})
(DEST/'cpu-summary.json').write_text(json.dumps(rows,indent=2)+'\n')
heaps=[]
for name in ['batch_gcvb_1t','long_pair']:
    text=(DEST/(name+'_heap_summary.stdout')).read_text()
    totals=re.findall(r'^Total bytes\s+(\d+)',text,re.M)
    heaps.append({'profile':name,'peak_tracked_heap_bytes':int(re.search(r'Heap size bytes\s+(\d+)',text).group(1)),'allocation_calls':int(re.search(r'Total allocations\s+(\d+)',text).group(1)),'cumulative_allocated_bytes':int(totals[0]),'outstanding_at_exit_bytes':int(totals[1]),'note':'Tracked heap metrics differ from RSS. Outstanding at exit does not by itself establish an application leak.'})
(DEST/'heap-summary.json').write_text(json.dumps(heaps,indent=2)+'\n')
print(json.dumps({'cpu':[{'profile':r['profile'],'samples':r['samples'],'sampled_minus_process_cpu_s':r['sampled_minus_process_cpu_s']} for r in rows],'heap':heaps},indent=2))
