#!/usr/bin/env python3
"""Untimed correctness checks including exact base-pair lists."""
import json,subprocess,hashlib
from benchmark import ROOT,BIN,ENV,canonical
from cache_match import main as run_matched_cache_cases
# This stage follows supplement.py and precedes profiling; add the newly discovered
# effective-cap control while clean timing is still isolated from instrumentation.
run_matched_cache_cases()
cases={c['name']:c for c in json.loads((ROOT/'cases.json').read_text())}
dest=ROOT/'results/structural';dest.mkdir(exist_ok=True)
rows=[]
for name in ['throughput_1t','throughput_12t','cache_computed_control','cache_reused','cache_matched_computed','cache_matched_reused']:
    args=[('--outCsvCols=id1,id2,start1,end1,start2,end2,E,bpList') if a.startswith('--outCsvCols=') else a for a in cases[name]['args']]
    cmd=[str(BIN)]+args
    p=subprocess.run(cmd,env=ENV,capture_output=True,timeout=120)
    (dest/(name+'.stdout')).write_bytes(p.stdout);(dest/(name+'.stderr')).write_bytes(p.stderr)
    rows.append(dict(case=name,command=cmd,exit_code=p.returncode,canonical_sha256=hashlib.sha256(canonical(p.stdout)).hexdigest()))
lookup={r['case']:r for r in rows}
comparisons={label:lookup[a]['exit_code']==lookup[b]['exit_code']==0 and lookup[a]['canonical_sha256']==lookup[b]['canonical_sha256']
    for label,a,b in [('threads_1_vs_12','throughput_1t','throughput_12t'),('cached_vs_computed','cache_computed_control','cache_reused'),
                     ('cache_matched_caps','cache_matched_computed','cache_matched_reused')]}
(dest/'checks.json').write_text(json.dumps(dict(runs=rows,comparisons=comparisons),indent=2)+'\n')
print(json.dumps(comparisons))
