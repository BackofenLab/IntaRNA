#!/usr/bin/env python3
"""Post-profile diagnostic repeats for the observed ChiX energy instability."""
import csv,hashlib,io,json,random,subprocess,time
from pathlib import Path
from benchmark import ROOT,BIN,ENV,dump,sha
DEST=ROOT/'diagnostics/parallel';DEST.mkdir(exist_ok=True)
cases={x['name']:x for x in json.loads((ROOT/'cases.json').read_text())}
def rows(data):
    text='\n'.join(x for x in data.decode().splitlines() if x and not x.startswith('#'))
    return {(r['id1'],r['id2']):r for r in csv.DictReader(io.StringIO(text),delimiter=';')}
reference=rows((ROOT/'results/batch_chix_1t.measure1.stdout').read_bytes())
schedule=[(t,n) for t,count in [(1,3),(6,5),(12,20)] for n in range(1,count+1)]
random.Random(202609212).shuffle(schedule);checks=[]
for threads,n in schedule:
    name=f't{threads}_repeat{n:02d}';command=[str(BIN)]+cases[f'batch_chix_{threads}t']['args']
    if (DEST/(name+'.stdout')).exists():raise RuntimeError('Refusing to overwrite '+name)
    p=subprocess.run(command,env=ENV,capture_output=True,timeout=60)
    (DEST/(name+'.stdout')).write_bytes(p.stdout);(DEST/(name+'.stderr')).write_bytes(p.stderr)
    actual=rows(p.stdout);differences=[]
    for key in reference.keys()|actual.keys():
        a,b=reference.get(key),actual.get(key)
        if a!=b:differences.append({'id':key,'fields':{k:{'reference':a.get(k),'observed':b.get(k)} for k in a.keys()|b.keys() if a.get(k)!=b.get(k)}} if a and b else {'id':key,'reference':a,'observed':b})
    row={'name':name,'threads':threads,'command':command,'exit_code':p.returncode,'binary_sha256':sha(BIN),'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'output_equal':not differences and p.returncode==0,'differences':differences}
    checks.append(row)
    with (DEST/'checks.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(name,'equal' if row['output_equal'] else json.dumps(differences),flush=True)
summary={'scope':'Untimed post-profile repeat check, not additional performance measurements','reference':'results/batch_chix_1t.measure1.stdout','random_seed':202609212,'counts':{str(t):{'runs':sum(r['threads']==t for r in checks),'different':sum(r['threads']==t and not r['output_equal'] for r in checks)} for t in [1,6,12]},'observed_differences':[r for r in checks if not r['output_equal']]}
dump(DEST/'summary.json',summary);print(json.dumps(summary,indent=2))
