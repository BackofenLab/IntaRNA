#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
from benchmark import ROOT,canonical
from profile import profile_output
dest=ROOT/'profiles'
rows=[]
for name in ['biological_batch','hybrid_noaccess','exact_x','ensemble_m']:
    a=profile_output((dest/(name+'_reference.stdout')).read_bytes())
    b=profile_output((dest/(name+'_cpu.stdout')).read_bytes())
    table=list(csv.reader(b.decode().splitlines(),delimiter=';'))
    valid=bool(table) and table[0]==['id1','id2','start1','end1','start2','end2','E'] and all(len(r)==7 for r in table)
    rows.append(dict(name=name,kind='cpu',schema_valid=valid,output_equal=a==b,
        rows=len(table)-1,sha256=hashlib.sha256(b).hexdigest()))
for name in ['bio_phob_gcvb','algorithm_exact_x','memory_t3000_q800','memory_s']:
    a=canonical((ROOT/'results'/(name+'.measure1.stdout')).read_bytes())
    b=profile_output((dest/(name+'_heap.stdout')).read_bytes())
    rows.append(dict(name=name,kind='heap',output_equal=a==b,sha256=hashlib.sha256(b).hexdigest()))
(dest/'validated-outputs.json').write_text(json.dumps(dict(note='Only recognized gprofng creation banner and IntaRNA # log lines excluded; original files/validation retained.',checks=rows),indent=2)+'\n')
print(json.dumps(rows,indent=2))
if not all(r['output_equal'] and r.get('schema_valid',True) for r in rows):raise SystemExit(1)
