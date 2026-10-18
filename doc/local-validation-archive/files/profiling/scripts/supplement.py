#!/usr/bin/env python3
"""Run after benchmark.py: matched cache/BLAS and full-matrix memory checks."""
import json, subprocess, time, hashlib
from pathlib import Path
from benchmark import ROOT, BIN, ENV, INPUT, run, summarize

def main():
    cases=json.loads((ROOT/'cases.json').read_text())
    by_name={c['name']:c for c in cases}
    cache=ROOT/'inputs/cache';cache.mkdir(exist_ok=True)
    source=by_name['accessibility_computed']
    cmd=[str(BIN)]+source['args']+['--out=tAcc:'+str(cache/'target.ed'),'--out=qAcc:'+str(cache/'query.ed')]
    start=time.perf_counter()
    p=subprocess.run(cmd,env=ENV,capture_output=True)
    (ROOT/'logs/cache-generation.stdout').write_bytes(p.stdout)
    (ROOT/'logs/cache-generation.stderr').write_bytes(p.stderr)
    (ROOT/'logs/cache-generation.json').write_text(json.dumps(dict(command=cmd,wall_s=time.perf_counter()-start,exit_code=p.returncode,
        note='One-time ED generation includes a complete prediction and text serialization.'),indent=2)+'\n')
    p.check_returncode()
    new=[]
    def variant(base,name,extra=(),**kw):
        c=dict(by_name[base],name=name,args=by_name[base]['args']+list(extra),**kw)
        new.append(c)
    variant('accessibility_computed','cache_computed_control')
    variant('accessibility_computed','cache_reused',['--tAcc=E','--tAccFile='+str(cache/'target.ed'),
        '--qAcc=E','--qAccFile='+str(cache/'query.ed')])
    variant('bio_fhla_oxys','blas_one_control')
    variant('bio_fhla_oxys','blas_unset',env_overrides={'OPENBLAS_NUM_THREADS':None})
    variant('bio_fhla_oxys','blas_twelve',env_overrides={'OPENBLAS_NUM_THREADS':'12'})
    variant('memory_t3000_q800','memory_s',['--model=S'])
    variant('memory_window300','memory_s_window300',['--model=S'])
    cases+=new
    (ROOT/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    eligible=[c for c in new if run(c,'warmup',0,120)['exit_code']==0]
    import random
    rng=random.Random(7312027)
    for n in range(1,4):
        order=eligible[:];rng.shuffle(order)
        for c in order:run(c,'measure',n,120)
    summarize()
    (cache/'manifest.json').write_text(json.dumps({p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in cache.glob('*.ed')},indent=2)+'\n')

if __name__=='__main__':main()
