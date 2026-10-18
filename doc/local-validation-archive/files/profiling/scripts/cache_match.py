#!/usr/bin/env python3
"""Match effective caps after discovering the ED reload's one-column reduction."""
import json,random
from benchmark import ROOT,run,summarize

def main():
    cases=json.loads((ROOT/'cases.json').read_text())
    by_name={c['name']:c for c in cases}
    new=[]
    for base,name in [('cache_computed_control','cache_matched_computed'),('cache_reused','cache_matched_reused')]:
        new.append(dict(by_name[base],name=name,args=[a for a in by_name[base]['args'] if not a.startswith('--intLenMax=')]
                        +['--qIntLenMax=99','--tIntLenMax=149']))
    (ROOT/'cases.json').write_text(json.dumps(cases+new,indent=2)+'\n')
    eligible=[c for c in new if run(c,'warmup',0,90)['exit_code']==0]
    rng=random.Random(7312028)
    for rep in range(1,4):
        order=eligible[:];rng.shuffle(order)
        for c in order:run(c,'measure',rep,90)
    summarize()

if __name__=='__main__':main()
