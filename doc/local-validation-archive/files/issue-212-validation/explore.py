import subprocess, csv, io, random, json
from pathlib import Path
binary='/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-184/src/bin/IntaRNA'
base=[binary,'--energy=B','--acc=N','--model=S','--mode=M','--noSeed','--outMode=C','--outCsvCols=start1,end1,start2,end2,E,hybridDB','--outNumber=1000','--outDeltaE=100','--threads=1','--outNoLP=false']
def run(t,q,ov,extra=()):
    p=subprocess.run(base+['-t',t,'-q',q,'--outOverlap='+ov]+list(extra),capture_output=True,text=True)
    if p.returncode: raise Exception(p.stderr)
    return list(csv.DictReader(io.StringIO(p.stdout),delimiter=';'))
def overlap(a,b,axis):
    return max(int(a['start'+axis]),int(b['start'+axis]))<=min(int(a['end'+axis]),int(b['end'+axis]))
def allowed(a,b,ov): return all(not overlap(a,b,axis) for axis in {'N':'12','T':'2','Q':'1','B':''}[ov])
rng=random.Random(212)
found={}
for length in range(3,10):
    for n in range(200):
        t=''.join(rng.choices('ACGU',k=length)); q=''.join(rng.choices('ACGU',k=length))
        b=run(t,q,'B')
        if len(b)<2: continue
        for ov in 'NTQ':
            if ov in found: continue
            a=run(t,q,ov)
            missing=[x for x in b if all(allowed(x,y,ov) for y in a)]
            if missing:
                found[ov]={'t':t,'q':q,'actual':a,'missing':missing[:5]}
                print(ov,json.dumps(found[ov]),flush=True)
        if len(found)==3: break
    if len(found)==3: break
Path('/home/alex/CodexProjects/IntaRNA/issue-212-validation/explore-results.json').write_text(json.dumps(found,indent=2))
