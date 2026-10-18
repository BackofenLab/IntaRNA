import subprocess,csv,io,random,json
binary='/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-212/src/bin/IntaRNA'
base=[binary,'--energy=B','--acc=N','--model=S','--mode=H','--noSeed','--outMode=C','--outCsvCols=start1,end1,start2,end2,E,hybridDB','--outNumber=10','--outDeltaE=100','--threads=1','--outNoLP=false','--outNoGUend']
rng=random.Random(212)
for length in range(3,10):
    found=False
    for n in range(300):
        t=''.join(rng.choices('ACGU',k=length));q=''.join(rng.choices('ACGU',k=length))
        for ov in 'NTQ':
            cmd=base+['-t',t,'-q',q,'--outOverlap='+ov]
            p=subprocess.run(cmd,capture_output=True,text=True)
            if p.returncode: raise Exception(p.stderr)
            rows=list(csv.DictReader(io.StringIO(p.stdout),delimiter=';'))
            invalid=[r for r in rows if any(t[int(r[a])-1]+q[int(r[b])-1] in ('GU','UG') for a,b in [('start1','end2'),('end1','start2')])]
            if invalid:
                print(json.dumps({'t':t,'q':q,'ov':ov,'rows':rows,'invalid':invalid,'cmd':cmd}),flush=True)
                found=True;break
        if found:break
    if found:break
