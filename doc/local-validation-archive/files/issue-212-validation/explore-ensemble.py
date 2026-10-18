import subprocess,csv,io,random,json
b='/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-212/src/bin/IntaRNA'
base=[b,'--energy=B','--acc=C','--model=P','--mode=H','--noSeed','--outMode=C','--outCsvCols=start1,end1,start2,end2,E,ED1,ED2','--outNumber=1000','--outDeltaE=100','--outMaxE=100','--threads=1','--outNoLP=false']
def run(t,q,ov):
 p=subprocess.run(base+['-t',t,'-q',q,'--outOverlap='+ov],capture_output=True,text=True)
 if p.returncode: raise Exception(p.stderr)
 return list(csv.DictReader(io.StringIO(p.stdout),delimiter=';'))
rng=random.Random(213)
found={}
for length in range(3,9):
 for n in range(200):
  t=''.join(rng.choices('ACGU',k=length));q=''.join(rng.choices('ACGU',k=length)); bb=run(t,q,'B')
  match={tuple(r[x] for x in ('start1','end1','start2','end2')):r for r in bb}
  for ov in 'NTQ':
   rows=run(t,q,ov)
   mismatches=[(r,match.get(tuple(r[x] for x in ('start1','end1','start2','end2')))) for r in rows if tuple(r[x] for x in ('start1','end1','start2','end2')) in match and float(r['E'])>-1000 and r['E']!=match[tuple(r[x] for x in ('start1','end1','start2','end2'))]['E']]
   if mismatches:
    print(json.dumps({'t':t,'q':q,'ov':ov,'rows':rows,'mismatch':mismatches}),flush=True);raise SystemExit
