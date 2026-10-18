import subprocess,random,json,csv,io
b='/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-212/src/bin/IntaRNA'
rng=random.Random(212)
for length in range(8,21):
 for n in range(100):
  t=''.join(rng.choices('ACGU',k=length));q='GG'
  cmd=[b,'-t',t,'-q',q,'--energy=B','--acc=C','--seedBP=2','--outOverlap=N','-n','10','--outMode=C','--outCsvCols=start1,end1,start2,end2,E','--tRegionLenMax=4','--outDeltaE=100','--threads=1','--outMaxE=100','-v']
  p=subprocess.run(cmd,capture_output=True,text=True)
  rows=list(csv.DictReader(io.StringIO('\n'.join(x for x in p.stdout.splitlines() if not x.startswith('#'))),delimiter=';'))
  if len(rows)>1:
   print(json.dumps({'t':t,'q':q,'rows':rows,'stderr':p.stderr,'cmd':cmd}),flush=True);raise SystemExit
