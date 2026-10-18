import json,os,pathlib,random,statistics,subprocess,time
R=pathlib.Path('/home/alex/CodexProjects/IntaRNA');P=R/'issue-258-validation/audit-performance';D=P/'serializer-control';D.mkdir(exist_ok=True)
versions={'before':P/'IntaRNA-before','fixed':R/'pr-worktrees/issue-258-bpsvg/src/bin/IntaRNA'}
cmd0=json.loads((P/'commands.json').read_text())['sparse400/head/svg']
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE');rng=random.Random(261);records=[];commands={}
for rep in range(-1,7):
 order=list(versions);rng.shuffle(order)
 for version in order:
  cmd=cmd0.copy();cmd[3]=str(versions[version]);cmd=[('--out=bpsvg:/dev/null' if x.startswith('--out=bpsvg:') else '--out=STDOUT' if x.startswith('--out=') else x) for x in cmd];commands[version]=cmd
  start=time.perf_counter();r=subprocess.run(['/usr/bin/time','-f','%e,%U,%S,%M','-o',str(D/'time.txt')]+cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,env=env);wall=time.perf_counter()-start
  if r.returncode:raise RuntimeError(r.stderr.decode())
  elapsed,user,system,rss=map(float,(D/'time.txt').read_text().strip().split(','))
  if rep>=0:records.append(dict(version=version,rep=rep,wall=wall,user=user,system=system,cpu=user+system,rss=int(rss)))
summary=[]
for version in versions:
 a=[x for x in records if x['version']==version]
 summary.append(dict(version=version,median=statistics.median(x['wall'] for x in a),minimum=min(x['wall'] for x in a),maximum=max(x['wall'] for x in a),cpu_median=statistics.median(x['cpu'] for x in a),cpu_min=min(x['cpu'] for x in a),cpu_max=max(x['cpu'] for x in a),user_median=statistics.median(x['user'] for x in a),system_median=statistics.median(x['system'] for x in a),rss_max=max(x['rss'] for x in a)))
for name,data in [('raw',records),('summary',summary),('commands',commands)]: (D/f'{name}.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(summary,indent=2))
