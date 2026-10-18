#!/usr/bin/env python3
import hashlib,json,os,pathlib,random,statistics,subprocess,time
ROOT=pathlib.Path('/home/alex/CodexProjects/IntaRNA')
OUT=ROOT/'issue-258-validation/audit-performance';D=OUT/'fixed'
commands=json.loads((OUT/'commands.json').read_text())
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
versions={'before':OUT/'IntaRNA-before','fixed':ROOT/'pr-worktrees/issue-258-bpsvg/src/bin/IntaRNA'}
records=[];rng=random.Random(260);configs=[(v,m) for v in versions for m in ['none','csv']];hashes={}
for rep in range(-1,5):
    order=configs.copy();rng.shuffle(order)
    for version,mode in order:
        cmd=commands[f'dense40/head/{mode}'].copy();cmd[3]=str(versions[version])
        cmd=[arg.replace(str(OUT),str(D)) if arg.startswith('--out=') else arg for arg in cmd]
        start=time.perf_counter();p=subprocess.run(['/usr/bin/time','-f','%e,%U,%S,%M','-o',str(D/'time.txt')]+cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=ENV);duration=time.perf_counter()-start
        if p.returncode:raise RuntimeError(p.stderr.decode())
        elapsed,user,system,rss=map(float,(D/'time.txt').read_text().strip().split(','))
        if rep<0:
            hashes[(version,mode)]=(hashlib.sha256((D/'ensemble.txt').read_bytes()).hexdigest(),hashlib.sha256((D/'pairs.csv').read_bytes()).hexdigest() if mode=='csv' else None)
        else:
            records.append(dict(case='dense40',version=version,mode=mode,rep=rep,wall_seconds=duration,rss_kib=int(rss)))
for mode in ['none','csv']:assert hashes['before',mode]==hashes['fixed',mode]
summary=[]
for version,mode in configs:
    r=[x for x in records if (x['version'],x['mode'])==(version,mode)]
    x=dict(case='dense40',version=version,mode=mode,wall_median=statistics.median(x['wall_seconds'] for x in r),rss_max_kib=max(x['rss_kib'] for x in r));summary.append(x);print(x)
(D/'control-raw.json').write_text(json.dumps(records,indent=2)+'\n');(D/'control-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Dense controls matched byte-for-byte. DONE')
