#!/usr/bin/env python3
import collections, hashlib, json, os, pathlib, random, statistics, subprocess, time, xml.etree.ElementTree as ET
ROOT=pathlib.Path('/home/alex/CodexProjects/IntaRNA')
OUT=ROOT/'issue-258-validation/audit-performance'
BASE=ROOT/'pr-worktrees/issue-257-plan'
HEAD=ROOT/'pr-worktrees/issue-258-bpsvg'
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
CPU=min(os.sched_getaffinity(0)); NS={'s':'http://www.w3.org/2000/svg'}
common=['--model=P','--mode=M','--seedBP=4','--intLenMax=30','--qIntLoopMax=3','--tIntLoopMax=3','--threads=1','--outNoLP=false','--outNoGUend=false','--outNumber=0','--outMode=E','--default-log-file=/dev/null']
cases={}
cases['OxyS-fhlA']=[f'--target={HEAD}/doc/handson/fhlA.fasta',f'--query={HEAD}/doc/handson/OxyS.fasta']
for name,n in [('dense40',40),('sparse400',400),('sparse1000',1000)]:
    t='G'*n if name.startswith('dense') else 'A'*(n//2)+'G'*8+'A'*(n-n//2-8)
    q='C'*n if name.startswith('dense') else 'A'*(n//2)+'C'*8+'A'*(n-n//2-8)
    for label,seq in [('target',t),('query',q)]:
        (OUT/f'{name}-{label}.fa').write_text(f'>{label}\n{seq}\n')
    cases[name]=[f'--target={OUT}/{name}-target.fa',f'--query={OUT}/{name}-query.fa','--tAcc=N','--qAcc=N']
configs=[('base','none'),('head','none'),('base','csv'),('head','csv'),('head','svg'),('head','both')]
commands={}; records=[]; checks=[]; rng=random.Random(258)
def run(name,version,mode,rep):
    tree=BASE if version=='base' else HEAD
    cmd=['taskset','-c',str(CPU),str(tree/'src/bin/IntaRNA')]+common+cases[name]+[f'--out={OUT}/ensemble.txt']
    if mode in ('csv','both'):cmd += [f'--out=bpProb:{OUT}/pairs.csv']
    if mode in ('svg','both'):cmd += [f'--out=bpsvg:{OUT}/plot.svg']
    commands[f'{name}/{version}/{mode}']=cmd
    begin=time.perf_counter()
    proc=subprocess.run(['/usr/bin/time','-f','%e,%U,%S,%M','-o',str(OUT/'time.txt')]+cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=ENV)
    duration=time.perf_counter()-begin
    if proc.returncode:raise RuntimeError(f'{name}/{version}/{mode}: {proc.stderr.decode()}')
    elapsed,user,system,rss=map(float,(OUT/'time.txt').read_text().strip().split(','))
    data=dict(case=name,version=version,mode=mode,rep=rep,wall_seconds=duration,time_seconds=elapsed,user_seconds=user,system_seconds=system,rss_kib=int(rss),csv_bytes=(OUT/'pairs.csv').stat().st_size if mode in ('csv','both') else 0,svg_bytes=(OUT/'plot.svg').stat().st_size if mode in ('svg','both') else 0)
    if rep>=0:
        records.append(data)
        (OUT/'raw.json').write_text(json.dumps(records,indent=2)+'\n')
    return (OUT/'ensemble.txt').read_bytes(),(OUT/'pairs.csv').read_bytes() if mode in ('csv','both') else None
for name in cases:
    print('START',name,flush=True)
    baseline_ens=None;baseline_csv=None
    # Warm each configuration and independently compare full old/new outputs.
    for version,mode in configs:
        ensemble,csv=run(name,version,mode,-1)
        if baseline_ens is None:baseline_ens=ensemble
        assert ensemble==baseline_ens,(name,version,mode,'ensemble differs')
        if csv is not None:
            if baseline_csv is None:baseline_csv=csv
            assert csv==baseline_csv,(name,version,mode,'CSV differs')
        if mode=='both':
            # Stream parsing bounds verifier memory separately from timing.
            rows=[row.split(';')[1:] for row in csv.decode().strip().splitlines()[1:]]
            count=seed_count=0
            for event,el in ET.iterparse(OUT/'plot.svg',events=['end']):
                if el.get('data-type')=='base-pair':
                    i=count//len(rows[0]);j=count%len(rows[0]);assert float(el.get('data-probability'))==float(rows[i][j]);count+=1
                if el.tag.endswith('rect') and el.get('class')=='seed' and el.get('data-target') is not None:seed_count+=1
                el.clear()
            assert count==len(rows)*len(rows[0])
            checks.append(dict(case=name,cells=count,seed_cells=seed_count,ensemble=ensemble.decode(),csv_sha256=hashlib.sha256(csv).hexdigest(),csv_bytes=len(csv),svg_bytes=(OUT/'plot.svg').stat().st_size))
    # Interleaved randomized order limits drift between base and candidate.
    for rep in range(7):
        order=configs.copy();rng.shuffle(order)
        for version,mode in order:run(name,version,mode,rep)
    for version,mode in configs:
        rs=[r for r in records if (r['case'],r['version'],r['mode'])==(name,version,mode)]
        print(name,version,mode,'median seconds',round(statistics.median(r['wall_seconds'] for r in rs),4),'max RSS KiB',max(r['rss_kib'] for r in rs),flush=True)
    (OUT/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    (OUT/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
summary=[]
for name in cases:
    for version,mode in configs:
        rs=[r for r in records if (r['case'],r['version'],r['mode'])==(name,version,mode)]
        summary.append(dict(case=name,version=version,mode=mode,wall_median=statistics.median(r['wall_seconds'] for r in rs),wall_min=min(r['wall_seconds'] for r in rs),wall_max=max(r['wall_seconds'] for r in rs),rss_max_kib=max(r['rss_kib'] for r in rs),rss_median_kib=statistics.median(r['rss_kib'] for r in rs),csv_bytes=rs[0]['csv_bytes'],svg_bytes=rs[0]['svg_bytes']))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('DONE',flush=True)
