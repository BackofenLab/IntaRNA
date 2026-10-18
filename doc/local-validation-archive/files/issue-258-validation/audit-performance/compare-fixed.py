#!/usr/bin/env python3
import hashlib,json,os,pathlib,random,statistics,subprocess,time
ROOT=pathlib.Path('/home/alex/CodexProjects/IntaRNA')
OUT=ROOT/'issue-258-validation/audit-performance'
FIXED=OUT/'fixed';FIXED.mkdir(exist_ok=True)
original_commands=json.loads((OUT/'commands.json').read_text())
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
versions={'before':OUT/'IntaRNA-before','fixed':ROOT/'pr-worktrees/issue-258-bpsvg/src/bin/IntaRNA'}
records=[];checks=[];commands={};rng=random.Random(259)
def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        while data:=stream.read(1048576):h.update(data)
    return h.hexdigest()
def run(case,version,mode,rep):
    cmd=original_commands[f'{case}/head/{mode}'].copy();cmd[3]=str(versions[version])
    cmd=[arg.replace(str(OUT),str(FIXED)) if arg.startswith('--out=') else arg for arg in cmd]
    commands[f'{case}/{version}/{mode}']=cmd
    start=time.perf_counter()
    p=subprocess.run(['/usr/bin/time','-f','%e,%U,%S,%M','-o',str(FIXED/'time.txt')]+cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=ENV)
    duration=time.perf_counter()-start
    if p.returncode:raise RuntimeError(p.stderr.decode())
    elapsed,user,system,rss=map(float,(FIXED/'time.txt').read_text().strip().split(','))
    if rep>=0:
        records.append(dict(case=case,version=version,mode=mode,rep=rep,wall_seconds=duration,time_seconds=elapsed,user_seconds=user,system_seconds=system,rss_kib=int(rss),svg_bytes=(FIXED/'plot.svg').stat().st_size))
        (FIXED/'raw.json').write_text(json.dumps(records,indent=2)+'\n')
    return dict(ensemble=digest(FIXED/'ensemble.txt'),svg=digest(FIXED/'plot.svg'),csv=digest(FIXED/'pairs.csv') if mode=='both' else None)
cases=['OxyS-fhlA','sparse400','sparse1000'];configs=[(version,mode) for version in versions for mode in ['svg','both']]
for case in cases:
    print('START',case,flush=True)
    warm={}
    for version,mode in configs:
        warm[(version,mode)]=run(case,version,mode,-1)
    for mode in ['svg','both']:
        assert warm['before',mode]==warm['fixed',mode],(case,mode,warm)
    checks.append({'case':case,'identical_before_fixed':True,'hashes':warm['fixed','both']})
    for rep in range(5):
        order=configs.copy();rng.shuffle(order)
        for version,mode in order:run(case,version,mode,rep)
    for version,mode in configs:
        rs=[r for r in records if (r['case'],r['version'],r['mode'])==(case,version,mode)]
        print(case,version,mode,'median seconds',round(statistics.median(r['wall_seconds'] for r in rs),4),'max RSS KiB',max(r['rss_kib'] for r in rs),flush=True)
summary=[]
for case in cases:
    for version,mode in configs:
        rs=[r for r in records if (r['case'],r['version'],r['mode'])==(case,version,mode)]
        summary.append(dict(case=case,version=version,mode=mode,wall_median=statistics.median(r['wall_seconds'] for r in rs),wall_min=min(r['wall_seconds'] for r in rs),wall_max=max(r['wall_seconds'] for r in rs),rss_max_kib=max(r['rss_kib'] for r in rs),rss_median_kib=statistics.median(r['rss_kib'] for r in rs),svg_bytes=rs[0]['svg_bytes']))
(FIXED/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(FIXED/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
(FIXED/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(FIXED/'binaries.json').write_text(json.dumps({v:dict(path=str(p),sha256=digest(p)) for v,p in versions.items()},indent=2)+'\n')
print('DONE',flush=True)
