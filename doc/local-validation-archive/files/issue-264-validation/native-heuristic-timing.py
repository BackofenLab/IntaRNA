from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import time

root=Path(__file__).resolve().parent
repo=root.parent/'pr-worktrees'/'issue-264'
old=root.parent/'pr-worktrees'/'issue-257-plan'
directory=root/'native-heuristic-timing'
directory.mkdir(exist_ok=True)
for name,source in [('base',old),('candidate',repo)]:
    shutil.copy2(source/'src/bin/IntaRNA',directory/name)
common=['--model=P','--mode=H','--seedBP=4','--intLenMax=30','--qIntLoopMax=3','--tIntLoopMax=3',
        '--threads=1','--outNoLP=false','--outNoGUend=false','--outNumber=0','--outMode=E',
        '--default-log-file=/dev/null','--out=STDOUT']
cases={'OxyS-fhlA':['--target='+str(repo/'doc/handson/fhlA.fasta'),'--query='+str(repo/'doc/handson/OxyS.fasta')],
       'dense40':['--target='+'G'*40,'--query='+'C'*40,'--tAcc=N','--qAcc=N']}
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
records=[]
for case,arguments in cases.items():
    for repeat in range(3):
        for variant in (('base','candidate') if repeat%2==0 else ('candidate','base')):
            command=[str(directory/variant),*arguments,*common]
            rss=directory/'rss.txt'
            start=time.perf_counter()
            result=subprocess.run(['/usr/bin/time','-f','%M','-o',str(rss),*command],env=env,text=True,capture_output=True,timeout=60)
            elapsed=time.perf_counter()-start
            result.check_returncode()
            records.append({'case':case,'variant':variant,'repeat':repeat,'seconds':elapsed,
                            'rss_kib':int(rss.read_text()),'output':result.stdout.strip()})
            print(records[-1],flush=True)
manifest={'baseline':'7a522358274f8ae636634bd89e51e2d8d7ee86b6','notes':'Exploratory native timing; older baseline, intentionally different heuristic ensembles, no pair output, build/test activity may overlap.',
          'commands':{case:args+common for case,args in cases.items()},'records':records,
          'binary_sha256':{variant:hashlib.sha256((directory/variant).read_bytes()).hexdigest() for variant in ('base','candidate')}}
(directory/'results.json').write_text(json.dumps(manifest,indent=2)+'\n')
