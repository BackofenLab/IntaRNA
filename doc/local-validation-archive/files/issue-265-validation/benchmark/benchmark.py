#!/usr/bin/env python3
"""Reproducible end-to-end accessibility benchmark; unsafe base runs are compared to base serial."""
import argparse, hashlib, json, os, pathlib, random, shutil, statistics, subprocess, time

p=argparse.ArgumentParser()
p.add_argument('--base',type=pathlib.Path,required=True)
p.add_argument('--candidate',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--repetitions',type=int,default=3)
p.add_argument('--cpus',default='8-11')
a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
for label, binary in [('base',a.base),('candidate',a.candidate)]:
    shutil.copy2(binary,a.output/('IntaRNA-'+label))
rng=random.Random(265)
queries=a.output/'queries.fa'
queries.write_text(''.join('>query%02d\n%s\n' % (i+1,''.join(rng.choice('ACGU') for j in range(400))) for i in range(32)))
target=a.output/'target.fa'
target.write_text('>target\nGCGCGC\n')
options=['--query='+str(queries.resolve()),'--target='+str(target.resolve()),'--model=X','--mode=K','--seedBP=3','--qIntLenMax=8','--tIntLenMax=6','--qAccW=120','--qAccL=120','--tAccW=6','--tAccL=6','--outMode=C','--outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,ED1,ED2']
env=dict(os.environ,OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES='cores',OPENBLAS_NUM_THREADS='1')

def hashes(directory):
    lines=(directory/'predictions.csv').read_text().splitlines()
    norm=lines[0]+'\n'+'\n'.join(sorted(lines[1:]))+'\n'
    result={'predictions.csv':hashlib.sha256(norm.encode()).hexdigest()}
    for f in sorted(directory.glob('query*.ed')):
        result[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
    if len(result)!=33:
        raise RuntimeError('Expected 32 query ED files and CSV: '+str(directory))
    return result

def run(label,threads,iteration,warmup=False):
    name=('%s-t%d-%s%d' % (label,threads,'warmup' if warmup else 'run',iteration))
    directory=a.output/name
    directory.mkdir()
    command=['/usr/bin/time','-f','%e %U %S %M','-o',str((directory/'time.txt').resolve()),'taskset','-c',a.cpus,str((a.output/('IntaRNA-'+label)).resolve()),*options,'--threads='+str(threads),'--out='+str((directory/'predictions.csv').resolve()),'--out=qAcc:'+str((directory/'query.ed').resolve())]
    with (directory/'stdout.txt').open('w') as stdout,(directory/'stderr.txt').open('w') as stderr:
        start=time.perf_counter()
        completed=subprocess.run(command,env=env,stdout=stdout,stderr=stderr)
        elapsed=time.perf_counter()-start
    if completed.returncode:
        raise RuntimeError('Failed '+name+': '+(directory/'stderr.txt').read_text())
    wall,user,system,rss=(directory/'time.txt').read_text().split()
    result={'name':name,'label':label,'threads':threads,'wall_s':float(wall),'monotonic_wall_s':elapsed,'user_s':float(user),'system_s':float(system),'peak_rss_kib':int(rss),'command':command,'hashes':hashes(directory)}
    print(name,wall,rss,flush=True)
    return result

warmups=[run(label,threads,0,True) for label,threads in [('base',1),('candidate',1),('base',4),('candidate',4)]]
reference=warmups[0]['hashes']
results=[]
for i in range(a.repetitions):
    configs=[('base',1),('candidate',1),('base',4),('candidate',4)]
    if i%2: configs.reverse()
    for label,threads in configs: results.append(run(label,threads,i+1))
for result in warmups+results:
    result['matches_base_serial']=result['hashes']==reference
    result['different_files']=[k for k,v in result['hashes'].items() if reference.get(k)!=v]
summary=[]
for label,threads in [('base',1),('candidate',1),('base',4),('candidate',4)]:
    selected=[r for r in results if r['label']==label and r['threads']==threads]
    summary.append({'label':label,'threads':threads,'median_wall_s':statistics.median(r['wall_s'] for r in selected),'min_wall_s':min(r['wall_s'] for r in selected),'max_wall_s':max(r['wall_s'] for r in selected),'median_peak_rss_kib':statistics.median(r['peak_rss_kib'] for r in selected),'matches_base_serial':all(r['matches_base_serial'] for r in selected)})
meta={'cpu_affinity':a.cpus,'environment':{key:env[key] for key in ['OMP_DYNAMIC','OMP_PROC_BIND','OMP_PLACES','OPENBLAS_NUM_THREADS']},'binary_hashes':{label:hashlib.sha256((a.output/('IntaRNA-'+label)).read_bytes()).hexdigest() for label in ['base','candidate']},'input_hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [queries,target]},'base_version':subprocess.check_output([str((a.output/'IntaRNA-base').resolve()),'--version'],text=True),'candidate_version':subprocess.check_output([str((a.output/'IntaRNA-candidate').resolve()),'--version'],text=True),'summary':summary,'warmups':warmups,'results':results}
(a.output/'results.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
if not all(r['matches_base_serial'] for r in warmups+results):
    raise SystemExit('Scientific output mismatch found; inspect results.json')
