#!/usr/bin/env python3
"""Default-mode study: fixed scientific defaults, sequential randomized timing."""
import argparse,csv,hashlib,io,json,os,platform,random,signal,statistics,subprocess,time,threading
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'source/src/bin/IntaRNA'
INPUT=ROOT/'inputs';RESULT=ROOT/'results'
REVISION='0a4568ad6e3da52f935221ff0da98a4a171e9833'
CONTROL={'LC_ALL':'C','OPENBLAS_NUM_THREADS':'1','OMP_DYNAMIC':'FALSE','OMP_PROC_BIND':'close','OMP_PLACES':'cores'}
ENV=dict(os.environ,**CONTROL)
DEFAULTS=['--model=X','--mode=H','--acc=C','--accW=150','--accL=100','--windowWidth=0','--seedBP=7','--seedMaxUP=0','--intLenMax=0']
HEADER=['id1','id2','start1','end1','start2','end2','E','bpList']
TIME_KEYS=['wall_s','user_s','system_s','max_rss_kib','major_faults','minor_faults','voluntary_context_switches','involuntary_context_switches','fs_inputs','fs_outputs','time_exit_code']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
def fasta(p):
    records=[]
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):records.append([line[1:],''])
        elif line.strip():records[-1][1]+=line.strip().upper().replace('T','U')
    return records

def write_fasta(name,records):
    p=INPUT/(name+'.fa');p.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in records));return p

def canonical(data):return b'\n'.join(sorted(line for line in data.splitlines() if line and not line.startswith(b'#') and not line.startswith(b'Creating experiment directory ')))+b'\n'
def validate(data):
    lines=[l for l in data.decode().splitlines() if l and not l.startswith('#') and not l.startswith('Creating experiment directory ')]
    rows=list(csv.reader(lines,delimiter=';'))
    return bool(rows) and rows[0]==HEADER and all(len(r)==len(HEADER) for r in rows),max(0,len(rows)-1)

def structure_hash(data):
    text='\n'.join(x for x in data.decode().splitlines() if x and not x.startswith('#'))
    rows=csv.DictReader(io.StringIO(text),delimiter=';')
    values=sorted(tuple(r[k] for k in HEADER if k!='E') for r in rows)
    return hashlib.sha256(json.dumps(values).encode()).hexdigest()

def prepare():
    INPUT.mkdir(exist_ok=True);RESULT.mkdir(exist_ok=True)
    docs=ROOT/'source/doc/handson'
    for name in ['fhlA','OxyS','phoB','GcvB','GcvB.ST','ilvE','ChiX_NC_000913']:
        write_fasta(name,fasta(docs/(name+'.fasta')))
    records=fasta(docs/'NC_000913.fa');indices=[i*(len(records)-1)//47 for i in range(48)]
    write_fasta('promoters48_systematic',[records[i] for i in indices])
    dump(INPUT/'selection.json',{'source':'doc/handson/NC_000913.fa','source_sha256':sha(docs/'NC_000913.fa'),'total_records':len(records),'indices_zero_based':indices,'rule':'48 evenly spaced record indices including endpoints; no outcome-based selection','note':'Bounded sample, not the full target collection; differs from the first-48 sample used on 2026-09-19.'})
    rng=random.Random(20260921)
    target=''.join(rng.choices('ACGU',k=10000));query=''.join(rng.choices('ACGU',k=800))
    for n in [100,300,1000,3000,10000]:write_fasta(f't{n}',[(f'synthetic_target_{n}',target[:n])])
    for n in [100,200,400,800]:write_fasta(f'q{n}',[(f'synthetic_query_{n}',query[:n])])
    cases=[]
    def add(name,category,t,q,threads=1):
        args=['--target='+str(INPUT/(t+'.fa')),'--query='+str(INPUT/(q+'.fa'))]+DEFAULTS+['--threads='+str(threads),'--outMode=C','--outCsvCols='+','.join(HEADER)]
        cases.append({'name':name,'category':category,'target':t,'query':q,'threads':threads,'predictions':len(fasta(INPUT/(t+'.fa')))*len(fasta(INPUT/(q+'.fa'))),'args':args})
    for t,q in [('fhlA','OxyS'),('phoB','GcvB'),('ilvE','GcvB.ST')]:add('bio_'+t.lower(),'single_biological',t,q)
    for n in [100,300,1000,3000,10000]:add(f'length_t{n}_q100','target_length',f't{n}','q100')
    for n in [200,400,800]:add(f'length_t300_q{n}','query_length','t300',f'q{n}')
    add('long_pair','long_pair','t3000','q800')
    for query,name in [('GcvB','gcvb'),('ChiX_NC_000913','chix')]:
        for threads in [1,6,12]:add(f'batch_{name}_{threads}t','biological_batch','promoters48_systematic',query,threads)
    dump(ROOT/'cases.json',cases)
    dump(INPUT/'manifest.json',{p.name:{'sha256':sha(p),'records':[{'id':n,'length':len(s),'gc_fraction':(s.count('G')+s.count('C'))/len(s)} for n,s in fasta(p)]} for p in sorted(INPUT.glob('*.fa'))})
    dump(ROOT/'study-plan.json',{'revision':REVISION,'cases':len(cases),'measured_repetitions':5,'excluded_warmups_per_case':1,'random_seed':20260921,'scientific_defaults':DEFAULTS,'controlled_environment':CONTROL,'output':'CSV including full bpList, retaining traceback; output formatting differs from CLI default','profiles':['batch_gcvb_1t','batch_chix_1t repeated 8x as one batch','long_pair'],'limitations':['One host; bounded biological sample','No comparison with historical revision','No algorithm/window/accessibility-off sweep','Diagnostics and sampling run after clean timings']})
    print('Prepared',len(cases),'cases',flush=True)

def capture(cmd):
    p=subprocess.run(cmd,capture_output=True,text=True,env=ENV);return {'returncode':p.returncode,'output':p.stdout+p.stderr}

def environment():
    deps=Path(os.environ.get('INTARNA_DEPS','/home/alex/CursorProjects/IntaRNA_optimization/intaRNA_legacy/.conda-env'))
    limits={}
    cg=Path('/sys/fs/cgroup'+Path('/proc/self/cgroup').read_text().split('::',1)[1].strip())
    for parent in [cg,*cg.parents]:
        if not parent.is_relative_to('/sys/fs/cgroup'):continue
        for n in ['cpu.max','memory.max','cpuset.cpus.effective']:
            p=parent/n
            if p.exists():limits[str(p)]=p.read_text().strip()
    dump(ROOT/'environment.json',{'recorded_utc':datetime.now(timezone.utc).isoformat(),'revision':REVISION,'binary_sha256':sha(BIN),'platform':platform.platform(),'version':capture([str(BIN),'--version']),'cpu':capture(['lscpu']),'compiler':capture([str(deps/'bin/x86_64-conda-linux-gnu-g++'),'--version']),'linked_libraries':capture(['ldd',str(BIN)]),'environment_controlled':CONTROL,'affinity':sorted(os.sched_getaffinity(0)),'load_average':os.getloadavg(),'memory':Path('/proc/meminfo').read_text(),'governors':{str(p):p.read_text().strip() for p in Path('/sys/devices/system/cpu').glob('cpu*/cpufreq/scaling_governor')},'perf_event_paranoid':Path('/proc/sys/kernel/perf_event_paranoid').read_text().strip(),'cgroup_limits':limits})

def run(case,phase,repetition,timeout=240):
    assert all(x in case['args'] for x in DEFAULTS)
    dest=RESULT/f"{case['name']}.{phase}{repetition}"
    if dest.with_suffix(dest.suffix+'.stdout').exists():raise RuntimeError('Refusing to overwrite '+str(dest))
    output=Path(str(dest)+'.stdout');err=Path(str(dest)+'.stderr');timing=Path(str(dest)+'.time')
    cmd=['/usr/bin/time','-f','%e\t%U\t%S\t%M\t%F\t%R\t%w\t%c\t%I\t%O\t%x','-o',str(timing),str(BIN)]+case['args']
    start=time.perf_counter();stamp=datetime.now(timezone.utc).isoformat();load=os.getloadavg();timed_out=False
    with output.open('wb') as out,err.open('wb') as stderr:
        p=subprocess.Popen(cmd,env=ENV,stdout=out,stderr=stderr,start_new_session=True)
        def expire():
            nonlocal timed_out
            if p.poll() is None:
                timed_out=True
                try:os.killpg(p.pid,signal.SIGKILL)
                except ProcessLookupError:pass
        watchdog=threading.Timer(timeout,expire);watchdog.daemon=True;watchdog.start()
        p.wait();watchdog.cancel()
    elapsed=time.perf_counter()-start
    row={'case':case['name'],'category':case['category'],'phase':phase,'repetition':repetition,'command':cmd,'started_utc':stamp,'ended_utc':datetime.now(timezone.utc).isoformat(),'elapsed_s':elapsed,'exit_code':p.returncode,'timeout':timed_out,'load_before':load,'load_after':os.getloadavg(),'binary_sha256':sha(BIN),'output_sha256':sha(output),'canonical_sha256':hashlib.sha256(canonical(output.read_bytes())).hexdigest()}
    values=timing.read_text().strip().splitlines()
    if values and len(values[-1].split('\t'))==len(TIME_KEYS):row.update(dict(zip(TIME_KEYS,map(float,values[-1].split('\t')))))
    row['csv_valid'],row['output_rows']=validate(output.read_bytes())
    with (RESULT/'runs.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(f"{phase}{repetition} {case['name']}: {elapsed:.3f}s rc={p.returncode} RSS={row.get('max_rss_kib',0)/1024:.2f}MiB",flush=True)
    if p.returncode or not row['csv_valid']:raise RuntimeError('Run failed '+str(dest))
    return row

def summarize():
    records=[json.loads(line) for line in (RESULT/'runs.jsonl').read_text().splitlines()];summary=[]
    for case in json.loads((ROOT/'cases.json').read_text()):
        rows=[r for r in records if r['case']==case['name'] and r['phase']=='measure'];values=[r['elapsed_s'] for r in rows]
        if not rows:continue
        structural={structure_hash((RESULT/f"{r['case']}.{r['phase']}{r['repetition']}.stdout").read_bytes()) for r in rows}
        summary.append({'case':case['name'],'category':case['category'],'threads':case['threads'],'predictions':case['predictions'],'n':len(rows),'median_s':statistics.median(values),'min_s':min(values),'max_s':max(values),'range_over_median_pct':100*(max(values)-min(values))/statistics.median(values),'median_user_s':statistics.median(r['user_s'] for r in rows),'median_system_s':statistics.median(r['system_s'] for r in rows),'median_rss_mib':statistics.median(r['max_rss_kib']/1024 for r in rows),'peak_rss_mib':max(r['max_rss_kib']/1024 for r in rows),'stable_output':len({r['canonical_sha256'] for r in rows})==1,'stable_bpList':len(structural)==1,'reported_rows':rows[0]['output_rows']})
    dump(RESULT/'summary.json',summary)
    with (RESULT/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=summary[0].keys(),lineterminator='\n');writer.writeheader();writer.writerows(summary)
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','summarize','environment']);a=p.parse_args()
    if a.action=='prepare':return prepare()
    if a.action=='environment':return environment()
    if a.action=='summarize':return summarize()
    cases=json.loads((ROOT/'cases.json').read_text());environment()
    for case in cases:run(case,'warmup',0)
    rng=random.Random(20260921)
    for repetition in range(1,6):
        order=cases[:];rng.shuffle(order)
        for case in order:run(case,'measure',repetition)
    summarize()
if __name__=='__main__':main()
