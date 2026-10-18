#!/usr/bin/env python3
"""Separate intrusive profiles from benchmark runs. Execute only with host idle."""
import hashlib, json, os, subprocess, time
from pathlib import Path
from benchmark import ROOT, BIN, ENV, INPUT, canonical, write_fasta, fasta

DEST=ROOT/'profiles'

def profile_output(data):
    return canonical(b'\n'.join(line for line in data.splitlines()
        if not line.startswith(b'Creating experiment directory ')))

def execute(name, cmd, timeout=180):
    print(name, flush=True)
    start=time.perf_counter()
    wrapped=['/usr/bin/time','-f','%e\t%U\t%S\t%M','-o',str(DEST/(name+'.time'))]+cmd
    with (DEST/(name+'.stdout')).open('wb') as out, (DEST/(name+'.stderr')).open('wb') as err:
        p=subprocess.run(wrapped, stdout=out, stderr=err, env=ENV, timeout=timeout)
    row=dict(name=name,command=cmd,exit_code=p.returncode,wall_s=time.perf_counter()-start)
    vals=(DEST/(name+'.time')).read_text().strip().splitlines()
    if vals and len(vals[-1].split('\t'))==4:
        row.update(dict(zip(['gnu_wall_s','user_s','system_s','max_rss_kib'],map(float,vals[-1].split('\t')))))
    with (DEST/'commands.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    return row

def main():
    DEST.mkdir(exist_ok=True)
    cases={c['name']:c for c in json.loads((ROOT/'cases.json').read_text())}
    exact_targets=write_fasta('profile_exact_targets',[(f'exact_repeat_{i}',fasta(INPUT/'t100.fa')[0][1]) for i in range(128)])
    ensemble_targets=write_fasta('profile_ensemble_targets',[(f'ensemble_repeat_{i}',fasta(INPUT/'t60.fa')[0][1]) for i in range(32)])
    (DEST/'input-manifest.json').write_text(json.dumps({str(p):dict(sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest(),
        records=len(fasta(p))) for p in [exact_targets,ensemble_targets]},indent=2)+'\n')
    specs=[
      ('biological_batch', ['--target='+str(INPUT/'promoters48.fa'),'--query='+str(INPUT/'GcvB.fa'),
          '--threads=1','--outMode=C','--outCsvCols=id1,id2,start1,end1,start2,end2,E']),
      ('hybrid_noaccess', [('--query='+str(INPUT/'q400.fa')) if a.startswith('--query=') else a
          for a in cases['memory_t3000_q800']['args']]+['--acc=N','--intLenMax=150']),
      ('exact_x', [('--target='+exact_targets) if a.startswith('--target=') else a for a in cases['algorithm_exact_x']['args']]),
      ('ensemble_m', ['--target='+ensemble_targets,'--query='+str(INPUT/'q60.fa'),
          '--threads=1','--model=P','--mode=M','--noSeed','--outMode=C',
          '--outCsvCols=id1,id2,start1,end1,start2,end2,E']),
    ]
    for name,args in specs:
        clean=execute(name+'_reference',[str(BIN)]+args)
        exp=str(DEST/(name+'.er'))
        prof=execute(name+'_cpu',['gprofng','collect','app','-p','100','-a','off','-o',exp,str(BIN)]+args)
        for kind,sort in [('self','e.totalcpu'),('inclusive','i.totalcpu')]:
            execute(name+'_'+kind,['gprofng','display','text','-metrics','e.%totalcpu:i.%totalcpu:name',
                '-sort',sort,'-limit','40','-functions',exp])
        execute(name+'_statistics',['gprofng','display','text','-header','-statistics',exp])
        execute(name+'_phases',[str(BIN)]+args+['--verbose'])
        ref=(DEST/(name+'_reference.stdout')).read_bytes()
        sampled=(DEST/(name+'_cpu.stdout')).read_bytes()
        with (DEST/'validation.jsonl').open('a') as f:
            f.write(json.dumps(dict(name=name,output_equal=profile_output(ref)==profile_output(sampled),
                clean_s=clean['wall_s'],sampled_s=prof['wall_s'],
                note='collector startup/postprocessing included in sampled wall time'))+'\n')
    for name in ['bio_phob_gcvb','algorithm_exact_x','memory_t3000_q800','memory_s']:
        args=cases[name]['args']; exp=str(DEST/(name+'_heap.er'))
        execute(name+'_heap',['gprofng','collect','app','-p','off','-H','on','-a','off','-o',exp,str(BIN)]+args)
        execute(name+'_heap_summary',['gprofng','display','text','-limit','30','-heapstat','-allocs','-leaks',exp])
    for name in ['bio_phob_gcvb','accessibility_computed','accessibility_none','accessibility_global',
                 'length_t3000_q100','window_300','throughput_1t','throughput_6t']:
        execute(name+'_phases',[str(BIN)]+cases[name]['args']+['--verbose'])

if __name__=='__main__':main()
