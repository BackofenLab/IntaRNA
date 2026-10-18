#!/usr/bin/env python3
"""Separate CPU, diagnostic-phase, and heap runs of the requested default mode."""
import json,subprocess,time,hashlib
from datetime import datetime,timezone
from benchmark import ROOT,BIN,ENV,INPUT,fasta,write_fasta,canonical,validate,dump,sha
DEST=ROOT/'profiles'
def execute(name,command,timeout=420):
    if (DEST/(name+'.stdout')).exists():raise RuntimeError('Refusing to overwrite '+name)
    wrapped=['/usr/bin/time','-f','%e\t%U\t%S\t%M','-o',str(DEST/(name+'.time'))]+command
    started=datetime.now(timezone.utc).isoformat()
    with (DEST/(name+'.stdout')).open('wb') as out,(DEST/(name+'.stderr')).open('wb') as err:
        p=subprocess.run(wrapped,env=ENV,stdout=out,stderr=err,timeout=timeout)
    row={'name':name,'command':command,'exit_code':p.returncode,'started_utc':started,'ended_utc':datetime.now(timezone.utc).isoformat(),'resources':(DEST/(name+'.time')).read_text().strip()}
    with (DEST/'commands.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(name,p.returncode,flush=True);p.check_returncode()
    return (DEST/(name+'.stdout')).read_bytes()

def main():
    DEST.mkdir(exist_ok=True)
    cases={x['name']:x for x in json.loads((ROOT/'cases.json').read_text())};checks=[]
    records=fasta(INPUT/'promoters48_systematic.fa')
    repeated=write_fasta('profile_chix384',[(f'repeat{k}_{n}',seq) for k in range(8) for n,seq in records])
    dump(DEST/'input-manifest.json',{'profile_chix384.fa':{'sha256':sha(repeated),'records':384,'construction':'Eight copies of the same systematically selected 48 targets, with unique IDs; no genome-wide claim.'}})
    cpu=[('gcvb',cases['batch_gcvb_1t']['args']),('gcvb_repeat',cases['batch_gcvb_1t']['args']),('chix384',[('--target='+str(repeated)) if a.startswith('--target=') else a for a in cases['batch_chix_1t']['args']]),('long_pair',cases['long_pair']['args'])]
    for name,args in cpu:
        if name=='chix384':reference=execute(name+'_reference',[str(BIN)]+args)
        else:reference=(ROOT/'results'/(({'gcvb':'batch_gcvb_1t','gcvb_repeat':'batch_gcvb_1t','long_pair':'long_pair'}[name])+'.measure1.stdout')).read_bytes()
        exp=DEST/(name+'.er')
        actual=execute(name+'_cpu',['gprofng','collect','app','-p','100','-a','off','-o',str(exp),str(BIN)]+args)
        valid,rows=validate(actual);equal=canonical(actual)==canonical(reference)
        checks.append({'profile':name,'kind':'cpu','csv_valid':valid,'output_equal_bpList':equal,'reported_rows':rows});assert valid and equal,name
        for kind,sort in [('self','e.totalcpu'),('inclusive','i.totalcpu')]:
            execute(name+'_'+kind,['gprofng','display','text','-metrics','e.%totalcpu:i.%totalcpu:name','-sort',sort,'-limit','50','-functions',str(exp)])
        execute(name+'_statistics',['gprofng','display','text','-header','-statistics',str(exp)])
        execute(name+'_lines',['gprofng','display','text','-metrics','e.%totalcpu:i.%totalcpu:name','-sort','e.totalcpu','-limit','30','-lines',str(exp)])
        if name in ['gcvb','long_pair']:
            execute(name+'_loop_source',['gprofng','display','text','-metrics','e.%totalcpu:i.%totalcpu:name','-source','IntaRNA::InteractionEnergyVrna::getE_interLeft(unsigned long, unsigned long, unsigned long, unsigned long) const',str(exp)])
    for name in ['batch_gcvb_1t','batch_chix_1t','long_pair']:
        actual=execute(name+'_phases',[str(BIN)]+cases[name]['args']+['--verbose'])
        equal=canonical(actual)==canonical((ROOT/'results'/(name+'.measure1.stdout')).read_bytes())
        checks.append({'profile':name,'kind':'phases','output_equal_bpList':equal});assert equal,name
    for name in ['batch_gcvb_1t','long_pair']:
        exp=DEST/(name+'_heap.er')
        actual=execute(name+'_heap',['gprofng','collect','app','-p','off','-H','on','-a','off','-o',str(exp),str(BIN)]+cases[name]['args'])
        equal=canonical(actual)==canonical((ROOT/'results'/(name+'.measure1.stdout')).read_bytes())
        checks.append({'profile':name,'kind':'heap','output_equal_bpList':equal});assert equal,name
        execute(name+'_heap_summary',['gprofng','display','text','-limit','30','-heapstat','-allocs','-leaks',str(exp)])
    dump(DEST/'checks.json',checks)
if __name__=='__main__':main()
