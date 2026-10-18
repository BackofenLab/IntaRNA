#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import os
import statistics
import subprocess

root=Path(__file__).resolve().parent
env=dict(os.environ,OMP_NUM_THREADS='1')
cases=[('dense','many',40),('dense','few',40),('band','many',60),('band','few',60)]
samples=[]
checks=[]

def call(variant,case,pairs,dump='-'):
    topology,density,n=case
    command=[str(root/variant/'benchmark'),topology,density,str(n),str(pairs),dump]
    rss=root/'last-rss.txt'
    result=subprocess.check_output(['/usr/bin/time','-f','%M','-o',str(rss),*command],text=True,env=env)
    seconds,z,moment,boundaries,checksum,seeds=result.strip().split('\t')
    return {'variant':variant,'case':f'{topology}{n}-{density}','pairs':pairs,
            'seconds':float(seconds),'rss_kib':int(rss.read_text().strip()),'z':float(z),
            'pair_moment':float(moment),'boundaries':int(boundaries),'checksum':float(checksum),'seeds':int(seeds)}

def read_dump(path):
    values={}
    for line in path.read_text().splitlines():
        fields=line.split()
        if fields[0]=='b':
            values[tuple(fields[:5])]=tuple(map(float,fields[5:]))
        else:
            values[tuple(fields[:-1])]=(float(fields[-1]),)
    return values

for case in cases:
    name=f'{case[0]}{case[2]}-{case[1]}'
    paths={variant:root/f'{name}-{variant}.txt' for variant in ('base','candidate')}
    for variant,path in paths.items():
        call(variant,case,1,str(path))
    before,after=(read_dump(paths[variant]) for variant in ('base','candidate'))
    assert before.keys()==after.keys(),name
    maximum=0.0
    for key,values in before.items():
        for a,b in zip(values,after[key]):
            relative=abs(a-b)/max(abs(a),abs(b)) if a or b else 0.0
            maximum=max(maximum,relative)
            assert relative<1e-12,(name,key,a,b,relative)
    checks.append({'case':name,'compared_entries':len(before),'max_relative_error':maximum})
    for pairs in (0,1):
        for repetition in range(3):
            # Alternate order to reduce a consistent order bias.
            variants=('base','candidate') if repetition%2==0 else ('candidate','base')
            for variant in variants:
                sample=call(variant,case,pairs)
                sample['repetition']=repetition
                samples.append(sample)
    print(name,checks[-1],flush=True)

with (root/'samples.tsv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(samples[0]),delimiter='\t')
    writer.writeheader();writer.writerows(samples)
summary=[]
for case in cases:
    name=f'{case[0]}{case[2]}-{case[1]}'
    for pairs in (0,1):
        groups={variant:[row for row in samples if row['case']==name and row['pairs']==pairs and row['variant']==variant]
                for variant in ('base','candidate')}
        result={'case':name,'pairs':pairs,'seeds':groups['base'][0]['seeds']}
        for variant,rows in groups.items():
            result[variant+'_median_seconds']=statistics.median(row['seconds'] for row in rows)
            result[variant+'_min_seconds']=min(row['seconds'] for row in rows)
            result[variant+'_max_seconds']=max(row['seconds'] for row in rows)
            result[variant+'_median_rss_kib']=statistics.median(row['rss_kib'] for row in rows)
        result['runtime_ratio']=result['candidate_median_seconds']/result['base_median_seconds']
        result['rss_ratio']=result['candidate_median_rss_kib']/result['base_median_rss_kib']
        summary.append(result)
(root/'results.json').write_text(json.dumps({'checks':checks,'summary':summary},indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
