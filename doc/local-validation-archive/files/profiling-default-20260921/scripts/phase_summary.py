#!/usr/bin/env python3
"""Parse diagnostic phase logs, retaining integer-duration truncation bounds."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for path in sorted((ROOT/'profiles').glob('*_phases.stdout')):
    data=path.read_text()
    groups={k:dict(calls=0,lower_s=0,upper_s=0) for k in ['accessibility','prediction_including_seed','seed_nested']}
    detailed=[]
    for fun,n,unit in re.findall(r'Executed \[(.*?)\] in \[(\d+) ([^\]]+)\]',data):
        scale={'ms':.001,'seconds':1,'minutes':60,'hours':3600}.get(unit)
        if scale is None:raise ValueError(unit)
        lo=int(n)*scale;hi=lo+scale
        if 'AccessibilityVrna::' in fun:group='accessibility'
        elif '::predict(' in fun:group='prediction_including_seed'
        elif '::fillSeed(' in fun:group='seed_nested'
        else:group=None
        if group:
            groups[group]['calls']+=1;groups[group]['lower_s']+=lo;groups[group]['upper_s']+=hi
        detailed.append(dict(function=fun,displayed=f'{n} {unit}',lower_s=lo,upper_s=hi,group=group))
    rows.append(dict(profile=path.stem.removesuffix('_phases'),groups=groups,phases=detailed,
        valid_seed_count_sum=sum(map(int,re.findall(r'valid seeds = (\d+)',data)))))
(ROOT/'profiles/phase-summary.json').write_text(json.dumps(rows,indent=2)+'\n')
for row in rows:
    print(row['profile'],row['groups'])
