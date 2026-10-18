#!/usr/bin/env python3
"""Static scientific figure, generated only from saved summaries."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
S={x['case']:x for x in json.loads((ROOT/'results/summary.json').read_text())}
CPU=json.loads((ROOT/'profiles/cpu-summary.json').read_text())
PH=json.loads((ROOT/'profiles/phase-summary.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.titlelocation':'left','axes.grid':True,'grid.alpha':.18})
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
blue='#175577';green='#237c70';orange='#be6a2e';gray='#cbd4d9'
a=axs[0,0];ns=[100,300,1000,3000,10000];rows=[S[f'length_t{n}_q100'] for n in ns]
med=np.array([r['median_s'] for r in rows]);lo=med-np.array([r['min_s'] for r in rows]);hi=np.array([r['max_s'] for r in rows])-med
a.errorbar(ns,med,yerr=[lo,hi],color=blue,marker='o',capsize=4)
a.set(xscale='log',yscale='log',xlabel='Target length (nt); query fixed at 100 nt',ylabel='Wall time (s)',title='A  Length scaling within the measured range')
a.set_xticks(ns,[f'{n:,}' for n in ns]);a.grid(True,which='minor',alpha=.06)
a=axs[0,1]
for rna,color,label in [('gcvb',blue,'GcvB (201 nt)'),('chix',green,'ChiX (84 nt)')]:
    ts=[1,6,12] if rna=='gcvb' else [1,6];base=S[f'batch_{rna}_1t']['median_s'];a.plot(ts,[base/S[f'batch_{rna}_{t}t']['median_s'] for t in ts],marker='o',color=color,label=label)
    if rna=='chix':a.scatter([12],[base/S['batch_chix_12t']['median_s']],marker='x',color=orange,s=60,label='ChiX 12t: energy mismatch')
a.set(xlabel='OpenMP threads (6 physical cores / 12 hardware threads)',ylabel='Observed time ratio relative to one thread',xticks=[1,6,12],title='B  48-target biological batch throughput');a.legend(frameon=False)
a=axs[1,0];left=np.zeros(4);short=['GcvB run 1','GcvB run 2','ChiX (384 pairs)','Long pair']
for text,needle,color in [('Loop-energy evaluation','::getE_interLeft(',blue),('Left extension','PredictorMfe2dHeuristicSeedExtension::fillHybridE_left(',green),('Right extension','PredictorMfe2dHeuristicSeedExtension::fillHybridE_right(',orange)]:
    values=np.array([sum(f['self_pct'] for f in r['self_functions'] if needle in f['function']) for r in CPU]);a.barh(short,values,left=left,color=color,label=text);left+=values
a.barh(short,100-left,left=left,color=gray,label='Other');a.invert_yaxis();a.set(xlim=(0,100),xlabel='Exclusive sampled CPU (%)',title='C  Default-mode hotspots (100 ms sampling)');a.legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(.5,-.45),ncol=2)
a=axs[1,1];short=[];acc=[];pred=[]
for name,label in [('batch_gcvb_1t','GcvB 48'),('batch_chix_1t','ChiX 48'),('long_pair','Long pair')]:
    r=next(x for x in PH if x['profile']==name);g=r['groups'];aa=g['accessibility']['lower_s'];pp=g['prediction_including_seed']['lower_s'];short.append(label);acc.append(aa);pred.append(pp)
a.bar(short,acc,color=green,label='Local accessibility');a.bar(short,pred,bottom=acc,color=blue,label='Predictor, including seed')
a.set(ylabel='Logged phase time (s; lower truncation bound)',title='D  Phase attribution, measured separately');a.legend(frameon=False,fontsize=8)
fig.suptitle('IntaRNA default mode: seed-extension heuristic + local accessibility\nPrediction windowing disabled | master 0a4568a | 21 September 2026',fontsize=14,fontweight='bold')
fig.savefig(ROOT/'results/overview.png',dpi=180);fig.savefig(ROOT/'results/overview.svg')
