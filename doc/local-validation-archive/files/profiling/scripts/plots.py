#!/usr/bin/env python3
import json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/intarna-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
ROOT=Path(__file__).resolve().parents[1]
D={r['case']:r for r in json.loads((ROOT/'results/summary.json').read_text())}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.titleweight':'bold','axes.titlesize':12,'figure.facecolor':'#fafaf8','axes.facecolor':'#fafaf8'})
blue='#226e8b'; orange='#b86132'; green='#448065'
fig,axes=plt.subplots(2,3,figsize=(17,10),layout='constrained')

def line(ax,xx,names,color,label):
    yy=[D[n]['median_s'] for n in names]
    lo=[D[n]['median_s']-D[n]['min_s'] for n in names]
    hi=[D[n]['max_s']-D[n]['median_s'] for n in names]
    ax.errorbar(xx,yy,yerr=[lo,hi],fmt='o-',color=color,capsize=3,label=label)

ax=axes[0,0]; ns=[100,300,1000,3000,10000]
line(ax,ns,[f'length_t{n}_q100' for n in ns],blue,'Query = 100 nt')
ax.set(xscale='log',yscale='log',xlabel='Target length (nt)',ylabel='Wall time (seconds)',title='Target length: fixed query and folding window')
ax.set_xticks(ns,labels=[str(n) for n in ns]);ax.grid(alpha=.15)

ax=axes[0,1]; ts=[1,2,4,6,12]; base=D['throughput_1t']['median_s']
speed=[base/D[f'throughput_{n}t']['median_s'] for n in ts]
ax.plot(ts,speed,'o-',color=blue,label='Measured speedup')
ax.plot(ts,ts,'--',color='#aaaaaa',label='Ideal');ax.set(xlabel='IntaRNA threads',ylabel='Speedup vs one thread',title='48 biological targets × ChiX');ax.set_xticks(ts);ax.legend(frameon=False)

ax=axes[0,2]; names=['algorithm_default','algorithm_exact_x','algorithm_heuristic_s','algorithm_exact_s','algorithm_helix','algorithm_no_seed']
labels=['X heuristic','X exact','S heuristic','S exact','B heuristic','S / no seed']
values=[D[n]['median_s'] for n in names]; bars=ax.barh(labels[::-1],values[::-1],color=blue)
ax.set(xscale='log',xlabel='Wall time (seconds; log scale)',title='Algorithm choices: same 100 × 100 nt input')
for bar,v in zip(bars,values[::-1]):ax.text(v*1.06,bar.get_y()+bar.get_height()/2,f'{v:.3f}',va='center',fontsize=9)
ax.set_xlim(right=max(values)*2.5)

ax=axes[1,0]; names=['no_pair','composition_gc20','composition_gc50','composition_gc80','au_repeat','gc_repeat']
labels=['All A (no pairing)','20% GC distribution','50% GC distribution','80% GC distribution','AU repeat','GC repeat']
v=[D[n]['median_s'] for n in names]; bars=ax.barh(labels[::-1],v[::-1],color=green)
ax.set(xlabel='Wall time (seconds)',title='Sequence composition: 300 × 100 nt')
for bar,x in zip(bars,v[::-1]):ax.text(x+.035,bar.get_y()+bar.get_height()/2,f'{x:.3f}',va='center',fontsize=9)
ax.set_xlim(right=max(v)*1.2)

ax=axes[1,1]; names=['accessibility_computed','accessibility_global','accessibility_none','cache_matched_reused'];labels=['Local (default)','Global','Disabled','Reused ED*']
v=[D[n]['median_s'] for n in names];bars=ax.bar(labels,v,color=[blue,orange,'#829ca6',green])
ax.set(ylabel='Wall time (seconds)',title='Accessibility: 1000 × 100 nt; interaction cap 150')
for bar,x in zip(bars,v):ax.text(bar.get_x()+bar.get_width()/2,x+.06,f'{x:.3f}',ha='center')
ax.text(.02,.97,'Global/disabled alter the model. *ED caps = 149/99.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')
ax.set_ylim(top=max(v)*1.25)

ax=axes[1,2]; names=['memory_t3000_q800','memory_window300','memory_s','memory_s_window300'];labels=['X whole','X window 300','S whole','S window 300']
v=[D[n]['median_rss_mib'] for n in names if n in D]
if len(v)==4:
    bars=ax.bar(labels,v,color=[blue,'#88b5c6',orange,'#d5a38a'])
    for bar,x in zip(bars,v):ax.text(bar.get_x()+bar.get_width()/2,x+max(v)*.015,f'{x:.1f}',ha='center')
ax.set(ylabel='Peak resident memory (MiB)',title='Windowing memory effect: 3000 × 800 nt');ax.tick_params(axis='x',labelrotation=12)
fig.suptitle('IntaRNA profiling • commit 1b83779 • Ryzen 5 7530U\nMedians of 3 runs; clean timings, separate from instrumentation',fontsize=18,fontweight='bold')
out=ROOT/'results';fig.savefig(out/'overview.png',dpi=150);fig.savefig(out/'overview.svg');plt.close(fig)
print(out/'overview.png')
