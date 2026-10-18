from pathlib import Path
import subprocess, shlex, os
p=Path(__file__).resolve().parent
repo=Path('/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-265')
deps=Path('/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain')
cxx='/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
flags=[cxx,'-std=c++23','-O2','-DNDEBUG','-DINTARNA_MULITHREADING=0','-pthread','-DELPP_FEATURE_PERFORMANCE_TRACKING','-DELPP_NO_DEFAULT_LOG_FILE','-Wno-deprecated','-I'+str(repo/'src'),'-I'+str(deps/'include'),'-I'+str(deps/'include/ViennaRNA')]
log=(p/'commands.log').open('w')
def run(args):
 print(shlex.join(args),file=log,flush=True)
 subprocess.run(['taskset','-c','0-7']+args,stdout=log,stderr=subprocess.STDOUT,check=True)
for name in ['VrnaHandler','InteractionEnergyVrna']:
 run(flags+['-c',str(repo/'src/IntaRNA'/ (name+'.cpp')),'-o',str(p/(name+'.o'))])
for mode in ['fallback','control']:
 f=flags[:]
 if mode=='fallback': f[1:1]=['-I'+str(p/'override')]
 run(f+['-c',str(repo/'src/IntaRNA/AccessibilityVrna.cpp'),'-o',str(p/('AccessibilityVrna-'+mode+'.o'))])
 run(f+[str(p/'probe.cpp'),str(p/'VrnaHandler.o'),str(p/'InteractionEnergyVrna.o'),str(p/('AccessibilityVrna-'+mode+'.o')),str(repo/'src/IntaRNA/libIntaRNA.a'),str(repo/'src/libeasylogging.a'),'-L'+str(deps/'lib'),'-Wl,-rpath,'+str(deps/'lib'),'-Wl,--wrap=vrna_probs_window','-lboost_regex','-lboost_iostreams','-lboost_serialization','-lRNA','-lgomp','-lgsl','-lgslcblas','-lmpfr','-lgmp','-lz','-o',str(p/('probe-'+mode))])
 run([str(p/('probe-'+mode))])
print('Both branch probes passed')
