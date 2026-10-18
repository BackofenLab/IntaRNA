#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import subprocess

root=Path(__file__).resolve().parent
repo=root.parents[1]/'pr-worktrees'/'issue-264'
toolchain=root.parents[1]/'issue-246-validation'/'toolchain-gcc14'
dependencies=root.parents[1]/'issue-246-validation'/'toolchain'
compiler=toolchain/'bin'/'x86_64-conda-linux-gnu-g++'
commands=[]
for variant in ('base','candidate'):
    command=[str(compiler),'-std=c++23','-O3','-DNDEBUG','-ffunction-sections','-fdata-sections',
             '-I'+str(root/variant),'-I'+str(repo/'src'),'-I'+str(dependencies/'include'),
             str(root/'driver.cpp'),str(root/variant/'IntaRNA'/'SeededPartitionFunction.cpp'),str(root/'StackSeedDomain.cpp'),
             '-Wl,--gc-sections','-Wl,-rpath,'+str(toolchain/'lib'),'-o',str(root/variant/'benchmark')]
    subprocess.run(command,check=True)
    commands.append(command)
manifest={'baseline':'cb1c9d419f96c7d2a10bae864d922dcb0234f4d4',
          'compiler':subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],
          'commands':commands,'sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [root/'driver.cpp',root/'StackSeedDomain.cpp',*(root/variant/'IntaRNA'/name for variant in ('base','candidate') for name in ('SeededPartitionFunction.cpp','SeededPartitionFunction.h'))]}}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
