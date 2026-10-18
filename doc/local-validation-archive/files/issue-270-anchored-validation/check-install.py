from pathlib import Path
import os,shlex,subprocess
root=Path(__file__).resolve().parent
repo=root.parent/'pr-worktrees/issue-270-anchored'
prefix=root/'install-release'
deps=root.parent/'issue-246-validation/toolchain'
cxx=root.parent/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
subprocess.run(['make','install',f'prefix={prefix}'],cwd=root/'release',check=True)
for name in ('SeededPartitionFunction.h','StackSeedDomain.h','PredictorMfeEns2dHeuristicSeedExtension.h'):
    probe=root/'header-probe.cpp';probe.write_text(f'#include <IntaRNA/{name}>\nint main() {{ return 0; }}\n')
    subprocess.run([str(cxx),'-std=c++23',f'-I{prefix}/include',f'-I{deps}/include','-fsyntax-only',str(probe)],check=True)
    print('Standalone installed header:',name,flush=True)
env=dict(os.environ,PKG_CONFIG_PATH=str(prefix/'lib/pkgconfig'))
flags=shlex.split(subprocess.check_output(['pkg-config',f'--define-variable=prefix={prefix}','--cflags','--libs','IntaRNA'],env=env,text=True))
consumer=root/'installed-probabilities'
subprocess.run([str(cxx),'-std=c++23',str(repo/'doc/seeded-probabilities-example.cpp'),f'-Wl,-rpath,{prefix}/lib',f'-Wl,-rpath,{deps}/lib',*flags,'-o',str(consumer)],check=True)
result=subprocess.check_output([str(consumer)])
assert b'bpProb;' in result
(root/'installed-probabilities.csv').write_bytes(result)
print('Installed pkg-config probability consumer passed',flush=True)
