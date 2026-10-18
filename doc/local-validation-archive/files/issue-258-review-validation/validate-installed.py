from pathlib import Path
import os, resource, shlex, subprocess

root=Path('/home/alex/CodexProjects/IntaRNA')
validation=root/'issue-258-review-validation'
deps=root/'issue-246-validation/toolchain'
compiler=root/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
source=root/'pr-worktrees/issue-257-plan'

def run(args, **kwargs):
    return subprocess.run([str(x) for x in args],check=True,**kwargs)

def no_core():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))

for mode,build in [('release',source),('debug',validation/'debug')]:
    prefix=root/f'issue-258-validation/install-review-{mode}'
    with (validation/f'logs/{mode}-install.log').open('w') as log:
        run(['make','install'],cwd=build,stdout=log,stderr=subprocess.STDOUT)
    env=os.environ.copy()
    env['PKG_CONFIG_PATH']=str(prefix/'lib/pkgconfig')
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','IntaRNA'],env=env,text=True))
    link=[f'-L{deps}/lib',f'-Wl,-rpath,{deps}/lib']
    if mode=='release':
        for header in ['StackSeedDomain.h','SeededPartitionFunction.h']:
            probe=validation/'header-probe.cpp'
            probe.write_text(f'#include <IntaRNA/{header}>\nint main() {{}}\n')
            run([compiler,'-std=c++23',f'-I{prefix}/include',f'-I{deps}/include','-fsyntax-only',probe])
            print(f'PASS installed standalone header {header}',flush=True)
        consumer=validation/'installed-probabilities-example'
        run([compiler,'-std=c++23',source/'doc/seeded-probabilities-example.cpp',*link,*flags,'-o',consumer])
        result=subprocess.check_output([consumer],text=True)
        assert 'bpProb;' in result
        (validation/'installed-example.csv').write_text(result)
        print('PASS installed pkg-config probability consumer',flush=True)
    executable=validation/f'check-invariant-{mode}'
    run([compiler,'-std=c++23','-DHAVE_CONFIG_H',f'-I{build}/src',validation/'check-invariant.cpp',*link,*flags,'-o',executable])
    valid=subprocess.check_output([executable],text=True)
    assert 'secondLengthCalls=0' in valid if mode=='release' else 'secondLengthCalls=4' in valid
    invalid=subprocess.run([executable,'mismatch'],text=True,capture_output=True,preexec_fn=no_core)
    if mode=='release':
        assert invalid.returncode==0 and 'secondLengthCalls=0' in invalid.stdout
    else:
        assert invalid.returncode==-6 and 'len == handler.getSeedLength2(i,j)' in invalid.stderr
    print(f'PASS {mode} invariant: {valid.strip()}; mismatch exit={invalid.returncode}',flush=True)
