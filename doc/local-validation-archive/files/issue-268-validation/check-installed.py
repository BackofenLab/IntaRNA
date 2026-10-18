from pathlib import Path
import os, shlex, subprocess
root=Path(__file__).resolve().parent
build=root/'release';install=root/'install-release'
deps=root.parent/'issue-246-validation/toolchain'
cxx=root.parent/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
log=root/'logs/installed-api.log'
with log.open('w') as stream:
    def run(command,**kwargs):
        print(shlex.join(map(str,command)),file=stream,flush=True)
        return subprocess.run(command,check=True,stdout=stream,stderr=subprocess.STDOUT,**kwargs)
    run(['make','-C',str(build),'install'])
    for name in ('StackSeedDomain.h','SeededPartitionFunction.h'):
        probe=root/('installed-'+name+'.cpp')
        probe.write_text('#include <IntaRNA/'+name+'>\nint main() { return 0; }\n')
        run([str(cxx),'-std=c++23','-I'+str(install/'include'),'-I'+str(deps/'include'),'-fsyntax-only',str(probe)])
    env=dict(os.environ,PKG_CONFIG_PATH=str(install/'lib/pkgconfig'))
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','IntaRNA'],env=env,text=True))
    exe=root/'installed-probabilities-example'
    run([str(cxx),'-std=c++23',str(build/'doc/seeded-probabilities-example.cpp'),
         '-L'+str(deps/'lib'),'-Wl,-rpath,'+str(deps/'lib'),*flags,'-o',str(exe)])
    result=subprocess.check_output([str(exe)],text=True,env=dict(env,OMP_NUM_THREADS='1'))
    (root/'installed-probabilities-example.csv').write_text(result)
    assert result.startswith('bpProb;'),result
    print('Standalone headers and installed pkg-config probability consumer passed',file=stream,flush=True)
print('Installed API checks passed:',log)
