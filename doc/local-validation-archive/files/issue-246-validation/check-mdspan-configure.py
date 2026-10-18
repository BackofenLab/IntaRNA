from pathlib import Path
import os
import subprocess

root=Path(__file__).resolve().parent
source=root/'mdspan-probe-source'
deps=root/'toolchain'
env=dict(os.environ,PKG_CONFIG_PATH=str(deps/'lib/pkgconfig'),LD_LIBRARY_PATH=str(deps/'lib'),OPENBLAS_NUM_THREADS='1')
cases=[('gcc14-auto','toolchain-gcc14','auto',0),('gcc14-std','toolchain-gcc14','std',None),('gcc16-auto','toolchain','auto',1),('gcc16-kokkos','toolchain','kokkos',0),('invalid-backend','toolchain','invalid',None)]
for label,compiler,backend,expected in cases:
    build=root/('probe-'+label)
    build.mkdir()
    args=[str(source/'configure'),f'CXX={root/compiler}/bin/x86_64-conda-linux-gnu-g++',f'CC={root/compiler}/bin/x86_64-conda-linux-gnu-gcc',f'CPPFLAGS=-I{deps}/include',f'LDFLAGS=-L{deps}/lib -Wl,-rpath,{deps}/lib',f'--with-boost={deps}',f'--with-boost-libdir={deps}/lib',f'--with-vrna={deps}',f'--with-mdspan={backend}']
    log=root/'logs'/f'probe-{label}.log'
    with log.open('w') as output:
        result=subprocess.run(args,cwd=build,env=env,stdout=output,stderr=subprocess.STDOUT)
    if expected is None:
        assert result.returncode!=0, label
        assert ('--with-mdspan=std requires' if backend=='std' else '--with-mdspan must be') in log.read_text(), label
    else:
        assert result.returncode==0,(label,log)
        for header in ('src/config.h','src/IntaRNA/intarna_config.h'):
            assert f'#define INTARNA_USE_STD_MDSPAN {expected}' in (build/header).read_text(),(label,header)
    print(f'{label}: expected outcome confirmed',flush=True)
# Confirm the fallback probe rejects an unusable bundled tree with its own diagnostic.
headers=source/'src/mdspan'
hidden=source/'src/mdspan-unavailable'
headers.rename(hidden)
try:
    build=root/'probe-missing-bundled'
    build.mkdir()
    args[-1]='--with-mdspan=kokkos'
    log=root/'logs/probe-missing-bundled.log'
    with log.open('w') as output:
        result=subprocess.run(args,cwd=build,env=env,stdout=output,stderr=subprocess.STDOUT)
    assert result.returncode!=0 and 'bundled Kokkos mdspan headers do not compile' in log.read_text()
    print('missing-bundled: expected diagnostic confirmed',flush=True)
finally:
    hidden.rename(headers)
