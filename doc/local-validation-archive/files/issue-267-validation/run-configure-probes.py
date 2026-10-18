from pathlib import Path
import json, os, re, shutil, subprocess
root=Path('/home/alex/CodexProjects/IntaRNA')
source=root/'pr-worktrees/issue-267'
validation=root/'issue-267-validation'
probe_source=validation/'probe-source'
deps=root/'issue-246-validation/toolchain'
gcc14=root/'issue-246-validation/toolchain-gcc14/bin'
logs=validation/'logs'
files=subprocess.check_output(['git','ls-files','-z','--cached','--others','--exclude-standard'],cwd=source).decode().split('\0')
for name in files:
    if not name: continue
    target=probe_source/name
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source/name,target)
with (logs/'bootstrap-probes.log').open('w') as log:
    subprocess.run(['bash','autotools-init.sh'],cwd=probe_source,stdout=log,stderr=subprocess.STDOUT,check=True)
fixtures=validation/'configure-fixtures'
old=fixtures/'old-headers/include/ViennaRNA'
old.mkdir(parents=True,exist_ok=True)
oldconfig=(deps/'include/ViennaRNA/vrna_config.h').read_text()
oldconfig=re.sub(r'(#define VRNA_VERSION\s+)"[^"]+"',r'\1"2.6.4"',oldconfig)
oldconfig=re.sub(r'(#define VRNA_VERSION_MINOR\s+)\d+',r'\g<1>6',oldconfig)
oldconfig=re.sub(r'(#define VRNA_VERSION_PATCH\s+)\d+',r'\g<1>4',oldconfig)
(old/'vrna_config.h').write_text(oldconfig)
oldpc=fixtures/'old-metadata/lib/pkgconfig'
oldpc.mkdir(parents=True,exist_ok=True)
metadata=(deps/'lib/pkgconfig/RNAlib2.pc').read_text().replace('Version: 2.7.2','Version: 2.4.14')
(oldpc/'RNAlib2.pc').write_text(metadata)
missing_decl=fixtures/'missing-logger-declaration/include/ViennaRNA/utils'
missing_decl.mkdir(parents=True,exist_ok=True)
(missing_decl/'log.h').write_text('/* Test fixture: headers omit the required logger declaration. */\n')
missing_symbol=fixtures/'missing-logger-symbol/include/ViennaRNA/utils'
missing_symbol.mkdir(parents=True,exist_ok=True)
(missing_symbol/'log.h').write_text('#define vrna_log_lock_set intarna_test_unavailable_logger_lock_set\n#include "'+str(deps/'include/ViennaRNA/utils/log.h')+'"\n')
common=[str(probe_source/'configure'),
    'CC='+str(gcc14/'x86_64-conda-linux-gnu-gcc'),
    'CXX='+str(gcc14/'x86_64-conda-linux-gnu-g++'),
    'LDFLAGS=-L'+str(deps/'lib')+' -Wl,-rpath,'+str(deps/'lib'),
    '--with-boost='+str(deps),'--with-boost-libdir='+str(deps/'lib'),
    '--with-zlib='+str(deps),'--disable-debug']
static_libs='LIBS=-lgomp -lgsl -lgslcblas -lpthread -lmpfr -lgmp -lstdc++ -lm'
cases=[
    ('no-pkg-config',0,[],['--disable-pkg-config',static_libs],deps,None),
    ('old-pkg-metadata',1,[],[],fixtures/'old-metadata','Package dependency requirement'),
    ('old-headers-no-pkg',1,[old.parent],['--disable-pkg-config',static_libs],deps,'ViennaRNA 2.7.0 or newer is required'),
    ('old-headers-current-pkg',1,[old.parent],[],deps,'ViennaRNA 2.7.0 or newer is required'),
    ('missing-logger-declaration',1,[missing_decl.parents[1]],[],deps,'vrna_log_lock_set'),
    ('missing-logger-symbol',1,[missing_symbol.parents[1]],[],deps,'undefined reference to `intarna_test_unavailable_logger_lock_set'),
]
results=[]
for name,expected,inc,extra,vrna,diagnostic in cases:
    build=validation/'configure-probes'/name
    build.mkdir(parents=True,exist_ok=True)
    cpp='CPPFLAGS='+' '.join('-I'+str(path) for path in inc+[deps/'include'])
    args=common+[cpp,'--with-vrna='+str(vrna)]+extra
    env=os.environ.copy()
    for variable in ('VRNA_CFLAGS','VRNA_LIBS','PKG_CONFIG_PATH','PKG_CONFIG_LIBDIR'):
        env.pop(variable,None)
    with (logs/('configure-'+name+'.log')).open('w') as log:
        completed=subprocess.run(args,cwd=build,env=env,stdout=log,stderr=subprocess.STDOUT)
    config=(build/'config.log').read_text()
    output=(logs/('configure-'+name+'.log')).read_text()
    passed=(completed.returncode==0 if expected==0 else completed.returncode!=0)
    if expected:
        passed=passed and 'Vienna RNA C library version >= 2.7.0 is required.' in output
    if diagnostic:
        passed=passed and diagnostic in config
    result={'case':name,'expected':'success' if expected==0 else 'rejection','exit_code':completed.returncode,'passed':passed,'command':args,'log':str(logs/('configure-'+name+'.log'))}
    results.append(result)
    print(json.dumps({key:result[key] for key in ('case','expected','exit_code','passed')}),flush=True)
(validation/'configure-probes/results.json').write_text(json.dumps(results,indent=2)+'\n')
if not all(r['passed'] for r in results):
    raise SystemExit(1)
