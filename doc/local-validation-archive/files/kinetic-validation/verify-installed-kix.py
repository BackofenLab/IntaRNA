import os
from pathlib import Path
import shlex
import subprocess
root = Path('/home/alex/CodexProjects/IntaRNA')
validation = root/'kinetic-validation'
prefix = validation/'kix-stage'/str(validation/'install-release').lstrip('/')
toolchain = root/'issue-246-validation/toolchain'
cxx = root/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
env = dict(os.environ, PKG_CONFIG_PATH=str(prefix/'lib/pkgconfig')+':'+str(toolchain/'lib/pkgconfig'))
pc = ['pkg-config', '--define-variable=prefix='+str(prefix)]
cflags = shlex.split(subprocess.check_output([*pc,'--cflags','IntaRNA'],env=env,text=True))
libs = shlex.split(subprocess.check_output([*pc,'--libs','IntaRNA'],env=env,text=True))
assert not (prefix/'include/IntaRNA/PredictorSeedExtensionKineticPruned.h').exists()
assert (prefix/'bin/IntaRNAkix').is_symlink()
for header in ['PredictorSeedExtensionKinetic.h','InteractionEnergyVrna.h']:
    source = validation/'header-probe-kix.cpp'
    source.write_text('#include <IntaRNA/'+header+'>\nint main() { return 0; }\n')
    subprocess.run([str(cxx),'-std=c++23','-fsyntax-only',*cflags,str(source)],check=True)
subprocess.run([str(cxx),'-std=c++23',*cflags,str(validation/'consumer.cpp'),'-L'+str(toolchain/'lib'),
                '-Wl,-rpath,'+str(toolchain/'lib'),*libs,'-leasylogging','-o',str(validation/'consumer-kix')],check=True)
subprocess.run([str(validation/'consumer-kix')],check=True)
args = ['--target=GGGGGG','--query=CCCCCC','--energy=B','--acc=N',
        '--seedTQ=3||&3||','--outNumber=10','--outMode=C','--outCsvCols=hybridDB,E','--default-log-file=/dev/null']
actual = subprocess.check_output([str(prefix/'bin/IntaRNAkix'),*args],text=True)
assert actual == 'hybridDB;E\n1||||||&1||||||;-6\n1||||&3||||;-4\n3||&3||;-2\n', actual
helptext = subprocess.check_output([str(prefix/'bin/IntaRNAkix'),'--fullhelp'],text=True)
assert 'IntaRNAkix' in helptext
print('Installed K/Vrna headers compile independently; installed K consumer links and predicts -400; installed IntaRNAkix link reports the expected prefixes; removed header absent.')
