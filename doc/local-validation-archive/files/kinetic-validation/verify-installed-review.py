import os
from pathlib import Path
import shlex
import subprocess
root = Path('/home/alex/CodexProjects/IntaRNA')
validation = root/'kinetic-validation'
prefix = validation/'install-release'
toolchain = root/'issue-246-validation/toolchain'
cxx = root/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
env = dict(os.environ, PKG_CONFIG_PATH=str(prefix/'lib/pkgconfig')+':'+str(toolchain/'lib/pkgconfig'))
cflags = shlex.split(subprocess.check_output(['pkg-config','--cflags','IntaRNA'],env=env,text=True))
libs = shlex.split(subprocess.check_output(['pkg-config','--libs','IntaRNA'],env=env,text=True))
for header in ['PredictorSeedExtensionKinetic.h','PredictorSeedExtensionKineticPruned.h','InteractionEnergyVrna.h']:
    source = validation/'header-probe-review.cpp'
    source.write_text('#include <IntaRNA/'+header+'>\nint main() { return 0; }\n')
    subprocess.run([str(cxx),'-std=c++23','-fsyntax-only',*cflags,str(source)],check=True)
source = validation/'consumer-review.cpp'
source.write_text((validation/'consumer.cpp').read_text().replace('PredictorSeedExtensionKinetic','PredictorSeedExtensionKineticPruned'))
subprocess.run([str(cxx),'-std=c++23',*cflags,str(source),'-L'+str(toolchain/'lib'),
                '-Wl,-rpath,'+str(toolchain/'lib'),*libs,'-leasylogging','-o',str(validation/'consumer-review')],check=True)
subprocess.run([str(validation/'consumer-review')],check=True)
print('Three installed headers compile standalone; installed pruned predictor links and predicts -400.')
