from pathlib import Path
import shlex
import shutil
import subprocess
root = Path('/home/alex/CodexProjects/IntaRNA')
source = root/'pr-worktrees/kinetic-seed-extension/src/IntaRNA/PredictorSeedExtensionKinetic.cpp'
build = root/'kinetic-validation/release'
temp = root/'kinetic-validation/uncached-comparison'
temp.mkdir(exist_ok=True)
s = source.read_text()
old = 'buildCandidates(sides[best.left ? 0 : 1], bounds, best.left, last1, last2);'
assert s.count(old) == 1
s = s.replace(old, 'buildCandidates(sides[0], bounds, true, last1, last2);\n\t\tbuildCandidates(sides[1], bounds, false, last1, last2);')
(temp/'PredictorSeedExtensionKinetic.cpp').write_text(s)
recipe = subprocess.check_output(['make','-n','V=1','-W',str(source),'libIntaRNA_a-PredictorSeedExtensionKinetic.o'], cwd=build/'src/IntaRNA',text=True).splitlines()[0]
compile_command = shlex.split(recipe.split(' -MT ')[0]) + ['-c',str(temp/'PredictorSeedExtensionKinetic.cpp'),'-o',str(temp/'libIntaRNA_a-PredictorSeedExtensionKinetic.o')]
subprocess.run(compile_command,cwd=build/'src/IntaRNA',check=True)
shutil.copy2(build/'src/IntaRNA/libIntaRNA.a',temp/'libIntaRNA.a')
subprocess.run(['ar','r',str(temp/'libIntaRNA.a'),str(temp/'libIntaRNA_a-PredictorSeedExtensionKinetic.o')],check=True)
subprocess.run(['ranlib',str(temp/'libIntaRNA.a')],check=True)
recipe = subprocess.check_output(['make','-n','V=1','-W','IntaRNA-IntaRNA.o','IntaRNA'],cwd=build/'src/bin',text=True).splitlines()[-1]
link = shlex.split(recipe)
link[link.index('-o')+1] = str(temp/'IntaRNA')
link[link.index('../../src/IntaRNA/libIntaRNA.a')] = str(temp/'libIntaRNA.a')
subprocess.run(link,cwd=build/'src/bin',check=True)
print('Built uncached comparison with the same revised move/seed semantics and toolchain.')
