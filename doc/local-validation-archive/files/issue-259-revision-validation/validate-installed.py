from pathlib import Path
import os,shlex,subprocess
root=Path('/home/alex/CodexProjects/IntaRNA')
repo=root/'pr-worktrees/issue-258-bpsvg'
work=root/'issue-259-revision-validation'
prefix=root/'issue-258-validation/install-release'
deps=root/'issue-246-validation/toolchain'
cxx=root/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'
common=[str(cxx),'-std=c++23','-fopenmp','-I'+str(prefix/'include'),'-I'+str(deps/'include')]
probe=work/'header-probe.cpp'
headers=sorted((prefix/'include/IntaRNA').glob('*.h'))
for header in ([] if os.environ.get("SKIP_HEADERS") else headers):
    probe.write_text('#include <IntaRNA/'+header.name+'>\nint main() {return 0;}\n')
    subprocess.run(common+['-fsyntax-only',str(probe)],check=True)
print(f'{len(headers)} installed headers passed',flush=True)
env=dict(os.environ,PKG_CONFIG_PATH=str(prefix/'lib/pkgconfig'))
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','IntaRNA'],env=env,text=True))
original=(repo/'doc/seeded-probabilities-example.cpp').read_text()
for name in ['PredictorMfeEns2dSeedExtension','PredictorMfeEns2d','PredictorMfeEns2dHeuristic','PredictorMfeEns2dHeuristicSeedExtension']:
    text=original.replace('PredictorMfeEns2dSeedExtension',name)
    if not name.endswith('SeedExtension'):
        text=text.replace('new SeedHandlerNoBulge(energy, seeds), &probabilities','&probabilities')
    text=text.replace('probabilities.finalize();','probabilities.finalize();\n\tif(probabilities.isApproximate() != '+str('Heuristic' in name).lower()+') return 2;')
    text=text.replace('BasePairProbabilityWriter::write(std::cout, probabilities, target, query);','BasePairProbabilityWriter::write(std::cout, probabilities, target, query);\n\tBasePairProbabilityWriter::writeSvg(std::cerr, probabilities, energy, 1, 1);')
    source=work/(name+'-consumer.cpp'); source.write_text(text)
    binary=work/(name+'-consumer')
    subprocess.run(common+[str(source),'-L'+str(deps/'lib'),'-Wl,-rpath,'+str(deps/'lib'),*flags,'-o',str(binary)],check=True)
    result=subprocess.run([str(binary)],capture_output=True,check=True)
    assert result.stdout.startswith(b'bpProbApproximate;' if 'Heuristic' in name else b'bpProb;')
    assert b'<svg' in result.stderr
    print(name+' installed consumer passed',flush=True)
