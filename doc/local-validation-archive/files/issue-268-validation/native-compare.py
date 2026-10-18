from pathlib import Path
import json, os, subprocess
root=Path('/home/alex/CodexProjects/IntaRNA')
out=root/'issue-268-validation/native-comparison'
out.mkdir(exist_ok=True)
old=json.loads((root/'issue-266-validation/native-comparison/commands.json').read_text())
commands=[]
for command in old:
    variant='base' if 'issue-264' in command[0] else 'candidate'
    binary=root/'issue-268-validation'/('baseline/IntaRNA' if variant=='base' else 'release/src/bin/IntaRNA')
    command=[str(binary),*[a.replace('issue-266-validation/native-comparison','issue-268-validation/native-comparison') for a in command[1:]]]
    commands.append(command)
(out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
for command in commands:
    variant='base' if 'baseline/IntaRNA' in command[0] else 'candidate'
    if os.environ.get('VARIANT') and os.environ['VARIANT']!=variant: continue
    name=next(Path(a.split(':',1)[1]).stem for a in command if a.startswith('--out=bpProb:'))
    result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
    (out/(name+'.stdout')).write_text(result.stdout)
    (out/(name+'.stderr')).write_text(result.stderr)
    print(name, 'completed',flush=True)
if not os.environ.get('VARIANT'):
    for before in sorted(out.glob('*-base.csv')):
        after=before.with_name(before.name.replace('-base.csv','-candidate.csv'))
        assert before.read_bytes()==after.read_bytes(),(before,after)
        stdout=before.with_suffix('.stdout')
        assert stdout.read_bytes()==after.with_suffix('.stdout').read_bytes(),stdout
    print('All 12 native CLI stdout and pair CSV comparisons byte-identical')
