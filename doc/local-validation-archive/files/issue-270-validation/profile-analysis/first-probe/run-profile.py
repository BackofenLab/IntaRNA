from pathlib import Path
import subprocess
import shlex
import json
import hashlib
import time
import argparse
parser=argparse.ArgumentParser(description='Run one seeded M/H attribution sample on the OxyS/fhlA issue reproduction.')
parser.add_argument('--binary',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--reference',type=Path,help='optional directory with M-seed.stdout and H-seed.stdout')
args=parser.parse_args()
out=args.output.resolve()
binary=args.binary.resolve()
out.mkdir(parents=True,exist_ok=True)
query='GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGC'
target='AGUUAGUCAAUGACCUUUUGCACCGCUUUGCGGUGCUUUCCUGGAACAACAAAAUGUCAUAUACACCGAUGAGUGAUCUCGGACAACAAGGGUUGUUCGACAUCACUCGGACA'
commands=[]
for mode in ['M','H']:
    command=[str(binary),'--qId=OxyS','--query='+query,'--tId=fhlA','--target='+target,'--model=P','--mode='+mode,'--threads=1']
    commands.append(shlex.join(command))
    before=time.perf_counter()
    with (out/(mode+'-seed.stdout')).open('w') as stdout, (out/(mode+'-seed.stderr')).open('w') as stderr:
        subprocess.run(['/usr/bin/time','-v','-o',str(out/(mode+'-seed.time')),*command],stdout=stdout,stderr=stderr,check=True)
    print(mode, time.perf_counter()-before, flush=True)
    if args.reference:
        reference=(args.reference/(mode+'-seed.stdout')).read_bytes()
        actual=(out/(mode+'-seed.stdout')).read_bytes()
        print(mode,'stdout_exact_match',actual==reference,flush=True)
        if actual!=reference:
            raise RuntimeError(mode+' instrumented output differs from reference')
(out/'commands.txt').write_text('\n'.join(commands)+'\n')
(out/'binary-sha256.txt').write_text(hashlib.sha256(binary.read_bytes()).hexdigest()+'  '+str(binary)+'\n')
