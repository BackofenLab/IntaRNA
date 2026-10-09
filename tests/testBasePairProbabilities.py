"""Native CLI checks with an independent exhaustive base-pair-energy oracle."""
import csv
import gzip
import math
from pathlib import Path
import subprocess
import sys
import tempfile

BIN = sys.argv[1]
COMMON = ['--energy=B', '--acc=N', '--model=P', '--mode=M', '--seedBP=2',
          '--threads=1', '--outOverlap=B', '--outNoLP=false', '--outNoGUend=false',
          '--outMode=C', '--outCsvCols=Zall,Eall', '--default-log-file=/dev/null']

def run(t='GGAGGG', q='CCCC', options=(), success=True):
    keys = {o.split('=')[0] for o in options}
    if any(o.startswith('--outMode=') and o!='--outMode=C' for o in options):
        keys.add('--outCsvCols')
    common = [o for o in COMMON if o.split('=')[0] not in keys]
    p = subprocess.run([BIN, f'--target={t}', f'--query={q}', *common, *options],
                       text=True, capture_output=True)
    assert (p.returncode == 0) == success, (p.args, p.returncode, p.stdout, p.stderr)
    if not success:
        assert p.returncode>0 and 'ERROR' in p.stdout+p.stderr, (p.returncode,p.stdout,p.stderr)
    return p

def matrix(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        rows = list(csv.reader(f, delimiter=';'))
    assert rows[0][0] == 'bpProb', rows
    return rows[0][1:], [r[0] for r in rows[1:]], [[None if v=='NA' else float(v) for v in r[1:]] for r in rows[1:]]

def oracle(t, q, regions=None, no_lp=False, seeds=None, unseeded=False):
    # Both internal coordinates increase; only printing reverses query labels.
    n, m = len(t), len(q)
    qr = q[::-1]
    valid = {(i,j) for i in range(n) for j in range(m) if t[i]+qr[j] in ('GC','CG','AU','UA','GU','UG')}
    z = 0.0
    mass = [[0.0]*m for _ in t]
    if regions is None: regions = [(0,n-1,0,m-1)]
    def stack(a,b): return b==(a[0]+1,a[1]+1)
    def visit(path, end1, end2):
        nonlocal z
        admitted = any(stack(a,b) for a,b in zip(path,path[1:])) if seeds is None else any(
            any(path[k:k+len(s)]==s for k in range(len(path))) for s in seeds)
        lp_ok = all((k and stack(path[k-1],p)) or (k+1<len(path) and stack(p,path[k+1])) for k,p in enumerate(path))
        if (unseeded or admitted) and (not no_lp or lp_ok):
            weight = math.exp(len(path))
            z += weight
            for i,j in path: mass[i][m-1-j] += weight
        i,j = path[-1]
        for a in range(i+1,end1+1):
            for b in range(j+1,end2+1):
                if (a,b) in valid: visit(path+[(a,b)],end1,end2)
    for lo,hi,lq,hq in regions:
        for i in range(lo,hi+1):
            for j in range(lq,hq+1):
                if (i,j) in valid: visit([(i,j)],hi,hq)
    return z, [[x/z if z else None for x in row] for row in mass]

def same(a,b):
    assert len(a)==len(b)
    for x,y in zip(a,b):
        for v,w in zip(x,y):
            assert (v is None and w is None) or (v is not None and w is not None and math.isclose(v,w,rel_tol=2e-12,abs_tol=2e-14)), (v,w)

with tempfile.TemporaryDirectory(prefix='intarna-bpp-') as tmp:
    d=Path(tmp); path=d/'pairs.csv'
    def predict(t='GGAGGG',q='CCCC',extra=()):
        p=run(t,q,[f'--out=bpProb:{path}',*extra]);return p,matrix(path)[2]
    base, values=predict()
    z,expected=oracle('GGAGGG','CCCC');same(values,expected)
    assert math.isclose(float(base.stdout.splitlines()[1].split(';')[0]),z,rel_tol=1e-5)
    no_output=run();assert no_output.stdout==base.stdout
    for extra in (['--outNumber=0'],['--outNumber=10'],['--outMaxE=-100'],['--outDeltaE=0'],['--outPerRegion']):
        same(predict(extra=extra)[1],values)
    same(predict(extra=['--outNoLP'])[1],oracle('GGAGGG','CCCC',no_lp=True)[1])
    # Disjoint regions define a union, excluding chains crossing any boundary.
    extra=['--tRegion=1-3,4-6','--qRegion=1-2,3-4']
    regions=[(a,a+2,b,b+1) for a in (0,3) for b in (0,2)]
    regional=predict(extra=extra)[1];same(regional,oracle('GGAGGG','CCCC',regions)[1])
    same(predict(extra=extra+['--outPerRegion'])[1],regional)
    # Exact unseeded output uses the same original coordinates and region union.
    for flags, ranges in (([], None), (extra, regions)):
        for no_lp in (False, True):
            opts=['--noSeed', f'--outNoLP={str(no_lp).lower()}', *flags]
            predicted, probabilities=predict(extra=opts)
            same(probabilities, oracle('GGAGGG','CCCC',ranges,no_lp=no_lp,unseeded=True)[1])
            assert predicted.stdout==run(options=opts).stdout
    _,single_unseeded=predict('G','C',['--noSeed'])
    assert single_unseeded==[[1.0]]
    # Explicit singleton and mixed patterns; coordinates written in original order.
    seedopt=['--seedTQ=1|&4|,2||&2||']
    _,single=predict('GGGG','CCCC',seedopt)
    same(single,oracle('GGGG','CCCC',seeds=[[(0,0)],[(1,1),(2,2)]])[1])
    # Fixed admitted family makes energyAdd a representable common factor.
    tiny,tiny_values=predict('GGGG','CCCC',seedopt+['--energyAdd=100','--outMaxE=999'])
    same(tiny_values,single)
    ztext,etext=tiny.stdout.splitlines()[1].split(';');assert 0<float(ztext)<1e-19 and etext!='NA'
    diagnostics=run('GGGG','CCCC',seedopt+['--energyAdd=100','--outMaxE=999','--outCsvCols=Zall,Eall,EallTotal,P_E'])
    assert 'NA' not in diagnostics.stdout and float(diagnostics.stdout.splitlines()[1].split(';')[-1])>0
    ensemble=run('GGGG','CCCC',seedopt+['--energyAdd=100','--outMode=E'])
    assert 'Eall 0.00\n' not in ensemble.stdout
    # Explicit GU seed admission is separate from complete-interaction GU ends.
    _,gu=predict('GGG','CUC',['--seedTQ=2|&2|','--seedNoGU','--seedNoGUend','--outNoGUend'])
    same(gu,[[0,0,1],[0,1,0],[1,0,0]])
    # Labels use original sequence order and signed display shifts that skip zero.
    _,normal=predict('GAGGG','CCUCC')
    _,shifted=predict('GAGGG','CCUCC',['--tIdxPos0=-2','--qIdxPos0=5'])
    cols,rows,_=matrix(path);assert cols==['C_5','C_6','U_7','C_8','C_9'];assert rows==['G_-2','A_-1','G_1','G_2','G_3'];same(normal,shifted)
    _,empty=predict('AAAA','AAAA');assert all(v is None for row in empty for v in row)
    _,short=predict('G','C');assert short==[[None]]
    # Coverage and actual pairing differ at the unpaired A bulge.
    spot=d/'spot.csv';_,with_spot=predict(extra=[f'--out=spotProb:{spot}'])
    same(with_spot,values)
    covered=list(csv.reader(spot.open(),delimiter=';'))
    assert any(float(v)>0 for v in covered[3][1:]);assert all(v==0 for v in with_spot[2])
    compressed=d/'pairs.csv.gz';run(options=[f'--out=bpProb:{compressed}']);same(matrix(compressed)[2],values)
    # Unsupported requests fail explicitly and never create a probability file.
    for opts in (['--mode=H'],['--mode=S'],['--model=X','--mode=K'],['--windowWidth=3','--windowOverlap=2'],['--seedMaxUP=1'],['--seedTQ=1|.|&2||']):
        fail=d/'unsupported.csv';run(options=[f'--out=bpProb:{fail}',*opts],success=False);assert not fail.exists()
    fail=d/'range.csv';run('GGGG','CCCC',[f'--out=bpProb:{fail}',*seedopt,'--energyAdd=900'],success=False);assert not fail.exists()
    run(options=['--out=bpProb:/dev/full'],success=False)
    # Independent pairs and complete shared-stream blocks at both thread counts.
    fasta=d/'targets.fa';fasta.write_text('>a\nGGGG\n>b\nGAGG\n')
    for threads in (1,2):
        prefix=d/f'multi{threads}.csv'
        run(str(fasta),'CCCC',[f'--threads={threads}',f'--out=bpProb:{prefix}','--outNumber=0'])
        files=sorted(d.glob(f'multi{threads}*.csv'));assert len(files)==2,files
        for file,t in zip(files,('GGGG','GAGG')):same(matrix(file)[2],oracle(t,'CCCC')[1])
        shared=run(str(fasta),'CCCC',[f'--threads={threads}','--out=bpProb:STDOUT','--outNumber=0'])
        lines=shared.stdout.splitlines();starts=[i for i,l in enumerate(lines) if l.startswith('bpProb;')]
        assert len(starts)==2 and starts[1]-starts[0]==5,lines
print('base-pair probability CLI regressions passed')
