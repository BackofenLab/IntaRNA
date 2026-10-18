#!/usr/bin/env python3
"""Reproducible sequential IntaRNA benchmarks; no third-party Python dependencies."""
import argparse, csv, hashlib, json, os, platform, random, shutil, signal
import statistics, subprocess, time, threading
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / 'source/src/bin/IntaRNA'
INPUT = ROOT / 'inputs'
RESULT = ROOT / 'results'
ENV = dict(os.environ, LC_ALL='C', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
           OMP_PLACES='cores', OPENBLAS_NUM_THREADS='1')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            records.append([line[1:], ''])
        elif line.strip():
            records[-1][1] += line.strip().upper().replace('T', 'U')
    return records

def write_fasta(name, records):
    path = INPUT / (name + '.fa')
    path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in records))
    return str(path)

def prepare():
    INPUT.mkdir(exist_ok=True)
    docs = ROOT / 'source/doc/handson'
    for name in ['fhlA', 'OxyS', 'phoB', 'GcvB', 'GcvB.ST', 'ilvE', 'ChiX_NC_000913']:
        write_fasta(name, fasta(docs / (name + '.fasta')))
    write_fasta('promoters48', fasta(docs / 'NC_000913.fa')[:48])
    rng = random.Random(20260919)
    target = ''.join(rng.choices('ACGU', k=10000))
    query = ''.join(rng.choices('ACGU', k=800))
    for n in [40, 60, 100, 300, 1000, 3000, 10000]:
        write_fasta(f't{n}', [(f'synthetic_target_{n}', target[:n])])
    for n in [40, 50, 60, 100, 200, 400, 800]:
        write_fasta(f'q{n}', [(f'synthetic_query_{n}', query[:n])])
    for gc in [20, 50, 80]:
        weights = [(100-gc)/2, gc/2, gc/2, (100-gc)/2]
        write_fasta(f'gc{gc}t', [(f'gc{gc}target', ''.join(rng.choices('ACGU', weights=weights, k=300)))])
        write_fasta(f'gc{gc}q', [(f'gc{gc}query', ''.join(rng.choices('ACGU', weights=weights, k=100)))])
    for alphabet, name in [('AU', 'au_repeat'), ('GC', 'gc_repeat'), ('A', 'no_pair')]:
        write_fasta(name+'t', [(name+'target', (alphabet*300)[:300])])
        write_fasta(name+'q', [(name+'query', (alphabet*100)[:100])])
    cases = []
    def add(name, category, t, q, extra=(), threads=1):
        args = ['--target='+str(INPUT/(t+'.fa')), '--query='+str(INPUT/(q+'.fa')),
                '--outMode=C', '--outCsvCols=id1,id2,start1,end1,start2,end2,E',
                '--threads='+str(threads)] + list(extra)
        cases.append(dict(name=name, category=category, args=args, threads=threads))
    add('bio_fhla_oxys', 'biological', 'fhlA', 'OxyS')
    add('bio_phob_gcvb', 'biological', 'phoB', 'GcvB')
    add('bio_ilve_gcvb', 'biological', 'ilvE', 'GcvB.ST')
    for n in [100, 300, 1000, 3000, 10000]:
        add(f'length_t{n}_q100', 'target_length', f't{n}', 'q100')
    for n in [50, 200, 400, 800]:
        add(f'length_t300_q{n}', 'query_length', 't300', f'q{n}')
    for gc in [20, 50, 80]:
        add(f'composition_gc{gc}', 'composition', f'gc{gc}t', f'gc{gc}q')
    for name in ['au_repeat', 'gc_repeat', 'no_pair']:
        add(name, 'composition', name+'t', name+'q')
    for name, args in [
        ('default', []), ('exact_x', ['--mode=M']), ('heuristic_s', ['--model=S']),
        ('exact_s', ['--mode=M', '--model=S']), ('helix', ['--model=B']),
        ('no_seed', ['--noSeed']), ('gapped_seed', ['--seedMaxUP=2']),
        ('seed_only', ['--mode=S']), ('no_lonely_pairs', ['--outNoLP']),
        ('loop0', ['--intLoopMax=0']), ('loop30', ['--intLoopMax=30']),
        ('out10', ['--outNumber=10']), ('out100', ['--outNumber=100'])]:
        add('algorithm_'+name, 'algorithm', 't100', 'q100', args)
    for mode in ['H', 'M']:
        add('ensemble_'+mode, 'ensemble', 't40', 'q40', ['--model=P', '--mode='+mode, '--noSeed'])
    for name, args in [
        ('computed', ['--intLenMax=150']), ('none', ['--intLenMax=150', '--acc=N']),
        ('global', ['--intLenMax=150', '--accW=0', '--accL=0']),
        ('window75', ['--intLenMax=75', '--accW=75', '--accL=75']),
        ('cap75_control', ['--intLenMax=75'])]:
        add('accessibility_'+name, 'accessibility', 't1000', 'q100', args)
    for width in [300, 600]:
        add(f'window_{width}', 'prediction_windows', 't3000', 'q100',
            ['--intLenMax=150', '--windowWidth='+str(width), '--windowOverlap=150'])
    for threads in [1, 2, 4, 6, 12]:
        add(f'throughput_{threads}t', 'thread_scaling', 'promoters48', 'ChiX_NC_000913', threads=threads)
    add('memory_t3000_q800', 'memory_stress', 't3000', 'q800')
    add('memory_window300', 'memory_stress', 't3000', 'q800',
        ['--intLenMax=150', '--windowWidth=300', '--windowOverlap=150'])
    (ROOT/'cases.json').write_text(json.dumps(cases, indent=2)+'\n')
    (ROOT/'inputs/manifest.json').write_text(json.dumps({p.name: dict(sha256=sha(p),
        records=[dict(id=i, length=len(s), gc_fraction=round((s.count('G')+s.count('C'))/len(s),5)) for i,s in fasta(p)])
        for p in sorted(INPUT.glob('*.fa'))}, indent=2)+'\n')
    print(f'Prepared {len(cases)} cases', flush=True)

def canonical(data):
    lines = [line for line in data.splitlines() if not line.startswith(b'#')]
    return b'\n'.join(lines[:1]+sorted(lines[1:])) + b'\n'

def run(case, phase, repetition, timeout):
    dest = RESULT / f"{case['name']}.{phase}{repetition}"
    if Path(str(dest)+'.stdout').exists():
        raise RuntimeError('Refusing to overwrite prior evidence: '+str(dest))
    cmd = ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M\t%F\t%R\t%w\t%c\t%I\t%O\t%x',
           '-o', str(dest)+'.time', str(BIN)] + case['args']
    started_utc=datetime.now(timezone.utc).isoformat()
    binary_hash=sha(BIN)
    run_env=ENV.copy()
    for key,value in case.get('env_overrides',{}).items():
        if value is None:
            run_env.pop(key,None)
        else:
            run_env[key]=value
    start = time.perf_counter()
    with open(str(dest)+'.stdout', 'wb') as out, open(str(dest)+'.stderr', 'wb') as err:
        proc = subprocess.Popen(cmd, stdout=out, stderr=err, env=run_env, start_new_session=True)
        timed_out = False
        def expire():
            nonlocal timed_out
            if proc.poll() is None:
                timed_out = True
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        watchdog = threading.Timer(timeout, expire)
        watchdog.start()
        rc = proc.wait()  # blocking wait avoids timeout polling's ~50 ms timing bias
        watchdog.cancel()
        if timed_out:
            rc = 124
    elapsed = time.perf_counter()-start
    csv_rows=list(csv.reader([line for line in Path(str(dest)+'.stdout').read_text().splitlines()
                             if not line.startswith('#')],delimiter=';'))
    csv_valid=bool(csv_rows) and csv_rows[0]==['id1','id2','start1','end1','start2','end2','E'] and all(len(r)==7 for r in csv_rows)
    vals = Path(str(dest)+'.time').read_text().strip().splitlines()
    row = dict(case=case['name'], category=case['category'], phase=phase,
               repetition=repetition, elapsed_s=elapsed, exit_code=rc, timeout=timed_out,
               started_utc=started_utc, ended_utc=datetime.now(timezone.utc).isoformat(),
               binary_sha256=binary_hash, csv_valid=csv_valid,
               env_overrides=case.get('env_overrides',{}),
               command=cmd, load_average=os.getloadavg(),
               output_sha256=sha(str(dest)+'.stdout'),
               canonical_sha256=hashlib.sha256(canonical(Path(str(dest)+'.stdout').read_bytes())).hexdigest(),
               output_rows=max(0,len(csv_rows)-1))
    if vals and len(vals[-1].split('\t')) == 11:
        keys=['wall_s','user_s','system_s','max_rss_kib','major_faults','minor_faults',
              'voluntary_context_switches','involuntary_context_switches','fs_inputs','fs_outputs','time_exit_code']
        row.update(dict(zip(keys, map(float, vals[-1].split('\t')))))
    with (RESULT/'runs.jsonl').open('a') as f:
        f.write(json.dumps(row)+'\n')
    print(f"{phase} {repetition} {case['name']}: {elapsed:.3f}s rc={rc} RSS={row.get('max_rss_kib', 0)/1024:.1f}MiB", flush=True)
    return row

def summarize():
    rows = [json.loads(line) for line in (RESULT/'runs.jsonl').read_text().splitlines()]
    # Re-derive schema and equality from preserved stdout, ignoring IntaRNA's # INFO lines.
    # Original run records remain unchanged as evidence of collection-time parsing.
    for row in rows:
        path=RESULT/f"{row['case']}.{row['phase']}{row['repetition']}.stdout"
        data=canonical(path.read_bytes())
        cells=list(csv.reader(data.decode().splitlines(),delimiter=';'))
        row['canonical_sha256']=hashlib.sha256(data).hexdigest()
        row['csv_valid']=bool(cells) and cells[0]==['id1','id2','start1','end1','start2','end2','E'] and all(len(r)==7 for r in cells)
        row['output_rows']=max(0,len(cells)-1)
    (RESULT/'validated-runs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    summary = []
    for name in dict.fromkeys(r['case'] for r in rows):
        rr = [r for r in rows if r['case']==name and r['phase']=='measure' and r['exit_code']==0]
        if not rr:
            continue
        times=[r['elapsed_s'] for r in rr]
        summary.append(dict(case=name, category=rr[0]['category'], n=len(rr),
            median_s=statistics.median(times), min_s=min(times), max_s=max(times),
            median_user_s=statistics.median(r['user_s'] for r in rr),
            median_system_s=statistics.median(r['system_s'] for r in rr),
            median_rss_mib=statistics.median(r['max_rss_kib']/1024 for r in rr),
            peak_rss_mib=max(r['max_rss_kib']/1024 for r in rr),
            stable_output=len(set(r['canonical_sha256'] for r in rr))==1,
            rows=rr[0]['output_rows'], canonical_sha256=rr[0]['canonical_sha256']))
    (RESULT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (RESULT/'failures.json').write_text(json.dumps([r for r in rows if r['exit_code'] or not r['csv_valid']],indent=2)+'\n')
    with (RESULT/'summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare','run','summarize'])
    parser.add_argument('--only', default='')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=90)
    args=parser.parse_args()
    RESULT.mkdir(exist_ok=True)
    if args.action=='prepare':
        prepare(); return
    if args.action=='summarize':
        summarize(); return
    cases=json.loads((ROOT/'cases.json').read_text())
    if args.only:
        names=set(args.only.split(',')); cases=[c for c in cases if c['name'] in names]
    eligible=[]
    for case in cases:
        if run(case,'warmup',0,args.timeout)['exit_code']==0:
            eligible.append(case)
    rng=random.Random(7312026)
    for repeat in range(1,args.repetitions+1):
        order=eligible[:]; rng.shuffle(order)
        for case in order:
            run(case,'measure',repeat,args.timeout)
    summarize()

if __name__=='__main__':
    main()
