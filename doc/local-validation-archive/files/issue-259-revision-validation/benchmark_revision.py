#!/usr/bin/env python3
"""Compare existing prediction modes and measure new outside/output overhead."""
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import time
import xml.etree.ElementTree as ET

ROOT = Path('/home/alex/CodexProjects/IntaRNA')
REPO = ROOT/'pr-worktrees/issue-258-bpsvg'
HERE = ROOT/'issue-259-revision-validation'
OUT = HERE/'performance'
OUT.mkdir(exist_ok=True)
shutil.copy2(REPO/'src/bin/IntaRNA', OUT/'IntaRNA-final')
BINARIES = {'before': HERE/'IntaRNA-before', 'final': OUT/'IntaRNA-final'}
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_DYNAMIC='FALSE')
CPU = min(os.sched_getaffinity(0))
COMMON = ['--model=P', '--seedBP=4', '--intLenMax=30', '--qIntLoopMax=3', '--tIntLoopMax=3', '--threads=1', '--outNoGUend=false', '--outNumber=0', '--outMode=E', '--default-log-file=/dev/null', '--out=STDOUT']
PREDICTORS = {'noSeed-M': ['--mode=M', '--noSeed'], 'noSeed-H': ['--mode=H', '--noSeed'], 'seeded-H': ['--mode=H']}
CASES = {'OxyS-fhlA': [f'--target={REPO}/doc/handson/fhlA.fasta', f'--query={REPO}/doc/handson/OxyS.fasta', '--outNoLP=false']}
for name, length, no_lp in [('dense40', 40, False), ('dense32-noLP', 32, True)]:
    for strand, nucleotide in [('target', 'G'), ('query', 'C')]:
        (OUT/f'{name}-{strand}.fa').write_text(f'>{strand}\n{nucleotide*length}\n')
    CASES[name] = [f'--target={OUT}/{name}-target.fa', f'--query={OUT}/{name}-query.fa', '--tAcc=N', '--qAcc=N', f'--outNoLP={str(no_lp).lower()}']
CONFIGS = [('before', 'none'), ('final', 'none'), ('final', 'csv'), ('final', 'svg')]
records, checks, commands = [], [], {}
randomizer = random.Random(259)


def command(case, predictor, version, output, check=False):
    cmd = ['taskset', '-c', str(CPU), str(BINARIES[version])] + COMMON + PREDICTORS[predictor] + CASES[case]
    if output != 'none':
        prefix, filename = ('bpProb', 'pairs.csv') if output == 'csv' else ('bpsvg', 'plot.svg')
        destination = str(OUT/filename) if check else '/dev/null'
        cmd.append(f'--out={prefix}:{destination}')
    return cmd


def run(case, predictor, version, output, repetition):
    check = repetition < 0
    cmd = command(case, predictor, version, output, check)
    if not check:
        commands[f'{case}/{predictor}/{version}/{output}'] = cmd
    before = time.perf_counter()
    completed = subprocess.run(['/usr/bin/time', '-f', '%e,%U,%S,%M', '-o', str(OUT/'time.txt')] + cmd,
                               stdout=subprocess.PIPE if check else subprocess.DEVNULL, stderr=subprocess.PIPE, env=ENV)
    wall = time.perf_counter()-before
    if completed.returncode:
        raise RuntimeError(f'{case}/{predictor}/{version}/{output}: {completed.stderr.decode()}')
    elapsed, user, system, rss = map(float, (OUT/'time.txt').read_text().strip().split(','))
    if not check:
        records.append(dict(case=case, predictor=predictor, version=version, output=output, repetition=repetition,
                            wall_seconds=wall, user_seconds=user, system_seconds=system, cpu_seconds=user+system, rss_kib=int(rss)))
        (OUT/'raw.json').write_text(json.dumps(records, indent=2)+'\n')
    return completed.stdout


def check_plot(csv_data):
    rows = [line.split(';') for line in csv_data.decode().splitlines()]
    expected = {(row[0].rsplit('_', 1)[1], query.rsplit('_', 1)[1]): value for row in rows[1:] for query, value in zip(rows[0][1:], row[1:])}
    count = 0
    for _, element in ET.iterparse(OUT/'plot.svg', events=['end']):
        if element.get('data-type') == 'base-pair':
            pair = (element.get('data-target'), element.get('data-query'))
            actual, value = element.get('data-probability'), expected[pair]
            assert (actual == value == 'NA') or float(actual) == float(value), (pair, actual, value)
            count += 1
        element.clear()
    assert count == len(expected), (count, len(expected))
    return count


for case in CASES:
    for predictor in PREDICTORS:
        print(f'Checking {case}/{predictor}', flush=True)
        reference = None
        csv_data = None
        for version, output in CONFIGS:
            ensemble = run(case, predictor, version, output, -1)
            if reference is None:
                reference = ensemble
            assert ensemble == reference, (case, predictor, version, output, ensemble, reference)
            if output == 'csv':
                csv_data = (OUT/'pairs.csv').read_bytes()
                assert csv_data.startswith(b'bpProbApproximate' if predictor.endswith('H') else b'bpProb;')
            if output == 'svg':
                count = check_plot(csv_data)
                checks.append(dict(case=case, predictor=predictor, ensemble_sha256=hashlib.sha256(reference).hexdigest(),
                                   ensemble=reference.decode(), cells=count, csv_bytes=len(csv_data), svg_bytes=(OUT/'plot.svg').stat().st_size))
        for repetition in range(7):
            order = CONFIGS.copy()
            randomizer.shuffle(order)
            for version, output in order:
                run(case, predictor, version, output, repetition)
        (OUT/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
        (OUT/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
        print(f'Completed {case}/{predictor}', flush=True)

summary = []
for case in CASES:
    for predictor in PREDICTORS:
        for version, output in CONFIGS:
            group = [row for row in records if (row['case'], row['predictor'], row['version'], row['output']) == (case, predictor, version, output)]
            summary.append(dict(case=case, predictor=predictor, version=version, output=output,
                                wall_median=statistics.median(row['wall_seconds'] for row in group),
                                wall_min=min(row['wall_seconds'] for row in group), wall_max=max(row['wall_seconds'] for row in group),
                                cpu_median=statistics.median(row['cpu_seconds'] for row in group), rss_max_kib=max(row['rss_kib'] for row in group)))
(OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
metadata = dict(cpu=CPU, repeats=7, binaries={name: dict(path=str(binary), sha256=hashlib.sha256(binary.read_bytes()).hexdigest()) for name, binary in BINARIES.items()},
                compiler='Same configured GCC 14 release toolchain; retained pre-revision executable versus final worktree executable.',
                timing='Fresh processes; seven deterministically interleaved repetitions; one thread; CSV/SVG and ensemble sent to /dev/null for timings. Separate file-output warmups compare ensemble bytes and every CSV/SVG cell.')
(OUT/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
print(json.dumps(summary, indent=2), flush=True)
