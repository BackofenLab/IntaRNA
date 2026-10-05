#!/usr/bin/env python3
"""Reproducible end-to-end kinetic benchmark; emits timings and output hashes.

Use an otherwise identical binary with both ends rebuilt after every move as
--uncached to isolate candidate reuse. All runs are single-threaded. Timings
include startup, accessibility, seed generation, table setup and prediction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cached', required=True, type=Path)
    parser.add_argument('--uncached', type=Path)
    parser.add_argument('--repeat', type=int, default=5)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error('--repeat must be positive')
    rng = random.Random(254)
    target = ''.join(rng.choices('ACGU', k=600))
    query = ''.join(rng.choices('ACGU', k=80))
    cases = [('random-no-ED', target, query, ['--acc=N']),
             ('random-folded', target, query, ['--acc=C', '--accW=150', '--accL=100']),
             ('stack-rich', 'G'*100, 'C'*30, ['--acc=N'])]
    variants = [('cached-K', str(args.cached.resolve()), 'K'),
                ('pruned-L', str(args.cached.resolve()), 'L')]
    if args.uncached:
        variants.insert(0, ('uncached-K', str(args.uncached.resolve()), 'K'))
    results = []
    with tempfile.TemporaryDirectory(prefix='intarna-kinetic-bench-') as temp:
        for case, t, q, extra in cases:
            common = ['--target='+t, '--query='+q, '--model=X', '--seedBP=7',
                      '--intLenMax=60', '--intLoopMax=10', '--threads=1',
                      '--outNoLP', '--outNumber=10', '--outMode=C',
                      '--outCsvCols=hybridDB,E', '--default-log-file=/dev/null'] + extra
            measurements = {label: [] for label, _, _ in variants}
            outputs = {}
            # Warm up every variant, then rotate run order across repetitions.
            for iteration in range(args.repeat+1):
                order = variants[iteration % len(variants):] + variants[:iteration % len(variants)]
                for label, binary, mode in order:
                    stats = Path(temp)/'time.txt'
                    command = [binary, '--mode='+mode] + common
                    start = time.perf_counter()
                    run = subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(stats)] + command,
                                         check=True, capture_output=True)
                    elapsed = time.perf_counter()-start
                    print(f"{case}/{label} run {iteration}: {elapsed:.3f}s", file=sys.stderr, flush=True)
                    if label in outputs and outputs[label] != run.stdout:
                        raise RuntimeError(f'Non-deterministic output: {case}/{label}')
                    outputs[label] = run.stdout
                    if iteration:
                        measurements[label].append((elapsed, int(stats.read_text().strip())))
            if 'uncached-K' in outputs and outputs['uncached-K'] != outputs['cached-K']:
                raise RuntimeError(f'Candidate caching changed the result: {case}')
            for label, _, _ in variants:
                samples = measurements[label]
                results.append(dict(case=case, variant=label,
                                    seconds_median=statistics.median(x[0] for x in samples),
                                    rss_KiB_max=max(x[1] for x in samples),
                                    samples_seconds=[x[0] for x in samples],
                                    reported_rows=max(0, outputs[label].count(b'\n')-1),
                                    sha256=hashlib.sha256(outputs[label]).hexdigest(),
                                    equals_K=outputs[label] == outputs['cached-K']))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
