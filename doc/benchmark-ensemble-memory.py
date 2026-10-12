#!/usr/bin/env python3
"""Compare fresh-process peak RSS/runtime and prediction bytes for issue #260."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import random
import statistics
import subprocess


def canonical(output):
    rows = list(csv.reader(io.StringIO(output), delimiter=';'))
    return [rows[0], *sorted(rows[1:])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'feat-before', 'feat-after'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--benchmark-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--cases', nargs='+', help='Optional subset of named fixtures')
    parser.add_argument('--models', nargs='+', choices=('X', 'P'), default=['X', 'P'])
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    fixtures = args.outdir/'fixtures'
    fixtures.mkdir(exist_ok=True)
    records = (args.benchmark_dir/'input/NC_000913/target/NC_000913.fa').read_text().split('>')[1:]
    targets = {}
    for count in (1, 8, 32):
        targets[count] = fixtures/f'ecoli-{count}.fa'
        targets[count].write_text('>'+'>'.join(records[:count]))
    target, query = fixtures/'dense-target.fa', fixtures/'dense-query.fa'
    target.write_text('>dense-target\n'+'G'*256+'\n')
    query.write_text('>dense-query\n'+'C'*24+'\n')
    query_dir = args.benchmark_dir/'input/NC_000913/query'
    cases = [
        ('OxyS-one', targets[1], query_dir/'OxyS_NC_000913.fasta', 1, False, []),
        ('OxyS-eight', targets[8], query_dir/'OxyS_NC_000913.fasta', 4, False, []),
        ('GcvB-32', targets[32], query_dir/'GcvB_NC_000913.fasta', 4, False, []),
        ('GcvB-32-noLP', targets[32], query_dir/'GcvB_NC_000913.fasta', 4, True, []),
        ('dense256x24', target, query, 1, False,
         ['--energy=B', '--acc=N', '--seedBP=4', '--intLoopMax=0']),
    ]
    binaries = {name: getattr(args, name.replace('-', '_')).resolve()
                for name in ('before', 'after', 'feat-before', 'feat-after')}
    raw, summary, commands = [], [], {}
    for case, target, query, threads, noLP, extra in cases:
        if args.cases and case not in args.cases:
            continue
        for model in args.models:
            options = ['--target='+str(target.resolve()), '--query='+str(query.resolve()),
                       '--threads='+str(threads), '--model='+model, '--mode=H',
                       '--outNoLP='+str(noLP).lower(), '--tIntLenMax=80',
                       '--outMode=C', '--outNumber=1', '--default-log-file=/dev/null',
                       '--outCsvCols=id1,id2,start1,end1,start2,end2,E,Zall', *extra]
            expected = None
            for version, binary in binaries.items():
                command = [str(binary), *options]
                commands[f'{case}/{model}/{version}'] = command
                output = subprocess.run(command, check=True, capture_output=True, text=True, timeout=300).stdout
                value = canonical(output)
                if expected is None:
                    expected = value
                assert value == expected, (case, model, version, value, expected)
                (args.outdir/f'{case}-{model}-{version}.csv').write_text(output)
            for repeat in range(args.repeats):
                order = list(binaries)
                random.Random(260+repeat).shuffle(order)
                for version in order:
                    timing = args.outdir/'time.txt'
                    result = subprocess.run(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(timing),
                                             *commands[f'{case}/{model}/{version}']],
                                            check=True, capture_output=True, text=True, timeout=300)
                    if canonical(result.stdout) != expected:
                        (args.outdir/'mismatch.csv').write_text(result.stdout)
                        (args.outdir/'mismatch.stderr').write_text(result.stderr)
                        (args.outdir/'mismatch.json').write_text(json.dumps(dict(case=case, model=model,
                            version=version, repeat=repeat, expected=expected,
                            observed=canonical(result.stdout)), indent=2)+'\n')
                        raise AssertionError((case, model, version, repeat))
                    wall, user, system, rss = map(float, timing.read_text().split())
                    raw.append(dict(case=case, model=model, version=version, repeat=repeat,
                                    wall=wall, user=user, system=system, rss_kib=rss))
            for version in binaries:
                runs = [r for r in raw if (r['case'], r['model'], r['version']) == (case, model, version)]
                summary.append(dict(case=case, model=model, version=version,
                                    wall_median=statistics.median(r['wall'] for r in runs),
                                    wall_min=min(r['wall'] for r in runs), wall_max=max(r['wall'] for r in runs),
                                    rss_median_kib=statistics.median(r['rss_kib'] for r in runs)))
            for filename, data in [('raw.json', raw), ('summary.json', summary), ('commands.json', commands)]:
                (args.outdir/filename).write_text(json.dumps(data, indent=2)+'\n')
            print('verified', case, model, flush=True)
    metadata = dict(repeats=args.repeats, binaries={name: dict(path=str(path),
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name, path in binaries.items()},
                    methodology='Fresh processes, one correctness warmup; deterministically interleaved repetitions. '
                    'GNU time wall/user/system seconds and peak RSS KiB. Compare all CSV rows after sorting by contents '
                    'to account for thread completion order. Both model X and P are compared across all four binaries.')
    (args.outdir/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
