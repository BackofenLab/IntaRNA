#!/usr/bin/env bash
# Per-seed selection across all dispatched modes and CLI merge/validation paths.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
python3 - "$bin" "$tmp" <<'PY'
import csv
import io
from pathlib import Path
import subprocess
import sys

binary, tmp = sys.argv[1:]
common = [binary, '--energy=B', '--acc=N', '--seedBP=2', '--threads=1',
          '--outMode=C', '--outCsvCols=start1,end1,start2,end2,E',
          '--outNumber=1000', '--outMaxE=100', '--outDeltaE=100',
          '--default-log-file=' + tmp + '/info.log']

def run(*args, ok=True):
    overrides = {arg.split("=", 1)[0] for arg in args if arg.startswith("--")}
    command = [arg for arg in common if arg.split("=", 1)[0] not in overrides] + list(args)
    if any(arg.startswith("--outMode=") and arg != "--outMode=C" for arg in args):
        command = [arg for arg in command if not arg.startswith("--outCsvCols=")]
    process = subprocess.run(command, text=True, capture_output=True)
    assert (process.returncode == 0) == ok, (args, process.stdout, process.stderr)
    return process.stdout

def rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter=';'))

for model, modes in [('S', 'HMS'), ('X', 'HMRKS'), ('B', 'H'), ('P', 'HMS')]:
    for mode in modes:
        args = ['-t', 'GGCGGG', '-q', 'CCCGCC', '--model=' + model, '--mode=' + mode,
                '--helixMinBP=2', '--helixMaxBP=4']
        default = run(*args)
        assert run(*args, '--outNumberPerSeed=0') == default, (model, mode, 'zero changed output')
        one = rows(run(*args, '--outNumberPerSeed=1'))
        two = rows(run(*args, '--outNumberPerSeed=2'))
        assert one and two, (model, mode, one, two)
        assert len(one) <= len(two), (model, mode)
        assert len(rows(run(*args, '--outNumberPerSeed=2', '--outNumber=1'))) == 1
        assert not rows(run(*args, '--outNumberPerSeed=1', '--outNumber=0'))
        # Quotas never change the ordinary partition/ensemble statistics.
        if mode != 'K':
            assert run(*args, '--outMode=E') == run(*args, '--outMode=E', '--outNumberPerSeed=1'), (model, mode, 'partition changed')

args = ['-t', 'GGCGGG', '-q', 'CCCGCC', '--model=X', '--mode=M']
for invalid in ['-1', '1001', '2147483648', 'abc']:
    run(*args, '--outNumberPerSeed=' + invalid, ok=False)
run(*args, '--outNumberPerSeed=1', '--noSeed', ok=False)
assert run(*args, '--outNumberPerSeed=0', '--noSeed')
# Evaluation ignores prediction options before numeric validation.
run('-t', 'GG', '-q', 'CC', '--rri=1||&1||', '--outNumberPerSeed=-1')
parameter_file = Path(tmp) / 'parameters'
parameter_file.write_text('outNumberPerSeed=2\n')
assert run(*args, '--parameterFile=' + str(parameter_file)) == run(*args, '--outNumberPerSeed=2')
assert not rows(run('-t', 'AAAA', '-q', 'AAAA', '--outNumberPerSeed=1'))

regions = ['-t', 'GGGGAAAAGGGG', '-q', 'CCCCAAAACCCC', '--model=X', '--mode=M',
           '--tRegion=1-4,9-12', '--qRegion=1-4,9-12', '--outNumberPerSeed=1']
for overlap in 'NTQB':
    merged = rows(run(*regions, '--outOverlap=' + overlap))
    separate = rows(run(*regions, '--outOverlap=' + overlap, '--outPerRegion=true'))
    assert merged and separate
    if overlap == 'N':
        assert len(merged) == 2 and len(separate) == 4, (merged, separate)

# Overlapping windows share seed identities in original coordinates, including
# reversed query offsets. Stable results must not depend on worker count.
windows = ['-t', 'GGCGGGCGGGCGGGCG', '-q', 'CGCCCGCCCGCCCGCC', '--model=X', '--mode=M',
           '--tIntLenMax=4', '--qIntLenMax=4', '--outNumberPerSeed=2']
full = run(*windows)
for threads in ['1', '2']:
    assert run(*windows, '--windowWidth=10', '--windowOverlap=4', '--threads=' + threads) == full

# Prediction trackers observe the same ordinary candidates exactly once.
for model in 'SXP':
    tracker = ['-t', 'GGCGGG', '-q', 'CCCGCC', '--model=' + model, '--mode=M',
               '--out=tMinE:' + tmp + '/min.csv']
    run(*tracker)
    expected = (Path(tmp) / 'min.csv').read_bytes()
    run(*tracker, '--outNumberPerSeed=1')
    assert (Path(tmp) / 'min.csv').read_bytes() == expected

# Probability outputs retain the same raw ensemble when reporting is limited.
for mode in 'MH':
    probability = ['-t', 'GGCGGG', '-q', 'CCCGCC', '--model=P', '--mode=' + mode,
                   '--out=bpProb:' + tmp + '/pairs.csv']
    run(*probability)
    expected = (Path(tmp) / 'pairs.csv').read_bytes()
    run(*probability, '--outNumberPerSeed=1')
    assert (Path(tmp) / 'pairs.csv').read_bytes() == expected
print('Per-seed mode matrix, validation, region/window merges and probabilities passed')
PY
