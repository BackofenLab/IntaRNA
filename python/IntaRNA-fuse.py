#!/usr/bin/env python3
"""Fuse compatible IntaRNA interaction sites using joint accessibility costs.

Python replacement for IntaRNA-fuse.pl by Martin Raden. Requires IntaRNA and
ViennaRNA's RNAfold; all energies are in kcal/mol and coordinates are 1-based.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import io
from pathlib import Path
import re
import shutil
import subprocess
import sys


ESSENTIAL_COLUMNS = ('E', 'start1', 'end1', 'start2', 'end2', 'ED1', 'ED2', 'E_init')
DEFAULT_COLUMNS = ESSENTIAL_COLUMNS + ('hybridDB',)


def sequence(value):
    """Accept one literal sequence, as required for RNAfold constraints."""
    if not re.fullmatch('[ACGUTNacgutn]+', value):
        raise argparse.ArgumentTypeError('expected a nonempty sequence of A, C, G, U, T or N')
    return value


def positive_integer(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('expected a positive integer')
    return number


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Report non-overlapping IntaRNA sites and compatible fused interactions.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', required=True, type=sequence, help='target sequence (=seq1)')
    parser.add_argument('-q', required=True, type=sequence, help='query sequence (=seq2)')
    parser.add_argument('-n', type=positive_integer, default=10,
                        help='maximum number of predictions, including the MFE')
    parser.add_argument('-l', type=positive_integer, default=100,
                        help='maximum intramolecular base-pair span')
    parser.add_argument('-o', default=','.join(DEFAULT_COLUMNS), help='CSV columns to report')
    args = parser.parse_args(argv)
    args.columns = [column.strip() for column in args.o.split(',')]
    if not all(args.columns) or len(set(args.columns)) != len(args.columns):
        parser.error('-o must contain distinct, nonempty CSV column names')
    missing = [column for column in ESSENTIAL_COLUMNS if column not in args.columns]
    if missing:
        parser.error('-o is missing essential CSV columns: ' + ','.join(missing))
    return args


def find_program(name):
    """Prefer installed siblings, then the source-tree binary, then PATH."""
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir / name]
    if name == 'IntaRNA':
        candidates.append(script_dir.parent / 'src' / 'bin' / name)
    for candidate in candidates:
        found = shutil.which(str(candidate))
        if found:
            return found
    found = shutil.which(name)
    if not found:
        raise RuntimeError('{} executable not found beside the script or on PATH'.format(name))
    return found


def run_command(command, input_text=None):
    """Pass sequences and constraints without shell interpretation; check failures."""
    result = subprocess.run(command, input=input_text, universal_newlines=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError('{} exited with status {}: {}'.format(
            Path(command[0]).name, result.returncode, detail))
    if result.stderr:
        sys.stderr.write(result.stderr)
    return result.stdout


def energy(value):
    """Read a finite energy without introducing binary rounding in the sum."""
    try:
        number = Decimal(value)
    except InvalidOperation:
        raise ValueError('invalid energy: ' + value)
    if not number.is_finite():
        raise ValueError('nonfinite energy: ' + value)
    return number


def ensemble_energy(program, seq, span, constraint=None):
    """Compute the full or constrained single-strand ensemble free energy."""
    command = [program, '--maxBPspan={}'.format(span), '--noPS', '-p0']
    input_text = seq + '\n'
    if constraint is not None:
        command.append('-C')
        input_text += constraint + '\n'
    output = run_command(command, input_text)
    match = re.search(r'energy of ensemble\s*=\s*(\S+)\s+kcal', output)
    if match is None:
        raise ValueError('RNAfold did not report an ensemble energy')
    return energy(match.group(1))


def predictions(command, columns, lengths):
    """Read CSV and validate 1-based inclusive sites before constructing masks."""
    records = list(csv.reader(io.StringIO(run_command(command)), delimiter=';'))
    if not records or [col.strip() for col in records[0]] != columns:
        raise ValueError('IntaRNA returned an unexpected CSV header')
    rows = []
    for record in records[1:]:
        if not record:
            continue
        if len(record) != len(columns):
            raise ValueError('IntaRNA returned an incomplete CSV row')
        row = dict(zip(columns, record))
        for strand, length in enumerate(lengths, 1):
            start, end = int(row['start' + str(strand)]), int(row['end' + str(strand)])
            if not 1 <= start <= end <= length:
                raise ValueError('IntaRNA returned an out-of-range interaction site')
        for column in ('E', 'ED1', 'ED2', 'E_init'):
            energy(row[column])
        rows.append(row)
    return rows


def overlaps_mfe(mfe, constrained):
    """Preserve the Perl criterion: overlap spans at least two bases per strand."""
    return all(min(int(mfe['end' + s]), int(constrained['end' + s])) >
               max(int(mfe['start' + s]), int(constrained['start' + s]))
               for s in ('1', '2'))


def unpaired_constraint(length, strand, sites):
    """Block both inclusive interaction sites in the original strand orientation."""
    mask = ['.'] * length
    for site in sites:
        start, end = int(site['start' + strand]), int(site['end' + strand])
        mask[start - 1:end] = ['x'] * (end - start + 1)
    return ''.join(mask)


def format_energy(value):
    text = format(value, 'f')
    return text.rstrip('0').rstrip('.') if '.' in text else text


def fuse(args):
    """Return original predictions followed by each compatible fused row."""
    intarna, rnafold = find_program('IntaRNA'), find_program('RNAfold')
    command = [intarna, '-t', args.t, '--tAccW=0', '--tAccL={}'.format(args.l),
               '-q', args.q, '--qAccW=0', '--qAccL={}'.format(args.l),
               '--outOverlap=N', '-n', str(args.n), '--outMode=C',
               '--outCsvCols=' + ','.join(args.columns)]
    lengths = (len(args.t), len(args.q))
    subopts = predictions(command, args.columns, lengths)
    if len(subopts) < 2:
        return subopts

    # Reuse each unconstrained single-strand ensemble for all candidate pairs.
    target_full = ensemble_energy(rnafold, args.t, args.l)
    query_full = ensemble_energy(rnafold, args.q, args.l)
    output = [subopts[0]]
    for subopt in subopts[1:]:
        output.append(subopt)
        constraints = ['--tAccConstr=b:{}-{}'.format(subopt['start1'], subopt['end1']),
                       '--qAccConstr=b:{}-{}'.format(subopt['start2'], subopt['end2'])]
        constrained = predictions(command + constraints, args.columns, lengths)
        if not constrained or not overlaps_mfe(subopts[0], constrained[0]):
            continue

        # Recompute the cost of opening the union of sites on each strand.
        partner = constrained[0]
        ed1 = ensemble_energy(rnafold, args.t, args.l,
                              unpaired_constraint(lengths[0], '1', (subopt, partner))) - target_full
        ed2 = ensemble_energy(rnafold, args.q, args.l,
                              unpaired_constraint(lengths[1], '2', (subopt, partner))) - query_full
        combined = sum(energy(site['E']) - energy(site['ED1']) - energy(site['ED2'])
                       for site in (subopt, partner)) + ed1 + ed2

        # Keep both initiation terms, as in the original fusion energy formula.
        fused = {col: subopt[col] + ':' + partner[col] for col in args.columns}
        fused.update(E=format_energy(combined), ED1=format_energy(ed1), ED2=format_energy(ed2))
        output.append(fused)
    return output


def main(argv=None):
    args = parse_args(argv)
    try:
        rows = fuse(args)
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        print('ERROR: {}'.format(error), file=sys.stderr)
        return 1
    # Publish only after every external call and energy calculation succeeded.
    writer = csv.DictWriter(sys.stdout, fieldnames=args.columns, delimiter=';', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return 0


if __name__ == '__main__':
    sys.exit(main())
