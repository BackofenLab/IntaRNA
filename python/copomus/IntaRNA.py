# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

"""Checked IntaRNA subprocess and CSV adapter used by CopomuS."""

import csv
import io
import math
import re
from intarna_tools import csv_header, csv_rows, find_program, run_command


class IntaRNA:
    """IntaRNA Python wrapper; no-interaction output is represented by False."""

    @staticmethod
    def version_tuple(versionstring):
        """Read numeric release components, accepting development suffixes."""
        match = re.match(r'(\d+)\.(\d+)\.(\d+)', versionstring)
        if not match:
            raise ValueError(f'Unrecognized IntaRNA version: {versionstring}')
        return tuple(map(int, match.groups()))

    @staticmethod
    def check_binary(required_version):
        """Check executable availability and the minimum supported release."""
        output = run_command([find_program('IntaRNA'), '--version'])
        match = re.search(r'IntaRNA\s+(\d+\.\d+\.\d+\S*)', output)
        if not match or IntaRNA.version_tuple(match[1]) < IntaRNA.version_tuple(required_version):
            raise RuntimeError(f'This tool requires IntaRNA version {required_version} or newer')

    def run(self, query, target, qidxpos0, tidxpos0, outcsvcols, threads, n=1,
            param_file='', extra_params=None, raw_stdout=False):
        """Return one typed CSV row, all suboptimal rows, or raw stdout."""
        command = [find_program('IntaRNA'), '-q', query, '-t', target, '--outMode=C',
                   f'--outCsvCols={outcsvcols}', f'--qIdxPos0={qidxpos0}',
                   f'--tIdxPos0={tidxpos0}', f'--outNumber={n}', f'--threads={threads}']
        if param_file:
            command.append(f'--parameterFile={param_file}')
        stdout = run_command(command + list(extra_params or []))
        if raw_stdout:
            return stdout
        rows = self.output_to_dict(stdout)
        return rows[0] if rows and n == 1 else rows

    @staticmethod
    def output_to_dict(data):
        r"""Parse quoted CSV, preserving IDs and numeric zero; only NA/empty is missing.

        >>> IntaRNA.output_to_dict('start1;end1;id1;id2;E\n34;42;target;query;-4.7\n')
        [{'start1': 34, 'end1': 42, 'id1': 'target', 'id2': 'query', 'E': -4.7}]
        """
        reader = csv.reader(io.StringIO(data), delimiter=';', strict=True)
        header = csv_header(reader)
        rows = []
        for fields in csv_rows(reader, header):
            row = {}
            for name, value in zip(header, fields):
                if name in ('id1', 'id2', 'seq1', 'seq2', 'subseq1', 'subseq2'):
                    row[name] = value
                    continue
                if value in ('', 'NA'):
                    row[name] = None
                    continue
                try:
                    number = int(value)
                except ValueError:
                    try:
                        number = float(value)
                    except ValueError:
                        number = value
                if isinstance(number, float) and not math.isfinite(number):
                    raise ValueError(f'IntaRNA returned a nonfinite {name}')
                row[name] = number
            rows.append(row)
        return rows or False

    @staticmethod
    def to_fasta(sequences):
        """Format a mapping of record names to sequences."""
        return ''.join(f'>{name}\n{seq}\n' for name, seq in sequences.items())
