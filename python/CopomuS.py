#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>


__intarna_min__ = '3.1.2'

import math
import sys


import csv
from tempfile import TemporaryDirectory
from pathlib import Path
from argparse import ArgumentParser, SUPPRESS
from typing import List
from copomus.mutation import Mutation
from copomus.IntaRNA import IntaRNA
from copomus.candidate_selectors import get_selector
from copomus.candidate_filters import get_filter
from copomus.mutation_generators import get_generator


from intarna_tools import atomic_text_output, nonnegative_integer, read_sequence


class CopomuS:
    def __init__(self):
        self.query = ''
        self.target = ''
        self.qidxpos0 = 1
        self.tidxpos0 = 1
        self.param_file = ''
        self.measures = []
        self.candidate_selection = ''
        self.candidate_filters = []
        self.output = ''
        self.delimiter = ''
        self.mutation_generator = ''
        self.mutation_encoding = ''
        self.alpha, self.beta = 0.0, 0.0
        self.threads = 0

    def main(self, argv=None):
        # PARSE COMMAND LINE ARGS
        self.parse_cmd_args(argv)

        # CHECK INTARNA BIN AND VERSION
        IntaRNA.check_binary(__intarna_min__)

        # SELECT CANDIDATES
        bp_func = get_selector(self.candidate_selection) if not self.mutation_encoding else get_selector('mutEnc')
        bp = bp_func(self.query, self.target, self.qidxpos0, self.tidxpos0, self.param_file, self.threads,
                     enc=self.mutation_encoding)

        # FILTER OUT CANDIDATES
        if not self.mutation_encoding:
            for f in self.candidate_filters:
                func = get_filter(f)
                bp = func(bp, self.query, self.target, self.qidxpos0, self.tidxpos0)

        if not bp:
            raise ValueError('There are no base pair candidates left after filtering')

        # CREATE MUTATIONS
        if not self.mutation_encoding:
            gen_func = get_generator(self.mutation_generator)
            mutations = gen_func(self.query, self.target, self.qidxpos0, self.tidxpos0, bp, self.measures,
                                 self.alpha, self.beta)
        else:
            gen_func = get_generator('specific')
            mutations = gen_func(self.query, self.target, self.qidxpos0, self.tidxpos0, self.measures,
                                 self.alpha, self.beta, self.mutation_encoding)

        # RUN MEASUREMENTS
        self.run_tests(mutations)

        # SORT MUTATIONS
        mutations.sort()

        # CALCULATE RANKS
        self.calculate_ranks(mutations)

        # OUTPUT
        self.output_csv(mutations)

    def parse_cmd_args(self, argv=None):
        defaults = {
            'measures': ['mfeCover', 'E', 'minDeltaE'],
            'selector': 'mfe',
            'filters': [],
            'generator': 'flip',
            'alpha': 1.0,
            'beta': 1.0,
            'threads': 0
        }

        parser = ArgumentParser(description='Checks different measures for rating mutations')
        parser.add_argument('-q', '--query', dest='query', required=True, help='Query RNA/DNA sequence or single-record FASTA file.')
        parser.add_argument('-t', '--target', dest='target', required=True, help='Target RNA/DNA sequence or single-record FASTA file.')
        parser.add_argument('--qIdxPos0', dest='qidxpos0', default=1,
                            type=int, help='The starting index for the query. (Default: 1)')
        parser.add_argument('--tIdxPos0', dest='tidxpos0', default=1, type=int,
                            help='The starting index for the target. (Default: 1)')
        parser.add_argument('-m', '--measure', dest='measures', action='append', default=[],
                            help='Which measure to add to the output, can be used multiple times. '
                                 'Output will be sorted in order of measures specified. '
                                 f"(Default: {defaults['measures']})",
                            choices=['E', 'minDeltaE', 'mfeCover', 'EDqi', 'cEDqi', 'Eqi', 'Eti'])
        parser.add_argument('-s', '--candidateSelection', dest='candidate_selection', default=defaults['selector'],
                            help='Defines the method used to select candidate base pairs. '
                                 f"(Default: {defaults['selector']})",
                            choices=['mfe', 'mfeSO'])
        parser.add_argument('-f', '--candidateFilter', dest='candidate_filters', action='append',
                            default=[],
                            help='Filters candidate base pairs, can be used multiple times. '
                                 f"(Default: {defaults['filters']})",
                            choices=['GU', 'AU', 'CG', 'GC', 'lp', 'lpMfe', 'he', 'heMfe'])
        parser.add_argument('-g', '--generator', dest='mutation_generator', default=defaults['generator'],
                            choices=['flip', 'any'], help='Defines the method used for generating mutated sequences. '
                                                          f"(Default: {defaults['generator']})")
        parser.add_argument('--mutationEncoding', dest='mutation_encoding',
                            help='Allows direct candidate selection by specifying a mutation encoding. '
                                 'Overwrites options -s, -f, and -g')
        parser.add_argument('-o', '--output', dest='output', default=None,
                            help='Which file the output should be written to. (Default: STDOUT)')
        parser.add_argument('-d', '--delimiter', dest='delimiter', default='\t',
                            help='Defines the delimiter used to separate columns in the output, default tab. '
                                 '(Default: \\t)')
        parser.add_argument('-p', '--parameterFile', dest='param_file', default='',
                            help='Optional parameter file for IntaRNA to provide further parameters and '
                                 'prediction constraints.')
        parser.add_argument('--threads', dest='threads', default=defaults['threads'], type=nonnegative_integer,
                            help='Threads used for IntaRNA call')
        parser.add_argument('--alpha', dest='alpha', type=float, default=defaults['alpha'], help=SUPPRESS)
        parser.add_argument('--beta', dest='beta', type=float, default=defaults['beta'], help=SUPPRESS)

        args = parser.parse_args(argv)
        args.measures = list(dict.fromkeys(args.measures or defaults['measures']))
        args.candidate_filters = list(dict.fromkeys(args.candidate_filters))
        try:
            args.query = read_sequence(args.query)
            args.target = read_sequence(args.target)
            if not all(math.isfinite(value) for value in (args.alpha, args.beta)):
                raise ValueError('alpha and beta must be finite')
            if args.param_file and not Path(args.param_file).is_file():
                raise ValueError(f'Cannot find provided parameter file: {args.param_file}')
            if len(args.delimiter) != 1 or args.delimiter in '\r\n\0':
                raise ValueError('delimiter must be one non-newline character')
        except (OSError, ValueError) as error:
            parser.error(str(error))
        for key, value in vars(args).items():
            setattr(self, key, value)

    def run_tests(self, mutations: List[Mutation]):
        # Evaluate wildtype once: its sequences and requested outputs are identical
        # for every candidate. Profile contents are immutable and shared by measures.
        wildtype = None
        with TemporaryDirectory(prefix='CopomuS-') as directory:
            for mut in mutations:
                columns = list(dict.fromkeys(c for m in mut.measures for c in m.outcsvcols))
                pipes = list(dict.fromkeys(p for m in mut.measures for p in m.out_pipes))
                paths = {pipe: Path(directory) / pipe for pipe in pipes}
                out_params = [f'--out={pipe}:{path}' for pipe, path in paths.items()]
                for mode in ('ww', 'wm', 'mw', 'mm'):
                    if mode == 'ww' and wildtype is not None:
                        data, profiles = wildtype
                    else:
                        query = mut.qw if mode[0] == 'w' else mut.qm
                        target = mut.tw if mode[1] == 'w' else mut.tm
                        # Remove prior outputs so a missing profile cannot reuse stale data.
                        for path in paths.values():
                            path.unlink(missing_ok=True)
                        data = IntaRNA().run(query, target, mut.qidxpos0, mut.tidxpos0,
                                             ','.join(columns) or 'id1', self.threads,
                                             1, self.param_file, out_params)
                        profiles = {pipe: path.read_text(encoding='utf-8')
                                    for pipe, path in paths.items()}
                        if mode == 'ww':
                            wildtype = data, profiles
                    for measure in mut.measures:
                        measure.profiles = profiles
                        measure.calculate_score(data, mode, mut, self.param_file)
                for measure in mut.measures:
                    measure.post_compute()

    def calculate_ranks(self, mutations: List[Mutation]):
        # Calculate total ranking order
        rank = 1
        mutations[0].ranks['total'] = rank
        for i in range(1, len(mutations)):
            if mutations[i - 1] < mutations[i]:
                rank += 1
            mutations[i].ranks['total'] = rank

        # Calculate ranking for each measurement
        for i in range(len(self.measures)):
            rank = 1
            m_id = self.measures[i]

            measures = [(x.measures[i], j) for j, x in enumerate(mutations)]
            measures.sort(key=lambda x: x[0])

            mutations[measures[0][1]].ranks[m_id] = rank
            for j in range(1, len(measures)):
                measure, mut_id = measures[j]
                if measures[j - 1][0] < measures[j][0]:
                    rank += 1
                mutations[mut_id].ranks[m_id] = rank

    def output_csv(self, mutations: List[Mutation]):
        # Open a destination only after all predictions and ranking succeeded.
        if self.output:
            with atomic_text_output(self.output) as handle:
                self.write_csv(mutations, handle)
        else:
            self.write_csv(mutations, sys.stdout)

    def write_csv(self, mutations: List[Mutation], handle):
        # CREATE CSV WRITER
        writer = csv.writer(handle, delimiter=self.delimiter, lineterminator="\n")

        # WRITE CSV
        headers = [j for i in
                   [[f'{m.id}_rank'] + [f'{m.id}_{k}' for k in m.score.keys()] for m in mutations[0].measures]
                   for j in i]
        writer.writerow(['mutation', 'rank', 'qIndex', 'tIndex', 'bpWildtype', 'bpMutated'] + headers)

        for m in mutations:
            values = []
            for x in m.measures:
                values.append(m.ranks[x.id])
                scores = list(x.score.values())
                for i in range(len(scores)):  # replace None with 'NA'
                    if scores[i] is None:
                        scores[i] = 'NA'
                values += scores
            # values = [j for i in [list(x.score.values()) for x in r.measures] for j in i]
            writer.writerow([str(m), m.ranks['total'], m.query_n, m.target_n, m.bp_ww, m.bp_mm] + values)  # write csv


def main(argv=None):
    try:
        CopomuS().main(argv)
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
