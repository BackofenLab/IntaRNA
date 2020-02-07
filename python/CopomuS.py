#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>


__intarna_min__ = '3.1.2'
__python_min__ = '3.7.0'

import os
import sys
from platform import python_version

if python_version() < __python_min__:
    sys.stderr.write(f'CopomuS requires python >= {__python_min__}\n')
    sys.exit(1)

import csv
from tempfile import gettempdir
from argparse import ArgumentParser, FileType, SUPPRESS
from typing import List
from copomus.mutation import Mutation
from copomus.IntaRNA import IntaRNA
from copomus.candidate_selectors import get_selector
from copomus.candidate_filters import get_filter
from copomus.mutation_generators import get_generator



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

    def main(self):
        # PARSE COMMAND LINE ARGS
        self.parse_cmd_args()

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
            sys.stderr.write('There are no base pair candidates left after filtering!\n')
            sys.exit(1)

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

    def parse_cmd_args(self):
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
        parser.add_argument('-q', '--query', dest='query', required=True, help='The query sequence.')
        parser.add_argument('-t', '--target', dest='target', required=True, help='The target sequence.')
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
                            choices=['GU', 'AU', 'CG', 'lp', 'lpMfe', 'he', 'heMfe'])
        parser.add_argument('-g', '--generator', dest='mutation_generator', default=defaults['generator'],
                            choices=['flip', 'any'], help='Defines the method used for generating mutated sequences. '
                                                          f"(Default: {defaults['generator']})")
        parser.add_argument('--mutationEncoding', dest='mutation_encoding',
                            help='Allows direct candidate selection by specifying a mutation encoding. '
                                 'Overwrites options -s, -f, and -g')
        parser.add_argument('-o', '--output', dest='output', nargs='?', type=FileType('w'), default=sys.stdout,
                            help='Which file the output should be written to. (Default: STDOUT)')
        parser.add_argument('-d', '--delimiter', dest='delimiter', default='\t',
                            help='Defines the delimiter used to separate columns in the output, default tab. '
                                 '(Default: \\t)')
        parser.add_argument('-p', '--parameterFile', dest='param_file', default='',
                            help='Optional parameter file for IntaRNA to provide further parameters and '
                                 'prediction constraints.')
        parser.add_argument('--threads', dest='threads', default=defaults['threads'], type=int,
                            help='Threads used for IntaRNA call')
        parser.add_argument('--alpha', dest='alpha', type=float, default=defaults['alpha'], help=SUPPRESS)
        parser.add_argument('--beta', dest='beta', type=float, default=defaults['beta'], help=SUPPRESS)

        args = parser.parse_args()

        mm = []  # filter out any duplicate entries for measures
        args.measures = defaults['measures'] if not args.measures else args.measures
        for m in args.measures:
            if m not in mm:
                mm.append(m)
        args.measures = mm

        ff = []  # filter out any duplicate entries for filters
        args.candidate_filters = defaults['filters'] if not args.candidate_filters else args.candidate_filters
        for f in args.candidate_filters:
            if f not in ff:
                ff.append(f)
        args.candidate_filters = ff

        args.query = args.query.upper().replace('T', 'U')  # capitalize sequences and replace T -> U
        args.target = args.target.upper().replace('T', 'U')

        if args.param_file and not os.path.exists(args.param_file):
            sys.stderr.write(f'Cannot find provided parameter file: {args.param_file}\n')
            sys.exit(1)

        for key, value in args.__dict__.items():
            self.__setattr__(key, value)

    def run_tests(self, mutations: List[Mutation]):
        for mut in mutations:
            outcsvcols = list(set([j for i in [mm.outcsvcols for mm in mut.measures] for j in i]))  # combine lists
            outcsvcols = ','.join(outcsvcols)  # join columns
            outcsvcols = 'id1' if not outcsvcols else outcsvcols  # if no col -> select id1

            out_pipes = [j for i in [m.out_pipes for m in mut.measures] for j in i]  # combine lists

            # create args for IntaRNA call
            out_param = [f'--out={x}:{os.path.join(gettempdir(), f"CopomuS_{x}")}_{os.getpid()}.temp' for x in out_pipes]

            for mode in ['ww', 'wm', 'mw', 'mm']:  # go through each mutation mode
                q = mut.qw if mode[0] == 'w' else mut.qm  # get actual sequences
                t = mut.tw if mode[1] == 'w' else mut.tm

                i = IntaRNA()  # call IntaRNA, we have 4 calls per base pair, one for each mutation mode
                data = i.run(q, t, mut.qidxpos0, mut.tidxpos0, outcsvcols, self.threads, 1, self.param_file, out_param)

                for m in mut.measures:  # add results to each measure
                    m.calculate_score(data, mode, mut, self.param_file)

            for m in mut.measures:
                m.post_compute()

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
        # CREATE CSV WRITER
        writer = csv.writer(self.output, delimiter=self.delimiter)

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


if __name__ == '__main__':
    c = CopomuS()
    c.main()
