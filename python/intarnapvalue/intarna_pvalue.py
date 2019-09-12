#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from subprocess import PIPE, Popen, run
import random
import os
import argparse
import sys
from typing import List, Tuple
import numpy as np
from scipy.integrate import quad as integ
from scipy.stats import norm as gauss
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_l as gum
from intarnapvalue import dinucl_shuffle


class IntaRNApvalue:
    def __init__(self, test_args=None):
        self.bin = self.find_binary()
        self.query = ''
        self.target = ''
        self.n = 0
        self.shuffle_query = False
        self.shuffle_target = False
        self.threads = ''
        self.dist = ''
        self.output = ''
        self.process_cmd_args(test_args)
        if not test_args:
            self.main()

    def main(self) -> None:
        """The main function"""
        if self.output == 'scores':  # output scores and exit the process
            scores, non_interactions = self.get_scores()
            print('\n'.join(iter([str(x) for x in scores])))
            sys.exit(0)

        original_score = self.get_single_score(self.query, self.target)
        if not original_score:  # exit if given seq has no interaction and not scores output only mode
            print('The query/target combination you specified has no favorable interaction')
            sys.exit(1)

        scores, non_interactions = self.get_scores()
        pvalue = 0.0
        if self.dist == 'gauss':
            pvalue = self.calculate_pvalue_gauss(original_score, scores)
        elif self.dist == 'none':
            pvalue = self.calculate_pvalue_empirical(original_score, scores)
        elif self.dist == 'gumbel':
            pvalue = self.calculate_pvalue_gumbel(original_score, scores)
        elif self.dist == 'gev':
            pvalue = self.calculate_pvalue_gev(original_score, scores)
        print(pvalue)

    @staticmethod
    def find_binary() -> str:
        """Tries to find the IntaRNA executable and returns path to binary or exits with error code 1 if not found"""
        if not run('IntaRNA --version', shell=True, stdout=PIPE, stderr=PIPE).returncode:
            return 'IntaRNA'

        # if binary not found in path, search in parent of this script recursively
        bin_name = 'IntaRNA.exe' if os.name == 'nt' else 'IntaRNA'
        for dir_path, dir_name, file_names in os.walk(os.path.abspath(os.path.join(os.curdir, '..'))):
            if bin_name in file_names:
                return os.path.join(dir_path, bin_name)

        print('Error: Cannot find IntaRNA binary executable, please add it to your PATH')
        sys.exit(1)

    def process_cmd_args(self, test_args=None) -> None:
        """Processes all commandline args and sets them as instance variables.

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '-n', '10', '-m', 'b', '-d', 'gauss',
        ... '--threads', '3'])
        >>> i.query
        'AGGAUG'
        >>> i.target
        'UUUAUCGUU'
        >>> i.n
        10
        >>> i.shuffle_query
        True
        >>> i.shuffle_target
        True
        >>> i.threads
        '3'
        >>> i.dist
        'gauss'
        >>> i.output
        'pvalue'
        >>> i = IntaRNApvalue(['-q', 'Z', '-t', 'UUUAUCGUU', '-n', '10', '-m', 'b', '--threads', '3'])
        Traceback (most recent call last):
        ...
        SystemExit: 1
        >>> open('test.fasta', 'w').write('>someseq\\nGACU')
        13
        >>> i = IntaRNApvalue(['-q', 'test.fasta', '-t', 'UUUAUCGUU', '-n', '10', '-m', 'b', '--threads', '3'])
        >>> i.query
        'GACU'
        >>> os.remove('test.fasta')
        """
        parser = argparse.ArgumentParser(description='Calculates p-values to IntaRNA scores.')
        parser.add_argument('-q', '--query', dest='query', type=str, help='Query sequence', required=True)
        parser.add_argument('-t', '--target', dest='target', type=str, help='Target sequence', required=True)
        parser.add_argument('-n', '--scores', dest='n', type=int, required=True,
                            help='How many randomly generated scores are used to calculate the p-value.')
        parser.add_argument('-m', '--shuffle-mode', dest='sm', required=True, choices=['q', 't', 'b'],
                            help='Which sequences are going to be shuffled: both, query only or target only.')
        parser.add_argument('-d', '--distribution', dest='dist', choices=['gev', 'gumbel', 'gauss'], default='gev',
                            help='Which distribution is fitted and used to calculate the pvalue.')
        parser.add_argument('-o', '--output', dest='output', choices=['pvalue', 'scores'], default='pvalue',
                            help='If set to scores, outputs all IntaRNA scores from random sequences to STDOUT. '
                                 'This is useful for pipeing the scores. Otherwise outputs just the pvalue.')
        parser.add_argument('--threads', type=str, default='0', help='Sets the amount of threads used for IntaRNA.')
        parser.add_argument('--seed', type=str, default=None,
                            help='Random seed to make sequence generation deterministic.')

        args = parser.parse_args(test_args)

        allowed = ['G', 'A', 'C', 'U', 'T']
        if os.path.exists(args.query):
            args.query = self.read_fasta_file(args.query)

        if os.path.exists(args.target):
            args.target = self.read_fasta_file(args.target)

        if False in [n in allowed for n in args.query.upper()] or False in [n in allowed for n in args.target.upper()]:
            print('A sequence you specified contains illegal characters, allowed: G, A, C, U (T)')
            sys.exit(1)

        self.shuffle_query = True if args.sm in ['b', 'q'] else False
        self.shuffle_target = True if args.sm in ['b', 't'] else False
        for key, value in args.__dict__.items():
            self.__setattr__(key, value)

        random.seed(a=args.seed)

    @staticmethod
    def read_fasta_file(filename: str) -> str:
        """Reads a FASTA file and returns the sequence. Exits the program if something goes wrong.

        >>> with open('test.fasta', 'w') as f:
        ...     f.write('>somename\\nGACUGGAGUGC')
        21
        >>> IntaRNApvalue.read_fasta_file('test.fasta')
        'GACUGGAGUGC'
        >>> IntaRNApvalue.read_fasta_file('test2.fasta')
        Traceback (most recent call last):
        ...
        SystemExit: 1
        >>> os.remove('test.fasta')
        """
        try:
            with open(filename) as f:
                data = f.read().split('\n')
                if not data[0].startswith('>'):
                    print('Error: {} is not in FASTA format!'.format(filename))
                    sys.exit(1)
                f.close()
                return data[1].upper().replace('T', 'U')
        except FileNotFoundError as e:  # In theory this can never happen, lets prevent it anyways
            print(e)
            sys.exit(1)
        except IndexError:
            print('Error: {} is not in FASTA format!'.format(filename))
            sys.exit(1)

    @staticmethod
    def shuffle_sequence(seq: str, n: int) -> List[str]:
        """Shuffles a sequence n times and returns a list of sequences, duplicate entries are possible

        >>> random.seed('IntaRNA')
        >>> IntaRNApvalue.shuffle_sequence('AGGAUGGGGGA', 5)
        ['AUGGAGGGGGA', 'AUGGGGAGGGA', 'AGGGGAUGGGA', 'AUGGGGGAGGA', 'AUGGGAGGGGA']
        """
        return [dinucl_shuffle.dinucl_shuffle(seq) for _ in range(n)]

    @staticmethod
    def to_fasta(sequences: List[str]) -> str:
        """Combines a list of sequences into a string in FASTA format

        >>> IntaRNApvalue.to_fasta(['AUGGAGGGGGA', 'AUGGGGAGGGA', 'AGGGGAUGGGA', 'AUGGGGGAGGA', 'AUGGGAGGGGA'])
        '>0\\nAUGGAGGGGGA\\n>1\\nAUGGGGAGGGA\\n>2\\nAGGGGAUGGGA\\n>3\\nAUGGGGGAGGA\\n>4\\nAUGGGAGGGGA\\n'
        """
        fasta_str = ''
        n = 0
        for seq in sequences:
            fasta_str += '>{}\n{}\n'.format(n, seq)
            n += 1
        return fasta_str

    def get_scores(self) -> Tuple[List[float], int]:
        """Calculates n IntaRNA scores from random sequences with given parameters as class variables"""
        scores = []
        missing = self.n
        non_interactions = 0

        while missing > 0:
            if self.shuffle_query and self.shuffle_target:  # shuffle both
                query = self.shuffle_sequence(self.query, 1)[0]  # get a random query
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing))
            elif self.shuffle_query and not self.shuffle_target:  # only shuffle query
                query = 'STDIN'
                target = self.target  # target stays the same
                shuffles = self.to_fasta(self.shuffle_sequence(self.query, missing))
            else:  # only shuffle target
                query = self.query  # query stays the same
                target = 'STDIN'
                shuffles = self.to_fasta(self.shuffle_sequence(self.target, missing))

            p = Popen([self.bin, '-q', query, '-t', target, '--outMode=C', '--outCsvCols=E', '--threads', self.threads],
                      stdout=PIPE, stdin=PIPE, universal_newlines=True)
            stdout, stderr = p.communicate(input=shuffles)  # send shuffles as STDIN
            if p.returncode:  # If IntaRNA exits with a returncode != 0, skip this iteration
                continue
            stdout = stdout.split('\n')  # split on newline
            del stdout[0], stdout[-1]  # remove first element aka 'E' and trailing newline element
            scores.extend(stdout)  # add elements to scores
            missing = self.n - len(scores)
            non_interactions += missing  # count non-interactions

        # return list with all elements as float and amount of non-interactions
        return [float(x) for x in scores], non_interactions

    def get_single_score(self, query, target) -> float:
        """Gets an IntaRNA score to a single query/target combination"""
        o = run('{} -q {} -t {} --outMode=C --outCsvCols=E --threads {}'.format(self.bin, query, target, self.threads),
                stdout=PIPE, stdin=PIPE, shell=True, universal_newlines=True).stdout
        if o.startswith('E') and o != 'E\n':  # Check that we got a result
            return float(o.split('\n')[1])
        else:
            return 0  # no interaction

    @staticmethod
    def calculate_pvalue_empirical(original_score, scores: list = None) -> float:
        """Calculates a p-value to a target/query combination empirical with a given amount of shuffle iterations

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-m', 'b', '--threads', '3'])
        >>> i.calculate_pvalue_empirical(-10.0, [-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.25
        """
        return [score <= original_score for score in scores].count(True) / len(scores)

    @staticmethod
    def calculate_pvalue_gauss(original_score, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a gaussian distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-m', 'b', '--threads', '3'])
        >>> i.calculate_pvalue_gauss(-10.0, [-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.2429106747265256
        """
        loc, scale = gauss.fit(scores)

        def f(x):
            return gauss.pdf(x, loc=loc, scale=scale)
        return integ(f, -np.inf, original_score)[0]

    @staticmethod
    def calculate_pvalue_gumbel(original_score: float, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a gumbel distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-m', 'b', '--threads', '3'])
        >>> i.calculate_pvalue_gumbel(-10.0, [-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.19721934073203196
        """
        loc, scale = gum.fit(scores)

        def f(x):
            return gum.pdf(x, loc=loc, scale=scale)
        return integ(f, -np.inf, original_score)[0]

    @staticmethod
    def calculate_pvalue_gev(original_score: float, scores: list) -> float:
        """Calculates a p-value to a target/query combination by int. with a given amount of shuffle iterations by
        fitting a generalized extreme value distribution and integrating from -inf to the original score

        >>> i = IntaRNApvalue(['-q', 'AGGAUG', '-t', 'UUUAUCGUU', '--scores', '10', '-m', 'b', '--threads', '3'])
        >>> i.calculate_pvalue_gev(-10.0, [-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78])
        0.17611816922560236
        """
        shape, loc, scale = gev.fit(scores)

        def f(x):
            return gev.pdf(x, shape, loc=loc, scale=scale)
        return integ(f, -np.inf, original_score)[0]


if __name__ == '__main__':
    i = IntaRNApvalue()
