#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

"""Estimate interaction-energy significance against dinucleotide shuffles."""

import argparse
import csv
import io
import math
from pathlib import Path
import random
import sys
from tempfile import TemporaryDirectory

from intarna_tools import (csv_header, csv_rows, find_program, nonnegative_integer,
                           positive_integer, read_fasta, read_sequence, run_command)
from intarnapvalue.dinucl_shuffle import dinucl_shuffle


class IntaRNApvalue:
    """Configure from argv; explicit argument lists do not run predictions.

    ``get_scores()`` returns favorable shuffled energies and the number of
    rejected pairs. The fitted distribution is conditional on a reported
    interaction under the supplied IntaRNA constraints.
    """

    def __init__(self, test_args=None):
        self.process_cmd_args(test_args)
        self.bin = self.find_binary()
        if test_args is None:
            self.main()

    @staticmethod
    def find_binary():
        """Use the same executable lookup as the other companion tools."""
        return find_program('IntaRNA')

    def process_cmd_args(self, test_args=None):
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-q', '--query', required=True, help='RNA/DNA sequence or single-record FASTA')
        parser.add_argument('-t', '--target', required=True, help='RNA/DNA sequence or single-record FASTA')
        parser.add_argument('-s', '--samples', type=positive_integer, required=True,
                            help='number of favorable shuffled pairs to collect')
        parser.add_argument('-m', '--shuffle-mode', dest='shuffle_mode', required=True,
                            choices=['q', 't', 'b'], help='shuffle query, target, or both independently')
        parser.add_argument('-p', '--parameterFile', default='', help='IntaRNA parameter file for ALL predictions')
        parser.add_argument('-d', '--distribution', dest='dist',
                            choices=['gev', 'gumbel', 'gauss', 'empirical', 'none'], default='gev',
                            help='fitted distribution; empirical (alias none) counts scores <= observed')
        parser.add_argument('-o', '--output', choices=['pvalue', 'scores'], default='pvalue',
                            help='print one p-value or one shuffled energy per line')
        parser.add_argument('--threads', type=nonnegative_integer, default=0,
                            help='IntaRNA threads; zero selects all available CPUs')
        parser.add_argument('--randSeed', default=None, help='random seed string for reproducible shuffling')
        parser.add_argument('--max-attempts', type=positive_integer,
                            help='maximum shuffled pairs tried (default: 100 times --samples)')
        args = parser.parse_args(test_args)
        try:
            args.query = read_sequence(args.query)
            args.target = read_sequence(args.target)
            if args.parameterFile and not Path(args.parameterFile).is_file():
                raise ValueError(f'Cannot find provided parameter file: {args.parameterFile}')
            args.max_attempts = args.max_attempts or 100 * args.samples
            if args.max_attempts < args.samples:
                raise ValueError('--max-attempts must be at least --samples')
        except (OSError, ValueError) as error:
            parser.error(str(error))
        for key, value in vars(args).items():
            setattr(self, key, value)
        self.shuffle_query = args.shuffle_mode in ('q', 'b')
        self.shuffle_target = args.shuffle_mode in ('t', 'b')
        self.rng = random.Random(args.randSeed)

    read_fasta_file = staticmethod(read_fasta)

    @staticmethod
    def shuffle_sequence(seq, n, rng=None):
        """Generate n independent, composition-preserving RNA shuffles."""
        if n < 0:
            raise ValueError('shuffle count cannot be negative')
        return [dinucl_shuffle(seq, rng=rng) for _ in range(n)]

    @staticmethod
    def to_fasta(sequences):
        """Number the records so batch results can be reordered deterministically."""
        return ''.join(f'>{index}\n{seq}\n' for index, seq in enumerate(sequences))

    def _predict(self, query, target, input_text=None, *, pairwise=False):
        """Apply identical constraints to observed and shuffled pairs; fail once."""
        command = [self.bin, '-q', query, '-t', target, '--outMode=C',
                   '--outCsvCols=id1,id2,E', '--outSep=;', '--outNumber=1',
                   '--outPerRegion=0', f'--outPairwise={int(pairwise)}',
                   f'--threads={self.threads}']
        if self.parameterFile:
            command.append(f'--parameterFile={self.parameterFile}')
        reader = csv.reader(io.StringIO(run_command(command, input_text)), delimiter=';', strict=True)
        header = csv_header(reader)
        if header != ['id1', 'id2', 'E']:
            raise ValueError('IntaRNA returned an unexpected CSV header')
        rows = []
        for target_id, query_id, value in csv_rows(reader, header):
            score = float(value)
            if not math.isfinite(score):
                raise ValueError('IntaRNA returned a nonfinite energy')
            rows.append((target_id, query_id, score))
        return rows

    def get_single_score(self, query, target):
        """Return the observed energy, including zero, or None for no interaction."""
        rows = self._predict(query, target)
        if len(rows) > 1:
            raise ValueError('IntaRNA returned multiple energies for one pair')
        return rows[0][2] if rows else None

    def get_scores(self):
        """Collect favorable scores in bounded batches; reject invalid/duplicate IDs."""
        scores, attempts = [], 0
        # Two files allow IntaRNA's pairwise mode to pair independent query/target
        # shuffles by index, avoiding a Cartesian product and repeated query draws.
        with TemporaryDirectory(prefix='IntaRNApvalue-') as directory:
            query_file = Path(directory) / 'query.fa'
            target_file = Path(directory) / 'target.fa'
            while len(scores) < self.samples and attempts < self.max_attempts:
                count = min(1000, self.samples - len(scores), self.max_attempts - attempts)
                if self.shuffle_query and self.shuffle_target:
                    query_file.write_text(self.to_fasta(self.shuffle_sequence(self.query, count, self.rng)),
                                          encoding='utf-8')
                    target_file.write_text(self.to_fasta(self.shuffle_sequence(self.target, count, self.rng)),
                                           encoding='utf-8')
                    rows = self._predict(str(query_file), str(target_file), pairwise=True)
                elif self.shuffle_query:
                    fasta = self.to_fasta(self.shuffle_sequence(self.query, count, self.rng))
                    rows = self._predict('STDIN', self.target, fasta)
                else:
                    fasta = self.to_fasta(self.shuffle_sequence(self.target, count, self.rng))
                    rows = self._predict(self.query, 'STDIN', fasta)
                batch = {}
                for target_id, query_id, score in rows:
                    identifier = query_id if self.shuffle_query else target_id
                    index = int(identifier)
                    if not 0 <= index < count or index in batch:
                        raise ValueError('IntaRNA returned an invalid or duplicate shuffled-pair ID')
                    if self.shuffle_query and self.shuffle_target and target_id != query_id:
                        raise ValueError('IntaRNA returned an unpaired shuffle result')
                    batch[index] = score
                scores.extend(batch[index] for index in sorted(batch))
                attempts += count
        if len(scores) != self.samples:
            raise RuntimeError(f'Only {len(scores)} of {self.samples} favorable shuffled interactions '
                               f'found in {attempts} attempts; adjust constraints or --max-attempts')
        return scores, attempts - len(scores)

    @staticmethod
    def _sample(original_score, scores, *, fit=False):
        """Reject nonfinite, empty or degenerate fitting data before SciPy runs."""
        if not math.isfinite(original_score):
            raise ValueError('observed energy must be finite')
        values = list(scores)
        if not values or not all(math.isfinite(x) for x in values):
            raise ValueError('scores must contain finite energies')
        if fit and (len(values) < 3 or len(set(values)) < 3):
            raise ValueError('distribution fitting needs at least three distinct energies; use empirical otherwise')
        return values

    @staticmethod
    def calculate_pvalue_empirical(original_score, scores):
        """Fraction of favorable shuffled energies at least as low as observed."""
        values = IntaRNApvalue._sample(original_score, scores)
        return sum(score <= original_score for score in values) / len(values)

    @staticmethod
    def _fitted_pvalue(original_score, scores, name):
        from scipy import stats
        values = IntaRNApvalue._sample(original_score, scores, fit=True)
        distribution = getattr(stats, name)
        # Keep the historical fits on raw energies (including the GEV shape).
        # The genome-wide CSV tool fits negated energies: a different model.
        params = distribution.fit(values)
        if not all(math.isfinite(x) for x in params) or params[-1] <= 0:
            raise ValueError('distribution fit returned invalid parameters')
        if not all(math.isfinite(float(x)) for x in distribution.logpdf(values, *params)):
            raise ValueError('distribution fit does not cover all sampled energies')
        result = float(distribution.cdf(original_score, *params))
        if not math.isfinite(result) or not 0 <= result <= 1:
            raise ValueError('distribution fit returned an invalid p-value')
        return result

    @staticmethod
    def calculate_pvalue_gauss(original_score, scores):
        return IntaRNApvalue._fitted_pvalue(original_score, scores, 'norm')

    @staticmethod
    def calculate_pvalue_gumbel(original_score, scores):
        return IntaRNApvalue._fitted_pvalue(original_score, scores, 'gumbel_l')

    @staticmethod
    def calculate_pvalue_gev(original_score, scores):
        return IntaRNApvalue._fitted_pvalue(original_score, scores, 'genextreme')

    def main(self):
        """Publish output only after every required prediction and fit succeeds."""
        original_score = None
        if self.output != 'scores':
            original_score = self.get_single_score(self.query, self.target)
            if original_score is None:
                raise ValueError('The query/target combination has no favorable interaction')
        scores, non_interactions = self.get_scores()
        if non_interactions:
            print(f'Skipped {non_interactions} shuffled pairs without a reported interaction; '
                  'p-values are conditional on reported interactions.', file=sys.stderr)
        if self.output == 'scores':
            print('\n'.join(map(str, scores)))
        else:
            calculators = {'gauss': self.calculate_pvalue_gauss, 'gumbel': self.calculate_pvalue_gumbel,
                           'gev': self.calculate_pvalue_gev, 'empirical': self.calculate_pvalue_empirical,
                           'none': self.calculate_pvalue_empirical}
            print(calculators[self.dist](original_score, scores))


def main(argv=None):
    """CLI entry point; --help works without IntaRNA, NumPy or SciPy installed."""
    try:
        tool = IntaRNApvalue(sys.argv[1:] if argv is None else argv)
        tool.main()
    except ImportError as error:
        print(f'Missing Python dependency: {error}. Install numpy and scipy.', file=sys.stderr)
        return 1
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
