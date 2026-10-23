"""Sampling and statistical regressions, independent of an installed IntaRNA."""

from collections import Counter
import io
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from intarnapvalue.intarna_pvalue import IntaRNApvalue, main
from intarnapvalue.dinucl_shuffle import dinucl_shuffle
from intarna_tools import read_sequence

QUERY = 'AGGAUGGGGGA'
TARGET = 'CCGCUUCUCGUUCCU'
ARGS = ['-q', QUERY, '-t', TARGET, '-s', '5', '-m', 'b', '--randSeed', 'reproducible']


def fasta_sequences(text):
    return text.splitlines()[1::2]


class PvalueTests(unittest.TestCase):
    def tool(self, *args):
        with patch.object(IntaRNApvalue, 'find_binary', return_value='/tmp/IntaRNA with spaces'):
            return IntaRNApvalue([*ARGS, *args])

    def test_normalized_wrapped_single_fasta_and_invalid_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'seq.fa'
            path.write_text('\ufeff>sequence\nacgt\ntgca\n')
            self.assertEqual(read_sequence(str(path)), 'ACGUUGCA')
            self.assertEqual(read_sequence('a' * 10000), 'A' * 10000)
            for invalid in ('', '>empty\n', 'ACGU\n', '>one\nAC\n>two\nGU\n', '>one\nACN\n'):
                path.write_text(invalid)
                with self.assertRaises(ValueError):
                    read_sequence(str(path))
        for args in (['-q', ''], ['-q', 'AX'], ['-s', '0'], ['--threads=-1'], ['--max-attempts', '4']):
            with patch('sys.stderr', new_callable=io.StringIO), self.assertRaises(SystemExit):
                self.tool(*args)

    def test_shuffle_composition_endpoints_and_singletons(self):
        rng = random.Random(13)
        for seq in ('A', 'acgt', 'G' * 1000, QUERY, TARGET):
            normalized = seq.upper().replace('T', 'U')
            for _ in range(30):
                shuffled = dinucl_shuffle(seq, rng)
                self.assertEqual(len(shuffled), len(normalized))
                self.assertEqual((shuffled[0], shuffled[-1]), (normalized[0], normalized[-1]))
                self.assertEqual(Counter(shuffled), Counter(normalized))
                self.assertEqual(Counter(zip(shuffled, shuffled[1:])),
                                 Counter(zip(normalized, normalized[1:])))
        for seq in ('', 'N', 'AX'):
            with self.assertRaises(ValueError):
                dinucl_shuffle(seq)

    def test_local_rng_deterministic_and_no_global_reseed(self):
        state = random.getstate()
        a, b = self.tool(), self.tool()
        self.assertEqual(random.getstate(), state)
        self.assertEqual(a.shuffle_sequence(QUERY, 20, a.rng), b.shuffle_sequence(QUERY, 20, b.rng))
        self.assertEqual(random.getstate(), state)

    def test_both_strands_independent_pairwise_and_cleanup(self):
        tool = self.tool()
        seen = []
        def predict(query, target, input_text=None, *, pairwise=False):
            self.assertTrue(pairwise)
            self.assertIsNone(input_text)
            q, t = Path(query), Path(target)
            seen.extend((q, t))
            queries, targets = fasta_sequences(q.read_text()), fasta_sequences(t.read_text())
            self.assertEqual(len(queries), 5)
            self.assertEqual(len(targets), 5)
            self.assertGreater(len(set(queries)), 1)
            for sequences, original in ((queries, QUERY), (targets, TARGET)):
                for seq in sequences:
                    self.assertEqual(Counter(seq), Counter(original))
                    self.assertEqual(Counter(zip(seq, seq[1:])), Counter(zip(original, original[1:])))
            return [(str(i), str(i), -i) for i in reversed(range(5))]
        with patch.object(tool, '_predict', side_effect=predict):
            self.assertEqual(tool.get_scores(), ([0, -1, -2, -3, -4], 0))
        self.assertTrue(all(not p.exists() for p in seen))

    def test_each_single_strand_mode_and_retry_counts(self):
        for mode in ('q', 't'):
            tool = self.tool('-m', mode)
            attempts = []
            def predict(query, target, input_text=None):
                self.assertEqual((query, target), ('STDIN', TARGET) if mode == 'q' else (QUERY, 'STDIN'))
                seqs = fasta_sequences(input_text)
                attempts.append(len(seqs))
                # Reject one pair on the first pass, then fill the remainder.
                start = 1 if len(attempts) == 1 else 0
                return [('fixed' if mode == 'q' else str(i), str(i) if mode == 'q' else 'fixed', -i)
                        for i in range(start, len(seqs))]
            with patch.object(tool, '_predict', side_effect=predict):
                scores, skipped = tool.get_scores()
            self.assertEqual((len(scores), skipped), (5, 1))
            self.assertEqual(attempts, [5, 1])

    def test_bounded_no_interactions_and_fail_fast(self):
        tool = self.tool('--max-attempts', '7')
        with patch.object(tool, '_predict', return_value=[]) as run, self.assertRaisesRegex(RuntimeError, '7 attempts'):
            tool.get_scores()
        self.assertEqual(run.call_count, 2)
        with patch.object(tool, '_predict', side_effect=RuntimeError('command failed')) as run, \
                self.assertRaisesRegex(RuntimeError, 'command failed'):
            tool.get_scores()
        self.assertEqual(run.call_count, 1)
        for rows in ([('0', '0', -1), ('0', '0', -2)], [('0', '1', -1)], [('5', '5', -1)]):
            with patch.object(tool, '_predict', return_value=rows), self.assertRaises(ValueError):
                tool.get_scores()

    def test_parameter_file_all_predictions_and_zero_energy(self):
        tool = self.tool()
        tool.parameterFile = '/tmp/parameter file with spaces'
        with patch('intarnapvalue.intarna_pvalue.run_command', return_value='id1;id2;E\nt;q;0\n') as run:
            self.assertEqual(tool.get_single_score(QUERY, TARGET), 0)
            tool._predict('STDIN', TARGET, '>0\nAAA\n')
            for call in run.call_args_list:
                command = call.args[0]
                self.assertIn('--parameterFile=/tmp/parameter file with spaces', command)
                self.assertIn('--outNumber=1', command)
                self.assertNotIn('', command)
        for bad in ('E\n-1\n', 'id1;id2;E\na;b;nan\n', 'id1;id2;E\na;b\n'):
            with patch('intarnapvalue.intarna_pvalue.run_command', return_value=bad), self.assertRaises(ValueError):
                tool.get_single_score(QUERY, TARGET)
        with patch('intarnapvalue.intarna_pvalue.run_command', return_value='id1;id2;E\n'):
            self.assertIsNone(tool.get_single_score(QUERY, TARGET))

    def test_cdf_matches_independent_legacy_integration(self):
        from scipy import integrate, stats
        scores = [-1.235, -1.435645, -6.234234, -12.999, -15.23, -6.98, -6.23, -2.78]
        for name, calculate in (('norm', IntaRNApvalue.calculate_pvalue_gauss),
                                ('gumbel_l', IntaRNApvalue.calculate_pvalue_gumbel),
                                ('genextreme', IntaRNApvalue.calculate_pvalue_gev)):
            distribution = getattr(stats, name)
            params = distribution.fit(scores)
            expected = integrate.quad(lambda x: distribution.pdf(x, *params), -math.inf, -10)[0]
            self.assertAlmostEqual(calculate(-10, scores), expected, places=7)
            self.assertLessEqual(calculate(-20, scores), calculate(-10, scores))
        self.assertEqual(IntaRNApvalue.calculate_pvalue_empirical(-10, scores), 0.25)
        self.assertEqual(IntaRNApvalue.calculate_pvalue_empirical(0, [0]), 1)
        for scores in ([], [float('nan')], [float('inf')], [-1, -1, -1], [-1, -2]):
            for method in (IntaRNApvalue.calculate_pvalue_gev, IntaRNApvalue.calculate_pvalue_gauss,
                           IntaRNApvalue.calculate_pvalue_gumbel):
                with self.assertRaises(ValueError):
                    method(-1, scores)

    def test_cli_help_without_binary_or_dependencies(self):
        result = subprocess.run([sys.executable, '-S', '-m', 'intarnapvalue', '--help'],
                                env={**os.environ, 'PYTHONPATH': str(ROOT), 'PATH': ''},
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--samples', result.stdout)

    def test_cli_empirical_scores_and_errors(self):
        with patch.object(IntaRNApvalue, 'find_binary', return_value='IntaRNA'), \
                patch.object(IntaRNApvalue, 'get_single_score', return_value=0), \
                patch.object(IntaRNApvalue, 'get_scores', return_value=([-1, 0, -2, -3, -4], 2)), \
                patch('sys.stdout', new_callable=io.StringIO) as stdout, \
                patch('sys.stderr', new_callable=io.StringIO) as stderr:
            self.assertEqual(main([*ARGS, '-d', 'empirical']), 0)
            self.assertEqual(stdout.getvalue(), '1.0\n')
            self.assertIn('conditional', stderr.getvalue())
        with patch.object(IntaRNApvalue, 'find_binary', side_effect=RuntimeError('binary missing')), \
                patch('sys.stdout', new_callable=io.StringIO) as stdout, \
                patch('sys.stderr', new_callable=io.StringIO) as stderr:
            self.assertEqual(main(ARGS), 1)
            self.assertEqual(stdout.getvalue(), '')
            self.assertIn('binary missing', stderr.getvalue())


if __name__ == '__main__':
    unittest.main()
