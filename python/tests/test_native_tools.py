"""Optional integration tests against a native IntaRNA binary on PATH."""

import csv
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
QUERY = 'GCUGAAAAACAUAACCCAUAAAAUGCUAGCUGUACCAGGAACCA'
TARGET = 'GGUUUCUUCGCCUCUGCGUUCACCAAAGUGUUCACCC'


@unittest.skipUnless(shutil.which('IntaRNA'), 'native IntaRNA is not on PATH')
class NativeToolsTests(unittest.TestCase):
    def run_tool(self, command, *, success=True):
        result = subprocess.run([sys.executable, *command], cwd=ROOT.parent,
                                env={**os.environ, 'PYTHONPATH': str(ROOT)},
                                capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        return result

    def test_sampling_modes_reproducible_across_threads(self):
        for mode in ('q', 't', 'b'):
            args = ['-m', 'intarnapvalue', '-q', QUERY, '-t', TARGET, '-s', '5',
                    '-m', mode, '--randSeed', 'native', '-o', 'scores']
            a = self.run_tool([*args, '--threads', '1'])
            b = self.run_tool([*args, '--threads', '2'])
            self.assertEqual(a.stdout, b.stdout)
            self.assertEqual(len([float(x) for x in a.stdout.splitlines()]), 5)

    def test_native_no_interaction_attempt_limit(self):
        result = self.run_tool(['-m', 'intarnapvalue', '-q', 'AAAAAAAA', '-t', 'AAAAAAAA',
                               '-s', '1', '-m', 'b', '-o', 'scores', '--threads', '1',
                               '--max-attempts', '2'], success=False)
        self.assertEqual(result.stdout, '')
        self.assertIn('2 attempts', result.stderr)

    def test_native_parameter_file_and_empirical(self):
        with tempfile.TemporaryDirectory() as directory:
            params = Path(directory) / 'params with spaces.txt'
            params.write_text('seedBP=3\nqAcc=N\ntAcc=N\n')
            result = self.run_tool(['-m', 'intarnapvalue', '-q', QUERY, '-t', TARGET,
                                    '-s', '5', '-m', 'b', '--threads', '1', '-d', 'empirical',
                                    '--randSeed', 'test', '-p', str(params)])
            self.assertTrue(0 <= float(result.stdout) <= 1)

    def test_native_copomus_profiles_signed_indices_and_explicit_mutation(self):
        result = self.run_tool([str(ROOT / 'CopomuS.py'), '-q', QUERY, '-t', TARGET,
                               '--threads', '1', '-m', 'EDqi', '-m', 'cEDqi', '-m', 'Eqi', '-m', 'Eti',
                               '--qIdxPos0', '-20', '--tIdxPos0', '-15', '-f', 'lpMfe'])
        rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter='\t'))
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(row['EDqi_ww'], row['cEDqi_ww'])
            self.assertNotEqual(row['qIndex'], '0')
            self.assertNotEqual(row['tIndex'], '0')
        encoding = rows[0]['mutation']
        result = self.run_tool([str(ROOT / 'CopomuS.py'), '-q', QUERY, '-t', TARGET,
                               '--threads', '1', '--qIdxPos0', '-20', '--tIdxPos0', '-15',
                               '--mutationEncoding', encoding])
        explicit = list(csv.DictReader(io.StringIO(result.stdout), delimiter='\t'))
        self.assertEqual(len(explicit), 1)
        self.assertEqual(explicit[0]['mutation'], encoding)


if __name__ == '__main__':
    unittest.main()
