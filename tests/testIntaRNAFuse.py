"""Black-box fusion regressions, plus an optional native IntaRNA/RNAfold check."""
import csv
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'python' / 'IntaRNA-fuse.py'
COLUMNS = 'E,start1,end1,start2,end2,ED1,ED2,E_init,hybridDB'.split(',')
TARGET = 'AAAACCCCGGGGUUUU'
QUERY = 'CCCCAAAAGGGGUUUU'
MFE = ['-10', '2', '5', '3', '6', '1', '2', '4.1', '((((&))))']
SUBOPT = ['-6', '10', '13', '11', '14', '0.5', '1.5', '4.1', '((((&))))']
PARTNER = ['-9', '3', '6', '4', '7', '2', '3', '4.1', '((((&))))']

# Executable fixtures exercise actual argv/stdin handling, including paths with
# spaces. Values are deliberately distinct for the two strands and ensembles.
FAKE_TOOL = r'''
import csv, json, os, sys
from pathlib import Path
kind = Path(sys.argv[0]).name
data = json.loads(Path(os.environ['FUSE_FIXTURE']).read_text())
stdin = sys.stdin.read() if kind == 'RNAfold' else ''
with open(os.environ['FUSE_LOG'], 'a') as log:
    log.write(json.dumps([kind, sys.argv[1:], stdin]) + '\n')
if data.get('fail') == kind:
    print('fixture failure', file=sys.stderr)
    sys.exit(7)
if kind == 'RNAfold':
    if data.get('bad_energy'):
        print(data['bad_energy'])
    else:
        target = stdin.splitlines()[0] == data['target']
        value = (-1 if target else -2) if '-C' in sys.argv else (-4 if target else -7)
        print(' free energy of ensemble = {} kcal/mol'.format(value))
else:
    columns = next(a.split('=', 1)[1] for a in sys.argv if a.startswith('--outCsvCols=')).split(',')
    rows = data['constrained'] if any(a.startswith('--tAccConstr=') for a in sys.argv) else data['rows']
    writer = csv.writer(sys.stdout, delimiter=';', lineterminator='\n')
    writer.writerow(columns)
    for row in rows:
        writer.writerow([row[c] for c in columns])
'''


def csv_rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter=';'))


class FuseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='intarna fuse ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.script = self.root / 'scripts' / SCRIPT.name
        self.script.parent.mkdir()
        shutil.copyfile(SCRIPT, self.script)
        shutil.copyfile(SCRIPT.parent / "intarna_tools.py", self.script.parent / "intarna_tools.py")
        self.bin = self.root / 'tools with spaces'
        self.bin.mkdir()
        for name in ('IntaRNA', 'RNAfold'):
            path = self.bin / name
            path.write_text('#!' + sys.executable + '\n' + FAKE_TOOL)
            path.chmod(0o755)
        self.fixture = self.root / 'fixture.json'
        self.log = self.root / 'calls.jsonl'
        self.env = dict(os.environ, PATH=str(self.bin), FUSE_FIXTURE=str(self.fixture),
                        FUSE_LOG=str(self.log))
        self.data = dict(target=TARGET, rows=[dict(zip(COLUMNS, row)) for row in (MFE, SUBOPT)],
                         constrained=[dict(zip(COLUMNS, PARTNER))])

    def run_script(self, *options, sequences=True, success=True):
        self.fixture.write_text(json.dumps(self.data))
        command = [sys.executable, str(self.script)]
        if sequences:
            command += ['-t', TARGET, '-q', QUERY]
        result = subprocess.run(command + list(options), env=self.env, cwd=self.root,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        return result

    def calls(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()]

    def test_fused_energy_constraints_and_user_sequences(self):
        result = self.run_script()
        rows = csv_rows(result.stdout)
        self.assertEqual(rows[:2], self.data['rows'])
        self.assertEqual(len(rows), 3)
        # (-6 - .5 - 1.5) + (-9 - 2 - 3) + 3 + 5 = -14.
        self.assertEqual(rows[2], dict(zip(COLUMNS,
            ['-14', '10:3', '13:6', '11:4', '14:7', '3', '5', '4.1:4.1',
             '((((&)))):((((&))))'])))
        calls = self.calls()
        intarna = [args for name, args, stdin in calls if name == 'IntaRNA']
        self.assertEqual(len(intarna), 2)
        for args in intarna:
            self.assertEqual(args[args.index('-t') + 1], TARGET)
            self.assertEqual(args[args.index('-q') + 1], QUERY)
            self.assertEqual(args[args.index('-n') + 1], '10')
            for option in ('--tAccW=0', '--qAccW=0', '--tAccL=100', '--qAccL=100', '--outOverlap=N'):
                self.assertIn(option, args)
        self.assertIn('--tAccConstr=b:10-13', intarna[1])
        self.assertIn('--qAccConstr=b:11-14', intarna[1])
        fold = [(args, stdin) for name, args, stdin in calls if name == 'RNAfold']
        self.assertEqual([stdin for args, stdin in fold],
                         [TARGET + '\n', QUERY + '\n', TARGET + '\n..xxxx...xxxx...\n',
                          QUERY + '\n...xxxx...xxxx..\n'])
        for args, stdin in fold:
            self.assertIn('--maxBPspan=100', args)
            self.assertIn('--noPS', args)
            self.assertIn('-p0', args)
            self.assertEqual('-C' in args, len(stdin.splitlines()) == 2)

    def test_reordered_columns_and_option_overrides(self):
        columns = list(reversed(COLUMNS))
        result = self.run_script('-o', ', '.join(columns), '-n', '3', '-l', '8')
        self.assertEqual(result.stdout.splitlines()[0], ';'.join(columns))
        self.assertEqual(csv_rows(result.stdout)[2]['E'], '-14')
        for name, args, stdin in self.calls():
            if name == 'RNAfold':
                self.assertIn('--maxBPspan=8', args)
            else:
                self.assertEqual(args[args.index('-n') + 1], '3')
                self.assertIn('--tAccL=8', args)
                self.assertIn('--qAccL=8', args)

    def test_empty_or_single_initial_prediction(self):
        for rows in ([], [self.data['rows'][0]]):
            with self.subTest(rows=rows):
                self.data['rows'] = rows
                result = self.run_script()
                self.assertEqual(csv_rows(result.stdout), rows)
                self.assertEqual(result.stdout.splitlines()[0], ';'.join(COLUMNS))
        self.assertTrue(all(name == 'IntaRNA' for name, args, stdin in self.calls()))

    def test_no_constrained_prediction_or_insufficient_overlap(self):
        partner = self.data['constrained'][0]
        for constrained in ([], [dict(partner, start1='5')],
                            [dict(partner, start2='6')], [dict(partner, start1='7', end1='9')]):
            with self.subTest(constrained=constrained):
                self.data['constrained'] = constrained
                self.assertEqual(csv_rows(self.run_script().stdout), self.data['rows'])
        self.assertFalse(any('-C' in args for name, args, stdin in self.calls()))

    def test_help_and_invalid_arguments_need_no_tools(self):
        shutil.rmtree(self.bin)
        self.assertIn('maximum', self.run_script('-h', sequences=False).stdout)
        for options in ([], ['-t', TARGET], ['-t', 'ACGU;echo', '-q', QUERY],
                        ['-t', TARGET, '-q', QUERY, '-n', '0'],
                        ['-t', TARGET, '-q', QUERY, '-l', '-1'],
                        ['-t', TARGET, '-q', QUERY, '-o', 'E'],
                        ['-t', TARGET, '-q', QUERY, '-o', ','.join(COLUMNS + ['E'])]):
            with self.subTest(options=options):
                self.assertIn('error:', self.run_script(*options, sequences=False, success=False).stderr)
        self.assertFalse(self.log.exists())

    def test_missing_executable(self):
        for name in ('IntaRNA', 'RNAfold'):
            with self.subTest(name=name):
                path = self.bin / name
                hidden = path.with_suffix('.hidden')
                path.rename(hidden)
                self.assertIn(name + ' executable not found', self.run_script(success=False).stderr)
                hidden.rename(path)

    def test_external_errors_publish_no_csv(self):
        for name in ('IntaRNA', 'RNAfold'):
            with self.subTest(name=name):
                self.data['fail'] = name
                result = self.run_script(success=False)
                self.assertIn('status 7: fixture failure', result.stderr)
                self.assertEqual(result.stdout, '')

    def test_invalid_energy_or_coordinates_publish_no_csv(self):
        for text in ('no ensemble energy', 'energy of ensemble = NaN kcal/mol'):
            self.data['bad_energy'] = text
            result = self.run_script(success=False)
            self.assertEqual(result.stdout, '')
            self.assertIn('energy', result.stderr)
        del self.data['bad_energy']
        for field, value in (('start1', '0'), ('end2', '17'), ('E', 'NaN'), ('ED1', 'invalid')):
            original = self.data['rows'][0][field]
            self.data['rows'][0][field] = value
            self.assertEqual(self.run_script(success=False).stdout, '')
            self.data['rows'][0][field] = original

    def test_installed_and_source_tree_lookup(self):
        for destination in (self.script.parent, self.root / 'src' / 'bin'):
            with self.subTest(destination=destination):
                destination.mkdir(parents=True, exist_ok=True)
                local = destination / 'IntaRNA'
                (self.bin / 'IntaRNA').rename(local)
                self.assertEqual(len(csv_rows(self.run_script().stdout)), 3)
                local.rename(self.bin / 'IntaRNA')


class NativeFuseTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('INTARNABINPATH') and shutil.which('RNAfold'),
                         'native check requires INTARNABINPATH and RNAfold on PATH')
    def test_oxys_fhla_energy_against_independent_tool_calls(self):
        # The legacy demonstration pair has compatible, non-overlapping sites.
        target = 'aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca'
        query = 'AGAAUAGCAAUGAACGAUUAUCCCUAUCAAGCAUUCUGACUGAUAAUUGCUCACAGAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGCUUUUUUUUGCCAUAAAAAA'
        binary = Path(os.environ['INTARNABINPATH']).resolve() / 'src' / 'bin' / 'IntaRNA'
        common = [str(binary), '-t', target, '-q', query, '--tAccW=0', '--qAccW=0',
                  '--tAccL=100', '--qAccL=100', '--outOverlap=N', '-n', '3',
                  '--outMode=C', '--outCsvCols=' + ','.join(COLUMNS)]
        original = csv_rows(subprocess.check_output(common, universal_newlines=True))
        with tempfile.TemporaryDirectory(prefix='native fuse ') as tmp:
            script = Path(tmp) / SCRIPT.name
            shutil.copyfile(SCRIPT, script)
            shutil.copyfile(SCRIPT.parent / "intarna_tools.py", script.parent / "intarna_tools.py")
            env = dict(os.environ, PATH=str(binary.parent) + os.pathsep + os.environ['PATH'])
            output = csv_rows(subprocess.check_output(
                [sys.executable, str(script), '-t', target, '-q', query, '-n', '3'],
                env=env, universal_newlines=True))
        self.assertEqual([row for row in output if ':' not in row['start1']], original)
        fused = [row for row in output if ':' in row['start1']]
        self.assertTrue(fused, 'example must exercise a fused interaction')
        for row in fused:
            candidate = next(r for r in original if r['start1'] == row['start1'].split(':')[0])
            partner = csv_rows(subprocess.check_output(common + [
                '--tAccConstr=b:{start1}-{end1}'.format(**candidate),
                '--qAccConstr=b:{start2}-{end2}'.format(**candidate)], universal_newlines=True))[0]
            for column in COLUMNS:
                if column not in ('E', 'ED1', 'ED2'):
                    self.assertEqual(row[column], candidate[column] + ':' + partner[column])
            costs = []
            for strand, seq in (('1', target), ('2', query)):
                starts = [int(x) for x in row['start' + strand].split(':')]
                ends = [int(x) for x in row['end' + strand].split(':')]
                mask = ''.join('x' if any(a <= p <= b for a, b in zip(starts, ends)) else '.'
                               for p in range(1, len(seq) + 1))
                values = []
                for constraint in (None, mask):
                    command = ['RNAfold', '--noPS', '-p0', '--maxBPspan=100']
                    text = seq + '\n'
                    if constraint:
                        command.append('-C')
                        text += constraint + '\n'
                    output = subprocess.check_output(command, input=text, universal_newlines=True)
                    values.append(float(re.search(r'energy of ensemble = (\S+) kcal', output)[1]))
                costs.append(values[1] - values[0])
                self.assertAlmostEqual(float(row['ED' + strand]), costs[-1], places=8)
            expected = sum(float(site['E']) - float(site['ED1']) - float(site['ED2'])
                           for site in (candidate, partner)) + sum(costs)
            self.assertAlmostEqual(float(row['E']), expected, places=8)


if __name__ == '__main__':
    unittest.main()
