"""CopomuS regressions for coordinates, mutation selection, profiles and ranking."""

import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from CopomuS import CopomuS, main
from copomus.IntaRNA import IntaRNA
from copomus.candidate_filters import get_filter, filter_lp_mfe
from copomus.candidate_selectors import get_mutation_enc
from copomus.indexing import array_index_to_idx, idx_to_array_index, add_sub_idx
from copomus.measures import MinDeltaEnergy, MFECover, Measure
from copomus.mutation import Mutation, parse_mutation_encoding
from copomus.mutation_generators import gen_specific


class CopomuSTests(unittest.TestCase):
    def test_signed_coordinates_round_trip(self):
        for offset, expected in ((-3, [-3, -2, -1, 1, 2]), (0, [0, 1, 2, 3, 4]),
                                 (2, [2, 3, 4, 5, 6])):
            for i, value in enumerate(expected):
                self.assertEqual(array_index_to_idx(i, offset), value)
                self.assertEqual(idx_to_array_index(value, offset), i)
        self.assertEqual(add_sub_idx(-1, 1), 1)
        self.assertEqual(add_sub_idx(1, -1), -1)
        with self.assertRaises(ValueError):
            idx_to_array_index(0, -1)

    def test_zero_energies_ids_and_quoted_csv(self):
        rows = IntaRNA.output_to_dict('id1;id2;start1;E;ED1;other\n"a;b";001;0;0.0;NA;\n')
        self.assertEqual(rows, [{'id1': 'a;b', 'id2': '001', 'start1': 0, 'E': 0.0,
                                 'ED1': None, 'other': None}])
        for empty in ('id1\n', 'E\n', 'E;ED1\n\n'):
            self.assertIs(IntaRNA.output_to_dict(empty), False)
        for bad in ('E;E\n1;2\n', 'E;ED1\n1\n', 'E\nnan\n', 'E\ninf\n'):
            with self.assertRaises(ValueError):
                IntaRNA.output_to_dict(bad)

    def test_version_and_no_empty_subprocess_arguments(self):
        self.assertEqual(IntaRNA.version_tuple('4.0.0-dev'), (4, 0, 0))
        with patch('copomus.IntaRNA.find_program', return_value='/tmp/with spaces/IntaRNA'), \
                patch('copomus.IntaRNA.run_command', return_value='E\n0\n') as run:
            self.assertEqual(IntaRNA().run('GC', 'GC', 1, 1, 'E', 0), {'E': 0})
            self.assertNotIn('', run.call_args.args[0])
            IntaRNA().run('GC', 'GC', 1, 1, 'E', 0, param_file='with spaces.txt')
            self.assertIn('--parameterFile=with spaces.txt', run.call_args.args[0])

    def test_explicit_mutations_and_bounds(self):
        args = ('GC', 'GC', 1, 1)
        self.assertEqual(get_mutation_enc(*args, '', 2, enc='G1C&C2G'), [(1, 2)])
        mutations = gen_specific(*args, ['E'], 1, 1, 'G1C&C2G')
        self.assertEqual((mutations[0].qm, mutations[0].tm), ('CC', 'GG'))
        for bad in ('G0C&C2G', 'C-1G&C2G', 'G3C&C2G', 'G1N&C2G', 'G1G&C2G',
                    'A1C&C2G', 'G1C&C2G,', 'bad'):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                parse_mutation_encoding(bad, *args)
        with self.assertRaises(ValueError):
            Mutation([], 'GC', 'GC', 1, 1, -1, 2, 'CG')
        mut = gen_specific('AGC', 'GCU', -2, -2, [], 1, 1, 'A-2U&U1A')[0]
        self.assertEqual((mut.qm, mut.tm), ('UGC', 'GCA'))

    def test_filters_across_zero_and_gc_alias(self):
        self.assertIs(get_filter('GC'), get_filter('CG'))
        self.assertEqual(get_filter('CG')([(1, 2)], 'GC', 'GC', 1, 1), [])
        pairs = [(-1, 1), (1, -1)]
        self.assertEqual(filter_lp_mfe(pairs, 'GC', 'GC', -1, -1), pairs)
        zero_pairs = [(0, 1), (1, 0)]
        self.assertEqual(filter_lp_mfe(zero_pairs, 'GC', 'GC', 0, 0), zero_pairs)

    def test_zero_differences_total_order_and_missing_scores(self):
        measures = []
        for values in ((-3, 0, 0, -3), (-3, -3, -3, -3), (-3, -4, -4, -3)):
            measure = MinDeltaEnergy()
            measure.score.update(zip(('ww', 'wm', 'mw', 'mm'), values))
            measure.post_compute()
            measures.append(measure)
        self.assertEqual([m.score['min'] for m in measures], [3, 0, -1])
        self.assertEqual(sorted(measures[::-1]), measures)
        missing = MinDeltaEnergy()
        missing.post_compute()
        self.assertTrue(missing.non_i())
        cover = MFECover()
        cover.calculate_score(False, 'wm')
        self.assertIs(cover.score['wm'], False)
        measure = Measure()
        measure.score.update(ww=-5, mm=-2, wm=-1.5, mw=-1.5)
        self.assertTrue(measure.is_valid(1, 0))

    def test_profiles_shared_and_wildtype_reused(self):
        app = CopomuS()
        app.measures = ['EDqi', 'cEDqi', 'Eqi', 'Eti']
        muts = gen_specific('GC', 'GC', 1, 1, app.measures, 1, 1, 'G1C&C2G,C2G&G1C')
        paths_seen = []
        def prediction(*args):
            for option in args[-1]:
                pipe, filename = option.removeprefix('--out=').split(':', 1)
                path = Path(filename)
                paths_seen.append(path)
                path.write_text('# header\n# columns\n1\t0\n2\t1\n' if pipe == 'qAcc'
                                else 'idx;nt;E\n1;G;-2\n2;C;-3\n')
            return {'id1': 'target'}
        with patch.object(IntaRNA, 'run', side_effect=prediction) as run:
            app.run_tests(muts)
        self.assertEqual(run.call_count, 7)  # one wildtype + three per mutation
        self.assertTrue(all(not path.exists() for path in paths_seen))
        for mut in muts:
            self.assertEqual(mut.measures[0].score, mut.measures[1].score | {'is_valid': 0})
            self.assertIsNotNone(mut.measures[0].score['ww'])

    def test_cli_validation_and_failed_output_preservation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'results.tsv'
            output.write_text('previous results')
            for extra in (['-q', '', '-t', 'GC'], ['-q', 'GC', '-t', 'GC', '--threads=-1'],
                          ['-q', 'GC', '-t', 'GC', '-d', ';;']):
                with patch('sys.stderr', new_callable=io.StringIO), self.assertRaises(SystemExit):
                    CopomuS().parse_cmd_args([*extra, '-o', str(output)])
                self.assertEqual(output.read_text(), 'previous results')
            with patch.object(IntaRNA, 'check_binary'), \
                    patch.object(IntaRNA, 'run', side_effect=RuntimeError('prediction failed')), \
                    patch('sys.stderr', new_callable=io.StringIO):
                self.assertEqual(main(['-q', 'GC', '-t', 'GC', '--mutationEncoding', 'G1C&C2G',
                                       '-o', str(output)]), 1)
            self.assertEqual(output.read_text(), 'previous results')

    def test_cli_explicit_mutation_zero_energy_csv(self):
        stdout = io.StringIO()
        with patch.object(IntaRNA, 'check_binary'), patch.object(IntaRNA, 'run', return_value={'E': 0}), \
                patch('sys.stdout', stdout):
            self.assertEqual(main(['-q', 'GC', '-t', 'GC', '--mutationEncoding', 'G1C&C2G',
                                   '-m', 'minDeltaE']), 0)
        row = list(csv.DictReader(io.StringIO(stdout.getvalue()), delimiter='\t'))[0]
        self.assertEqual(row['mutation'], 'G1C&C2G')
        self.assertEqual(row['minDeltaE_min'], '0')
        self.assertNotIn('\r', stdout.getvalue())


if __name__ == '__main__':
    unittest.main()
