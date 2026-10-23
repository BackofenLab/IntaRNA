"""Contracts shared by standalone and installed Python companion commands."""

import importlib.util
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from intarna_tools import atomic_text_output, find_program, run_command


class CommonToolsTests(unittest.TestCase):
    def test_atomic_publication_and_rollback(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.csv'
            output.write_text('old')
            output.chmod(0o640)
            with self.assertRaises(RuntimeError):
                with atomic_text_output(output) as handle:
                    handle.write('partial')
                    raise RuntimeError('failed computation')
            self.assertEqual(output.read_text(), 'old')
            with atomic_text_output(output) as handle:
                handle.write('α\n')
            self.assertEqual(output.read_bytes(), 'α\n'.encode())
            self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o640)
            self.assertEqual(list(Path(directory).iterdir()), [output])

    def test_subprocess_arguments_stdin_and_failure(self):
        payload = '$(touch SHOULD_NOT_EXIST); quoted spaces'
        result = run_command([sys.executable, '-c', 'import sys; print(sys.argv[1]); print(sys.stdin.read())', payload], 'ACGU')
        self.assertEqual(result, payload + '\nACGU\n')
        with self.assertRaisesRegex(RuntimeError, 'status 7: specific failure'):
            run_command([sys.executable, '-c', 'import sys; print("specific failure", file=sys.stderr); sys.exit(7)'])

    def test_lookup_siblings_source_tree_path_and_frozen_executable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = root / 'python' / 'script.py'
            source = root / 'src' / 'bin' / 'IntaRNA'
            sibling = script.parent / 'IntaRNA'
            source.parent.mkdir(parents=True)
            sibling.parent.mkdir()
            for path in (source, sibling):
                path.write_text('#!/bin/sh\nexit 0\n')
                path.chmod(0o755)
            with patch.dict('os.environ', {'PATH': ''}):
                self.assertEqual(find_program('IntaRNA', script), str(sibling))
                sibling.unlink()
                self.assertEqual(find_program('IntaRNA', script), str(source))
                source.unlink()
                with self.assertRaises(RuntimeError):
                    find_program('IntaRNA', script)
            sibling.write_text('#!/bin/sh\nexit 0\n')
            sibling.chmod(0o755)
            with patch.object(sys, 'frozen', True, create=True), \
                    patch.object(sys, 'executable', str(sibling.parent / 'IntaRNApvalue')):
                self.assertEqual(find_program('IntaRNA'), str(sibling))

    def test_optional_builder_uses_paths_and_propagates_failures(self):
        spec = importlib.util.spec_from_file_location('companion_builder', ROOT / 'bin' / 'build.py')
        builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)
        with patch.object(builder.subprocess, 'run') as run:
            self.assertEqual(builder.main(), 0)
            command = run.call_args.args[0]
            self.assertEqual(command[:3], [sys.executable, '-m', 'PyInstaller'])
            self.assertEqual(command[command.index('--paths') + 1], str(ROOT))
            self.assertTrue(run.call_args.kwargs['check'])
            temporary = Path(command[command.index('--workpath') + 1])
        self.assertFalse(temporary.exists())
        with patch.object(builder.subprocess, 'run', side_effect=subprocess.CalledProcessError(7, 'PyInstaller')), \
                self.assertRaises(subprocess.CalledProcessError):
            builder.main()

    def test_all_standalone_help_without_optional_dependencies(self):
        for name in ('IntaRNA-fuse.py', 'CopomuS.py', 'IntaRNA_CSV_p-value.py', 'IntaRNA_plotRegions.py'):
            result = subprocess.run([sys.executable, '-S', str(ROOT / name), '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('usage:', result.stdout)


if __name__ == '__main__':
    unittest.main()
