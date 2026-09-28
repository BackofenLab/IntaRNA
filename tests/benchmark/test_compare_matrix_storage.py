#!/usr/bin/env python3
"""Linux regression checks for the benchmark's process-group timeout."""
import importlib.util
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location(
    "benchmark", Path(__file__).with_name("compare-matrix-storage.py"))
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class BenchmarkProcessTests(unittest.TestCase):
    def assert_child_stopped(self, pid_file):
        self.assertTrue(pid_file.exists(), "child must start before testing cleanup")
        pid = int(pid_file.read_text())
        try:
            state = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0]
        except FileNotFoundError:
            state = "gone"
        self.assertIn(state, ("gone", "Z", "X"), "child is still running")

    def test_capture_exit_status_and_output(self):
        result = benchmark.run_timed(
            [sys.executable, "-c", "import sys; print('output'); print('error', file=sys.stderr); sys.exit(7)"],
            os.environ, timeout=5)
        self.assertEqual(result.returncode, 7)
        self.assertEqual(result.stdout, b"output\n")
        self.assertEqual(result.stderr, b"error\n")

    @unittest.skipUnless(sys.platform.startswith("linux"), "uses Linux /proc and GNU time")
    def test_timeout_kills_child_of_time_wrapper(self):
        with tempfile.TemporaryDirectory() as tmp:
            pid_file = Path(tmp) / "child.pid"
            child = ("import os,pathlib,time; pathlib.Path(" + repr(str(pid_file)) +
                     ").write_text(str(os.getpid())); print('started', flush=True); time.sleep(60)")
            started = time.monotonic()
            try:
                with self.assertRaises(subprocess.TimeoutExpired) as caught:
                    benchmark.run_timed(["/usr/bin/time", sys.executable, "-c", child],
                                        os.environ, timeout=2)
                self.assertLess(time.monotonic() - started, 10)
                self.assertIn(b"started", caught.exception.output)
                self.assert_child_stopped(pid_file)
            finally:
                # Keep the regression test safe even against a broken runner.
                if pid_file.exists():
                    try:
                        os.kill(int(pid_file.read_text()), signal.SIGKILL)
                    except ProcessLookupError:
                        pass

    @unittest.skipUnless(sys.platform.startswith("linux"), "uses Linux /proc and GNU time")
    def test_interrupt_kills_child_of_time_wrapper(self):
        with tempfile.TemporaryDirectory() as tmp:
            pid_file = Path(tmp) / "child.pid"
            child = ("import os,pathlib,time; pathlib.Path(" + repr(str(pid_file)) +
                     ").write_text(str(os.getpid())); time.sleep(60)")
            communicate = subprocess.Popen.communicate
            interrupted = False

            def interrupt_once(process, *args, **kwargs):
                nonlocal interrupted
                if not interrupted:
                    interrupted = True
                    try:
                        communicate(process, timeout=2)
                    except subprocess.TimeoutExpired:
                        raise KeyboardInterrupt from None
                return communicate(process, *args, **kwargs)

            started = time.monotonic()
            try:
                with patch.object(subprocess.Popen, "communicate", interrupt_once):
                    with self.assertRaises(KeyboardInterrupt):
                        benchmark.run_timed(["/usr/bin/time", sys.executable, "-c", child],
                                            os.environ, timeout=5)
                self.assertLess(time.monotonic() - started, 10)
                self.assert_child_stopped(pid_file)
            finally:
                if pid_file.exists():
                    try:
                        os.kill(int(pid_file.read_text()), signal.SIGKILL)
                    except ProcessLookupError:
                        pass


if __name__ == "__main__":
    unittest.main()
