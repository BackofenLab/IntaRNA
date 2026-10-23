#!/usr/bin/env python3
"""Exercise archive downloads against a local HTTP server, without Emscripten."""
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from threading import Thread
import unittest


ARCHIVE = b"pinned dependency archive\n"
CHECKSUM = hashlib.sha256(ARCHIVE).hexdigest()
DOWNLOADER = Path(__file__).with_name("download.sh")


class DownloadHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.server.requests.append(self.path)
        if self.path == "/missing":
            self.send_error(404)
            return
        body = ARCHIVE if self.path == "/archive" else b"not an archive\n"
        self.send_response(200)
        self.send_header("Content-Length", str(len(body) + (10 if self.path == "/partial" else 0)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *_args):
        pass


class DownloadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), DownloadHandler)
        cls.thread = Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()

    def setUp(self):
        self.directory = TemporaryDirectory(prefix="intarna-download-")
        self.addCleanup(self.directory.cleanup)
        # Exercise quoted paths and the same .part/cache names used by builds.
        self.destination = Path(self.directory.name) / "dependency archive.tar.gz"
        self.server.requests = []

    def download(self, path="/archive"):
        url = f"http://127.0.0.1:{self.server.server_port}{path}"
        return subprocess.run(
            ["bash", str(DOWNLOADER), url, str(self.destination), CHECKSUM],
            capture_output=True, text=True, timeout=30,
        )

    def test_download_and_reuse_verified_cache(self):
        result = self.download()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.destination.read_bytes(), ARCHIVE)
        self.assertFalse(Path(str(self.destination) + ".part").exists())
        # A valid cache must work even when the remote archive is unavailable.
        result = self.download("/missing")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.server.requests, ["/archive"])

    def test_replace_invalid_cache_and_stale_partial(self):
        self.destination.write_bytes(b"old archive or cached HTTP error")
        Path(str(self.destination) + ".part").write_bytes(b"interrupted download")
        result = self.download()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.destination.read_bytes(), ARCHIVE)
        self.assertFalse(Path(str(self.destination) + ".part").exists())
        self.assertEqual(self.server.requests, ["/archive"])

    def test_failed_responses_are_not_cached_and_next_run_recovers(self):
        for path in ("/invalid", "/missing", "/partial"):
            with self.subTest(path=path):
                result = self.download(path)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.destination.exists())
                self.assertFalse(Path(str(self.destination) + ".part").exists())
        result = self.download()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.destination.read_bytes(), ARCHIVE)


if __name__ == "__main__":
    unittest.main()
