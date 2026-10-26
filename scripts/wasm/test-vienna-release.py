#!/usr/bin/env python3
"""Check stable-major selection, pagination and archive identity without GitHub."""
from copy import deepcopy
import importlib.util
from io import BytesIO
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from urllib.error import HTTPError


SCRIPT = Path(__file__).with_name("vienna-release.py")
SPEC = importlib.util.spec_from_file_location("vienna_release", SCRIPT)
resolver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resolver)
DIGEST = "a" * 64
PYODIDE_VERSION = "314.0.7"
PYODIDE_FILES = ["pyodide.mjs", "pyodide.asm.mjs", "pyodide.asm.wasm",
                  "pyodide-lock.json", "python_stdlib.zip"]


def release(version, **fields):
    return {
        "tag_name": f"v{version}", "draft": False, "prerelease": False,
        "assets": [{
            "name": f"ViennaRNA-{version}.tar.gz", "digest": f"sha256:{DIGEST}",
            "browser_download_url": f"{resolver.REPOSITORY}/releases/download/v{version}/ViennaRNA-{version}.tar.gz",
        }], **fields,
    }


class ReleaseTests(unittest.TestCase):
    def resolve(self, pages):
        with patch.object(resolver, "urlopen", side_effect=[BytesIO(json.dumps(page).encode()) for page in pages]) as fetch:
            result = resolver.latest_release()
            return result, fetch

    def test_numeric_selection_excludes_major3_drafts_and_prereleases(self):
        result, _ = self.resolve([[
            release("3.0.0"), release("2.9.0"), release("2.10.0"),
            release("2.11.0", prerelease=True), release("2.12.0", draft=True),
            release("2.13.0-rc.1"), release("2.7.2"), release("2.7.1"),
        ]])
        self.assertEqual(result["version"], "2.10.0")
        self.assertEqual(result["sha256"], DIGEST)

    def test_pagination_finds_version2_after_one_hundred_newer_releases(self):
        result, fetch = self.resolve([[release(f"3.0.{n}") for n in range(100)], [release("2.7.2")]])
        self.assertEqual(result["version"], "2.7.2")
        self.assertEqual(fetch.call_count, 2)
        self.assertTrue(fetch.call_args.args[0].full_url.endswith("per_page=100&page=2"))

    def test_newest_release_must_have_exact_asset_url_and_digest(self):
        latest = release("2.8.0")
        for field, value in [("name", "source.tar.gz"), ("digest", None), ("digest", "sha256:bad"), ("digest", "md5:" + DIGEST), ("browser_download_url", "https://example.com/archive.tar.gz")]:
            with self.subTest(field=field, value=value):
                invalid = deepcopy(latest)
                invalid["assets"][0][field] = value
                with self.assertRaises(ValueError):
                    self.resolve([[release("2.7.2"), invalid]])
        with self.assertRaises(ValueError):
            self.resolve([[release("2.8.0", assets=latest["assets"] * 2)]])

    def test_no_supported_release_or_api_failure_is_an_error(self):
        with self.assertRaises(ValueError):
            self.resolve([[release("3.0.0"), release("2.7.1"), release("2.8.0", prerelease=True)]])
        with patch.object(resolver, "urlopen", side_effect=HTTPError(resolver.API, 403, "rate limited", {}, None)):
            with self.assertRaises(HTTPError):
                resolver.latest_release()

    def test_checked_in_pin_and_shell_fields(self):
        manifest = json.loads(subprocess.check_output(["python3", str(SCRIPT)], text=True))
        self.assertEqual(manifest["version"], "2.7.2")
        fields = subprocess.check_output(["python3", str(SCRIPT), "--fields"], text=True).split()
        self.assertEqual(fields, list(manifest.values()))

    def test_invalid_pin_fails_without_emitting_shell_fields(self):
        manifest = json.loads(SCRIPT.with_suffix(".json").read_text())
        with TemporaryDirectory(prefix="vienna manifest ") as directory:
            path = Path(directory) / "release.json"
            for field, value in [("version", "3.0.0"), ("version", "2.7.0"), ("url", "https://example.com"), ("sha256", "bad")]:
                with self.subTest(field=field):
                    path.write_text(json.dumps({**manifest, field: value}))
                    result = subprocess.run(["python3", str(SCRIPT), "--manifest", str(path), "--fields"], capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertEqual(result.stdout, "")

    def test_provenance_rejects_stale_engine_and_records_selected_archive(self):
        manifest = json.loads(SCRIPT.with_suffix(".json").read_text())
        with TemporaryDirectory(prefix="vienna engine ") as directory:
            site = Path(directory)
            assets = site / "wasm"
            assets.mkdir()
            pyodide = site / "vendor" / "pyodide"
            pyodide.mkdir(parents=True)
            for name in PYODIDE_FILES:
                (pyodide / name).write_bytes(f"test {name}".encode())
            (assets / "package.json").write_text('{"type":"module"}')
            (assets / "intarna.wasm").write_bytes(b"test engine")
            for version in ["2.7.0", manifest["version"]]:
                (assets / "intarna.js").write_text(
                    "export default async ({print}) => ({callMain() {"
                    f"print('IntaRNA 3.4.1 using Vienna RNA package {version} and boost 1.87.0');"
                    "return 0;}});"
                )
                result = subprocess.run([
                    "node", str(SCRIPT.with_name("write-build-info.mjs")),
                    str(SCRIPT.parents[2]), str(site), "4.0.23", "test-commit", "",
                    manifest["version"], manifest["url"], manifest["sha256"],
                    PYODIDE_VERSION,
                ], capture_output=True, text=True)
                if version == "2.7.0":
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("differs from selected", result.stderr)
                    self.assertFalse((site / "build-info.json").exists())
                else:
                    self.assertEqual(result.returncode, 0, result.stderr)
                    metadata = json.loads((site / "build-info.json").read_text())
                    self.assertEqual(metadata["viennarna"], version)
                    self.assertEqual(metadata["viennarnaArchive"], manifest)
                    self.assertEqual(metadata["pyodide"], PYODIDE_VERSION)
                    self.assertEqual(set(metadata["pyodideSha256"]), set(PYODIDE_FILES))
                    self.assertIn(manifest["sha256"], (site / "build-info.txt").read_text())


if __name__ == "__main__":
    unittest.main()
