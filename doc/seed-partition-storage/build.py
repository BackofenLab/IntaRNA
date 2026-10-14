#!/usr/bin/env python3
"""Build the benchmark against a chosen real SeededPartitionFunction.cpp.

The source root must already be configured and built by Autotools. Only the
driver and chosen kernel translation unit need compilation; the remaining
dependencies come from the built IntaRNA library and logging object.
Example:

  python3 build.py --root ../.. --kernel ../../src/IntaRNA/SeededPartitionFunction.cpp \\
    --output /tmp/dense

Repeat with the patched kernel for /tmp/hash. The compiler, compile options and
dependency flags are read from the configured Makefile, ensuring both variants
match the library build. The command is printed and saved beside the executable.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess


def configured_flags(root):
    names = ("CXX", "DEFS", "CPPFLAGS", "AM_CXXFLAGS", "CXXFLAGS", "AM_LDFLAGS", "LDFLAGS", "LIBS")
    # Make expands its own recursive variables. The additional phony target only
    # prints configuration; it does not build or modify any repository target.
    recipe = ".PHONY: seed-storage-build-flags\nseed-storage-build-flags:\n"
    recipe += "".join("\t$(info " + name + "=$(" + name + "))\n" for name in names)
    recipe += "\t@:\n"
    result = subprocess.run(["make", "-s", "--no-print-directory", "-C", str(root / "src/bin"),
                             "-f", "Makefile", "-f", "-", "seed-storage-build-flags"],
                            input=recipe, text=True, capture_output=True, check=True)
    values = dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line)
    return {name: shlex.split(values[name]) for name in names}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--kernel", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, kernel, output = args.root.resolve(), args.kernel.resolve(), args.output.resolve()
    if not (root / "src/IntaRNA/intarna_config.h").is_file():
        parser.error("--root must contain a configured IntaRNA source tree")
    for dependency in ("src/IntaRNA/libIntaRNA.a", "src/easylogging++.o"):
        if not (root / dependency).is_file():
            parser.error("build IntaRNA first: missing " + str(root / dependency))
    output.parent.mkdir(parents=True, exist_ok=True)
    configured = configured_flags(root)
    command = (configured["CXX"] + configured["DEFS"] + configured["CPPFLAGS"]
               + configured["AM_CXXFLAGS"] + configured["CXXFLAGS"]
               + ["-I" + str(root / "src"), "-I" + str(root),
                  str(Path(__file__).with_name("benchmark.cpp").resolve()),
                  str(kernel), str(root / "src/IntaRNA/libIntaRNA.a"),
                  str(root / "src/easylogging++.o"),
                  "-o", str(output)] + configured["AM_LDFLAGS"]
               + configured["LDFLAGS"] + configured["LIBS"])
    print(shlex.join(command), flush=True)
    subprocess.run(command, check=True)
    output.with_suffix(output.suffix + ".build.json").write_text(json.dumps({
        "command": command, "compiler_version": subprocess.check_output(
            configured["CXX"] + ["--version"], text=True),
        "kernel": str(kernel), "kernel_sha256": hashlib.sha256(kernel.read_bytes()).hexdigest(),
        "benchmark_sha256": hashlib.sha256(Path(__file__).with_name("benchmark.cpp").read_bytes()).hexdigest(),
        "configured_source_root": str(root)
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
