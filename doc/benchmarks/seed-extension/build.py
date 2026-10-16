#!/usr/bin/env python3
"""Build isolated baseline/current seeded kernels using a configured source tree."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="configured IntaRNA source checkout")
    parser.add_argument("--output", type=Path, required=True,
                        help="directory for source snapshots and executables")
    parser.add_argument("--cxx", default="g++", help="C++23 compiler executable")
    parser.add_argument("--baseline", default="985ba44bdc16e92162e954b7baffa3f3923d4776")
    parser.add_argument("--include", type=Path, action="append", default=[],
                        help="extra dependency include directory; repeat as needed")
    parser.add_argument("--runtime-library-dir", type=Path, action="append", default=[],
                        help="optional runtime library directory to embed as rpath")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    driver = Path(__file__).with_name("driver.cpp").resolve()
    (output / "driver.cpp").write_bytes(driver.read_bytes())
    for name in ("build.py", "run.py"):
        (output / name).write_bytes(Path(__file__).with_name(name).read_bytes())
    names = ("SeededPartitionFunction.cpp", "SeededPartitionFunction.h",
             "StackSeedDomain.cpp", "StackSeedDomain.h")
    for variant in ("base", "candidate"):
        (output / variant / "IntaRNA").mkdir(parents=True, exist_ok=True)
        for name in names:
            content = (subprocess.check_output(
                ["git", "-C", str(source), "show", args.baseline + ":src/IntaRNA/" + name])
                if variant == "base" else (source / "src" / "IntaRNA" / name).read_bytes())
            (output / variant / "IntaRNA" / name).write_bytes(content)
    commands = []
    for variant in ("base", "candidate"):
        command = [args.cxx, "-std=c++23", "-O3", "-DNDEBUG",
                   "-ffunction-sections", "-fdata-sections",
                   "-I" + str(output / variant), "-I" + str(source / "src"),
                   *("-I" + str(p.resolve()) for p in args.include),
                   str(output / "driver.cpp"),
                   str(output / variant / "IntaRNA/SeededPartitionFunction.cpp"),
                   str(output / variant / "IntaRNA/StackSeedDomain.cpp"), "-Wl,--gc-sections",
                   *("-Wl,-rpath," + str(p.resolve()) for p in args.runtime_library_dir),
                   "-o", str(output / variant / "benchmark")]
        subprocess.run(command, check=True)
        commands.append(command)
    files = [output / "driver.cpp", output / "build.py", output / "run.py",
             *(output / variant / "IntaRNA" / name
               for variant in ("base", "candidate") for name in names)]
    manifest = {
        "baseline": subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", args.baseline], text=True).strip(),
        "compiler": subprocess.check_output([args.cxx, "--version"], text=True).splitlines()[0],
        "commands": commands,
        "sha256": {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in files},
        "shared_header_sha256": {
            name: hashlib.sha256((source / "src" / "IntaRNA" / name).read_bytes()).hexdigest()
            for name in ("PartitionArithmetic.h", "Matrix.h", "Interaction.h",
                         "general.h", "intarna_config.h")},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
