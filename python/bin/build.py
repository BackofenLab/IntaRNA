#!/usr/bin/env python3
# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

"""Build IntaRNApvalue with PyInstaller without relying on the working directory."""

from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def main():
    directory = Path(__file__).resolve().parent
    source = directory.parent
    # Keep build metadata private; leave only the requested binary in bin/build.
    with TemporaryDirectory(prefix='IntaRNApvalue-build-') as temporary:
        subprocess.run([sys.executable, '-m', 'PyInstaller', '--noconfirm', '--onefile',
                        '--name', 'IntaRNApvalue', '--paths', str(source),
                        '--distpath', str(directory / 'build'), '--workpath', temporary,
                        '--specpath', temporary, str(source / 'intarnapvalue' / '__main__.py')],
                       check=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
