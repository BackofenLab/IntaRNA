#!/usr/bin/env python3

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

import os
import sys
import subprocess
from shutil import rmtree, move


if __name__ == '__main__':
    build_dir = os.path.abspath(os.path.join(__file__, os.pardir))
    main = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir, 'intarnapvalue', '__main__.py'))

    subprocess.call( f'{sys.executable} -m PyInstaller -y -F "{main}"', shell=True )

    os.remove(os.path.join(build_dir, '__main__.spec'))
    rmtree(os.path.join(build_dir, 'build'))
    os.mkdir(os.path.join(build_dir, 'build'))
    if os.name == 'posix':
        move(os.path.join(build_dir, 'dist', '__main__'), os.path.join(build_dir, 'build', 'IntaRNApvalue'))
    else:
        move(os.path.join(build_dir, 'dist', '__main__.exe'), os.path.join(build_dir, 'build', 'IntaRNApvalue.exe'))
    rmtree(os.path.join(build_dir, 'dist'))
