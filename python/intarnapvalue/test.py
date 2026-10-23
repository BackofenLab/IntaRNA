#!/usr/bin/env python3
# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

"""Compatibility entry point for the source-tree p-value regression tests."""

from pathlib import Path
import unittest

if __name__ == '__main__':
    tests = Path(__file__).resolve().parents[1] / 'tests'
    suite = unittest.defaultTestLoader.discover(str(tests), pattern='test_intarnapvalue.py')
    raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
