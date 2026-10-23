#!/usr/bin/env python3
"""Install IntaRNA's companion modules and command-line tools."""

from setuptools import setup

setup(
    name='intarnapvalue',
    version='0.1',
    description='IntaRNA energy significance, mutation selection and output tools',
    author='Fabio Gutmann',
    python_requires='>=3.10',
    url='https://github.com/BackofenLab/IntaRNA/',
    packages=['intarnapvalue', 'copomus'],
    py_modules=['intarna_tools'],
    scripts=['CopomuS.py', 'IntaRNA-fuse.py', 'IntaRNA_CSV_p-value.py', 'IntaRNA_plotRegions.py'],
    entry_points={'console_scripts': ['IntaRNApvalue=intarnapvalue.intarna_pvalue:main']},
    install_requires=['numpy>=1.23', 'scipy>=1.9'],
    extras_require={'plot': ['matplotlib>=3.6']},
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
    ],
)
