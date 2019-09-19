#!/usr/bin/env python3


from setuptools import setup

setup(
    name='intarnapvalue',
    version='0.1',
    description='Calculates p-values of IntaRNA energy scores',
    author='Fabio Gutmann',
    python_requires='>=3.6.0',
    url='https://github.com/BackofenLab/IntaRNA/',
    packages=['intarnapvalue'],
#    install_requires=['scipy', 'numpy'],
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ]
)
