#!/usr/bin/python3

from setuptools import setup


description='BLS12-381 and Signatures in python'

long_description=('Implements the BLS12 curve and optimal ate pairing, '
                  'as well as BLS signatures and aggregation. Use for '
                  'reference / educational purposes only. '
                  'Based on reference implementation from Chia BLS Signatures')


setup(
    name='python-bls',
    version='0.1.5',
    url='https://github.com/zebra-lucky/python-bls',
    license='Apache License 2.0',
    install_requires=[
        'typing; python_version < "3.5"',
    ],
    python_requires='>=3.4',
    packages=['bls_py'],
    platforms='any',
    description=description,
    long_description=long_description,
)
