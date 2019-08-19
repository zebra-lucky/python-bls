#!/usr/bin/python3

import sys
from setuptools import setup
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False


description='BLS12-381 and Signatures in python'


long_description=('Implements the BLS12 curve and optimal ate pairing, '
                  'as well as BLS signatures and aggregation. Use for '
                  'reference / educational purposes only. '
                  'Based on reference implementation from Chia BLS Signatures')


if sys.platform in ('windows', 'win32'):
    MP_LIB = 'mpir'
else:
    MP_LIB = 'gmp'


ext = '.pyx' if USE_CYTHON else '.c'


extensions = [
    Extension('bls_py.fields_t_c', [f'extmod/bls_py/fields_t_c{ext}'],
              libraries=[MP_LIB]),
]


if USE_CYTHON:
    extensions = cythonize(extensions)


setup(
    name='python-bls',
    version='0.1.8',
    url='https://github.com/zebra-lucky/python-bls',
    license='Apache License 2.0',
    python_requires='>=3.6',
    packages=['bls_py'],
    ext_modules=extensions,
    platforms='any',
    description=description,
    long_description=long_description,
)
