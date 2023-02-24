#!/bin/bash
set -ev

if [[ $PYARCH1 == 'win32' ]]; then
  echo win32 build
  export MPIR_DIR=mpir/lib/win32/Release
else
  echo win64 build
  export MPIR_DIR=mpir/lib/x64/Release
fi

python -m pip install --upgrade pip wheel cython==0.29.33
python setup.py build_ext -I $MPIR_DIR -L $MPIR_DIR
python setup.py bdist_wheel

mkdir tests && cp bls_py/tests.py tests/
cd tests
python -m pip install pytest
python -m pip install ../dist/python_bls-$PKG_VERSION-$PYVER1-$PYARCH1.whl
pytest -vvv tests.py
