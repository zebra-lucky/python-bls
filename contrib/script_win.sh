#!/bin/bash
set -ev

PATH_BACKUP=$PATH

mkdir tests
cp bls_py/tests.py tests/

# build x64 wheels
export MPIR_DIR=mpir/lib/x64/Release

# build python3.6 x64 wheel
export PATH="/c/Python36:/c/Python36/Scripts:$PATH_BACKUP"
python -m pip install --upgrade pip wheel cython==0.29.13

rm -rf build
python setup.py build_ext -I $MPIR_DIR -L $MPIR_DIR
python setup.py bdist_wheel

# Tests
cd tests
python -m pip install pytest
python -m pip install ../dist/python_bls-$PKG_VERSION-cp36-cp36m-win_amd64.whl
pytest -vvv tests.py
cd ..

# build python3.7 x64 wheel
choco uninstall -y python3 --version 3.6.8
rm -rf /c/Python36/
choco install -y python3 --version 3.7.4

export PATH="/c/Python37:/c/Python37/Scripts:$PATH_BACKUP"
python -m pip install --upgrade pip wheel cython==0.29.13

rm -rf build
python setup.py build_ext -I $MPIR_DIR -L $MPIR_DIR
python setup.py bdist_wheel

# Tests
cd tests
python -m pip install pytest
python -m pip install ../dist/python_bls-$PKG_VERSION-cp37-cp37m-win_amd64.whl
pytest -vvv tests.py
cd ..

# build win32 wheels
export MPIR_DIR=mpir/lib/win32/Release

# build python3.6 x64 wheel
choco uninstall -y python3 --version 3.7.4
rm -rf /c/Python37/
choco install -y python3 --version 3.6.8 --x86

export PATH="/c/Python36:/c/Python36/Scripts:$PATH_BACKUP"
python -m pip install --upgrade pip wheel cython==0.29.13

rm -rf build
python setup.py build_ext -I $MPIR_DIR -L $MPIR_DIR
python setup.py bdist_wheel

# Tests
cd tests
python -m pip install pytest
python -m pip install ../dist/python_bls-$PKG_VERSION-cp36-cp36m-win32.whl
pytest -vvv tests.py
cd ..

# build python3.7 x64 wheel
choco uninstall -y python3 --version 3.6.8
rm -rf /c/Python36/
choco install -y python3 --version 3.7.4 --x86

export PATH="/c/Python37:/c/Python37/Scripts:$PATH_BACKUP"
python -m pip install --upgrade pip wheel cython==0.29.13

rm -rf build
python setup.py build_ext -I $MPIR_DIR -L $MPIR_DIR
python setup.py bdist_wheel

# Tests
cd tests
python -m pip install pytest
python -m pip install ../dist/python_bls-$PKG_VERSION-cp37-cp37m-win32.whl
pytest -vvv tests.py
cd ..
