#!/bin/bash
set -ev

cd /io/gmp-6.2.1
./configure --enable-fat
make
make check
make install

cd /io/
mkdir tests
cp bls_py/tests.py tests/

for PYV in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
    cd /io
    /opt/python/$PYV/bin/pip install Cython==0.29.33
    /opt/python/$PYV/bin/pip wheel . -w dist/
    auditwheel repair dist/python_bls-$PKG_VERSION-$PYV-linux*.whl --plat $PLAT -w dist/
    rm dist/python_bls-$PKG_VERSION-$PYV-linux*.whl

    # Tests
    /opt/python/$PYV/bin/pip install pytest
    /opt/python/$PYV/bin/pip install dist/python_bls-$PKG_VERSION-$PYV-$PLAT_TAGS.whl
    cd tests
    /opt/python/$PYV/bin/pytest -vvv tests.py
done

if [[ $PLAT == manylinux2014_x86_64 ]]; then
    cd /io
    /opt/python/cp37-cp37m/bin/python setup.py sdist
fi
