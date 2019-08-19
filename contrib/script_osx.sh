#!/bin/bash
set -ev

sudo pip3 install --upgrade pip wheel delocate==0.7.4 cython==0.29.13
python3 setup.py build_ext
python3 setup.py bdist_wheel

# Delocate and Tests
sudo pip3 install pytest
case $PYTHON_VERSION in
  3.6.8)
    export PATH=/Library/Frameworks/Python.framework/Versions/3.6/bin:$PATH
    delocate-wheel -v dist/python_bls-$PKG_VERSION-cp36-cp36m-macosx_10_6_intel.whl
    sudo pip3 install dist/python_bls-$PKG_VERSION-cp36-cp36m-macosx_10_6_intel.whl
    ;;
  3.7.4)
    export PATH=/Library/Frameworks/Python.framework/Versions/3.7/bin:$PATH
    delocate-wheel -v dist/python_bls-$PKG_VERSION-cp37-cp37m-macosx_10_6_intel.whl
    sudo pip3 install dist/python_bls-$PKG_VERSION-cp37-cp37m-macosx_10_6_intel.whl
    ;;
esac
mkdir tests
cp bls_py/tests.py tests/
cd tests
python3 -m pytest -vvv tests.py
