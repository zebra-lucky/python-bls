#!/bin/bash
set -ev


PYPKG_NAME=python-$PYTHON_VERSION-macosx10.6.pkg
case $PYTHON_VERSION in
  3.6.8)
    PY_SHA256=3c5fd87a231eca3ee138b0cdc2be6517a7ca428304d41901a86b51c6a22b910c
    ;;
  3.7.4)
    PY_SHA256=061c4006efd6374720613a7af1f4b66352900aedb7a374c63d300714ce1835ee
    ;;
esac
echo "$PY_SHA256  $PYPKG_NAME" > $PYPKG_NAME.sha256

PYFTP=https://www.python.org/ftp/python/$PYTHON_VERSION
curl -O $PYFTP/$PYPKG_NAME
shasum -a256 -s -c $PYPKG_NAME.sha256
sudo installer -pkg $PYPKG_NAME -target /
rm $PYPKG_NAME $PYPKG_NAME.sha256

mkdir gmp
cp /usr/local/Cellar/gmp/*/lib/libgmp.a gmp/
cp /usr/local/Cellar/gmp/*/include/gmp.h gmp/
