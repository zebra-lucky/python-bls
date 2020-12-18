#!/bin/bash
set -ev


GMP_BOTTLE=gmp-6.1.2_2.el_capitan.bottle.tar.gz
BINTRAY_URI=https://bintray.com/homebrew/bottles/
wget -O $GMP_BOTTLE ${BINTRAY_URI}download_file?file_path=${GMP_BOTTLE}
brew uninstall --ignore-dependencies gmp
brew install $GMP_BOTTLE


case $PYTHON_VERSION in
  3.6.8)
    PY_SHA256=3c5fd87a231eca3ee138b0cdc2be6517a7ca428304d41901a86b51c6a22b910c
    PYPKG_NAME=python-$PYTHON_VERSION-macosx10.6.pkg
    ;;
  3.7.6)
    PY_SHA256=1a595e7e6aca247b8229bd00a8e5fc8879fd85e18a2aa747ee7ecffb500cbfdd
    PYPKG_NAME=python-$PYTHON_VERSION-macosx10.6.pkg
    ;;
  3.8.6)
    PY_SHA256=c8195f8b056aff380f1e89dc7dd4c37372c37e5163b91c9391f4256bf5b44fe8
    PYPKG_NAME=python-$PYTHON_VERSION-macosx10.9.pkg
    ;;
esac
echo "$PY_SHA256  $PYPKG_NAME" > $PYPKG_NAME.sha256

PYFTP=https://www.python.org/ftp/python/$PYTHON_VERSION
echo "Downloading $PYFTP/$PYPKG_NAME"
curl -O $PYFTP/$PYPKG_NAME
shasum -a256 -s -c $PYPKG_NAME.sha256
sudo installer -pkg $PYPKG_NAME -target /
rm $PYPKG_NAME $PYPKG_NAME.sha256

mkdir gmp
cp /usr/local/Cellar/gmp/*/lib/libgmp.a gmp/
cp /usr/local/Cellar/gmp/*/include/gmp.h gmp/
