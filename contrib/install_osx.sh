#!/bin/bash
set -ev


GMP_BOTTLE="gmp-6.2.1.high_sierra.bottle.tar.gz"
TOKEN="QQ=="
GHCR_AUTH="Authorization: Bearer $TOKEN"
GHCR_URI="https://ghcr.io/v2/homebrew/core/gmp/blobs/sha256"
GHCR_HASH="54191ce7fa888df64b9c52870531ac0ce2e8cbd40a7c4cdec74cb2c4a421af97"
curl -o $GMP_BOTTLE -L -H "$GHCR_AUTH" "$GHCR_URI:$GHCR_HASH"
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
  3.9.13)
    PY_SHA256=167c4e2d9f172a617ba6f3b08783cf376dec429386378066eb2f865c98030dd7
    PYPKG_NAME=python-$PYTHON_VERSION-macosx10.9.pkg
    ;;
  3.10.10)
    PY_SHA256=1c24eb452065f285249f94965f916f6a5422ec88c131e53597ded4507aa627f7
    PYPKG_NAME=python-$PYTHON_VERSION-macos11.pkg
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
