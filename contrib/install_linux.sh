#!/bin/bash
set -ev

docker pull quay.io/pypa/$PLAT

wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.lz
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.lz.sig
gpg --recv-keys 28C67298
gpg --verify gmp-6.1.2.tar.lz.sig
tar --lzip -xf gmp-6.1.2.tar.lz
