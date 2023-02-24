#!/bin/bash
set -ev

docker pull quay.io/pypa/$PLAT

wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz.sig
gpg --keyserver hkps://keyserver.ubuntu.com \
  --recv-keys 343C2FF0FBEE5EC2EDBEF399F3599FF828C67298
gpg --verify gmp-6.2.1.tar.lz.sig
tar --lzip -xf gmp-6.2.1.tar.lz
