#!/bin/bash
set -ev

docker run --rm -e PLAT=$PLAT -e PKG_VERSION=$PKG_VERSION \
    -v `pwd`:/io quay.io/pypa/$PLAT \
    $PRE_CMD /io/contrib/build_wheels_linux.sh
