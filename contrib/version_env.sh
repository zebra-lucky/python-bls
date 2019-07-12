#!/bin/bash


VERSION_LINE=(`grep 'version=' setup.py`)
if [[ $VERSION_LINE =~ version=\'([^\']+)\' ]]
then
    export PKG_VERSION=${BASH_REMATCH[1]};
    echo set PKG_VERSION to $PKG_VERSION
fi
