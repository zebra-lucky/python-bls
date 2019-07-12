#!/bin/bash
set -ev

wget http://www.tortall.net/projects/yasm/releases/vsyasm-1.3.0-win64.zip
mkdir "/c/Program Files/yasm"
unzip -d "/c/Program Files/yasm" vsyasm-1.3.0-win64.zip
rm vsyasm-1.3.0-win64.zip
choco install -y vcredist2010
export PATH="/c/Program Files/yasm:$PATH"

choco install -y vcbuildtools
choco install -y python3 --version 3.6.8

export PATH="/c/Python36:/c/Python36/Scripts:$PATH"

git clone https://github.com/BrianGladman/mpir/
cd mpir
echo "1" | python ./msvc/mpir_config.py 15
cd ./msvc/vs15
./msbuild.bat gc lib x64 Release
./msbuild.bat gc lib Win32 Release
