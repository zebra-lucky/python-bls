#!/bin/bash
set -ev

wget http://www.tortall.net/projects/yasm/releases/vsyasm-1.3.0-win64.zip
mkdir "/c/Program Files/yasm"
unzip -d "/c/Program Files/yasm" vsyasm-1.3.0-win64.zip
rm vsyasm-1.3.0-win64.zip
export PATH="/c/Program Files/yasm:$PATH"

#export PATH="/c/Program Files (x86)/Microsoft Visual Studio/2019/Enterprise/MSBuild/Current/Bin/:$PATH"
git clone https://github.com/BrianGladman/mpir/
cd mpir
echo "1" | python ./msvc/mpir_config.py 15
cd ./msvc/vs19
if [[ $PYARCH1 == 'win32' ]]; then
  ./msbuild.bat gc lib Win32 Release
else
  ./msbuild.bat gc lib x64 Release
fi
