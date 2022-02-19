#/usr/bin/env bash

if ! command -v vmd &> /dev/null
then
    echo "'vmd' command not be found - using hardcoded (macos) path that you will want to change"
    vmd='/Applications/VMD 1.9.4a51-x86_64-Rev9.app/Contents/MacOS/startup.command'
else
    vmd="vmd"
fi

python2 ../../scripts/psc2vmd.py -i movie.dat -o movie.pdb --psf movie.psf --length 30
"${vmd}" -e vmd.script
