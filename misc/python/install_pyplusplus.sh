#!/bin/bash

# This script will download pygccxml + Py++
# and install them into "pyplusplus" in the current
# working directory.

python_exe=$1
version="1.0.0"
mirror="http://downloads.sourceforge.net/pygccxml"
pgx_arc="pygccxml-${version}.zip"
ppp_arc="pyplusplus-${version}.zip"
prefix=pyplusplus

for file in $pgx_arc $ppp_arc
do
  # Download archives if missing
  if [ ! -e $file ]
  then
    wget ${mirror}/${file}
  fi
  # Unpack and install
  if [ -e $file ]
  then
    unzip -q $file
    if [ -e "Py++-1.0.0" ]
    then
      mv Py++-1.0.0 `basename ${ppp_arc} .zip`
    fi
    cd `basename ${file} .zip`
    ${python_exe} setup.py install --prefix=../$prefix
    cd ..
    rm -fR `basename ${file} .zip`
  fi
done

unset python_exe
unset version
unset mirror
unset pgx_arc
unset ppp_arc
unset prefix


