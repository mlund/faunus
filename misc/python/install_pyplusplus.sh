#!/bin/bash

# This script will download pygccxml + Py++
# and install them into "pyplusplus" in an
# auxiliary directory.

python_exe=$1
wget_exe=$2
version="1.0.0"
mirror="http://downloads.sourceforge.net/pygccxml"
pgx_arc="pygccxml-${version}.zip"
ppp_arc="pyplusplus-${version}.zip"
prefix=${PWD}/../auxiliary

for file in $pgx_arc $ppp_arc
do
  # Download archives if missing
  if [ ! -e $file ]
  then
    ${wget_exe} ${mirror}/${file}
  fi
  # Unpack and install
  if [ -e $file ]
  then
    unzip -q $file
    if [ -e "Py++-${version}" ]
    then
      mv Py++-${version} `basename ${ppp_arc} .zip`
    fi
    cd `basename ${file} .zip`
    ${python_exe} setup.py install --prefix=${prefix}
    cd ..
    rm -fR `basename ${file} .zip`
  fi
done
