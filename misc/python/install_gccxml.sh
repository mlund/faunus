#!/bin/bash

cvs=$1
prefix=${PWD}/../auxiliary

# checkout from cvs
if [ ! -e gccxml ]
then
  echo "Enter blank password when prompted."
  $cvs -d :pserver:anoncvs@www.gccxml.org:/cvsroot/GCC_XML login
  $cvs -d :pserver:anoncvs@www.gccxml.org:/cvsroot/GCC_XML co gccxml
fi

# build and install
if [ -e gccxml ]
then
  tmpdir=gccxml-build
  mkdir $tmpdir
  cd $tmpdir
  cmake ../gccxml -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
  make -j2
  make install
  cd ..
  rm -fR $tmpdir
fi
