#!/bin/bash
python_exe=$1
base=pyplusplus
export PATH=${base}/bin:$PATH
for version in 2.5 2.6
do
  if [ -e ${base}/lib/python${version} ]
  then
    ppath=${base}/lib/python${version}/site-packages
  fi
done
export PYTHONPATH=$ppath:$PYTHONPATH

if [ -x ${python_exe} ]
then
  rm generated/*.cpp generated/*.hpp generated/*.h generated/*.txt
  $python_exe generate.py
fi
