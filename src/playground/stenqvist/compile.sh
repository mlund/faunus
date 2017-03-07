#!/bin/sh

name=main
faunusV=faunus  # name of faunus-folder
path=$nb/faunus/src/playground/stenqvist/ # path to cpp-file

cd $path
mkdir CM
cp $name.cpp CM/.
cd CM
echo 'fau_example('$name' "./" '$name'.cpp)' > CMakeLists.txt

cd $nb/$faunusV
cmake . -DMYPLAYGROUND=$path/CM

cd $path/CM
make
cp $name ../.
rm $name $name.cpp
cd $path
