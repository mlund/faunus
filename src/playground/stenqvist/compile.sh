#!/bin/sh

faunusV=faunus  # name of faunus-folder
name=main # Name of cpp-file
path=$nb/$faunusV/src/playground/stenqvist/ # path to folder of cpp-file

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
