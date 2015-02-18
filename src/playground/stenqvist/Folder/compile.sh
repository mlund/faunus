mv ../nemo.cpp ../nemo.cpp.bak
cp nemo.cpp ../.
cd ../..
make -j2 stenqvist-nemo
cd stenqvist
mv nemo Folder/.
mv nemo.cpp.bak nemo.cpp



