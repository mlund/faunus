mv ../nemo.cpp ../nemo.cpp.bak.o987247
cp nemo.cpp ../.
cd ../..
make -j2 stenqvist-nemo
cd stenqvist
mv nemo Quad/.
mv nemo.cpp.bak.o987247 nemo.cpp



