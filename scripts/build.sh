export CC=clang
export CXX=clang++
cmake -DENABLE_MPI=on -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make faunus
make manual
make install

