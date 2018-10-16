CXX=clang++
CC=clang
env
CXX=clang++ CC=clang cmake -DENABLE_MPI=on -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make faunus
make faunus_nopbc
make faunus_sphere
make tests
ctest -V
#make manual
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make pyfaunus
make install

