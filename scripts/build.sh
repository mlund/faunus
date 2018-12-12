CXX=clang++
CC=clang
CXX=clang++ CC=clang cmake -DENABLE_MPI=on -DENABLE_OPENMP=on -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
alias make="make -s"
make faunus
#make faunus_nopbc
#make faunus_sphere
#make manual
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX # needed to build pyfaunus
make pyfaunus
make install
make tests
ctest -V
