CXX=clang++
CC=clang
CXX=clang++ CC=clang cmake -DENABLE_MPI=on -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX -DCMAKE_VERBOSE_MAKEFILE=on
make faunus
make tests
make test
#make manual
make install

