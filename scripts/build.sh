CXX=clang++
CC=clang
env
CXX=clang++ CC=clang cmake -DENABLE_MPI=on -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make faunus
make tests
make test
#make manual
make install

