cmake -DENABLE_MPI=off -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX
make faunus
make pyfaunus
make tests
make test
#make manual
make install

