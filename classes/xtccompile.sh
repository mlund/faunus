#!/bin/bash

incdir="/usr/local/gromacs/include/gromacs"
libdir="/usr/local/gromacs/lib/"

g++ -D CPLUSPLUS -I${incdir} -L${libdir} -lgmx xtc.C
