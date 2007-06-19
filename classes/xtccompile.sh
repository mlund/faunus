#!/bin/bash

incdir="/usr/local/gromacs/include/gromacs"
libdir="/usr/local/gromacs/lib/"

g++ -I${incdir} -L${libdir} xtc.C slump.C io.C point.C species.C -lgmx
