#!/bin/bash
#PBS -l nodes=1:toto7
#PBS -m ae
#PBS -M mikael.lund@teokem.lu.se
cd $PBS_O_WORKDIR

. use_modules
module add intel

