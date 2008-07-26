#!/bin/bash

#PBS -l nodes=1
#PBS -m ae
#PBS -M mikael.lund@teokem.lu.se
#PBS -l walltime=168:00:00
#cd $PBS_O_WORKDIR

biosim=../../
source ${biosim}/src/prot.sh

setjobid ".a"
cell_r=100
maxsep=99 
protein1="H-chain.aam"
protein2="lysozyme.aam"
hamaker=25000
vmdsupport="no"

adjust_dp="no"
prot_dp=3.0
prot_rot=0.5
ion_dp=100
titrate="yes"

for ph in 7.6 #10.95 7 9 10
do
  for salt in 42 #253 42 253
  do
    cleanup
    nion1=$[salt+39]
    nion2=$[salt]  
    micro=60
    smear="yes"
    simulate 
    exit

    output="lys-h.ph${ph}.nacl${salt}.vdw${hamaker}"
    micro=400000
    smear="yes"
    simulate > $output
    cleanup
  done
done
