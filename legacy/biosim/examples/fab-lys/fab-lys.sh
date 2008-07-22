#!/bin/bash

#PBS -l nodes=1
#PBS -m ae
#PBS -M mikael.lund@teokem.lu.se
#PBS -l walltime=168:00:00
#cd $PBS_O_WORKDIR

biosim=../../
source ${biosim}/src.bearnix/prot.sh

setjobid "fablys."
cell_r=100
minsep=30
maxsep=100 
protein1="../../struct/fab_double.aam"
protein2="../../struct/lysozyme.aam"
hamaker=25000
adjust_dp="no"
prot_dp=7.5
prot_rot=0.9
ion_dp=150
titrate="yes"
for ph in 7.6 #10.95 7 9 10
do
  for salt in 253
  do
    cleanup
    nion1=$[salt+72]
    nion2=$[salt]
    chion1=1
    chion2=-1 
#    micro=10000
#    save_coord="yes"
#    save_min=53.5
#    save_max=54.5
#    simulate  

    save_coord="no"
    adjust_dp="no"
    ion_dp=170
    prot_dp=7.5
    prot_rot=1.0
    minsnapshot="yes"
    output="../fablys/Wholefab-lys.100mM.vdw${hamaker}"
    micro=100000
    simulate > "${output}equlibration"
    micro=1500000
    simulate > $output
#    cleanup
  done
done
