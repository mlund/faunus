#!/bin/bash
#PBS -l walltime=48:00:00

biosim=/mnt/data/mikael/mc-temp/biosim
source ${biosim}/src/pka.sh
exe=${biosim}/src/pka

cell_r=40
setjobid .b
protein1="../struct/melittin.aam"
r_dielec=6.8
dielec_i=2
chion1=+1
chion2=-2

for ph in 6
do
  nion1=129
  nion2=64
  cleanup
  output="sletmig"
  micro=2000
  simulate > $output
  exit

  micro=200000
  simulate > $output
  cleanup
done

