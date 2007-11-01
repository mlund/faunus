#!/bin/bash

biosim=../../
source ${biosim}/src/pka.sh
#exe=${biosim}/src/pka_dielec

setjobid .cal1
protein1="3ICB-allatom.aam"
cell_r=100
r_dielec=7.9
dielec=80
dielec_i=2
salt=12
springconst=0.2

for springconst in 0 
do
  for ph in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 
  do
    cleanup
    output="out.ph${ph}"
    nion1=$[salt+19]
    nion2=$[salt]          
    micro=2000
    simulate #> .eq
  
    micro=60000
    simulate > $output
    cleanup
  done
done

