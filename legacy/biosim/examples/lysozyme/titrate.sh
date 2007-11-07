#!/bin/bash

biosim=../../
source ${biosim}/src/pka.sh

cell_r=80
setjobid .b
protein1="4LZT.aam"
r_dielec=6.8
dielec_i=2

for ph in 7.6
do
  for salt in 19
  do
    cleanup
    output="out.ph$ph"
    nion1=$[salt+13]
    nion2=$[salt]  
    micro=2000
    simulate
    exit
    
    micro=200000
    simulate > $output
    cleanup
  done
done

