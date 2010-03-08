#!/bin/bash

# Titration of calbindin
# in a 30 mM 1:1 salt solution

function mkinput() {
echo "
macrosteps 10
microsteps $microsteps
cellradius 90
pH         $pH
bjerrum    7.1
nion1      $nion1
nion2      $nion2
tion1      NA
tion2      CL
protein    calbindin.aam
atomfile   ../../../misc/faunatoms.dat
" > pka.conf
}

for salt in 50
do
  for pH in 7 6 5 4 3 2 8 9 10 11 12 
  do
    rm confout.aam
    nion1=$[salt+19]
    nion2=$[salt]
    suffix="calbindin.pH${pH}.salt${salt}"
    microsteps=1000
    mkinput
    ./pka
    microsteps=10000
    mkinput
    ./pka > ${suffix}.out
    mv smeared.aam > ${suffix}.aam
  done
done
