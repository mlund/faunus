#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps 10
microsteps $microsteps
dm_dp      $dm_dp
dm_minsep  $dm_minsep
dm_maxsep  $dm_maxsep
cellradius $cellradius
bjerrum    7.1
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
tion1      $tion1
tion2      $tion2
nprot1     1
protein1   $protein1
nprot2     1
protein2   $protein2
atomfile   ../../../misc/faunatoms.dat
" > twobody.conf
}

#--- Input parameters ---
cellradius=90
LJeps=0.3
nion1=10
nion2=28
tion1="NA"
tion2="I"
protein1="lysozyme-ph4.7.aam"
protein2="lysozyme-ph4.7.aam"
dm_dp=5
dm_minsep=0
dm_maxsep=40

for nion1 in 0
do
  nion2=$[nion1+18]
  for tion2 in "I"
  do
    rm confout.aam
    suffix="-${tion2}${nion2}-air" #-Sep-${dm_minsep}-${dm_maxsep}"
    microsteps=200
    mkinput
    ./twobody-hof #> sletmig

    rm rdfprot.dat cyl.dat
    microsteps=5000
    mkinput
    echo ${suffix}
    ./twobody-hof > out${suffix}

    mv rdfprot.dat rdfprot${suffix}.dat
    mv cyl.dat cyl${suffix}.dat
    #mv rdfsalt.dat rdfsalt${suffix}.dat
    #mv coord.xtc coord${suffix}.xtc
    #mv coord.gro coord${suffix}.gro
  done
done
