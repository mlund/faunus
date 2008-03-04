#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
dm_dp      $dm_dp
dm_minsep  $dm_minsep
dm_maxsep  $dm_maxsep
cellradius $cellradius
bjerrum    $bjerrum
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
tion1      $tion1
tion2      $tion2
nprot1     $nprot1
protein1   $protein1
nprot2     $nprot2
protein2   $protein2
" > twobody.conf
}

#--- Input parameters ---
macrosteps=10
cellradius=90
bjerrum=7.1
LJeps=0.3
nion1=10
nion2=28
tion1="NA"
tion2="I"
nprot1=1
protein1="lysozyme-ph4.7.aam"
nprot2=1
protein2="lysozyme-ph4.7.aam"
dm_dp=0
dm_minsep=0
dm_maxsep=40

for nion1 in 182
do
  nion2=$[nion1+18]
  for tion2 in "CL"
  do
    rm confout.aam
    suffix="-${tion2}${nion2}-air" #-Sep-${dm_minsep}-${dm_maxsep}"
    microsteps=500
    mkinput
    ./twobody-hof > sletmig
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
