#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
cellradius $cellradius
bjerrum    $bjerrum
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
tion1      $tion1
tion2      $tion2
atomfile   $atomfile
polymer    $polymer
springconstant  $springconstant
springeqdist    $springeqdist
dp_monomer      $dp_monomer
dp_salt         $dp_salt
pH              $pH
pKa_core        $pKa_core
" > pka.conf
}

#--- Input parameters ---
polymer="GLU3simptwostate.mol2"
macrosteps=10
microsteps=1000
cellradius=200
bjerrum=7.12
LJeps=0.0
pH=7.5
nion1=158
nion2=126
tion1="NA"
tion2="CL"
protein="calbindin.aam"
atomfile="faunatoms.dat"
springconstant=1.0
springeqdist=1.5
dp_monomer=2
dp_salt=50
pKa_core=10
#
for salt in 34 100 300 
do
for pH in 4.0 4.5 5.0 5.5 6.0 
do
  suffix="result/50mM/glu3b-salt${salt}-pH${pH}"
  microsteps=2000
  nion1=$[salt+32]
  nion2=$[salt]
  mkinput
  ./glu3b > ${suffix}.eq
  microsteps=10000
  mkinput
  ./glu3b > ${suffix}.out
  mv smeared.aam ${suffix}.aam
done
done
  rm confout.aam 

