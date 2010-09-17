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
polymer="GLU3coarse.mol2"
macrosteps=10
microsteps=1000
cellradius=100
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
springeqdist=3.0
dp_monomer=2
dp_salt=50
pKa_core=11.0
#
salt=139
for salt in 139 #100 300 900 
do
for pH in 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 
do
  suffix="result/55mM/glu3b-P11.0-C4.8b-pH${pH}"
  microsteps=3000
  nion1=$[salt+32]
  nion2=$[salt]
  mkinput
  ./glu3c > ${suffix}.eq
  microsteps=15000
  mkinput
  ./glu3c > ${suffix}.out
  mv smeared.aam ${suffix}.aam
done
done
  rm confout.aam 

