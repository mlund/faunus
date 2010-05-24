#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
cellradius $cellradius
boxlen     $boxlen
zboxlen    $zboxlen
bjerrum    $bjerrum
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
tion1      $tion1
tion2      $tion2
headarea   $headarea
scratio    $scratio
atomfile   $atomfile
polymer    $polymer
offset     $offset
springconstant  $springconstant
springeqdist    $springeqdist
dp_monomer      $dp_monomer
dp_salt         $dp_salt
pH              $pH
pKa_core        $pKa_core
MUs NA -12.075 CL -12.075
glu3_rotdp      $glu3_rotdp
glu3_zdp        $glu3_zdp
crankshaft_dp   $crankshaft_dp
crankshaft_min  $crankshaft_min
crankshaft_max  $crankshaft_max
crankshaft_num  $crankshaft_num
" > gcglu3pka.conf
}

#--- Input parameters ---
polymer="GLU3coarse.mol2"
macrosteps=10
microsteps=1000
cellradius=125.6
boxlen=250
zboxlen=250
bjerrum=7.12
LJeps=0.0
pH=7.5
nion1=158
nion2=126
tion1="NA"
tion2="CL"
protein="calbindin.aam"
atomfile="faunatoms.dat"
scratio=0
headarea=70
offset=45
springconstant=1.0
springeqdist=3.0
dp_monomer=2
dp_salt=50
pKa_core=10
glu3_rotdp=0.5
glu3_zdp=25
crankshaft_dp=1.0
crankshaft_min=1
crankshaft_max=1
crankshaft_num=40
#
for salt in 167 #100 300 900 
do
for pH in 6.5 #4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 
do
  rm confout.aam 
  rm gcgroup.conf
  suffix="Wall-Test" #result/50mM/glu3b-C4.8-coarse-eff60mM-pH${pH}"
#  microsteps=1000  
  nion1=$[salt+32]
  nion2=$[salt]
  mkinput
  ./glu3gc_wall > ${suffix}.eq
  exit
  microsteps=15000
  mkinput
  ./glu3gc_wall > ${suffix}.out
  mv smeared.aam ${suffix}.aam
done
done

