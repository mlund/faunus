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
nprot1     $nprot1
protein1   $protein1
nprot2     $nprot2
protein2   $protein2
atomfile   $atomfile
" > binding.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=10000
cellradius=90
bjerrum=7.1
LJeps=0.3
nion1=10
nion2=10
tion1="NA"
tion2="CL"
nprot1=1
protein1="substrate.aam"
nprot2=1
protein2="ligand.aam"
atomfile=../../../misc/faunatoms.dat
for LJeps in 0.3 
do
  microsteps=10000
  mkinput
  ./binding
  exit
  boxlen=$(less eq${suffix}| grep "length"|awk '{print $ 6}')
  microsteps=10000
#  mkinput
#  ./binding > prod${suffix}
  rm confout.aam 
done

