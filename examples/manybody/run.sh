#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
boxlen     $boxlen
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
" > manybody.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=1000
boxlen=80
bjerrum=7.1
LJeps=0.3
nion1=10
nion2=20
tion1="NA"
tion2="SO4"
nprot1=6
protein1="mrh4a.aam"
nprot2=0
protein2="mrh4a.aam"

for LJeps in 0.2 0.5
do
  suffix=".lj-${LJeps}"
  mkinput
  ./manybody > out${suffix}
  mv rdfprot.dat rdfprot${suffix}.dat
  mv rdfsalt.dat rdfsalt${suffix}.dat
done

