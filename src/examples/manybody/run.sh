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
voldp      $voldp
pressure   $pressure
e_r        $e_r
kappa      $kappa
atomfile   $atomfile
" > isobaric.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=10000
boxlen=350
bjerrum=7.12
LJeps=0.0
nion1=10
nion2=20
tion1="NA"
tion2="SO4"
nprot1=30
protein1="neu.aam"
nprot2=0
protein2="mrh4a.aam"
voldp=1
pressure=0.000001
kappa=0.050
e_r=78.67
atomfile="../../../misc/faunatoms.dat"
for LJeps in 0.0 
do
  suffix=".lj-${LJeps}"
  mkinput
  ./isobaric > eq${suffix}
  boxlen=$(less eq${suffix}| grep "length"|awk '{print $ 6}')
  microsteps=10000
#  mkinput
#  ./isobaric > prod${suffix}
  rm confout.aam 
done

