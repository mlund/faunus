#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
boxlen     $boxlen
bjerrum    $bjerrum
LJeps      $LJeps
nprot1     $nprot1
protein1   $protein1
nprot2     $nprot2
protein2   $protein2
voldp      $voldp
pressure   $pressure
e_r        $e_r
kappa      $kappa
max        $max
penalty    $penalty
" > isobaric.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=5000
boxlen=350
bjerrum=7.12
LJeps=0.00
nprot1=30
protein1="neu.aam"
nprot2=0
protein2="mrh4a.aam"
voldp=1.00
pressure=0.000001
kappa=0.050  #  0.147
e_r=78.67
max=1000
penalty=0.1
for LJeps in 0.20 
do
  rm confout.aam
  suffix=".lj-${LJeps}"
  mkinput
  ./isobaric > eq${suffix}
  boxlen=$(less eq${suffix}| grep "length"|awk '{print $ 6}')
  microsteps=50000
  mkinput
  ./isobaric > prod${suffix}
  rm confout.aam 
done

