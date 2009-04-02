#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps    $macrosteps
microsteps    $microsteps
boxlen        $boxlen
bjerrum       $bjerrum
LJeps         $LJeps
nprot1        $nprot1
protein1      $protein1
nprot2        $nprot2
protein2      $protein2
voldp         $voldp
mrdp          $mrdp
mtdp          $mtdp
pressure      $pressure
e_r           $e_r
debyelen      $debyelen
maxlen        $maxlen
minlen        $minlen
binlen        $binlen
penalize      $penalize
penalty       $penalty
atomfile      $atomfile
volr          $volr
tr            $tr
rr            $rr
penaltyupdate $penaltyupdate
scalepen      $scalepen
" > isobaric.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=5000
boxlen=350
LJeps=0.00
nprot1=30
protein1="neu.aam"
nprot2=0
protein2="mrh4a.aam"
#--- Markov Parameters ---
voldp=1.00        
mrdp=2.00
mtdp=10.00
#--- Interactions and potentials
bjerrum=7.12
pressure=0.000001
debyelen=20
e_r=78.67
atomfile=../../../misc/faunatoms.dat
LJeps=0.2
#--- Sampling constrictions and penaltyfunction
maxlen=1000
minlen=100
binlen=1
penalize="yes"
penalty=0.01
scalepen=0.5
penaltyupdate="yes"
#--- Variations and execution

for LJeps in 0.20 
do
#  rm confout.aam
#  suffix=".lj-${LJeps}"
  mkinput
#  ./isobaric > eq${suffix}
#  boxlen=$(less eq${suffix}| grep "length"|awk '{print $ 6}')
#  microsteps=50000
#  mkinput
#  ./isobaric > prod${suffix}
#  rm confout.aam 
done

