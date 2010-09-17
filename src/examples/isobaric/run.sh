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
nprot3        $nprot3
protein3      $protein3
voldp         $voldp
mrdp          $mrdp
mtdp          $mtdp
ctdp          $ctdp
crdp          $crdp
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
clt           $clt 
clr           $clr
cluster_sep   $cluster_sep
penaltyupdate $penaltyupdate
scalepen      $scalepen 
threshold_g2g 40
threshold_g2p 10000
movie         $movie"> isobaric.conf
}

#--- Input parameters ---
macrosteps=10
microsteps=1
boxlen=2500
LJeps=0.5
nprot1=32
protein1="p1.aam"
nprot2=12224
protein2="sm.aam"
nprot3=0
protein3="sp.aam"
#--- Markov Parameters ---
voldp=0.30        
mrdp=2.00
mtdp=300.00
ctdp=10.00
crdp=1
cluster_sep=8.0
volr=0
tr=0
rr=0
clt=1
clr=0
#--- Interactions and potentials
bjerrum=7.12
pressure=0.000024155
debyelen=42.99
e_r=78.67
atomfile=../../../misc/faunatoms.dat
LJeps=0.50
#--- Sampling constrictions and penaltyfunction
maxlen=1000
minlen=100
binlen=1
penalize="no"
penalty=0.005
scalepen=1.0
penaltyupdate="no"
movie="no"
#--- Variations and execution

#  export OMP_NUM_THREADS=2

#  cp ${suffix}-conf.aam confout.aam
  suffix="test-max"
  microsteps=10  
  mkinput
  ./isobaric > ${suffix}.out
  exit
  boxlen=$(less ${suffix}.out| grep "#   Final      side length"|awk '{print $ 6}')
  microsteps=10000
  mkinput
  ./isobaric > ${suffix}.out
  cp confout.aam ${suffix}-conf.aam
  mv confout.gro ${suffix}-conf.gro
  mv confout.xyz ${suffix}-conf.xyz
  mv coord.xtc ${suffix}-coord.xtc
  mv aggregates.dat ${suffix}agg.dat
  mv length-distribution.dat ${suffix}-Ldist.dat
  mv rdfprot.dat ${suffix}-rdf.dat
  mv rdfprot11.dat ${suffix}-rdf11.dat
  mv rdfprot12.dat ${suffix}-rdf12.dat
  mv rdfprot22.dat ${suffix}-rdf22.dat
  mv systemenergy.dat ${suffix}-energy.dat
  cp penalty.dat ${suffix}-pen.dat
  mv oldpenalty.dat ${suffix}-oldpen.dat
  mv confout.pqr ${suffix}-conf.pqr


