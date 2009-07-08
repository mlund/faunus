#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
macrosteps            $macrosteps
microsteps            $microsteps
boxlen                $boxlen
zboxlen               $zboxlen
bjerrum               $bjerrum
LJeps                 $LJeps
scratio               $scratio
headarea              $headarea
protein               $protein
latticetranslation_dp $latticetranslation_dp
monomer_dp            $monomer_dp
graftpoint_dp         $graftpoint_dp
springconstant        $springconstant
prot_zdp              $prot_zdp
prot_rotdp            $prot_rotdp
debyelen              $debyelen
maxlen                $maxlen
atomfile              $atomfile
mr                    $mr
tr                    $tr
rr                    $rr
place_dp              $place_dp
place                 $place
movie                 $movie"> memp.conf
}

#--- System ---
boxlen=150
zboxlen=200
protein="lysozyme100mM75PH.aam"
scratio=0.1
headarea=70
#--- Markov Parameters ---
macrosteps=10
microsteps=1000
prot_zdp=30.00        
prot_rotdp=3.00
monomer_dp=5.0
graftpoint_dp=5.0
latticetranslation_dp=0.0
mr=1
tr=9
rr=1
#--- Interactions and potentials
LJeps=0.2
bjerrum=7.13
debyelen=9.61
atomfile=../../../misc/faunatoms.dat
springconstant=0.7
springeqdist=5.0
#--- Sampling constrictions and penaltyfunction
maxlen=100
movie="no"
place="no"
place_dp=20
#--- Variations and execution

#  export OMP_NUM_THREADS=2
for debyelen in 9.61 #43.0
do
#  rm conf.aam
  suffix="test" #-K${debyelen}"
  mkinput
  ./memp > ${suffix}
#  microsteps=200000
#  mkinput
#  ./memp > ${suffix}.out
#  mv conf.aam ${suffix}-conf.aam
#  mv conf.pqr ${suffix}-conf.pqr
#  mv coord.xtc ${suffix}-coord.xtc
#  mv rdfw.dat ${suffix}-rdfwp.dat
#  mv end-dist.dat ${suffix}-enddist.dat
#  mv dist.dat ${suffix}-dist.dat
done

