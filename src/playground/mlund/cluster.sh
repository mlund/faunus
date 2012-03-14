#!/bin/bash

faunus=$HOME/faunus/branches/ny
exe=$faunus/src/playground/mlund/mlund-cluster

function mkatoms() {
echo "
Atom  MM    -20     10.0    0.1    1       no
Atom  NEU     0      1.0    0.1    1       no
Atom  La     +3      0.0    0.1    1       no
Atom  Na     +1      0.0    0.1    1       no
" > cluster.atoms
}

function mkstruct() {
echo "1
 MM  0   0.00   0.00   0.00   -20   1  10
" > cluster.aam
}


function mkinput() {
echo "
atomlist               cluster.atoms
loop_macrosteps        10
loop_microsteps        $micro
cuboid_len             232.45     # Box side length Angstrom
sphere_radius          100
cylinder_radius        100
cylinder_len           400
temperature            298     # Kelvin
epsilon_r              78.3    # Water dielectric const
lj_eps                 0.05    # kT
softrep_sigma          5
polymer_N              2
polymer_file           cluster.aam
transrot_transdp       $transdp
transrot_rotdp         0 #0
transrot_clustersize   0 #10
#dh_ionicstrength       0.0075
dh_debyelength         31.1432066

tion1                  La
nion1                  $nion1
dpion1                 100
 
#tion1                  Na
#nion1                  2
#dpion1                 40
" > cluster.input
}

mkatoms
mkstruct

for nion1 in 13
do
  rm -f rdf_p2p.dat
  rm -f state
  micro=10000
  transdp=400
  mkinput
  $exe

  micro=200000000
  transdp=15 #30
  mkinput
  $exe
done
