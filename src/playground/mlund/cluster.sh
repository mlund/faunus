#!/bin/bash

faunus=$HOME/faunus/branches/ny
exe=$faunus/src/playground/mlund/mlund-cluster

function mkatoms() {
echo "
Atom  MM    -15     20.0    0.1    1       no
Atom  NEU     0      1.0    0.1    1       no
Atom  La     +3      2.0    0.1    1       no
Atom  Na     +0      2.0    0.1    1       no
" > cluster.atoms
}

function mkstruct() {
echo "1
 MM  0   0.00   0.00   0.00   -15   1  20
" > cluster.aam
}


function mkinput() {
echo "
atomlist               cluster.atoms
loop_macrosteps        10
loop_microsteps        $micro
cuboid_len             231     # Box side length Angstrom
sphere_radius          100
cylinder_radius        100
cylinder_len           300
temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const
lj_eps                 0.05    # kT
softrep_sigma          5
polymer_N              2
polymer_file           cluster.aam
transrot_transdp       20
transrot_rotdp         0
transrot_clustersize   4
#dh_ionicstrength       0.0075
dh_debyelength         31.1432066

tion1                  La
nion1                  $nion1
dpion1                 30
 
#tion1                  Na
#nion1                  2
#dpion1                 40
" > cluster.input
}

mkatoms
mkstruct

for nion1 in 38
do
  #rm state
  micro=5000
  mkinput
  $exe

  micro=4000000
  mkinput
  $exe
done
