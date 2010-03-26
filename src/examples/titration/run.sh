#!/bin/bash

#--- Write config file ---
function mkinput() {
echo "
pH         $pH
macrosteps $macrosteps
microsteps $microsteps
pH         $pH
bjerrum    $bjerrum
cellradius $cellradius
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
tion1      $tion1
tion2      $tion2
<<<<<<< .mine
protein   $protein
=======
nprot1     $nprot1
protein   $protein
nprot2     $nprot2
protein2   $protein2
voldp      $voldp
pressure   $pressure
>>>>>>> .r544
e_r        $e_r
kappa      $kappa
atomfile   $atomfile
" > pka.conf
}

#--- Input parameters ---
macrosteps=10
<<<<<<< .mine
microsteps=3000        
cellradius=240
=======
microsteps=100
cellradius=90
>>>>>>> .r544
bjerrum=7.12
LJeps=0.2
pH=7.5
<<<<<<< .mine
nion1=61
nion2=36
=======
nion1=42
nion2=23
>>>>>>> .r544
tion1="NA"
tion2="CL"
#protein="lysozyme.aam"
protein="alphaLa1F6S.aam"
atomfile="../../../misc/faunatoms.dat"
e_r=80
pH=7.5
#

  suffix="aLac5mM"
  mkinput
  ./pka > out1
<<<<<<< .mine
  microsteps=10000           
=======
  microsteps=100
>>>>>>> .r544
  mkinput
  ./pka > ${suffix}.out
  mv smeared.aam > ${suffix}.aam
  rm confout.aam 
<<<<<<< .mine
=======

>>>>>>> .r544
