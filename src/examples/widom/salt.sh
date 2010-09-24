#!/bin/bash

# Script to run a electrolyte at different concentrations.
# Note: This is just a script that generates input and executes
#       the Monte Carlo program 

function geninput() {
echo "
macrosteps 10 
microsteps $microsteps
boxlen     $boxlen
cellradius 100
bjerrum    $bjerrum
atomfile   saltlabfaunatoms.dat
nion1      $nion1
nion2      $nion2
nion3      $nion3
nion4      $nion4
nion5      $nion5
nion6      $nion6
nion7      $nion7
nion8      $nion8
tion1      $tion1
tion2      $tion2
tion3      $tion3
tion4      $tion4
tion5      $tion5
tion6      $tion6
tion7      $tion7
tion8      $tion8
dp_salt    $dp_salt
" > saltsolution.conf
}

microsteps=100
tion1="NA"
tion2="CL"
tion3="CA"
tion4="MG"
tion5="SO4"
tion6="K"
tion7="PION"
tion8="NION"
dp_salt=100
bjerrum=
exe="./saltsolution"
for boxlen in 1020 
do
    rm -f widom.aam
    nion1=0
    nion2=0
    nion3=0
    nion4=0
    nion5=0
    nion6=0
    nion7=100
    nion8=100
    outfile="name-of-outputfile"
    geninput
    $exe > $outfile
    microsteps=10000
    geninput
    $exe > $outfile
done
