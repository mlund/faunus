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
bjerrum    7.1
atomfile   faunatoms.dat
nion1      $nion1
nion2      $nion2
nion3      $nion3
nion4      $nion4
nion5      $nion5
nion6      $nion6
tion1      $tion1
tion2      $tion2
tion3      $tion3
tion4      $tion4
tion5      $tion5
tion6      $tion6
dp_salt    $dp_salt
" > widom.conf
}

microsteps=10
tion1="NA"
tion2="CL"
tion3="K"
tion4="MG"
tion5="CA"
tion6="CL"
dp_salt=10
exe="./widom_cube"
for boxlen in 367.9 255 176.9 122.6 85 
do
    rm -f widom.aam
    nion1=0
    nion2=0
    nion3=0
    nion4=0
    nion5=50
    nion6=100
#    boxlen=129.8
    microsteps=1000
    outfile="anything1.out"
    geninput
    $exe > t1b 
#exit
    microsteps=1000
    outfile="${boxlen}-2-1.out"
    geninput
    $exe > $outfile
done
