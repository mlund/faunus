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
bjerrum    7.13
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
cluster_dp $cluster_dp
" > widom.conf
}

microsteps=10
tion2="NA"
tion3="CL"
tion1="COL"
tion4="MG"
tion5="CA"
tion6="CL"
dp_salt=100
cluster_dp=25 
exe="./saltcluster"
for boxlen in 1020 #998.7 1258.2
do
#    rm -f widom.aam
    nion3=12224
    nion2=0
    nion1=32
    nion4=0
    nion5=0
    nion6=0
#    boxlen=129.8
    microsteps=10  
    outfile="JANNE"
    geninput
    $exe > $outfile
exit
    microsteps=1000
    #outfile="${boxlen}.out"
    geninput
    $exe > $outfile
done
