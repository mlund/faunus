#!/bin/bash

# Script to run a 1:1 electrolyte at different concentrations.
# Note: This would be more convenient by using the widom.py script!

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
tion1      $tion1
tion2      $tion2
dp_salt    $dp_salt
testsuite_stable yes
" > widom.conf
}

microsteps=10
tion1="NA"
tion2="CL"
dp_salt=20
exe="./widom_cube"

for conc in 0.001 0.005 0.01 0.05 0.1 0.2 0.4 0.8 1.0 1.5;
do
  for salt in 50;
  do
    rm -f widom.aam
    nion1=$salt
    nion2=$salt
    boxlen=`python -c "print '%.4f' % (pow(${salt}/($conc)/6.022e23*1e27,1/3.))"`
    cellradius=`python -c "from math import pi; print '%.4f' % (pow( 3/4./pi,1/3.)*$boxlen)"`
    microsteps=1000
    geninput
    $exe > /dev/null

    outfile="$tion1$tion2-${conc}M.out"
    microsteps=10000
    geninput
    $exe > $outfile
    gamma=`cat $outfile | grep --binary-files=text "Mean activity" | gawk '{print $6}'`
    echo $conc $gamma
  done
done

