#!/bin/bash

# Script to automatically vary input parameters for "bulk".
# This example is set up to simulate a simple seawater
# solution at different concentrations and temperatures.

# P A R A M E T E R S :

# Ionic radii
r1=1.8 #Na+
r2=1.5 #K+
r3=3.0 #Mg++
r4=2.8 #Ca++
r5=1.7 #Cl-
r6=1.6 #SO4--

# Nr. of pass in the MC simulation
n1f=1000
n2f=200
n3f=10

# Subroutine to create input file 'bulk.inp'
function bulktext {
echo "number of species
6

species  number  hardcore   charge     dp
1           50       $r1      1       42.0 
2            1       $r2      1       42.0 
3            6       $r3      2       11.0
4            1       $r4      2       11.0
5           57       $r5     -1       42.0
6            4       $r6     -2       11.0

ny1  ny2  ny3
${n1}  ${n2}  ${n3}

ink  islu
 ${random}    13   

box 
${box}

dtemp  eps
${T}    ${dielec}

nwins  nwint  nfix
10      25    0"
}

# subroutine to adjust the dielectric constant
function setdielectric {
if     [ $T -eq 298 ] ; then dielec=78.4 ;
  elif [ $T -eq 293 ] ; then dielec=80.2 ;
  elif [ $T -eq 288 ] ; then dielec=82.1 ;
  elif [ $T -eq 283 ] ; then dielec=84.0 ;
  elif [ $T -eq 278 ] ; then dielec=85.9 ;
  elif [ $T -eq 273 ] ; then dielec=87.9 ;
  else "WARNING - DIELECTRIC CONSTANT UNKNOWN !" ;
fi
}

#  L O O P S

# Salinity |  35    25    15    5
# Box len  |  55.9  62.5  74.3  107.4

echo "SIMULATION:"
for T in 298 273                    #loop over temperatures
do
  for box in 55.9 62.5 74.3 107.4   #loop over box sizes (concentration/salinity)
  do
    setdielectric                   #Set the correct dielectric constant
    echo "T=$T box=$box e=$dielec"  #info on screen
    outfile="$box-$T"               #specify output file
    random=1 n1=1000 n2=10   n3=10   bulktext >bulk.inp ; ./bulk.run >.temp    #equilibration run (short)
    random=0 n1=$n1f n2=$n2f n3=$n3f bulktext >bulk.inp ; ./bulk.run >$outfile #final run (longer)
  done
done
