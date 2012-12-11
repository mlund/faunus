#!/bin/bash

#
# Interactive script to titrate a rigid molecular structures
# using an implicit salt description.
#

faunus=$HOME/faunus/
pkaexe=src/tools/titrate_implicit
base=pka

function mktit() {
echo "
Process  H3PO4 H2PO4   2.12    $pH
Process  H2PO4 HPO4    7.21    $pH
Process  HPO4  PO4     12.67   $pH
Process  HASP  ASP     4.0     $pH
Process  HCTR  CTR     2.6     $pH
Process  HGLU  GLU     4.4     $pH
Process  HHIS  HIS     6.3     $pH
Process  HNTR  NTR     7.5     $pH
Process  HTYR  TYR     9.6     $pH
Process  HLYS  LYS     10.4    $pH
Process  HCYS  CYS     10.8    $pH
Process  HARG  ARG     12.0    $pH
" > $base.titration
}

# ----------------------------------
#   GENERATE ATOM PARAMETERS
# ----------------------------------
function mkatoms() {
echo "
Atom  Na      +1     2.0    0.1    1       no
Atom  Cl      -1     2.0    0.1    1       no
Atom  I       -1     2.0    0.1    1       no
Atom  SCN     -1     2.0    0.1    1       no

Atom  ASP     -1     3.6    0.1    110     no
Atom  HASP     0     3.6    0.1    110     no
Atom  CTR     -1     2.0    0.1    16      no
Atom  HCTR     0     2.0    0.1    16      no
Atom  GLU     -1     3.8    0.1    122     no
Atom  HGLU     0     3.8    0.1    122     no
Atom  HIS      0     3.9    0.1    130     no
Atom  HHIS     1     3.9    0.1    130     no
Atom  NTR      0     2.0    0.1    14      no
Atom  HNTR     1     2.0    0.1    14      no
Atom  TYR     -1     4.1    0.1    154     no
Atom  HTYR     0     4.1    0.1    154     no
Atom  LYS      0     3.7    0.1    116     no
Atom  HLYS     1     3.7    0.1    116     no
Atom  CYS     -1     3.6    0.1    103     no 
Atom  HCYS     0     3.6    0.1    103     no 
Atom  ARG      0     4.0    0.1    144     no
Atom  HARG     1     4.0    0.1    144     no

Atom  ALA      0     3.1    0.1    66      yes
Atom  ILE      0     3.6    0.1    102     yes
Atom  LEU      0     3.6    0.1    102     yes
Atom  MET      0     3.8    0.1    122     yes
Atom  PHE      0     3.9    0.1    138     yes
Atom  PRO      0     3.4    0.1    90      yes
Atom  TRP      0     4.3    0.1    176     yes
Atom  VAL      0     3.4    0.1    90      yes
Atom  SER      0     3.3    0.1    82      no
Atom  THR      0     3.5    0.1    94      no
Atom  ASN      0     3.6    0.1    108     no
Atom  GLN      0     3.8    0.1    120     no
Atom  GLY      0     2.9    0.1    54      no
" > $base.atoms
}

function mkinput() {
echo "
atomlist               $base.atoms
eq_processfile         $base.titration
loop_macrosteps        10
loop_microsteps        $micro
temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const
lj_eps                 0.0     # kT
dh_ionicstrength       $ionicstr   # mol/l
sphere_radius          500     # Angstroms
molecule               $struct
" > input
}

function userinput() {
read -p "Path to faunus directory ["$faunus"]: " i
faunus=${i:-$faunus}
if [ ! -d $faunus/src/tools ]; then
  echo "Invalid faunus directory." ; exit 1
fi
if [ ! -f $faunus/$pkaexe ]; then
  echo "Could not find titration program. Did you remember to make it?" ; exit 1
fi

read -p "AAM structure file to titrate: " struct
if [ ! -f $struct ]; then
  echo "Structure file does not exist." ; exit 1
fi

ionicstr=0.020
read -p "Ionic strength ["$ionicstr" mol/l]: " i
ionicstr=${i:-$ionicstr}

pHrange="2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0"
read -p "pH range [$pHrange]: " i
pHrange=${i:-$pHrange}

echo
echo "Summary of simulation:"
echo "  Titration program        = $faunus/$pkaexe"
echo "  Molecular structure file = $struct"
echo "  Ionic strength           = $ionicstr mol/l"
echo "  pH range                 = $pHrange"
}

userinput

echo
read -p "Start simulations? "
echo

echo "pH    Z        C       mu"
for pH in $pHrange; do
  echo -n "$pH "
  mkatoms
  mktit

  # Equilibration run
  rm -f state
  micro=1000
  mkinput
  $faunus/$pkaexe > eq

  # Production run
  micro=10000
  prefix="pH$pH-salt$ionicstr"
  mkinput
  $faunus/$pkaexe > $prefix.out
  mv confout.pqr $prefix.pqr
  Z=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $2}'`
  C=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $3}'`
  mu=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $4}'`
  echo $Z $C $mu
done

